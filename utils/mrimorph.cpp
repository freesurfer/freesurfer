/**
 * @brief utilities for 3d morph of one volume into another
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

#define USE_ITERATIVE_AVERAGING 0
#define SCALE_INVARIANT 0
#define USE_ORIGINAL_PROPERTIES 0

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "histo.h"
#include "icosahedron.h"
#include "macros.h"
#include "matrix.h"
#include "mri.h"
#include "mri_circulars.h"
#include "mrimorph.h"
#include "mrinorm.h"
#include "mrishash.h"
#include "mrisurf.h"
#include "numerics.h"
#include "proto.h"
#include "region.h"
#include "voxlist.h"

extern const char *Progname;

#define MN_SUB(mns1, mns2, v) V3_LOAD(v, mns1->x - mns2->x, mns1->y - mns2->y, mns1->z - mns2->z)
#define MNP_SUB(mns1, mns2, v) V3_LOAD(v, mns1->ox - mns2->ox, mns1->oy - mns2->oy, mns1->oz - mns2->oz)

#define EXP_K 10.0
#define BACKGROUND_VAL 10
#define MAX_EXP 200

int MRIsampleReferenceWeightingGradient(MRI *mri, int x, int y, int z, double *pdx, double *pdy, double *pdz);
double MRIsampleReferenceWeighting(MRI *mri, int x, int y, int z);

static void computeRigidAlignmentGradient(float *p, float *g);
static float computeRigidAlignmentErrorFunctional(float *p);
static float computeEMAlignmentErrorFunctional(float *p);
static void computeEMAlignmentGradient(float *p, float *g);
static int m3dPositionBorderNodes(MORPH_3D *m3d);
static float m3dNodeAverageExpansion(MORPH_3D *m3d, int i, int j, int k);

static int m3dblurDx(MORPH_3D *m3d, float sigma);
static int m3dTranslate(MORPH_3D *m3d, float dx, float dy, float dz);
static int m3RecomputeTranslation(MORPH_3D *m3d, MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms);
static int m3dMorphSkull(MORPH_3D *m3d, MRI_SURFACE *mris_in_skull, MRI_SURFACE *mris_ref_skull, MRI *mri);
static int m3dcheck(MORPH_3D *m3d);
static MRI *find_midline(MRI *mri_src, MRI *mri_thresh, float *px);
static int find_spinal_fusion(MRI *mri_thresh, float *px, float *py, float *pz);
static int mriLinearAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MP *parms);
static int mriQuasiNewtonLinearAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MP *parms);
static int mriQuasiNewtonEMAlignPyramidLevel(MRI *mri_in, GCA *gca, MP *parms);
static int m3dAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MRI *mri_ref_blur, MP *parms, MORPH_3D *m3d);
static int mriOrthonormalizeTransform(MATRIX *m_L);
static double mriIntensityRMS(MRI *mri_in, MRI *mri_ref, LTA *lta, double l_intensity, NECK_PARMS *np);
static double mriIntensitySSE(MRI *mri_in, MRI *mri_ref, MATRIX *m_L);
static int mriWriteImageView(MRI *mri, const char *base_name, int target_size, int view, int slice);
static int writeSnapshot(MRI *mri, MORPH_PARMS *parms, int n);
static int write3DSnapshot(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MORPH_3D *m3d, int n);
static double ltaGradientStep(
    MRI *mri_in, MRI *mri_ref, LTA *lta, double dt, double momentum, NECK_PARMS *np, double tmul);
static int openLogFile(MORPH_PARMS *parms);
static int logIntegration(MORPH_PARMS *parms, int n, double rms);
static int log3DIntegration(MORPH_PARMS *parms,
                            MORPH_3D *m3d,
                            int n,
                            double dt,
                            double rms,
                            double intensity_rms,
                            double distance_rms,
                            double area_rms);
static int finitep(float f);
static int mriNormalizeStds(MRI *mri);

// the following two functions are added for gradient-descent type of linear alignment
static int ComputeStepTransform(
    VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, float cx, float cy, float cz, MATRIX *Minc);
static int farid_align(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L);

static int debug_x = -1;
static int debug_y = -1;
static int debug_z = -1;

#define MAX_LABELS 1000

void showTran(float p[])
{
  printf("after pass:transform: ( %.2f, %.2f, %.2f, %.2f)\n", p[1], p[2], p[3], p[4]);
  printf("                      ( %.2f, %.2f, %.2f, %.2f)\n", p[5], p[6], p[7], p[8]);
  printf("                      ( %.2f, %.2f, %.2f, %.2f)\n", p[9], p[10], p[11], p[12]);
}

int MatrixPrintHires(FILE *fp, MATRIX *mat)
{
  int row, col, rows, cols;

  if (fp == NULL) {
    fp = stdout;
    ErrorPrintf(ERROR_BADPARM, "MatrixPrint: fp = NULL!");
  }
  if (mat == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrintHires: mat = NULL!"));

  rows = mat->rows;
  cols = mat->cols;

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, "% 2.5f", mat->rptr[row][col]);
          break;
        case MATRIX_COMPLEX:
          fprintf(fp, "% 2.5f + % 2.5f i", MATRIX_CELT_REAL(mat, row, col), MATRIX_CELT_IMAG(mat, row, col));
          break;
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrintHires: unknown type %d\n", mat->type));
      }
      if (col < cols) fprintf(fp, "  ");
    }
    fprintf(fp, ";\n");
  }
  fflush(stdout);

  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
HISTOGRAM *MRIhorizontalHistogram(MRI *mri, int thresh_low, int thresh_hi)
{
  HISTOGRAM *histo;
  int x, y, z, width, height, depth, npix, val;

  histo = HISTOalloc(mri->height);

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  for (y = 0; y < height; y++) {
    if (y == Gdiag_no) DiagBreak();
    for (npix = z = 0; z < depth; z++) {
      for (x = 0; x < width; x++) {
        val = MRIvox(mri, x, y, z);
        if (val >= thresh_low && val <= thresh_hi) npix++;
      }
    }
    histo->counts[y] = npix;
    histo->bins[y] = y;
  }
  return (histo);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
HISTOGRAM *MRIhorizontalBoundingBoxHistogram(MRI *mri, int thresh)
{
  HISTOGRAM *histo;
  int x, y, z, width, height, depth, xmin, xmax, zmin, zmax;

  histo = HISTOalloc(mri->height);

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  histo->nbins = 0;
  for (y = 0; y < height; y++) {
    zmin = depth;
    xmin = width;
    zmax = 0;
    xmax = 0;
    for (z = 0; z < depth; z++) {
      for (x = 0; x < width; x++) {
        if (MRIvox(mri, x, y, z) >= thresh) {
          if (x < xmin) xmin = x;
          if (z < zmin) zmin = z;
          if (x > xmax) xmax = x;
          if (z > zmax) zmax = z;
        }
      }
    }
    if (zmin == depth) zmin = -1;
    histo->counts[y] = zmin;
    if (histo->counts[y] >= 0 && y >= histo->nbins) histo->nbins = y;
    histo->bins[y] = y;
  }
  return (histo);
}
/*-----------------------------------------------------
\fn int MRIcountAboveThreshold(MRI *mri, int thresh)
\brief Returns the number of voxels above the given threshold
------------------------------------------------------*/
int MRIcountAboveThreshold(MRI *mri, double thresh)
{
  int y,  width, height, depth, count;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  count = 0;
#ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : count)
#endif
  for (y = 0; y < height; y++) {
    int x,z;
    double val;
    for (z = 0; z < depth; z++) {
      for (x = 0; x < width; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (val >= thresh) count++;
      }
    }
  }
  return (count);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIlabel(MRI *mri_src, MRI *mri_dst, int *pnlabels)
{
  int l, x, y, z, width, height, depth, *pxi, *pyi, *pzi, xi, yi, zi, total_on, total_labeled, labeled, xk, yk, zk;
  BUFTYPE *psrc, *plabel;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  total_on = MRIcountAboveThreshold(mri_src, 1);
  /*  fprintf(stdout, "labeling %d voxels.\n", total_on) ;*/
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  total_labeled = 0;
  for (l = 1; total_labeled < total_on; l++) {
    /* first find an unlabeled 'on' voxel and label it to start a new seed */
    for (labeled = z = 0; !labeled && z < depth; z++) {
      for (y = 0; !labeled && y < height; y++) {
        psrc = &MRIvox(mri_src, 0, y, z);
        plabel = &MRIvox(mri_dst, 0, y, z);
        for (x = 0; !labeled && x < width; x++) {
          if (*psrc++ && !*plabel) /* unlabeled on voxel */
          {
            *plabel = l;
            labeled = 1;
            total_labeled++;
          }
          else
            plabel++;
        }
      }
    }
    do {
      labeled = 0;
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          plabel = &MRIvox(mri_dst, 0, y, z);
          for (x = 0; x < width; x++) {
            if (*plabel++ == l) /* current label */
            {
              for (zk = -1; zk <= 1; zk++) {
                zi = pzi[z + zk];
                for (yk = -1; yk <= 1; yk++) {
                  yi = pyi[y + yk];
                  for (xk = -1; xk <= 1; xk++) {
                    xi = pxi[x + xk];
                    if (MRIvox(mri_src, xi, yi, zi) && !MRIvox(mri_dst, xi, yi, zi)) {
                      MRIvox(mri_dst, xi, yi, zi) = l;
                      total_labeled++;
                      labeled++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } while (labeled > 0);
    /*    fprintf(stdout, "\r%d labels found     ", l) ;*/
  }
  /*  fprintf(stdout, "\n") ;*/
  *pnlabels = l;
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIlabelBoundingBoxes(MRI *mri, MRI_REGION *bboxes, int nlabels)
{
  int l, x, y, z, width, height, depth, x0, x1, y0, y1, z0, z1;
  BUFTYPE *plabel;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  for (l = 1; l < nlabels; l++) {
    x0 = width;
    y0 = height;
    z0 = depth;
    x1 = y1 = z1 = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        plabel = &MRIvox(mri, 0, y, z);
        for (x = 0; x < width; x++) {
          if (*plabel++ == l) {
            if (x < x0) x0 = x;
            if (x > x1) x1 = x;
            if (y < y0) y0 = y;
            if (y > y1) y1 = y;
            if (z < z0) z0 = z;
            if (z > z1) z1 = z;
          }
        }
      }
    }
    x0 *= mri->thick;
    y0 *= mri->thick;
    z0 *= mri->thick;
    x1 *= mri->thick;
    y1 *= mri->thick;
    z1 *= mri->thick;
    bboxes[l].x = x0;
    bboxes[l].dx = x1 - x0 + 1;
    bboxes[l].y = y0;
    bboxes[l].dy = y1 - y0 + 1;
    bboxes[l].z = z0;
    bboxes[l].dz = z1 - z0 + 1;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIeraseOtherLabels(MRI *mri_src, MRI *mri_dst, int label)
{
  int x, y, z, width, height, depth;
  BUFTYPE *psrc, *pdst;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        if (*psrc++ != label)
          *pdst++ = 0;
        else
          *pdst++ = label;
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
int MRIeraseLabel(MRI *mri_src, MRI *mri_dst, int label)
{
  int x, y, z, width, height, depth;
  BUFTYPE *psrc, *pdst, l;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        l = *psrc++;
        if (l == label)
          *pdst++ = 0;
        else
          *pdst++ = l;
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
int MRIfindHorizontalLabelLimits(MRI *mri, int label, int *xmins, int *xmaxs)
{
  int x, y, z, width, height, depth, ymin;
  BUFTYPE *psrc;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  for (ymin = -1, y = 0; y < height; y++) {
    if (ymin >= 0) {
      xmins[y - ymin] = width;
      xmaxs[y - ymin] = -1;
    }
    for (z = 0; z < depth; z++) {
      psrc = &MRIvox(mri, 0, y, z);
      for (x = 0; x < width; x++) {
        if (*psrc++ == label) {
          if (ymin < 0) {
            xmins[0] = xmaxs[0] = x;
            ymin = y;
          }
          if (x > xmaxs[y - ymin]) xmaxs[y - ymin] = x;
          if (x < xmins[y - ymin]) xmins[y - ymin] = x;
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
#define MIN_SPINAL_ASPECT .75 /* will eliminate corpus callosum */

/* provide some very broad bounds on where it can be and what size it is */
#define MIN_SPINAL_X 70
#define MAX_SPINAL_X 160
#define MIN_SPINAL_Y 85
#define MAX_SPINAL_Y 260 /* can extend to bottom of image */
#define MIN_SPINAL_AREA 190

#define M3D_VERSION 1

#define DISTANCE_BELOW_SPINAL_FUSION 40.0

MRI *MRIfindNeck(MRI *mri_src, MRI *mri_dst, int thresh_low, int thresh_hi, MP *parms, int dir, NECK_PARMS *np)
{
  MRI *mri_label, *mri_thresh, *mri_midline, *mri_rot;
  int nlabels, l, max_dy, spinal_label, xmins[256], xmaxs[256], y;
  MRI_REGION bboxes[MAX_LABELS], *spinal_bbox, in_bbox;
  float n, len, avg_dy, x0, y0, z0, xmid, avg_dz, aspect, areas[MAX_LABELS], thick;
  MATRIX *m_rot, *m_trans, *m_L, *m_inv;
  VECTOR *v;
  double xa;
  static int callno = 0;

  if (parms && !parms->lta) parms->lta = LTAalloc(1, NULL);
  v = VectorAlloc(4, MATRIX_REAL);
  thick = mri_src->thick;
  callno++; /* for diagnostics */
  MRIboundingBox(mri_src, 80, &in_bbox);

  /* first find a sagittal slice at the midline */
  mri_thresh = MRIthresholdRangeInto(mri_src, NULL, thresh_low, thresh_hi);
  MRIopen(mri_thresh, mri_thresh);
  mri_midline = find_midline(mri_src, mri_thresh, &xmid);
  if (!mri_dst) mri_dst = MRIcopy(mri_src, NULL);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "midline%d.mgh", callno);
    MRIwrite(mri_midline, fname);
  }
  mri_label = MRIlabel(mri_midline, NULL, &nlabels);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "label%d.mgh", callno);
    MRIwrite(mri_label, fname);
  }

  /* find the label with the biggest y extent, it will be spinal cord */
  /* note that bboxes are in mm which may NOT be the same as voxels */
  MRIlabelBoundingBoxes(mri_label, bboxes, nlabels);
  MRIlabelAreas(mri_label, areas, nlabels);
  spinal_label = -1;
  for (l = 1, max_dy = 0; l < nlabels; l++) {
    aspect = (float)bboxes[l].dy / (float)bboxes[l].dx;
    if (aspect > MIN_SPINAL_ASPECT && bboxes[l].x >= MIN_SPINAL_X && bboxes[l].x + bboxes[l].dx - 1 <= MAX_SPINAL_X &&
        bboxes[l].y >= MIN_SPINAL_Y && bboxes[l].y + bboxes[l].dy - 1 <= MAX_SPINAL_Y && areas[l] >= MIN_SPINAL_AREA &&
        bboxes[l].dy >= max_dy) {
      /* candidate for spinal chord */
      spinal_label = l;
      max_dy = bboxes[l].dy;
    }
  }

  if (spinal_label < 0) ErrorReturn(NULL, (ERROR_BADPARM, "MRIfindNeck: could not find spinal chord"));

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,
            "spinal label %d, area = %2.0f, "
            "box = (%d, %d) --> (%d, %d)\n",
            spinal_label,
            areas[spinal_label],
            bboxes[spinal_label].x,
            bboxes[spinal_label].y,
            bboxes[spinal_label].x + bboxes[spinal_label].dx - 1,
            bboxes[spinal_label].y + bboxes[spinal_label].dy - 1);
  MRIeraseOtherLabels(mri_label, mri_label, spinal_label);
  MRIfindHorizontalLabelLimits(mri_label, spinal_label, xmins, xmaxs);
  spinal_bbox = &bboxes[spinal_label];

  /* calculate avg dz, skip 1st 6 mm or so which typically contain the pons */
  for (y = 6 / thick, n = avg_dz = 0.0f; y < spinal_bbox->dy / thick - 1; y++) {
    avg_dz += xmins[y] - xmins[y - 1];
    n++;
  }
  avg_dz /= n;
  len = sqrt(1 + avg_dz * avg_dz);
  avg_dz /= len;
  avg_dy = 1.0f / len;

  /* now build a rotation matrix which will orient the brain so that
     the spinal cord is vertical.
  */
  x0 = (float)(spinal_bbox->x + spinal_bbox->dx - 1) / 2;
  y0 = (float)(spinal_bbox->y + spinal_bbox->dy - 1) / 2;
  z0 = (float)(spinal_bbox->z + spinal_bbox->dz - 1) / 2;
  xa = acos(avg_dy);
  if (avg_dz > 0.0f) xa *= -1.0f;
  m_trans = MatrixIdentity(4, NULL);
  *MATRIX_RELT(m_trans, 1, 4) = -x0 * mri_src->thick;
  *MATRIX_RELT(m_trans, 2, 4) = -y0 * mri_src->thick;
  *MATRIX_RELT(m_trans, 3, 4) = -z0 * mri_src->thick;
  m_rot = MatrixAllocRotation(4, xa, X_ROTATION);
  m_L = MatrixMultiply(m_rot, m_trans, NULL);
  *MATRIX_RELT(m_trans, 1, 4) = x0 * mri_src->thick;
  *MATRIX_RELT(m_trans, 2, 4) = y0 * mri_src->thick;
  *MATRIX_RELT(m_trans, 3, 4) = z0 * mri_src->thick;
  MatrixMultiply(m_trans, m_L, m_L);

  mri_rot = MRIlinearTransform(mri_src, NULL, m_L);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "rot%d.mgh", callno);
    MRIwrite(mri_rot, fname);
  }
  MRIthresholdRangeInto(mri_rot, mri_thresh, thresh_low, thresh_hi);
  MRIopen(mri_thresh, mri_thresh);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "thresh%d.mgh", callno);
    MRIwrite(mri_thresh, fname);
  }
  MRIboundingBox(mri_rot, 80, &in_bbox);

  /* find coordinates for cutting plane in rotated system */
  find_spinal_fusion(mri_thresh, &x0, &y0, &z0);

  /* transform them into original voxel coordinates */
  VECTOR_ELT(v, 4) = 1.0 / mri_thresh->thick;
  V3_X(v) = x0;
  V3_Y(v) = y0;
  V3_Z(v) = z0;
  m_inv = MatrixSVDInverse(m_L, NULL);
  MatrixMultiply(m_inv, v, v);
  x0 = V3_X(v);
  y0 = V3_Y(v);
  z0 = V3_Z(v);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "unerased%d.mgh", callno);
    fprintf(stdout, "writing volume before erasure to %s\n", fname);
    MRIwrite(mri_dst, fname);
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,
            "erasing neck around (%2.1f, %2.1f, %2.1f), "
            "n = (%2.2f,%2.2f, %2.2f)\n",
            thick * x0,
            thick * y0,
            thick * z0,
            0.0f,
            avg_dy,
            avg_dz);

  if (np) {
    x0 *= thick;
    y0 *= thick;
    z0 *= thick;
    np->neck_x0 = x0;
    np->neck_y0 = y0;
    np->neck_z0 = z0;
    np->neck_dx = 0.0;
    np->neck_dy = avg_dy;
    np->neck_dz = avg_dz;
  }
  if (parms != NULL) {
    MATRIX *parms_m_L = parms->lta->xforms[0].m_L;

    /* plane may be inclined, so have to continue off end of image */

    /*    MRIeraseNeck(mri_dst, np) ;*/
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "erased%d.mgh", callno);
      fprintf(stdout, "writing volume after erasure to %s\n", fname);
      MRIwrite(mri_dst, fname);
    }

    /*
      dir > 0 means that this is the template image, so that the
      transformation is forward while dir < 0 means that this is an
      input image, and the transformation is inverted.
    */
    if (dir > 0)
      MatrixMultiply(parms_m_L, m_L, parms_m_L);
    else
      MatrixMultiply(parms_m_L, m_inv, parms_m_L);
  }
  MatrixFree(&m_trans);
  MatrixFree(&m_rot);
  MatrixFree(&m_L);
  VectorFree(&v);
  MRIfree(&mri_label);
  MRIfree(&mri_thresh);
  MRIfree(&mri_rot);
  MatrixFree(&m_inv);
  fflush(stdout);
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#ifdef WSIZE
#undef WSIZE
#endif
#define WSIZE 21
#define WHALF ((WSIZE - 1) / 2)
static MRI *find_midline(MRI *mri_src, MRI *mri_thresh, float *px)
{
  MRI *mri_dst;
  MRI_REGION bbox;
  int count, min_count, x, x0, xmid;

  MRIboundingBox(mri_src, 80, &bbox);

  xmid = x0 = bbox.x + bbox.dx / 2;
  min_count = mri_thresh->width * mri_thresh->height;
  for (mri_dst = NULL, x = x0 - WHALF; x <= x0 + WHALF; x++) {
    mri_dst = MRIextractPlane(mri_thresh, mri_dst, MRI_SAGITTAL, x);
    count = MRIcountAboveThreshold(mri_dst, 2);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      fprintf(stdout, "%d: count = %d (min=%d @ %d)\n", x, count, min_count, xmid);
      sprintf(fname, "mid%d.mgh", x);
      MRIwrite(mri_dst, fname);
    }

    if (count < min_count) {
      min_count = count;
      xmid = x;
      /*      fprintf(stdout, "new min %d found at %d\n", min_count, xmid) ;*/
    }
  }

  MRIextractPlane(mri_thresh, mri_dst, MRI_SAGITTAL, xmid);
  *px = (float)xmid;
  fflush(stdout);
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Compute the centroid of a label (in mm coordinates)
------------------------------------------------------*/
int MRIlabelCentroid(MRI *mri, int target_l, float *px, float *py, float *pz)
{
  int l, x, y, z, width, height, depth, nvox;
  BUFTYPE *plabel;
  float thick, x0, y0, z0;

  thick = mri->thick;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  for (x0 = y0 = z0 = 0.0f, nvox = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      plabel = &MRIvox(mri, 0, y, z);
      for (x = 0; x < width; x++) {
        l = *plabel++;
        if (l == target_l) {
          x0 += (float)x * thick;
          y0 += (float)y * thick;
          z0 += (float)z * thick;
          nvox++;
        }
      }
    }
  }

  if (nvox) {
    *px = x0 / nvox;
    *py = y0 / nvox;
    *pz = z0 / nvox;
  }
  else
    *px = *py = *pz = 0.0f;

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the area of all labels (in mm^3)
------------------------------------------------------*/
int MRIlabelAreas(MRI *mri, float *areas, int nlabels)
{
  int l, x, y, z, width, height, depth;
  BUFTYPE *plabel;
  float voxel_area;

  voxel_area = mri->xsize * mri->ysize * mri->zsize;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  for (l = 0; l < nlabels; l++) areas[l] = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      plabel = &MRIvox(mri, 0, y, z);
      for (x = 0; x < width; x++) {
        l = *plabel++;
        areas[l] += voxel_area;
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
#undef MIN_SPINAL_X
#undef MAX_SPINAL_X
#undef MIN_SPINAL_Y
#undef MAX_SPINAL_Y
#undef MIN_SPINAL_ASPECT
#undef MAX_SPINAL_ASPECT
#undef MIN_SPINAL_AREA

#define MIN_SPINAL_AREA 600.0f
#define MIN_SPINAL_AREA_RATIO 3.0f
#define MIN_SPINAL_ASPECT .5
#define MAX_SPINAL_ASPECT 1.5

/* bounds on spatial extent */
#define MIN_SPINAL_X 60
#define MAX_SPINAL_X 180
#define MIN_SPINAL_Y 50
#define MAX_SPINAL_Y 170

/* bounds on centroid position */
#define MAX_SPINAL_WIDTH 80
#define MIN_SPINAL_CX (128 - MAX_SPINAL_WIDTH / 2)
#define MAX_SPINAL_CX (128 + MAX_SPINAL_WIDTH / 2)

static int find_spinal_fusion(MRI *mri_thresh, float *px, float *py, float *pz)
{
  MRI *mri_label, *mri_slice;
  float areas[MAX_LABELS], spinal_area, aspect, corner_dist, thick, x0, y0, z0;
  int y, found, nlabels, l, spinal_l = 0;
  MRI_REGION bboxes[MAX_LABELS], *bbox;

  thick = mri_thresh->thick;

  /*
    search from the bottom of the MRI upwards looking for the point at
    which the two branches of the spinal chord fuse together. This will
    be detectable using a bunch of hacky model-based parameters - where
    it is, when it appears (detected by increase in area), no other big stuff
    near it, in approximately the right place, with reasonable aspect
    ratio. Note that all these parameters are very loose.
  */
  for (found = 0, y = mri_thresh->height - 1; y >= 0 && !found; y--) {
    if (y < 160) DiagBreak();
    mri_slice = MRIextractPlane(mri_thresh, NULL, MRI_HORIZONTAL, y);
    mri_label = MRIlabel(mri_slice, NULL, &nlabels);

    /* first erase labels that are not in possible regions */
    MRIlabelBoundingBoxes(mri_label, bboxes, nlabels);
    for (l = 1; l < nlabels; l++) {
      if ((bboxes[l].x < MIN_SPINAL_X) || (bboxes[l].x + bboxes[l].dx - 1 > MAX_SPINAL_X) ||
          (bboxes[l].y < MIN_SPINAL_Y) || (bboxes[l].y + bboxes[l].dy - 1 > MAX_SPINAL_Y))
        MRIeraseLabel(mri_label, mri_label, l);
      else {
        MRIlabelCentroid(mri_label, l, &x0, &y0, &z0);
        if (x0 < MIN_SPINAL_CX || x0 > MAX_SPINAL_CX) MRIeraseLabel(mri_label, mri_label, l);
      }
    }
    MRIlabelAreas(mri_label, areas, nlabels); /* compute remaining areas */
    for (spinal_area = 0, l = 1; !found && l < nlabels; l++) {
      aspect = (float)bboxes[l].dy / bboxes[l].dx;
      if ((aspect < MIN_SPINAL_ASPECT) || (aspect > MAX_SPINAL_ASPECT)) continue;
      if (areas[l] >= spinal_area) /* find remaining label with biggest area */
      {
        spinal_l = l;
        spinal_area = areas[l];
      }
    }

    /* make sure it is much bigger than anything else close to it */
    if (spinal_area > MIN_SPINAL_AREA) {
      for (found = l = 1; found && l < nlabels; l++) {
        if ((l == spinal_l) || FZERO(areas[l])) continue;
        corner_dist = REGIONminCornerDistance(&bboxes[l], &bboxes[spinal_l]);
        if (corner_dist < 25 && (spinal_area / areas[l]) < MIN_SPINAL_AREA_RATIO) found = 0;
      }
    }
    if (found && Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "spinal_label%d.mgh", y);
      MRIwrite(mri_label, fname);
    }
    if (found) DiagBreak();
    MRIfree(&mri_label);
    MRIfree(&mri_slice);
    if (found) /* don't decrement y again */
      break;
  }

  /* convert bboxes to voxel coordinates */
  bbox = &bboxes[spinal_l];
  bbox->x /= thick;
  bbox->dx /= thick;
  bbox->y /= thick;
  bbox->dy /= thick;

  /* pick midpoint of label at appropriate slice */
  *px = (float)(bbox->x + bbox->dx / 2);
  *py = (float)y;
  *pz = (float)(bbox->y + bbox->dy / 2);
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "spinal fusion found at slice %d (%2.0f,%2.0f,%2.0f)\n", y, *px, *py, *pz);
  return (y);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIinitTranslation(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  MATRIX *m_translation;
  float in_means[4], ref_means[4];
  double dx, dy, dz;

  m_translation = MatrixIdentity(4, NULL);
  MRIfindCenterOfBrain(mri_in, in_means, in_means + 1, in_means + 2);
  MRIfindCenterOfBrain(mri_ref, ref_means, ref_means + 1, ref_means + 2);
  dx = (double)(ref_means[0] - in_means[0]) * mri_in->thick;
  dy = (double)(ref_means[1] - in_means[1]) * mri_in->thick;
  dz = (double)(ref_means[2] - in_means[2]) * mri_in->thick;
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout,
            "centering template around (%d,%d,%d) and input around"
            " (%d,%d,%d)\n",
            (int)ref_means[0],
            (int)ref_means[1],
            (int)ref_means[2],
            (int)in_means[0],
            (int)in_means[1],
            (int)in_means[2]);
    fprintf(stdout, "initial translation = (%2.0f, %2.0f, %2.0f).\n", dx, dy, dz);
  }
  *MATRIX_RELT(m_translation, 1, 4) = dx;
  *MATRIX_RELT(m_translation, 2, 4) = dy;
  *MATRIX_RELT(m_translation, 3, 4) = dz;
  MatrixMultiply(m_translation, m_L, m_L);
  MatrixFree(&m_translation);
  fflush(stdout);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_DX 1.1
#define MAX_DY 1.1
#define MAX_DZ 1.1
#define MIN_DX (1.0 / MAX_DX)
#define MIN_DY (1.0 / MAX_DY)
#define MIN_DZ (1.0 / MAX_DZ)
#define MAX_RATIO 1.1

int MRIinitScaling(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  MATRIX *m_scaling;
  float sx, sy, sz, dx, dy, dz;
  MRI_REGION in_bbox, ref_bbox;

  m_scaling = MatrixIdentity(4, NULL);

  MRIboundingBox(mri_in, 120, &in_bbox);
  MRIboundingBox(mri_ref, 120, &ref_bbox);
  sx = (float)ref_bbox.dx / (float)in_bbox.dx;
  sy = (float)ref_bbox.dy / (float)in_bbox.dy;
  sz = (float)ref_bbox.dz / (float)in_bbox.dz;
  dx = (ref_bbox.x + ref_bbox.dx - 1) / 2 - (in_bbox.x + in_bbox.dx - 1) / 2;
  dy = (ref_bbox.y + ref_bbox.dy - 1) / 2 - (in_bbox.y + in_bbox.dy - 1) / 2;
  dz = (ref_bbox.z + ref_bbox.dz - 1) / 2 - (in_bbox.z + in_bbox.dz - 1) / 2;

  if (sx > MAX_DX) sx = MAX_DX;
  if (sx < MIN_DX) sx = MIN_DX;
  if (sy > MAX_DY) sy = MAX_DY;
  if (sy < MIN_DY) sy = MIN_DY;
  if (sz > MAX_DZ) sz = MAX_DZ;
  if (sz < MIN_DZ) sz = MIN_DZ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,
            "initial scaling: (%2.2f, %2.2f, %2.2f) <-- "
            "(%d/%d,%d/%d,%d/%d)\n",
            sx,
            sy,
            sz,
            ref_bbox.dx,
            in_bbox.dx,
            ref_bbox.dy,
            in_bbox.dy,
            ref_bbox.dz,
            in_bbox.dz);
  *MATRIX_RELT(m_scaling, 1, 1) = sx;
  *MATRIX_RELT(m_scaling, 2, 2) = sy;
  *MATRIX_RELT(m_scaling, 3, 3) = sz;

  MatrixMultiply(m_scaling, m_L, m_L);
  fflush(stdout);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MIN_PYR_WIDTH 16
#define USE_PYRAMID 1
int MRIlinearAlign(MRI *mri_in, MRI *mri_ref, MP *parms)
{
#if USE_PYRAMID
  int nlevels, max_levels, i;
  MRI *mri_in_pyramid[MAX_LEVELS], *mri_ref_pyramid[MAX_LEVELS];
#endif
  double rms;
  char base_name[STRLEN];
  double tmul;

  parms->rigid = 0;
  mriNormalizeStds(mri_ref);
  if (DZERO(parms->dt)) parms->dt = 1e-6;

  if (parms->disable_neck) {
    if (Gdiag & DIAG_SHOW) fprintf(stdout, "disabling neck discounting...\n");
    parms->ref_np.neck_x0 = parms->ref_np.neck_y0 = parms->ref_np.neck_z0 = 1000;
    parms->ref_np.neck_dx = parms->ref_np.neck_dy = parms->ref_np.neck_dz = 1;
    parms->in_np.neck_x0 = parms->in_np.neck_y0 = parms->in_np.neck_z0 = 1000;
    parms->in_np.neck_dx = parms->in_np.neck_dy = parms->in_np.neck_dz = 1;
  }

  if (mri_ref->nframes > 1)
    mri_ref->mean = MRImeanFrame(mri_ref, 1);
  else
    mri_ref->mean = 1.0;
  if (FZERO(mri_ref->mean)) mri_ref->mean = 1.0f;
  if (parms->lta->num_xforms == 1) /* one (global) transformation */
  {

    if (Gdiag & DIAG_SHOW /*&& DIAG_VERBOSE_ON*/) {
      rms = mriIntensityRMS(mri_in, mri_ref, parms->lta, 1.0f, &parms->in_np);
      fprintf(stdout, "after initial alignment  rms = %2.3f\n", rms);
    }
  }

  strcpy(base_name, parms->base_name);
  openLogFile(parms);
#if !USE_PYRAMID
  {
    double dt, sigma;
    MRI *mri_in_blur = NULL, *mri_ref_blur = NULL, *mri_kernel;

    dt = parms->dt;
    for (sigma = 4.0f; sigma >= 1.0f; sigma /= 2.0f) {
      parms->dt = sigma > 1 ? dt * sigma : dt;
      if (sigma > 0.0f) {
        mri_kernel = MRIgaussian1d(sigma, 17);
        fprintf(stdout, "blurring volumes with sigma = %2.1f...", sigma);
        fflush(stdout);
        mri_in_blur = MRIconvolveGaussian(mri_in, mri_in_blur, mri_kernel);
        mri_ref_blur = MRIconvolveGaussianMeanAndStdByte(mri_ref, mri_ref_blur, mri_kernel);
        fprintf(stdout, "done.\n");
        fflush(stdout);
        MRIfree(&mri_kernel);
        mriLinearAlignPyramidLevel(mri_in_blur, mri_ref_blur, parms);
      }
      else
        mriLinearAlignPyramidLevel(mri_in, mri_ref, parms);
    }
    parms->dt = dt;
    MRIfree(&mri_in_blur);
    MRIfree(&mri_ref_blur);
  }
#else
  if (parms->lta->num_xforms == 1)
    max_levels = MAX_LEVELS;
  else
    max_levels = 2;

  /* build Gaussian pyramid */
  mri_in_pyramid[0] = mri_in;
  mri_ref_pyramid[0] = mri_ref;
  for (nlevels = 1; nlevels < max_levels; nlevels++) {
    if (mri_in_pyramid[nlevels - 1]->width <= MIN_PYR_WIDTH) break;
    mri_in_pyramid[nlevels] = MRIreduce(mri_in_pyramid[nlevels - 1], NULL);
    mri_ref_pyramid[nlevels] = MRIreduceMeanAndStd(mri_ref_pyramid[nlevels - 1], NULL);
  }

  for (i = nlevels - 1; i >= 1; i--) {
    for (tmul = 1; tmul >= .5; tmul /= 10) {
      fprintf(stdout, "aligning pyramid level %d.\n", i);
      fflush(stdout);
      if ((Gdiag & DIAG_WRITE) && parms->log_fp) fprintf(parms->log_fp, "aligning pyramid level %d.\n", i);
      parms->trans_mul = tmul;
      mriLinearAlignPyramidLevel(mri_in_pyramid[i], mri_ref_pyramid[i], parms);
    }
  }

  /* free Gaussian pyramid */
  for (i = 1; i < nlevels; i++) {
    MRIfree(&mri_in_pyramid[i]);
    MRIfree(&mri_ref_pyramid[i]);
  }
#endif
  strcpy(parms->base_name, base_name);
  if (parms->log_fp) {
    fclose(parms->log_fp);
    parms->log_fp = NULL;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRIfindMeans(MRI *mri, float *means)
{
  int width, height, depth, x, y, z, val, npoints, mx, my, mz;
  BUFTYPE *psrc;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  mx = my = mz = npoints = 0;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri, 0, y, z);
      for (x = 0; x < width; x++) {
        val = *psrc++;
        if (val > 30) {
          npoints++;
          mx += x;
          my += y;
          mz += z;
        }
      }
    }
  }

  if (npoints) {
    means[0] = (float)mx / npoints;
    means[1] = (float)my / npoints;
    means[2] = (float)mz / npoints;
  }
  else
    means[0] = means[1] = means[2] = 0.0f;

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRIfindCenterOfBrain(MRI *mri, float *px0, float *py0, float *pz0)
{
  int x, y, z, width, height, depth;
  float x0, y0, z0, total, val;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  x0 = y0 = z0 = 0.0f;
  total = 0.0f;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        val = (float)MRIvox(mri, x, y, z);
        if (val > 20) /* don't let differences in bg bias result */
        {
          x0 += (float)x * val;
          y0 += (float)y * val;
          z0 += (float)z * val;
          total += val;
        }
      }
    }
  }
  if (total > 0) {
    *px0 = x0 / (float)total;
    *py0 = y0 / (float)total;
    *pz0 = z0 / (float)total;
  }
  else {
    *px0 = width / 2;
    *py0 = height / 2;
    *pz0 = depth / 2;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/

static int mriLinearAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  double rms, old_rms, dt;
  int k, n, nsmall = 0;
  MATRIX *m_L;
  static int ncalls = 0;
  char fname[STRLEN], base_name[STRLEN];

  if (parms->mri_in == NULL) parms->mri_in = mri_in;
  if (parms->mri_ref == NULL) parms->mri_ref = mri_ref;

  /*  MRIeraseBorders(mri_in, 1) ; MRIeraseBorders(mri_ref, 1) ;*/
  if (mri_ref->nframes > 1)
    mri_ref->mean = MRImeanFrame(mri_ref, 1);
  else
    mri_ref->mean = 1.0;
  fprintf(stdout, "mean std = %2.2f\n", mri_ref->mean);
  if (DZERO(mri_ref->mean)) mri_ref->mean = 1;

  /* clear momentum from prior level */
  for (k = 0; k < parms->lta->num_xforms; k++) MatrixClear(parms->lta->xforms[k].m_last_dL);

  strcpy(base_name, parms->base_name);
  int req = snprintf(parms->base_name, 100, "level%d_%s", ncalls, base_name);
  if( req >= 100 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  ncalls++; /* for diagnostics */
  dt = parms->dt;

  m_L = parms->lta->xforms[0].m_L;
  old_rms = rms = mriIntensityRMS(mri_in, mri_ref, parms->lta, 1.0f, &parms->in_np);
  fprintf(stdout,
          "000: rms = %2.3f, t = [%2.2f, %2.2f, %2.2f], dt = %2.2e, "
          "thick=%2.0f, mul=%2.2f\n",
          rms,
          m_L->rptr[1][4],
          m_L->rptr[2][4],
          m_L->rptr[3][4],
          dt,
          mri_in->thick,
          parms->trans_mul);
  fflush(stdout);

  if ((Gdiag & DIAG_WRITE) && parms->log_fp) {
    fprintf(parms->log_fp,
            "\n000: rms = %2.3f, t = [%2.2f, %2.2f, %2.2f], dt = %2.2e\n",
            rms,
            m_L->rptr[1][4],
            m_L->rptr[2][4],
            m_L->rptr[3][4],
            dt);
  }
  if ((Gdiag & DIAG_WRITE) && (parms->write_iterations > 0) && !parms->start_t) {
    writeSnapshot(parms->mri_in, parms, 0);

    if (ncalls == 1) {
#if USE_INVERSE == 0
      int req = snprintf(fname, STRLEN, "%sref", base_name);  
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "writing reference views to %s...\n", fname);
      MRIwriteImageViews(parms->mri_ref, fname, IMAGE_SIZE);
#endif
    }
  }
  for (n = 0; n < parms->niterations; n++) {
    dt = ltaGradientStep(mri_in, mri_ref, parms->lta, parms->dt, parms->momentum, &parms->in_np, parms->trans_mul);
    if (parms->rigid) mriOrthonormalizeTransform(parms->lta->xforms[0].m_L);
    rms = mriIntensityRMS(mri_in, mri_ref, parms->lta, 1.0f, &parms->in_np);
    fprintf(stdout,
            "%3.3d: rms = %2.3f, t = [%2.2f, %2.2f, %2.2f], "
            "dt=%2.2e\n",
            n + 1,
            rms,
            m_L->rptr[1][4],
            m_L->rptr[2][4],
            m_L->rptr[3][4],
            dt);
    fflush(stdout);

    if ((Gdiag & DIAG_WRITE) && (parms->write_iterations > 0) && (((n + 1) % parms->write_iterations) == 0))
      writeSnapshot(parms->mri_in, parms, n + 1);
    logIntegration(parms, n, rms);
    if (FZERO(rms) || (((old_rms - rms) / rms) < parms->tol))
      nsmall++;
    else if (nsmall > 0)
      nsmall = 0;

    /* more than 2 small steps or error increased - integration asymptoted */
    if ((nsmall > 2 || old_rms < rms) && !parms->rigid) break;
    old_rms = rms;
  }
  fprintf(stdout, "\n");

  if ((Gdiag & DIAG_WRITE) && parms->log_fp) fprintf(parms->log_fp, "\n");
  strcpy(parms->base_name, base_name);
  fflush(stdout);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the rms intensity error between a reference image and
          an input image after a linear transformation of the input
          coordinate system.
------------------------------------------------------*/
static double mriIntensityRMS(MRI *mri_in, MRI *mri_ref, LTA *lta, double l_intensity, NECK_PARMS *np)
{
  int x1, x2, x3, width, height, depth;
  VECTOR *v_X, *v_Y; /* original and transformed coordinate systems */
  double val, y1, y2, y3, thick;
  double sse = 0.0f, delta = 0.0, std, mean_std, dot, nx, ny, nz, ndx, ndy, ndz;
  float in_val;

  mean_std = mri_ref->mean;
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;

  v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */

  nx = np->neck_x0;
  ny = np->neck_y0;
  nz = np->neck_z0;
  ndx = np->neck_dx;
  ndy = np->neck_dy;
  ndz = np->neck_dz;
  thick = mri_in->thick;
  for (x3 = 0; x3 < depth; x3++) {
    v_X->rptr[4][1] = 1.0f / thick;
    v_X->rptr[3][1] = x3;
    for (x2 = 0; x2 < height; x2++) {
      v_X->rptr[2][1] = x2;
      for (x1 = 0; x1 < width; x1++) {
        if (x1 == 3 * width / 4 && x2 == 3 * height / 4 && x3 == 3 * depth / 4) DiagBreak();

        v_X->rptr[1][1] = x1;
        LTAtransformPoint(lta, v_X, v_Y);

        dot = (x1 * thick - nx) * ndx + (x2 * thick - ny) * ndy + (x3 * thick - nz) * ndz;
        if (dot >= 0) DiagBreak();
        dot = 1 / (1 + exp(dot)); /* 'don't care' weighting */
        y1 = (double)v_Y->rptr[1][1];
        y2 = (double)v_Y->rptr[2][1];
        y3 = (double)v_Y->rptr[3][1];

        in_val = MRIgetVoxVal(mri_in, x1, x2, x3, 0);
        if (y1 > -1 && y1 < width && y2 > -1 && y2 < height && y3 > -1 && y3 < depth) {
          MRIsampleVolume(mri_ref, y1, y2, y3, &val);
          if (mri_ref->nframes > 1)
            MRIsampleVolumeFrame(mri_ref, y1, y2, y3, 1, &std);
          else
            std = mean_std;
          std /= mean_std;
          if (DZERO(std)) std = 1.0;
          delta = (val - (double)in_val) / std;
        }
        else /* out of bounds, assume in val is 0 */
          delta = (0.0 - in_val) / mean_std;
        if (fabs(delta) > 20) DiagBreak();
        delta *= dot;
        sse += delta * delta * mri_in->thick;
        if (sse > 2000) DiagBreak();
        if (!std::isfinite(sse)) DiagBreak();
      }
    }
  }

  MatrixFree(&v_X);
  MatrixFree(&v_Y);

  sse = sqrt(sse / ((double)(width * height * depth)));
  return (sse); /* rms */
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int writeSnapshot(MRI *mri, MORPH_PARMS *parms, int n)
{
  MRI *mri_tmp;
  char fname[STRLEN];
  TRANSFORM transform;
  GCA *gca;

  if (!(Gdiag & DIAG_WRITE)) return (NO_ERROR);

#define USE_INVERSE 0
#if USE_INVERSE
  if (!n) {
    sprintf(fname, "%s_input", parms->base_name);
    MRIwriteImageViews(parms->mri_in, fname, IMAGE_SIZE);
  }
  mri_tmp = MRIinverseLinearTransform(parms->mri_ref, NULL, parms->lta->xforms[0].m_L);
#else
  mri_tmp = LTAtransform(mri, NULL, parms->lta);
#endif
  gca = parms->gca_red;
  if (mri->xsize < gca->xsize || mri->ysize < gca->ysize || mri->zsize < gca->zsize) {
    MRI *mri_extracted;
    mri_extracted = MRIextract(mri_tmp, NULL, 0, 0, 0, gca->width, gca->height, gca->depth);
    GCAcopyDCToMRI(gca, mri_extracted);
    MRIfree(&mri_tmp);
    mri_tmp = mri_extracted;
    mri_tmp->xsize = gca->xsize;
    mri_tmp->ysize = gca->ysize;
    mri_tmp->zsize = gca->zsize;
  }

  if ((n % 50 * parms->write_iterations) == 0 && DIAG_VERBOSE_ON) {
    sprintf(fname, "%s%3.3d.mgh", parms->base_name, n);
    if (Gdiag & DIAG_SHOW) fprintf(stdout, "writing snapshot to %s...\n", fname);
    MRIwrite(mri_tmp, fname);
  }
  sprintf(fname, "%s%3.3d", parms->base_name, n);
  MRIwriteImageViews(mri_tmp, fname, IMAGE_SIZE);
  sprintf(fname, "%s%3.3d.mgz", parms->base_name, n);
  MRIwrite(mri_tmp, fname);
  MRIfree(&mri_tmp);

  sprintf(fname, "%s%3.3d_fsamples.mgz", parms->base_name, n);
  printf("writing transformed samples to %s....\n", fname);
  transform.xform = (void *)parms->lta;
  transform.type = LINEAR_VOX_TO_VOX;
  GCAtransformAndWriteSamples((GCA *)parms->vgca, mri, parms->gcas, parms->nsamples, fname, parms->transform);
  fflush(stdout);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double ltaGradientStep(
    MRI *mri_in, MRI *mri_ref, LTA *lta, double dt, double momentum, NECK_PARMS *np, double tmul)
{
  int width, height, depth;
  VECTOR *v_X, *v_Y, *v_Yk; /* original and transformed coordinate systems */
  MATRIX *m_tmp, *m_L, *m_dT_X_T;
  VECTOR *v_dT, *v_X_T; /* gradient of mri_ref */
  double val, x1, x2, x3, y1, y2, y3, thick /*, len*/, in_val;
  double delta = 0.0, n, std, sigma, dx, dy, dz, wtotal, dsq, w_k_p, mean_std, dot, nx, ny, nz, ndx, ndy, ndz;
  int only_translation, k;
  LT *lt;
  MRI_REGION bbox;

  mean_std = mri_ref->mean;
  MRIboundingBox(mri_in, 5, &bbox);

  only_translation = (getenv("ONLY_TRANSLATION") != NULL);

  m_L = lta->xforms[0].m_L;
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  m_tmp = MatrixAlloc(4, 4, MATRIX_REAL);
  m_dT_X_T = MatrixAlloc(4, 4, MATRIX_REAL);
  v_X = VectorAlloc(4, MATRIX_REAL);    /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL);    /* transformed (dst) coordinates */
  v_Yk = VectorAlloc(4, MATRIX_REAL);   /* transformed (dst) coordinates */
  v_X_T = RVectorAlloc(4, MATRIX_REAL); /* transpose of input coords */
  v_dT = VectorAlloc(4, MATRIX_REAL);   /* gradient of target image */

  v_X->rptr[4][1] = 1.0f / mri_in->thick;
  for (k = 0; k < lta->num_xforms; k++) MatrixClear(lta->xforms[k].m_dL);

  nx = np->neck_x0;
  ny = np->neck_y0;
  nz = np->neck_z0;
  ndx = np->neck_dx;
  ndy = np->neck_dy;
  ndz = np->neck_dz;
  thick = mri_in->thick;
  for (n = 0.0, x3 = bbox.z; x3 < bbox.z + bbox.dz; x3++) {
    V3_Z(v_X) = x3;
    for (x2 = bbox.y; x2 < bbox.y + bbox.dy; x2++) {
      V3_Y(v_X) = x2;
      for (x1 = bbox.x; x1 < bbox.x + bbox.dx; x1++) {
        dot = (x1 * thick - nx) * ndx + (x2 * thick - ny) * ndy + (x3 * thick - nz) * ndz;
        if (dot >= 0) DiagBreak();
        dot = 1 / (1 + exp(dot)); /* 'don't care' weighting */
        V3_X(v_X) = x1;
        MatrixTranspose(v_X, v_X_T);

        /* first compute weight normalization factor */

        wtotal = LTAtransformPointAndGetWtotal(lta, v_X, v_Y);

        y1 = V3_X(v_Y);
        y2 = V3_Y(v_Y);
        y3 = V3_Z(v_Y);

        if (y1 > -1 && y1 < width && y2 > -1 && y2 < height && y3 > -1 && y3 < depth) {
          MRIsampleVolume(mri_ref, y1, y2, y3, &val);
          if (mri_ref->nframes > 1)
            MRIsampleVolumeFrame(mri_ref, y1, y2, y3, 1, &std);
          else
            std = mean_std;
          MRIsampleVolumeGradient(mri_ref, y1, y2, y3, &dx, &dy, &dz);
        }
        else {
          std = mean_std;
          val = 0.0f;
          dx = dy = dz = 0.0f;
        }
        std /= mean_std;
        if (DZERO(std)) std = 1.0;
        in_val = MRIgetVoxVal(mri_in, nint(x1), nint(x2), nint(x3), 0);
        delta = (val - (double)in_val) / std;
        delta *= dot;

        if (only_translation) {
          RV3_X(v_X_T) = RV3_Y(v_X_T) = RV3_Z(v_X_T) = 0;
        }
        V3_X(v_dT) = dx;
        V3_Y(v_dT) = dy;
        V3_Z(v_dT) = dz;
        MatrixMultiply(v_dT, v_X_T, m_dT_X_T);

        for (k = 0; k < lta->num_xforms; k++) {
          lt = &lta->xforms[k];
          sigma = lt->sigma;
          dx = lt->x0 - x1 * mri_in->thick;
          dy = lt->y0 - x2 * mri_in->thick;
          dz = lt->z0 - x3 * mri_in->thick;
          dsq = dx * dx + dy * dy + dz * dz;
          w_k_p = exp(-dsq / (2 * sigma * sigma));
          w_k_p /= wtotal;
          MatrixScalarMul(m_dT_X_T, w_k_p * mri_in->thick * delta, m_tmp);
          m_tmp->rptr[1][4] *= tmul;
          m_tmp->rptr[2][4] *= tmul;
          m_tmp->rptr[3][4] *= tmul;
          MatrixCheck(m_tmp);
          MatrixAdd(m_tmp, lt->m_dL, lt->m_dL);
          MatrixCheck(lt->m_dL);
          n++;
        }
      }
    }
  }

  if (FZERO(n))
    DiagBreak();
  else {
    int r, c;
    double max_val = 0;

    dt *= (lta->num_xforms / n);
    for (k = 0; k < lta->num_xforms; k++) {
      /* take step in negative gradient direction */
      lt = &lta->xforms[k];

#define MAX_LINEAR_DEL (0.001 * mri_in->thick)
#define MAX_TRANSLATION_DEL 0.25
      for (r = 1; r <= 3; r++) {
        for (c = 1; c <= 3; c++) {
          if (fabs(*MATRIX_RELT(lt->m_dL, r, c)) > max_val) max_val = fabs(*MATRIX_RELT(lt->m_dL, r, c));
        }
      }

      for (max_val = 0.0, r = 1; r <= lt->m_dL->rows; r++) {
        if (fabs(*MATRIX_RELT(lt->m_dL, r, 4)) > max_val) max_val = fabs(*MATRIX_RELT(lt->m_dL, r, 4));
      }

      MatrixScalarMul(lt->m_dL, -dt, lt->m_dL);
      MatrixScalarMul(lt->m_last_dL, momentum, lt->m_last_dL);
      MatrixAdd(lt->m_last_dL, lt->m_dL, lt->m_dL);
      MatrixAdd(lt->m_dL, lt->m_L, lt->m_L);
      MatrixCopy(lt->m_dL, lt->m_last_dL);
      MatrixCheck(lt->m_L);
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(stdout, "dL:\n");
        MatrixPrintHires(stdout, lt->m_dL);
        fprintf(stdout, "L:\n");
        MatrixPrintHires(stdout, lt->m_L);
      }
    }
  }
  MatrixFree(&m_dT_X_T);
  MatrixFree(&v_X);
  MatrixFree(&v_Y);
  MatrixFree(&v_X_T);
  MatrixFree(&v_dT);
  MatrixFree(&m_tmp);
  fflush(stdout);
  return (dt);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int openLogFile(MORPH_PARMS *parms)
{
  char fname[STRLEN];

  if ((Gdiag & DIAG_WRITE) && (parms->log_fp == NULL)) {
    sprintf(fname, "%s.log", parms->base_name);
    if (!parms->start_t)
      parms->log_fp = fopen(fname, "w");
    else
      parms->log_fp = fopen(fname, "a");
    if (!parms->log_fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);

    fprintf(parms->log_fp,
            "dt=%2.2e, momentum=%2.2f, tol=%2.2e, niteratons=%d\n",
            parms->dt,
            parms->momentum,
            parms->tol,
            parms->niterations);
    if (!FZERO(parms->l_intensity)) fprintf(parms->log_fp, "l_intensity = %2.4f\n", parms->l_intensity);
    if (!FZERO(parms->l_priors)) fprintf(parms->log_fp, "l_priors = %2.4f\n", parms->l_priors);
    if (!FZERO(parms->l_dist)) fprintf(parms->log_fp, "l_dist = %2.4f\n", parms->l_dist);
    if (!FZERO(parms->l_compression)) fprintf(parms->log_fp, "l_compression = %2.4f\n", parms->l_compression);
    if (!FZERO(parms->l_area)) fprintf(parms->log_fp, "l_area = %2.4f\n", parms->l_area);
    if (!FZERO(parms->l_nlarea)) fprintf(parms->log_fp, "l_nlarea = %2.4f\n", parms->l_nlarea);
    if (!FZERO(parms->l_nldist)) fprintf(parms->log_fp, "l_nldist = %2.4f\n", parms->l_nldist);
    if (!FZERO(parms->sigma)) fprintf(parms->log_fp, "sigma   = %2.4f\n", parms->sigma);
    if (!FZERO(parms->exp_k)) fprintf(parms->log_fp, "exp k   = %2.1f\n", parms->exp_k);
    fflush(parms->log_fp);
  }
  fprintf(stdout, "dt = %2.2e, momentum=%2.2f, tol=%2.2e\n", parms->dt, parms->momentum, parms->tol);
  if (!FZERO(parms->l_intensity)) fprintf(stdout, "l_intensity = %2.4f\n", parms->l_intensity);
  if (!FZERO(parms->l_dist)) fprintf(stdout, "l_compression = %2.4f\n", parms->l_compression);
  if (!FZERO(parms->l_dist)) fprintf(stdout, "l_dist = %2.4f\n", parms->l_dist);
  if (!FZERO(parms->l_area)) fprintf(stdout, "l_area = %2.4f\n", parms->l_area);
  if (!FZERO(parms->l_nlarea)) fprintf(stdout, "l_nlarea = %2.4f\n", parms->l_nlarea);
  if (!FZERO(parms->l_nldist)) fprintf(stdout, "l_nldist = %2.4f\n", parms->l_nldist);
  if (!FZERO(parms->sigma)) fprintf(stdout, "sigma   = %2.4f\n", parms->sigma);
  fflush(stdout);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int logIntegration(MORPH_PARMS *parms, int n, double rms)
{
  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->log_fp, "%3.3d: rms = %2.3f\n", n + 1, rms);
    fflush(parms->log_fp);
  }
  return (NO_ERROR);
}
/*--------------------------------------------------------------
                                3D Morph
----------------------------------------------------------------*/

#if USE_ITERATIVE_AVERAGING
static int m3daverageGradient(MORPH_3D *m3d, int niter);
#else
static int m3dblurGradient(MORPH_3D *m3d, float sigma);
#endif
static int m3dCopyTempToGradient(MORPH_3D *m3d);
static double m3dSSE(MRI *mri_in,
                     MRI *mri_ref,
                     MORPH_PARMS *parms,
                     MORPH_3D *m3d,
                     double *pintensity_rms,
                     double *pdistance_rms,
                     double *parea_rms);
static double m3dCorrelationSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d, NECK_PARMS *np);
static double m3dDistanceSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d);
static int m3dCompressionTerm(
    MORPH_3D *m3d, double l_compression, int i, int j, int k, double *pdx, double *pdy, double *pdz);
static int m3dNonlinearDistanceTerm(
    MORPH_3D *m3d, double l_nldist, int i, int j, int k, double *pdx, double *pdy, double *pdz);
static double m3dCompressionSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d);
static double m3dAreaSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d);
static double m3dNonlinearAreaSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d, MORPH_PARMS *parms);
static int m3dCorrelationTerm(MRI *mri_in,
                              MRI *mri_ref,
                              MRI *mri_ref_blur,
                              double l_intensity,
                              int xin,
                              int yin,
                              int zin,
                              double xref,
                              double yref,
                              double zref,
                              double *pdx,
                              double *pdy,
                              double *pdz);
static int m3dDistanceTerm(
    MORPH_3D *m3d, double l_distance, int i, int j, int k, double *pdx, double *pdy, double *pdz);
static int m3dNonlinearAreaTerm(
    MORPH_3D *m3d, double l_nlarea, int i, int j, int k, double *pdx, double *pdy, double *pdz, MORPH_PARMS *parms);
static int m3dAreaTerm(MORPH_3D *m3d, double l_area, int i, int j, int k, double *pdx, double *pdy, double *pdz);
static MORPH_3D *m3dInit(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, float scale);
static MORPH_3D *m3dalloc(int width, int height, int depth, float spacing);
static int m3dAllocProperties(MORPH_3D *m3d);
static int m3dInitDistances(MORPH_3D *m3d);
static int m3dInitAreas(MORPH_3D *m3d);
static double m3dIntegrationStep(
    MRI *mri_in, MRI *mri_ref, MRI *mri_ref_blur, MORPH_PARMS *parms, MORPH_3D *m3d, double dt);
static int m3dStoreMetricProperties(MORPH_3D *m3d);
static int m3dComputeMetricProperties(MORPH_3D *m3d);
static int m3dclearGradient(MORPH_3D *m3d);
static int m3dapplyGradient(MORPH_3D *m3d, double dt);
static double m3dscaleDeltaT(MORPH_3D *m3d, double dt, double max_grad, MORPH_PARMS *parms);
static int MRIsample3Dmorph(MORPH_3D *m3d, float x, float y, float z, float *pxd, float *pyd, float *pzd);
static MORPH_3D *m3dscaleUp2(MORPH_3D *m3d_in, MORPH_3D *m3d_out, MORPH_PARMS *parms);

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute SSE for entire energy functional, and
          return rms error for the individual terms as well.
------------------------------------------------------*/
static double m3dSSE(MRI *mri_in,
                     MRI *mri_ref,
                     MORPH_PARMS *parms,
                     MORPH_3D *m3d,
                     double *pintensity_rms,
                     double *pdistance_rms,
                     double *parea_rms)
{
  double sse, intensity_sse, area_sse, distance_sse, nvox, nlarea_sse, compression_sse;

  intensity_sse = m3dCorrelationSSE(mri_in, mri_ref, m3d, &parms->in_np);
  area_sse = m3dAreaSSE(mri_in, mri_ref, m3d);
  nlarea_sse = m3dNonlinearAreaSSE(mri_in, mri_ref, m3d, parms);
  distance_sse = m3dDistanceSSE(mri_in, mri_ref, m3d);
  compression_sse = m3dCompressionSSE(mri_in, mri_ref, m3d);

  nvox = m3d->width * m3d->height * m3d->depth;
  *pintensity_rms = sqrt(intensity_sse / nvox);
  *pdistance_rms = sqrt(distance_sse / nvox);
  *parea_rms = sqrt(area_sse / nvox);
  sse = parms->l_intensity * intensity_sse + parms->l_area * area_sse + parms->l_nlarea * nlarea_sse +
        parms->l_compression * compression_sse + parms->l_dist * distance_sse;
  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute correlation SSE
------------------------------------------------------*/
static double m3dCorrelationSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d, NECK_PARMS *np)
{
  double sse;
  int width, height, depth, x, y, z, xsi, ysi, zsi, *pxi, *pyi, *pzi;
  MORPH_NODE *mn;
  float x1, x2, x3, scale, dot, nx, ny, nz, ndx, ndy, ndz, thick, pdx, pdy, pdz;
  double ref_val, in_val, delta, xd, yd, zd, ref_std, mean_std;

  mean_std = mri_ref->mean;
  thick = mri_in->thick;
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  scale = m3d->node_spacing / thick;
  pxi = mri_in->xi;
  pyi = mri_in->yi;
  pzi = mri_in->zi;
  nx = np->neck_x0;
  ny = np->neck_y0;
  nz = np->neck_z0;
  ndx = np->neck_dx;
  ndy = np->neck_dy;
  ndz = np->neck_dz;
  for (sse = 0.0f, x = 0; x < m3d->width; x++) {
    x1 = x * scale; /* convert to voxel space of mri_in */
    xsi = pxi[nint(x1)];
    for (y = 0; y < m3d->height; y++) {
      x2 = y * scale;
      ysi = pyi[nint(x2)]; /* voxel coords of mri_in */
      for (z = 0; z < m3d->depth; z++) {
        x3 = z * scale;
        zsi = pzi[nint(x3)];       /* voxel coords of mri_in */
        mn = &m3d->nodes[x][y][z]; /* find out where this voxel went */
        xd = mn->x / mri_ref->thick;
        yd = mn->y / mri_ref->thick;
        zd = mn->z / mri_ref->thick;
        if (xd > -1 && yd > -1 && zd > 0 && xd < width && yd < height && zd < depth) {
          MRIsampleVolume(mri_ref, xd, yd, zd, &ref_val);
          MRIsampleVolumeFrame(mri_ref, xd, yd, zd, 1, &ref_std);
          ref_std /= mean_std;
        }
        else {
          ref_std = mean_std;
          ref_val = 0.0;
        }
        if (x == 3 * width / 4 && y == 3 * height / 4 && z == 3 * depth / 4) DiagBreak();
        if (DZERO(ref_std)) ref_std = 1.0; /* don't let morph be driven by background */
        in_val = (double)MRIgetVoxVal(mri_in, xsi, ysi, zsi, 0);
        pdx = (x * thick - nx);
        pdy = (y * thick - ny);
        pdz = (z * thick - nz);
        dot = pdx * ndx + pdy * ndy + pdz * ndz;
        if (dot >= 0) DiagBreak();
        dot = 1 / (1 + exp(.25 * (dot - 10))); /* 'don't care' weighting */
        delta = (in_val - ref_val) / ref_std;
        delta *= dot;
        sse += delta * delta * mri_in->thick; /* delta^2 dx */
        if (!finitep(sse)) DiagBreak();
      }
    }
  }

  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute SSE for distance term.
------------------------------------------------------*/
static double m3dDistanceSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d)
{
  int i, j, k, width, height, depth, n;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;
  double dist, xd, yd, zd, delta, sse, node_spacing;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  node_spacing = m3d->node_spacing;
  for (sse = 0.0, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        for (n = 0; n < NEIGHBORS; n++) {
          switch (n) {
            default:
            case 0: /* i-1 */
              if (i > 0)
                mn_nbr = &m3d->nodes[i - 1][j][k];
              else
                mn_nbr = NULL;
              break;
            case 1: /* i+1 */
              if (i < width - 1)
                mn_nbr = &m3d->nodes[i + 1][j][k];
              else
                mn_nbr = NULL;
              break;
            case 2: /* j-1 */
              if (j > 0)
                mn_nbr = &m3d->nodes[i][j - 1][k];
              else
                mn_nbr = NULL;
              break;
            case 3: /* j+1 */
              if (j < height - 1)
                mn_nbr = &m3d->nodes[i][j + 1][k];
              else
                mn_nbr = NULL;
              break;
            case 4: /* k-1 */
              if (k > 0)
                mn_nbr = &m3d->nodes[i][j][k - 1];
              else
                mn_nbr = NULL;
              break;
            case 5: /* k+1 */
              if (k < depth - 1)
                mn_nbr = &m3d->nodes[i][j][k + 1];
              else
                mn_nbr = NULL;
              break;
          }
          if (mn_nbr) {
            xd = mn->x - mn_nbr->x;
            yd = mn->y - mn_nbr->y;
            zd = mn->z - mn_nbr->z;
            dist = sqrt(xd * xd + yd * yd + zd * zd);
            delta = (dist - mnp->orig_dist[n]);
#if SCALE_INVARIANT
            delta /= node_spacing;
#endif
            sse += delta * delta * mri_in->thick;
          }
        }
      }
    }
  }

  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute SSE for area term.
------------------------------------------------------*/
static double m3dNonlinearAreaSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d, MORPH_PARMS *parms)
{
  double sse, delta, n3, ratio, exponent;
  int i, j, k, width, height, depth;
  MORPH_NODE *mn;
  MNP *mnp;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
/*  m3d->neg = 0 ;*/
#if SCALE_INVARIANT
  n3 = m3d->node_spacing * m3d->node_spacing * m3d->node_spacing;
#else
  n3 = 1.0f;
#endif
  for (sse = 0.0, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        /* scale up the area coefficient if the area of the current node is
           close to 0 or already negative */
        if (!FZERO(mnp->orig_area))
          ratio = mnp->area / mnp->orig_area;
        else {
          ratio = 1;
        }
        exponent = -parms->exp_k * ratio;
        if (exponent > MAX_EXP)
          delta = 0.0;
        else
          delta = log(1 + exp(exponent)) / parms->exp_k;

        sse += delta * mri_in->thick;
        if (!finitep(delta) || !finitep(sse)) DiagBreak();
      }
    }
  }

  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute SSE for area term.
------------------------------------------------------*/
static double m3dAreaSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d)
{
  double sse, delta, n3;
  int i, j, k, width, height, depth;
  MORPH_NODE *mn;
  MNP *mnp;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
/*  m3d->neg = 0 ;*/
#if SCALE_INVARIANT
  n3 = m3d->node_spacing * m3d->node_spacing * m3d->node_spacing;
#else
  n3 = 1.0f;
#endif
  for (sse = 0.0, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        delta = (mnp->area - mnp->orig_area) / n3;
        sse += delta * delta * mri_in->thick;
        if (!finitep(delta) || !finitep(sse)) DiagBreak();
      }
    }
  }

  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate and initialize 3D properties
------------------------------------------------------*/
static int m3dAllocProperties(MORPH_3D *m3d)
{
  int x, y, width, height, depth;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  m3d->pnodes = (MORPH_NODE_PROPERTIES ***)calloc(width, sizeof(MORPH_NODE_PROPERTIES **));
  if (!m3d->pnodes) ErrorExit(ERROR_NOMEMORY, "m3dAllocProperties: could not allocate node width");
  for (x = 0; x < width; x++) {
    m3d->pnodes[x] = (MORPH_NODE_PROPERTIES **)calloc(height, sizeof(MORPH_NODE_PROPERTIES *));
    if (!m3d->pnodes) ErrorExit(ERROR_NOMEMORY, "m3dAllocProperties: could not allocate node height %d", x);
    for (y = 0; y < height; y++) {
      m3d->pnodes[x][y] = (MORPH_NODE_PROPERTIES *)calloc(depth, sizeof(MORPH_NODE_PROPERTIES));
      if (!m3d->pnodes) ErrorExit(ERROR_NOMEMORY, "m3dAllocProperties: could not allocate node depth %d,%d", x, y);
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate and initialize a 3D morph structure.
------------------------------------------------------*/
static MORPH_3D *m3dalloc(int width, int height, int depth, float node_spacing)
{
  MORPH_3D *m3d;
  int x, y;

  m3d = (MORPH_3D *)calloc(1, sizeof(MORPH_3D));
  if (!m3d) ErrorExit(ERROR_NOMEMORY, "m3dalloc: could not allocate m3d");

  m3d->node_spacing = node_spacing;
  m3d->width = width;
  m3d->height = height;
  m3d->depth = depth;

  m3d->nodes = (MORPH_NODE ***)calloc(width, sizeof(MORPH_NODE **));
  if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dalloc: could not allocate node width");
  for (x = 0; x < width; x++) {
    m3d->nodes[x] = (MORPH_NODE **)calloc(height, sizeof(MORPH_NODE *));
    if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dalloc: could not allocate node height %d", x);
    for (y = 0; y < height; y++) {
      m3d->nodes[x][y] = (MORPH_NODE *)calloc(depth, sizeof(MORPH_NODE));
      if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dalloc: could not allocate node depth %d,%d", x, y);
    }
  }

  return (m3d);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Allocate and initialize a 3D morph structure.
------------------------------------------------------*/
static MORPH_3D *m3dInit(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, float node_spacing)
{
  MORPH_3D *m3d;
  int width, height, depth, x, y, i, j, k;
  float scale;
  VECTOR *v_X, *v_Y;
  MORPH_NODE *mn;
  MNP *mnp;

  v_X = VectorAlloc(4, MATRIX_REAL);
  v_Y = VectorAlloc(4, MATRIX_REAL);
  scale = mri_in->thick / node_spacing;

  /* put one node at the end of each row so expansion always interpolates */
  width = mri_in->width * scale + 1;
  height = mri_in->height * scale + 1;
  depth = mri_in->depth * scale + 1;
  m3d = (MORPH_3D *)calloc(1, sizeof(MORPH_3D));
  if (!m3d) ErrorExit(ERROR_NOMEMORY, "m3dInit: could not allocate m3d");

  m3d->node_spacing = node_spacing;
  m3d->width = width;
  m3d->height = height;
  m3d->depth = depth;
  m3d->mri_in = mri_in;
  m3d->mri_ref = mri_ref;
  m3d->lta = parms->lta;

  m3d->nodes = (MORPH_NODE ***)calloc(width, sizeof(MORPH_NODE **));
  if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dInit: could not allocate node width");
  for (x = 0; x < width; x++) {
    m3d->nodes[x] = (MORPH_NODE **)calloc(height, sizeof(MORPH_NODE *));
    if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dInit: could not allocate node height %d", x);
    for (y = 0; y < height; y++) {
      m3d->nodes[x][y] = (MORPH_NODE *)calloc(depth, sizeof(MORPH_NODE));
      if (!m3d->nodes) ErrorExit(ERROR_NOMEMORY, "m3dInit: could not allocate node depth %d,%d", x, y);
    }
  }

  if (!m3d->pnodes) m3dAllocProperties(m3d);

  /*parms->lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;*/

  /* now initialize node positions */
  v_X->rptr[4][1] = 1.0f /*/ mri_in->thick*/;
  for (k = 0; k < depth; k++) {
    V3_Z(v_X) = (float)k / scale;
    for (j = 0; j < height; j++) {
      V3_Y(v_X) = (float)j / scale;
      for (i = 0; i < width; i++) {
        V3_X(v_X) = (float)i / scale;
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        LTAtransformPoint(parms->lta, v_X, v_Y);
#if !USE_ORIGINAL_PROPERTIES
        mnp->ox = mn->x = V3_X(v_Y) * node_spacing;
        mnp->oy = mn->y = V3_Y(v_Y) * node_spacing;
        mnp->oz = mn->z = V3_Z(v_Y) * node_spacing;
#else
        mnp->ox = V3_X(v_X);
        mnp->oy = V3_Y(v_X);
        mnp->oz = V3_Z(v_X);
        mn->x = V3_X(v_Y);
        mn->y = V3_Y(v_Y);
        mn->z = V3_Z(v_Y);
#endif
      }
    }
  }

  m3dInitDistances(m3d);
  m3dInitAreas(m3d);
  VectorFree(&v_X);
  VectorFree(&v_Y);
  return (m3d);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Store the current metric properties of the m3d structure
          as the reference values to be used for computing distortion.
------------------------------------------------------*/
static int m3dStoreMetricProperties(MORPH_3D *m3d)
{
  MORPH_NODE *mn;
  MNP *mnp;
  int i, j, k, width, height, depth;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        if (i == width / 2 && j == height / 2 && k == depth / 2) {
          DiagBreak();
          i = width / 2;
          j = height / 2;
          k = depth / 2;
        }

        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        mnp->ox = mn->x;
        mnp->oy = mn->y;
        mnp->oz = mn->z;
      }
    }
  }
  m3dInitDistances(m3d);
  m3dInitAreas(m3d);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MIN_3D_PYR_WIDTH 32
MORPH_3D *MRI3Dmorph(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  int nlevels, i, max_levels, navgs, base_avgs;
  MRI *mri_in_pyramid[MAX_LEVELS], *mri_ref_pyramid[MAX_LEVELS];
  MORPH_3D *m3d_tmp, *m3d;
  char base_name[STRLEN];
  double dt, base_intensity, base_sigma, sigma, base_tol, max_thick;
  MRI *mri_in_transformed;
  MRI_SURFACE *mris_ref_skull, *mris_in_skull;
  float dx, dy, dz;
  MATRIX *m_tmp;

  mriNormalizeStds(mri_ref);
  if (parms->levels >= 0)
    max_levels = parms->levels;
  else
    max_levels = MAX_LEVELS;
  if (MRIfindNeck(mri_in, mri_in, 90, 120, NULL, 0, &parms->in_np) != NO_ERROR)
    parms->disable_neck = 1;
  else if (MRIfindNeck(mri_ref, mri_ref, 90, 120, NULL, 0, &parms->ref_np) != NO_ERROR)
    parms->disable_neck = 1;
  if (parms->disable_neck)

  {
    if (Gdiag & DIAG_SHOW) fprintf(stdout, "disabling neck discounting...\n");
    parms->ref_np.neck_x0 = parms->ref_np.neck_y0 = parms->ref_np.neck_z0 = 1000;
    parms->ref_np.neck_dx = parms->ref_np.neck_dy = parms->ref_np.neck_dz = 1;
    parms->in_np.neck_x0 = parms->in_np.neck_y0 = parms->in_np.neck_z0 = 1000;
    parms->in_np.neck_dx = parms->in_np.neck_dy = parms->in_np.neck_dz = 1;
  }
  else {
    fprintf(stdout,
            "moving spinal estimate from (%2.0f,%2.0f,%2.0f) to ",
            parms->in_np.neck_x0,
            parms->in_np.neck_y0,
            parms->in_np.neck_z0);
    parms->in_np.neck_x0 += parms->in_np.neck_dx * DISTANCE_BELOW_SPINAL_FUSION;
    parms->in_np.neck_y0 += parms->in_np.neck_dy * DISTANCE_BELOW_SPINAL_FUSION;
    parms->in_np.neck_z0 += parms->in_np.neck_dz * DISTANCE_BELOW_SPINAL_FUSION;
    parms->ref_np.neck_x0 += parms->ref_np.neck_dx * DISTANCE_BELOW_SPINAL_FUSION;
    parms->ref_np.neck_y0 += parms->ref_np.neck_dy * DISTANCE_BELOW_SPINAL_FUSION;
    parms->ref_np.neck_z0 += parms->ref_np.neck_dz * DISTANCE_BELOW_SPINAL_FUSION;
    fprintf(stdout, "(%2.0f,%2.0f,%2.0f)\n", parms->in_np.neck_x0, parms->in_np.neck_y0, parms->in_np.neck_z0);
  }

  openLogFile(parms);

  /* build Gaussian pyramid */
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "building Gaussian pyramid...\n");
  mri_in_pyramid[0] = mri_in;
  mri_ref_pyramid[0] = mri_ref;
  for (nlevels = 1; nlevels < max_levels; nlevels++) {
    if (mri_in_pyramid[nlevels - 1]->width <= MIN_3D_PYR_WIDTH) break;
    mri_in_pyramid[nlevels] = MRIreduceByte(mri_in_pyramid[nlevels - 1], NULL);
    mri_ref_pyramid[nlevels] = MRIreduceMeanAndStdByte(mri_ref_pyramid[nlevels - 1], NULL);
  }
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "morphing pyramid levels from coarse to fine...\n");

/* now morph each level, and use the resulting morph as the input to
   the next (finer) level.
*/
#define MIN_LEVEL 0
#define MAX_LEVEL nlevels - 1
  strcpy(base_name, parms->base_name);
  max_thick = mri_in_pyramid[MAX_LEVEL]->thick;
  m_tmp = parms->lta->xforms[0].m_L;
  parms->lta->xforms[0].m_L =
      MRIrasXformToVoxelXform(mri_in_pyramid[MAX_LEVEL], mri_ref_pyramid[MAX_LEVEL], m_tmp, NULL);
  MatrixFree(&m_tmp);
  m3d = m3dInit(mri_in_pyramid[MAX_LEVEL], mri_ref_pyramid[MAX_LEVEL], parms, max_thick);

  if (parms->morph_skull) /* apply radial morph to account for skull shapes */
  {
    MRI *mri_tmp;

    /* build representation of mri_in skull */
    mri_tmp = MRIcopy(parms->mri_in, NULL);
    MRIeraseNeck(mri_tmp, &parms->in_np);
    mri_in_transformed = LTAtransform(mri_tmp, NULL, parms->lta);
    if (Gdiag & DIAG_WRITE) MRIwriteImageViews(mri_in_transformed, "in_erased", IMAGE_SIZE);
    MRIfree(&mri_tmp);
    mris_in_skull = MRISshrinkWrapSkull(mri_in_transformed, parms);
    MRIfree(&mri_in_transformed);

    /* build representation of mri_ref skull */
    mri_tmp = MRIcopy(parms->mri_ref, NULL);
    MRIeraseNeck(mri_tmp, &parms->ref_np);
    if (Gdiag & DIAG_WRITE) MRIwriteImageViews(mri_tmp, "ref_erased", IMAGE_SIZE);
    if (Gdiag & DIAG_SHOW) fprintf(stdout, "building representations of the inner skull...\n");
    mris_ref_skull = MRISshrinkWrapSkull(mri_tmp, parms);
    MRIfree(&mri_tmp);

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRISwrite(mris_in_skull, "in_skull.geo");
      MRISwrite(mris_ref_skull, "ref_skull.geo");
    }

    /* translate so centers of skulls align */
    dx = mris_ref_skull->xctr - mris_in_skull->xctr;
    dy = mris_ref_skull->yctr - mris_in_skull->yctr;
    dz = mris_ref_skull->zctr - mris_in_skull->zctr;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      fprintf(stdout, "writing out representations of inner scalp\n");
      MRISwrite(mris_ref_skull, "ref_skull.geo");
      MRISwrite(mris_in_skull, "in_skull_postlinear.geo");
    }

    fprintf(stdout, "translating surface and morph by (%2.2f, %2.2f, %2.2f)\n", dx, dy, dz);
    MRIStranslate(mris_in_skull, dx, dy, dz);
    m3dTranslate(m3d, dx, dy, dz);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      fprintf(stdout, "writing out in translated \n");
      MRISwrite(mris_in_skull, "in_skull_posttrans.geo");
    }

    fprintf(stdout, "morphing skull shapes into register...\n");
    if (Gdiag & DIAG_WRITE) {
      MRI *mri_tmp;

      fprintf(stdout, "writing pre skull-normalized volumes...\n");
      MRIwriteImageViews(mri_in, "input", IMAGE_SIZE);
      mri_tmp = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
      MRIwriteImageViews(mri_tmp, "pre_skull_morph", IMAGE_SIZE);
      MRIfree(&mri_tmp);
    }
    m3dMorphSkull(m3d, mris_in_skull, mris_ref_skull, mri_ref_pyramid[MAX_LEVEL]);
    if (Gdiag & DIAG_WRITE) {
      MRI *mri_tmp;

      fprintf(stdout, "writing post skull-normalized images...\n");
      mri_tmp = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
      MRIwriteImageViews(mri_tmp, "post_skull_morph", IMAGE_SIZE);
      MRIfree(&mri_tmp);
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        fprintf(stdout, "writing post skull morphing surface and image views...\n");
        MRISwrite(mris_in_skull, "post_skull_morph.geo");
      }
    }
    m3RecomputeTranslation(m3d, mri_in_pyramid[MIN_LEVEL], mri_ref_pyramid[MIN_LEVEL], parms);
    MRISfree(&mris_ref_skull);
    MRISfree(&mris_in_skull);
  } /* end of skull morphing */

  m3d->mri_in = mri_in;
  m3d->mri_ref = mri_ref;
  m3dStoreMetricProperties(m3d);
  MRIfree(&parms->mri_in);
  MRIfree(&parms->mri_ref);

  dt = parms->dt;
  navgs = base_avgs = parms->navgs;
  base_tol = parms->tol;
  base_intensity = parms->l_intensity;
  base_sigma = parms->sigma;
  for (i = MAX_LEVEL; i >= MIN_LEVEL; i--) {
#if !SCALE_INVARIANT
/*    parms->dt = dt / (mri_in_pyramid[i]->thick) ;*/
#else
    parms->dt = dt / (mri_in_pyramid[i]->thick);
#endif
    /*    parms->l_intensity = l_intensity * (sigma*sigma+1.0f) ;*/

    int req = snprintf(parms->base_name, 100, "%s_", base_name);
    if( req >= 100 ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    parms->tol = mri_in_pyramid[MIN_LEVEL]->thick * base_tol / mri_in_pyramid[i]->thick;
    for (sigma = base_sigma; sigma >= 0; sigma /= 4) {
      if (sigma < 1) sigma = 0;
#if SCALE_INVARIANT
      parms->dt = dt * (sigma + 1) / (mri_in_pyramid[i]->thick);
      parms->l_intensity = base_intensity * (sigma + 1) * (max_thick / mri_in_pyramid[i]->thick);
#endif
      /*
        if negative areas have been created (i.e. crossing of coordinate
        lines), then remove them by doing an integration epoch with
        no blurring, then continue with the normal integration.
        */
      if (m3d->neg > 0 && sigma > 0) {
        double tol = parms->tol;
        parms->tol = 0; /* only until negative nodes are fixed */
        fprintf(
            stdout, "unfolding lattice, dt=%2.3f, l_int=%2.2f, tol=%2.2f\n", parms->dt, parms->l_intensity, parms->tol);
        if ((Gdiag & DIAG_WRITE) && parms->log_fp)
          fprintf(parms->log_fp,
                  "unfolding lattice, dt=%2.2f, l_int=%2.2f, tol=%2.2f.\n",
                  parms->dt,
                  parms->l_intensity,
                  parms->tol);
        parms->sigma = 0;
        m3dAlignPyramidLevel(mri_in_pyramid[i], mri_ref_pyramid[i], mri_ref_pyramid[i], parms, m3d);
        parms->tol = tol;
      }
      fprintf(stdout,
              "aligning pyramid level %d, dt=%2.3f, sigma=%2.2f, l_int=%2.2f,"
              "tol=%2.2f, neg=%d\n",
              i,
              parms->dt,
              sigma,
              parms->l_intensity,
              parms->tol,
              m3d->neg);
      if ((Gdiag & DIAG_WRITE) && parms->log_fp)
        fprintf(parms->log_fp,
                "aligning pyramid level %d, dt=%2.2f, sigma=%2.2f, l_int=%2.2f,"
                "tol=%2.2f.\n",
                i,
                parms->dt,
                sigma,
                parms->l_intensity,
                parms->tol);
      parms->sigma = sigma;
      m3dAlignPyramidLevel(mri_in_pyramid[i], mri_ref_pyramid[i], mri_ref_pyramid[i], parms, m3d);
      if (FZERO(sigma)) break;
    }

    if (i > MIN_LEVEL) {
      MRIfree(&mri_in_pyramid[i]);
      MRIfree(&mri_ref_pyramid[i]);
      m3d_tmp = m3dscaleUp2(m3d, NULL, parms);
      m3d = m3d_tmp;
      m3d->mri_in = mri_in;
      m3d->mri_ref = mri_ref;
      m3dComputeMetricProperties(m3d);
    }
  }


  strcpy(parms->base_name, base_name);

/* free Gaussian pyramid */
  if (parms->log_fp) {
    fclose(parms->log_fp);
    parms->log_fp = NULL;
  }

  if (Gdiag & DIAG_WRITE) {
    MRI *mri_tmp;

    fprintf(stdout, "writing final morphed volume image views...\n");
    mri_tmp = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
    MRIwriteImageViews(mri_tmp, "morphed", IMAGE_SIZE);
    MRIfree(&mri_tmp);
  }
  return (m3d);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          For each point in the morph, compute the point to which it
          maps. Then use a linear interpolation scheme to build the
          morphed image.

          This should really use floating point volumes, or at
          the very least short volumes, but unfortunatly I
          can't fit them in memory, so it all has to be
          done in unsigned char volumes which forces things
          to be scaled up and then back down in various places.
------------------------------------------------------*/
MRI *MRIapply3DMorph(MRI *mri_in, MORPH_3D *m3d, MRI *mri_morphed)
{
  int width, height, depth, x, y, z, xm1, ym1, zm1, xp1, yp1, zp1;
  float xd, yd, zd, dx, dy, dz, thick;
  double weight, orig_val, val;
  MRI *mri_weights, *mri_ctrl, *mri_s_morphed;

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  thick = mri_in->thick;
#define SCALE_FACTOR 25.0f /* so we can use byte images for float values */

  mri_weights = MRIalloc(width, height, depth, MRI_UCHAR);
  mri_s_morphed = MRIalloc(width, height, depth, MRI_SHORT);
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        /* compute voxel coordinates of this morph point */
        MRIsample3Dmorph(m3d, (float)x * thick, (float)y * thick, (float)z * thick, &xd, &yd, &zd);
        xd /= thick;
        yd /= thick;
        zd /= thick; /* voxel coords */
        orig_val = MRIgetVoxVal(mri_in, x, y, z, 0);
        if (orig_val > 40) DiagBreak();

        /* now use trilinear interpolation */
        xm1 = (int)floor(xd);
        ym1 = (int)floor(yd);
        zm1 = (int)floor(zd);
        xp1 = xm1 + 1;
        yp1 = ym1 + 1;
        zp1 = zm1 + 1;

        /* make sure they are within bounds */
        xm1 = mri_in->xi[xm1];
        ym1 = mri_in->xi[ym1];
        zm1 = mri_in->xi[zm1];
        xp1 = mri_in->xi[xp1];
        yp1 = mri_in->xi[yp1];
        zp1 = mri_in->xi[zp1];

        dx = xd - xm1;
        dy = yd - ym1;
        dz = zd - zm1;

        /* now compute contribution to each of 8 nearest voxels */
        weight = (1 - dx) * (1 - dy) * (1 - dz);
        if (xm1 == 121 && ym1 == 172 && zm1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xm1, ym1, zm1) += weight * orig_val;
        MRIvox(mri_weights, xm1, ym1, zm1) += nint(weight * SCALE_FACTOR);

        weight = (dx) * (1 - dy) * (1 - dz);
        if (xm1 == 121 && ym1 == 172 && zm1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xp1, ym1, zm1) += weight * orig_val;
        MRIvox(mri_weights, xp1, ym1, zm1) += nint(weight * SCALE_FACTOR);

        weight = (1 - dx) * (dy) * (1 - dz);
        if (xm1 == 121 && yp1 == 172 && zm1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xm1, yp1, zm1) += weight * orig_val;
        MRIvox(mri_weights, xm1, yp1, zm1) += nint(weight * SCALE_FACTOR);

        weight = (1 - dx) * (1 - dy) * (dz);
        if (xm1 == 121 && ym1 == 172 && zp1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xm1, ym1, zp1) += weight * orig_val;
        MRIvox(mri_weights, xm1, ym1, zp1) += nint(weight * SCALE_FACTOR);

        weight = (dx) * (dy) * (1 - dz);
        if (xp1 == 121 && yp1 == 172 && zm1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xp1, yp1, zm1) += weight * orig_val;
        MRIvox(mri_weights, xp1, yp1, zm1) += nint(weight * SCALE_FACTOR);

        weight = (dx) * (1 - dy) * (dz);
        if (xp1 == 121 && ym1 == 172 && zp1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xp1, ym1, zp1) += weight * orig_val;
        MRIvox(mri_weights, xp1, ym1, zp1) += nint(weight * SCALE_FACTOR);

        weight = (1 - dx) * (dy) * (dz);
        if (xm1 == 121 && yp1 == 172 && zp1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xm1, yp1, zp1) += weight * orig_val;
        MRIvox(mri_weights, xm1, yp1, zp1) += nint(weight * SCALE_FACTOR);

        weight = (dx) * (dy) * (dz);
        if (xp1 == 121 && yp1 == 172 && zp1 == 127) DiagBreak();
        MRISvox(mri_s_morphed, xp1, yp1, zp1) += weight * orig_val;
        MRIvox(mri_weights, xp1, yp1, zp1) += nint(weight * SCALE_FACTOR);
      }
    }
  }

  /* now normalize weights and values */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == 121 && y == 172 && z == 127) DiagBreak();
        weight = (float)MRIvox(mri_weights, x, y, z) / SCALE_FACTOR;
        if (!FZERO(weight)) {
          val = (double)MRISvox(mri_s_morphed, x, y, z) / weight;
          if (val > 255.0) val = 255.0;
          MRISvox(mri_s_morphed, x, y, z) = (short)nint(val);
        }
      }
    }
  }

  /* copy from short image to BUFTYPE one */
  if (!mri_morphed)
    mri_morphed = MRIclone(mri_in, NULL);
  else
    MRIclear(mri_morphed);
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        MRIvox(mri_morphed, x, y, z) = (BUFTYPE)MRISvox(mri_s_morphed, x, y, z);
      }
    }
  }

  MRIfree(&mri_s_morphed);

/* run soap bubble to fill in remaining holes */
  mri_ctrl = mri_weights;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == 121 && y == 172 && z == 127) DiagBreak();
        weight = (float)MRIvox(mri_weights, x, y, z) / SCALE_FACTOR;
        if (weight > .1)
          MRIvox(mri_ctrl, x, y, z) = 1;
        else
          MRIvox(mri_ctrl, x, y, z) = 0;
      }
    }
  }

  MRIbuildVoronoiDiagram(mri_morphed, mri_ctrl, mri_morphed);
  MRIsoapBubble(mri_morphed, mri_ctrl, mri_morphed, 5, -1);

  MRIfree(&mri_ctrl);
  return (mri_morphed);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Apply the previously computed morph which takes
          mri_in into mri_ref space, and use it to generate a
          warping of mri_ref into mri_in space. This is much
          simpler than the inverse mapping as the sampling is
          uniform in this direction.
------------------------------------------------------*/
MRI *MRIapplyInverse3DMorph(MRI *mri_ref, MORPH_3D *m3d, MRI *mri_morphed)
{
  int width, height, depth, x, y, z;
  float xd, yd, zd;
  double val;

  width = mri_ref->width;
  height = mri_ref->height;
  depth = mri_ref->depth;
  if (!mri_morphed) mri_morphed = MRIclone(mri_ref, NULL);
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        MRIsample3Dmorph(
            m3d, (float)x * mri_ref->thick, (float)y * mri_ref->thick, (float)z * mri_ref->thick, &xd, &yd, &zd);
        xd /= mri_ref->thick;
        yd /= mri_ref->thick;
        zd /= mri_ref->thick;
        if (xd > -1 && yd > -1 && zd > 0 && xd < width && yd < height && zd < depth)
          MRIsampleVolume(mri_ref, xd, yd, zd, &val);
        else
          val = 0.0;
        switch (mri_morphed->type) {
          case MRI_UCHAR:
            MRIvox(mri_morphed, x, y, z) = val;
            break;
          case MRI_FLOAT:
            MRIFvox(mri_morphed, x, y, z) = val;
            break;
          default:
            ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIinverse3DMorph: unsupported volume type %d", mri_morphed->type));
            break;
        }
      }
    }
  }
  return (mri_morphed);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIsample3DmorphOrig(MORPH_3D *m3d, float x, float y, float z, float *pxd, float *pyd, float *pzd)
{
  int xm, xp, ym, yp, zm, zp, width, height, depth;
  float xmd, ymd, zmd, xpd, ypd, zpd; /* d's are distances */

  /* convert to 3d morph index space */
  x /= m3d->node_spacing;
  y /= m3d->node_spacing;
  z /= m3d->node_spacing;
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  if (x >= width) x = width - 1.0;
  if (y >= height) y = height - 1.0;
  if (z >= depth) z = depth - 1.0;
  if (x < 0.0) x = 0.0;
  if (y < 0.0) y = 0.0;
  if (z < 0.0) z = 0.0;

  xm = MAX((int)x, 0);
  xp = MIN(width - 1, xm + 1);
  ym = MAX((int)y, 0);
  yp = MIN(height - 1, ym + 1);
  zm = MAX((int)z, 0);
  zp = MIN(depth - 1, zm + 1);

  xmd = x - (float)xm;
  ymd = y - (float)ym;
  zmd = z - (float)zm;
  xpd = (1.0f - xmd);
  ypd = (1.0f - ymd);
  zpd = (1.0f - zmd);

  *pxd = xpd * ypd * zpd * m3d->pnodes[xm][ym][zm].ox + xpd * ypd * zmd * m3d->pnodes[xm][ym][zp].ox +
         xpd * ymd * zpd * m3d->pnodes[xm][yp][zm].ox + xpd * ymd * zmd * m3d->pnodes[xm][yp][zp].ox +
         xmd * ypd * zpd * m3d->pnodes[xp][ym][zm].ox + xmd * ypd * zmd * m3d->pnodes[xp][ym][zp].ox +
         xmd * ymd * zpd * m3d->pnodes[xp][yp][zm].ox + xmd * ymd * zmd * m3d->pnodes[xp][yp][zp].ox;
  *pyd = xpd * ypd * zpd * m3d->pnodes[xm][ym][zm].oy + xpd * ypd * zmd * m3d->pnodes[xm][ym][zp].oy +
         xpd * ymd * zpd * m3d->pnodes[xm][yp][zm].oy + xpd * ymd * zmd * m3d->pnodes[xm][yp][zp].oy +
         xmd * ypd * zpd * m3d->pnodes[xp][ym][zm].oy + xmd * ypd * zmd * m3d->pnodes[xp][ym][zp].oy +
         xmd * ymd * zpd * m3d->pnodes[xp][yp][zm].oy + xmd * ymd * zmd * m3d->pnodes[xp][yp][zp].oy;
  *pzd = xpd * ypd * zpd * m3d->pnodes[xm][ym][zm].oz + xpd * ypd * zmd * m3d->pnodes[xm][ym][zp].oz +
         xpd * ymd * zpd * m3d->pnodes[xm][yp][zm].oz + xpd * ymd * zmd * m3d->pnodes[xm][yp][zp].oz +
         xmd * ypd * zpd * m3d->pnodes[xp][ym][zm].oz + xmd * ypd * zmd * m3d->pnodes[xp][ym][zp].oz +
         xmd * ymd * zpd * m3d->pnodes[xp][yp][zm].oz + xmd * ymd * zmd * m3d->pnodes[xp][yp][zp].oz;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int MRIsample3Dmorph(MORPH_3D *m3d, float x, float y, float z, float *pxd, float *pyd, float *pzd)
{
  int xm, xp, ym, yp, zm, zp, width, height, depth;
  float xmd, ymd, zmd, xpd, ypd, zpd; /* d's are distances */

  /* convert to 3d morph index space */
  x /= m3d->node_spacing;
  y /= m3d->node_spacing;
  z /= m3d->node_spacing;
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  if (x >= width) x = width - 1.0;
  if (y >= height) y = height - 1.0;
  if (z >= depth) z = depth - 1.0;
  if (x < 0.0) x = 0.0;
  if (y < 0.0) y = 0.0;
  if (z < 0.0) z = 0.0;

  xm = MAX((int)x, 0);
  xp = MIN(width - 1, xm + 1);
  ym = MAX((int)y, 0);
  yp = MIN(height - 1, ym + 1);
  zm = MAX((int)z, 0);
  zp = MIN(depth - 1, zm + 1);

  xmd = x - (float)xm;
  ymd = y - (float)ym;
  zmd = z - (float)zm;
  xpd = (1.0f - xmd);
  ypd = (1.0f - ymd);
  zpd = (1.0f - zmd);

  *pxd = xpd * ypd * zpd * m3d->nodes[xm][ym][zm].x + xpd * ypd * zmd * m3d->nodes[xm][ym][zp].x +
         xpd * ymd * zpd * m3d->nodes[xm][yp][zm].x + xpd * ymd * zmd * m3d->nodes[xm][yp][zp].x +
         xmd * ypd * zpd * m3d->nodes[xp][ym][zm].x + xmd * ypd * zmd * m3d->nodes[xp][ym][zp].x +
         xmd * ymd * zpd * m3d->nodes[xp][yp][zm].x + xmd * ymd * zmd * m3d->nodes[xp][yp][zp].x;
  *pyd = xpd * ypd * zpd * m3d->nodes[xm][ym][zm].y + xpd * ypd * zmd * m3d->nodes[xm][ym][zp].y +
         xpd * ymd * zpd * m3d->nodes[xm][yp][zm].y + xpd * ymd * zmd * m3d->nodes[xm][yp][zp].y +
         xmd * ypd * zpd * m3d->nodes[xp][ym][zm].y + xmd * ypd * zmd * m3d->nodes[xp][ym][zp].y +
         xmd * ymd * zpd * m3d->nodes[xp][yp][zm].y + xmd * ymd * zmd * m3d->nodes[xp][yp][zp].y;
  *pzd = xpd * ypd * zpd * m3d->nodes[xm][ym][zm].z + xpd * ypd * zmd * m3d->nodes[xm][ym][zp].z +
         xpd * ymd * zpd * m3d->nodes[xm][yp][zm].z + xpd * ymd * zmd * m3d->nodes[xm][yp][zp].z +
         xmd * ypd * zpd * m3d->nodes[xp][ym][zm].z + xmd * ypd * zmd * m3d->nodes[xp][ym][zp].z +
         xmd * ymd * zpd * m3d->nodes[xp][yp][zm].z + xmd * ymd * zmd * m3d->nodes[xp][yp][zp].z;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Free the properties field of a 3d morph structure
------------------------------------------------------*/
int MRI3DmorphFreeProperties(MORPH_3D *m3d)
{
  int x, y;

  if (m3d->pnodes) {
    for (x = 0; x < m3d->width; x++) {
      for (y = 0; y < m3d->height; y++) free(m3d->pnodes[x][y]);
      free(m3d->pnodes[x]);
    }
    free(m3d->pnodes);
    m3d->pnodes = NULL;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Free a 3d morph structure and set the pointer to NULL.
------------------------------------------------------*/
int MRI3DmorphFree(MORPH_3D **pm3d)
{
  MORPH_3D *m3d = *pm3d;
  int x, y;

  *pm3d = NULL;
  if (!m3d) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRI3DmorphFree: NULL pointer!"));

  for (x = 0; x < m3d->width; x++) {
    for (y = 0; y < m3d->height; y++) free(m3d->nodes[x][y]);
    free(m3d->nodes[x]);
  }
  free(m3d->nodes);

  if (m3d->pnodes) {
    for (x = 0; x < m3d->width; x++) {
      for (y = 0; y < m3d->height; y++) free(m3d->pnodes[x][y]);
      free(m3d->pnodes[x]);
    }
    free(m3d->pnodes);
  }
  free(m3d);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the directional derivative of the correlation term
          at the input image position xin,yin,zin which currently
          corresponds to the reference image position xref,yref,zref.
------------------------------------------------------*/
static int m3dCorrelationTerm(MRI *mri_in,
                              MRI *mri_ref,
                              MRI *mri_ref_blur,
                              double l_intensity,
                              int xin,
                              int yin,
                              int zin,
                              double xref,
                              double yref,
                              double zref,
                              double *pdx,
                              double *pdy,
                              double *pdz)
{
  double ref_val, in_val, delta, Tdx, Tdy, Tdz, len, ref_std, mean_std;

  mean_std = mri_ref->mean;
  if (FZERO(l_intensity)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }
  if (xref > -1 && yref > -1 && zref > 0 && xref < mri_ref->width && yref < mri_ref->height && zref < mri_ref->depth) {
    MRIsampleVolume(mri_ref, xref, yref, zref, &ref_val);
    MRIsampleVolumeFrame(mri_ref, xref, yref, zref, 1, &ref_std);
    if (DZERO(ref_std)) ref_std = 1.0;
  }
  else {
    ref_std = mean_std;
    ref_val = 0.0;
  }
  ref_std /= mean_std;
  in_val = (double)MRIvox(mri_in, xin, yin, zin);
  delta = -l_intensity * (ref_val - in_val) / ref_std; /* -delta */
  MRIsampleVolumeGradient(mri_ref_blur, xref, yref, zref, &Tdx, &Tdy, &Tdz);
  len = sqrt(Tdx * Tdx + Tdy * Tdy + Tdz * Tdz);
  if (!FZERO(len)) {
    Tdx /= len;
    Tdy /= len;
    Tdz /= len;
  }
  *pdx = Tdx * delta;
  *pdy = Tdy * delta;
  *pdz = Tdz * delta;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the directional derivative of the distance term
          w.r.t the i,j,kth node.
------------------------------------------------------*/
static int m3dDistanceTerm(MORPH_3D *m3d, double l_distance, int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  int n, num;
  double dx, dy, dz, dist, delta, xgrad, ygrad, zgrad;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;

  if (FZERO(l_distance)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }

  mn = &m3d->nodes[i][j][k];
  mnp = &m3d->pnodes[i][j][k];
  for (xgrad = ygrad = zgrad = 0.0, num = n = 0; n < NEIGHBORS; n++) {
    switch (n) {
      default:
      case 0: /* i-1 */
        if (i > 0)
          mn_nbr = &m3d->nodes[i - 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 1: /* i+1 */
        if (i < m3d->width - 1)
          mn_nbr = &m3d->nodes[i + 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 2: /* j-1 */
        if (j > 0)
          mn_nbr = &m3d->nodes[i][j - 1][k];
        else
          mn_nbr = NULL;
        break;
      case 3: /* j+1 */
        if (j < m3d->height - 1)
          mn_nbr = &m3d->nodes[i][j + 1][k];
        else
          mn_nbr = NULL;
        break;
      case 4: /* k-1 */
        if (k > 0)
          mn_nbr = &m3d->nodes[i][j][k - 1];
        else
          mn_nbr = NULL;
        break;
      case 5: /* k+1 */
        if (k < m3d->depth - 1)
          mn_nbr = &m3d->nodes[i][j][k + 1];
        else
          mn_nbr = NULL;
        break;
    }
    if (mn_nbr) {
      num++;
      dx = mn_nbr->x - mn->x;
      dy = mn_nbr->y - mn->y;
      dz = mn_nbr->z - mn->z;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      delta = dist - mnp->orig_dist[n];
      if (!FZERO(dist)) /* make it a unit vector */
      {
        dx /= dist;
        dx /= dist;
        dz /= dist;
      }
#if SCALE_INVARIANT
      delta /= m3d->node_spacing;
#endif
      xgrad += delta * dx;
      ygrad += delta * dy;
      zgrad += delta * dz;
      if (!finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad) || (fabs(dx) > 1e5) || (fabs(dy) > 1e5) ||
          (fabs(dz) > 1e5))
        DiagBreak();
    }
  }

  if (num) l_distance /= (float)num;
  *pdx = l_distance * xgrad;
  *pdy = l_distance * ygrad;
  *pdz = l_distance * zgrad;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the directional derivative of the area term
          w.r.t the i,j,kth node.
------------------------------------------------------*/
#define AREA_NEIGHBORS 8
static int m3dAreaTerm(MORPH_3D *m3d, double l_area, int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  MORPH_NODE *mn, *mni, *mnj, *mnk;
  MNP *mnp, *mnpi, *mnpj, *mnpk;
  float delta, node_spacing, n3 /*, odx, ody, odz*/;
  int n, width, height, depth, num, invert;
  static VECTOR *v_i = NULL, *v_j, *v_k, *v_k_x_j, *v_j_x_i, *v_i_x_k, *v_grad, *v_tmp;


  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  if (FZERO(l_area)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }

  node_spacing = m3d->node_spacing;
#if SCALE_INVARIANT
  n3 = node_spacing * node_spacing * node_spacing;
#else
  n3 = 1.0f;
#endif
  if (!v_i) /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL);
    v_j = VectorAlloc(3, MATRIX_REAL);
    v_k = VectorAlloc(3, MATRIX_REAL);
    v_grad = VectorAlloc(3, MATRIX_REAL);
    v_k_x_j = VectorAlloc(3, MATRIX_REAL);
    v_j_x_i = VectorAlloc(3, MATRIX_REAL);
    v_i_x_k = VectorAlloc(3, MATRIX_REAL);
    v_tmp = VectorAlloc(3, MATRIX_REAL);
  }
  else {
    V3_CLEAR(v_grad);
  }

  for (num = n = 0; n < AREA_NEIGHBORS; n++) {
    /* assign mn pointers to appropriate nodes */
    invert = 1; /* assume non-inverted coordinate system */
    switch (n) {
      default:
      case 4:
        if ((i == 0) || (j == 0) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        mni = &m3d->nodes[i - 1][j][k];
        mnj = &m3d->nodes[i][j - 1][k];
        mnk = &m3d->nodes[i][j][k - 1];
        mnpi = &m3d->pnodes[i - 1][j][k];
        mnpj = &m3d->pnodes[i][j - 1][k];
        mnpk = &m3d->pnodes[i][j][k - 1];
        break;
      case 0: /* first do central node */
        if ((i + 1 >= width) || (j + 1 >= height) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        mni = &m3d->nodes[i + 1][j][k];
        mnj = &m3d->nodes[i][j + 1][k];
        mnk = &m3d->nodes[i][j][k + 1];
        mnpi = &m3d->pnodes[i + 1][j][k];
        mnpj = &m3d->pnodes[i][j + 1][k];
        mnpk = &m3d->pnodes[i][j][k + 1];
        break;
      case 5: /*  i+1 */
        if ((i + 1 >= width) || (j == 0) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i + 1][j][k];
        mnp = &m3d->pnodes[i + 1][j][k];
        mni = &m3d->nodes[i][j][k];
        mnj = &m3d->nodes[i + 1][j - 1][k];
        mnk = &m3d->nodes[i + 1][j][k - 1];
        mnpi = &m3d->pnodes[i][j][k];
        mnpj = &m3d->pnodes[i + 1][j - 1][k];
        mnpk = &m3d->pnodes[i + 1][j][k - 1];
        break;
      case 1: /*  i-1 */
        if ((i == 0) || (j + 1 >= height) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i - 1][j][k];
        mnp = &m3d->pnodes[i - 1][j][k];
        mni = &m3d->nodes[i][j][k];
        mnj = &m3d->nodes[i - 1][j + 1][k];
        mnk = &m3d->nodes[i - 1][j][k + 1];
        mnpi = &m3d->pnodes[i][j][k];
        mnpj = &m3d->pnodes[i - 1][j + 1][k];
        mnpk = &m3d->pnodes[i - 1][j][k + 1];
        break;
      case 6: /* j+1 */
        if ((i == 0) || (j + 1 >= height) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j + 1][k];
        mnp = &m3d->pnodes[i][j + 1][k];
        mni = &m3d->nodes[i - 1][j + 1][k];
        mnj = &m3d->nodes[i][j][k];
        mnk = &m3d->nodes[i][j + 1][k - 1];
        mnpi = &m3d->pnodes[i - 1][j + 1][k];
        mnpj = &m3d->pnodes[i][j][k];
        mnpk = &m3d->pnodes[i][j + 1][k - 1];
        break;
      case 2: /* j-1 */
        if ((i + 1 >= width) || (j == 0) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i][j - 1][k];
        mnp = &m3d->pnodes[i][j - 1][k];
        mni = &m3d->nodes[i + 1][j - 1][k];
        mnj = &m3d->nodes[i][j][k];
        mnk = &m3d->nodes[i][j - 1][k + 1];
        mnpi = &m3d->pnodes[i + 1][j - 1][k];
        mnpj = &m3d->pnodes[i][j][k];
        mnpk = &m3d->pnodes[i][j - 1][k + 1];
        break;
      case 7: /* k+1 */
        if ((i == 0) || (j == 0) || (k + 1 >= depth)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j][k + 1];
        mnp = &m3d->pnodes[i][j][k + 1];
        mni = &m3d->nodes[i - 1][j][k + 1];
        mnj = &m3d->nodes[i][j - 1][k + 1];
        mnk = &m3d->nodes[i][j][k];
        mnpi = &m3d->pnodes[i - 1][j][k + 1];
        mnpj = &m3d->pnodes[i][j - 1][k + 1];
        mnpk = &m3d->pnodes[i][j][k];
        break;
      case 3: /* k-1 */
        if ((i + 1 >= width) || (j + 1 >= height) || (k == 0)) continue;
        mn = &m3d->nodes[i][j][k - 1];
        mnp = &m3d->pnodes[i][j][k - 1];
        mni = &m3d->nodes[i + 1][j][k - 1];
        mnj = &m3d->nodes[i][j + 1][k - 1];
        mnk = &m3d->nodes[i][j][k];
        mnpi = &m3d->pnodes[i + 1][j][k - 1];
        mnpj = &m3d->pnodes[i][j + 1][k - 1];
        mnpk = &m3d->pnodes[i][j][k];
        break;
    }

    num++;

    /* compute cross products and area delta */
    MN_SUB(mni, mn, v_i);
    MN_SUB(mnj, mn, v_j);
    MN_SUB(mnk, mn, v_k);


    delta = invert * (mnp->orig_area - mnp->area) / n3;
    if (mnp->area < 0) DiagBreak();
    if (delta > 100000) DiagBreak();

    /* compute cross-products and add the appropriate
       (i.e. scaled by area difference) cross-products to the gradient */
    switch (n) {
      default:
      case 4: /* central node in inverted coordinate system */
      case 0: /* first do central node */
        V3_CROSS_PRODUCT(v_k, v_j, v_k_x_j);
        V3_CROSS_PRODUCT(v_j, v_i, v_j_x_i);
        V3_CROSS_PRODUCT(v_i, v_k, v_i_x_k);
        V3_ADD(v_k_x_j, v_j_x_i, v_tmp);
        V3_ADD(v_i_x_k, v_tmp, v_tmp);
        V3_SCALAR_MUL(v_tmp, delta, v_tmp);
        break;
      case 5: /*  i+1 */
      case 1: /*  i-1 */
        V3_CROSS_PRODUCT(v_k, v_j, v_k_x_j);
        V3_SCALAR_MUL(v_k_x_j, -delta, v_tmp);
        break;
      case 6: /* j+1 */
      case 2: /* j-1 */
        V3_CROSS_PRODUCT(v_i, v_k, v_i_x_k);
        V3_SCALAR_MUL(v_i_x_k, -delta, v_tmp);
        break;
      case 7: /* k+1 */
      case 3: /* k-1 */
        V3_CROSS_PRODUCT(v_j, v_i, v_j_x_i);
        V3_SCALAR_MUL(v_j_x_i, -delta, v_tmp);
        break;
    }
    V3_ADD(v_tmp, v_grad, v_grad);
  }
  if (num) {
    l_area /= (float)num;
    *pdx = l_area * V3_X(v_grad);
    *pdy = l_area * V3_Y(v_grad);
    *pdz = l_area * V3_Z(v_grad);
  }
  else {
    /* if on the border st no area term can be computed, move in a
       quasi-rigid fashion by applying strong distance term.
    */
    DiagBreak();
    fprintf(stdout, "node (%d,%d,%d) has no areal term!\n", i, j, k);
    /*    return(m3dDistanceTerm(m3d, 10*l_area, i, j, k, pdx, pdy, pdz)) ;*/
  }
  if (!std::isfinite(*pdx) || !std::isfinite(*pdy) || !std::isfinite(*pdz)) DiagBreak();
  return (NO_ERROR);
}
static int m3dNonlinearAreaTerm(
    MORPH_3D *m3d, double l_nlarea, int i, int j, int k, double *pdx, double *pdy, double *pdz, MORPH_PARMS *parms)
{
  MORPH_NODE *mn, *mni, *mnj, *mnk;
  MNP *mnp;
  float delta, delta_scale, ratio, node_spacing, n3;
  int n, width, height, depth, num, invert;
  static VECTOR *v_i = NULL, *v_j, *v_k, *v_k_x_j, *v_j_x_i, *v_i_x_k, *v_grad, *v_tmp;
  double exponent;

  if (FZERO(l_nlarea)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  node_spacing = m3d->node_spacing;
#if SCALE_INVARIANT
  n3 = node_spacing * node_spacing * node_spacing;
#else
  n3 = 1.0f;
#endif
  if (!v_i) /* initialize */
  {
    v_i = VectorAlloc(3, MATRIX_REAL);
    v_j = VectorAlloc(3, MATRIX_REAL);
    v_k = VectorAlloc(3, MATRIX_REAL);
    v_grad = VectorAlloc(3, MATRIX_REAL);
    v_k_x_j = VectorAlloc(3, MATRIX_REAL);
    v_j_x_i = VectorAlloc(3, MATRIX_REAL);
    v_i_x_k = VectorAlloc(3, MATRIX_REAL);
    v_tmp = VectorAlloc(3, MATRIX_REAL);
  }
  else {
    V3_CLEAR(v_grad);
  }

  for (num = n = 0; n < AREA_NEIGHBORS; n++) {
    /* assign mn pointers to appropriate nodes */
    invert = 1;
    switch (n) {
      default:
      case 4:
        if ((i == 0) || (j == 0) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        mni = &m3d->nodes[i - 1][j][k];
        mnj = &m3d->nodes[i][j - 1][k];
        mnk = &m3d->nodes[i][j][k - 1];
        break;
      case 0: /* first do central node */
        if ((i + 1 >= width) || (j + 1 >= height) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        mni = &m3d->nodes[i + 1][j][k];
        mnj = &m3d->nodes[i][j + 1][k];
        mnk = &m3d->nodes[i][j][k + 1];
        break;
      case 5: /*  i+1 */
        if ((i + 1 >= width) || (j == 0) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i + 1][j][k];
        mnp = &m3d->pnodes[i + 1][j][k];
        mni = &m3d->nodes[i][j][k];
        mnj = &m3d->nodes[i + 1][j - 1][k];
        mnk = &m3d->nodes[i + 1][j][k - 1];
        break;
      case 1: /*  i-1 */
        if ((i == 0) || (j + 1 >= height) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i - 1][j][k];
        mnp = &m3d->pnodes[i - 1][j][k];
        mni = &m3d->nodes[i][j][k];
        mnj = &m3d->nodes[i - 1][j + 1][k];
        mnk = &m3d->nodes[i - 1][j][k + 1];
        break;
      case 6: /* j+1 */
        if ((i == 0) || (j + 1 >= height) || (k == 0)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j + 1][k];
        mnp = &m3d->pnodes[i][j + 1][k];
        mni = &m3d->nodes[i - 1][j + 1][k];
        mnj = &m3d->nodes[i][j][k];
        mnk = &m3d->nodes[i][j + 1][k - 1];
        break;
      case 2: /* j-1 */
        if ((i + 1 >= width) || (j == 0) || (k + 1 >= depth)) continue;
        mn = &m3d->nodes[i][j - 1][k];
        mnp = &m3d->pnodes[i][j - 1][k];
        mni = &m3d->nodes[i + 1][j - 1][k];
        mnj = &m3d->nodes[i][j][k];
        mnk = &m3d->nodes[i][j - 1][k + 1];
        break;
      case 7: /* k+1 */
        if ((i == 0) || (j == 0) || (k + 1 >= depth)) continue;
        invert = -1;
        mn = &m3d->nodes[i][j][k + 1];
        mnp = &m3d->pnodes[i][j][k + 1];
        mni = &m3d->nodes[i - 1][j][k + 1];
        mnj = &m3d->nodes[i][j - 1][k + 1];
        mnk = &m3d->nodes[i][j][k];
        break;
      case 3: /* k-1 */
        if ((i + 1 >= width) || (j + 1 >= height) || (k == 0)) continue;
        mn = &m3d->nodes[i][j][k - 1];
        mnp = &m3d->pnodes[i][j][k - 1];
        mni = &m3d->nodes[i + 1][j][k - 1];
        mnj = &m3d->nodes[i][j + 1][k - 1];
        mnk = &m3d->nodes[i][j][k];
        break;
    }

    num++;

    /* compute cross products and area delta */
    MN_SUB(mni, mn, v_i);
    MN_SUB(mnj, mn, v_j);
    MN_SUB(mnk, mn, v_k);

    delta = invert * (mnp->orig_area - mnp->area) / n3;

    /* scale up the area coefficient if the area of the current node is
       close to 0 or already negative */
    if (!FZERO(mnp->orig_area))
      ratio = mnp->area / mnp->orig_area;
    else {
      ratio = 1;
      /*      fprintf(stdout, "orig area = 0 at (%d, %d, %d)!!!\n", i, j, k);*/
    }
    if (mnp->area < 0) DiagBreak();
    /*delta = mnp->area / mnp->orig_area ;delta = 100*pow(0.001, abs(10*x));*/
    exponent = parms->exp_k * ratio;
    if (exponent > MAX_EXP) exponent = MAX_EXP;
    delta_scale = l_nlarea / (1.0 + exp(exponent));
    if (delta_scale > 10000) DiagBreak();
    delta *= delta_scale; /* chain rule */
    if (fabs(delta) > 10000) DiagBreak();

    /* compute cross-products and add the appropriate
       (i.e. scaled by area difference) cross-products to the gradient */
    switch (n) {
      default:
      case 4:
      case 0: /* first do central node */
        V3_CROSS_PRODUCT(v_k, v_j, v_k_x_j);
        V3_CROSS_PRODUCT(v_j, v_i, v_j_x_i);
        V3_CROSS_PRODUCT(v_i, v_k, v_i_x_k);
        V3_ADD(v_k_x_j, v_j_x_i, v_tmp);
        V3_ADD(v_i_x_k, v_tmp, v_tmp);
        V3_SCALAR_MUL(v_tmp, delta, v_tmp);
        break;
      case 5: /*  i+1 */
      case 1: /*  i-1 */
        V3_CROSS_PRODUCT(v_k, v_j, v_k_x_j);
        V3_SCALAR_MUL(v_k_x_j, -delta, v_tmp);
        break;
      case 6: /* j+1 */
      case 2: /* j-1 */
        V3_CROSS_PRODUCT(v_i, v_k, v_i_x_k);
        V3_SCALAR_MUL(v_i_x_k, -delta, v_tmp);
        break;
      case 7: /* k+1 */
      case 3: /* k-1 */
        V3_CROSS_PRODUCT(v_j, v_i, v_j_x_i);
        V3_SCALAR_MUL(v_j_x_i, -delta, v_tmp);
        break;
    }
    V3_ADD(v_tmp, v_grad, v_grad);
  }
  *pdx = V3_X(v_grad);
  *pdy = V3_Y(v_grad);
  *pdz = V3_Z(v_grad);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define AVG_LEN 5
static int m3dAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MRI *mri_ref_blur, MP *parms, MORPH_3D *m3d)
{
  double dt, sse, rms, nvox, intensity_rms, distance_rms, area_rms, d, old_rms, last_avg, current_avg;
  int n, last_neg;

  if (mri_ref->nframes > 1)
    mri_ref->mean = MRImeanFrame(mri_ref, 1);
  else
    mri_ref->mean = 1.0;
  if (DZERO(mri_ref->mean)) {
    /*    ErrorPrintf(ERROR_BADPARM, "AlignPyramidLevel: mean std is 0!!") ;*/
    mri_ref->mean = 1.0;
  }

  dt = parms->dt;
  nvox = m3d->width * m3d->height * m3d->depth;
  m3dComputeMetricProperties(m3d);
  sse = m3dSSE(mri_in, mri_ref, parms, m3d, &intensity_rms, &distance_rms, &area_rms);
  old_rms = rms = sqrt(sse / nvox);
  old_rms += 0.5 * parms->tol * rms;
  log3DIntegration(parms, m3d, 0, dt, rms, intensity_rms, distance_rms, area_rms);

  last_neg = m3d->neg;
  if ((Gdiag & DIAG_WRITE) && (parms->write_iterations > 0) && !parms->start_t) {
    /*    char fname[STRLEN] ;*/
    write3DSnapshot(m3d->mri_in, m3d->mri_ref, parms, m3d, 0);

  }

  last_avg = old_rms;
  current_avg = rms;
  for (n = parms->start_t; n < parms->start_t + parms->niterations; n++) {

    m3dPositionBorderNodes(m3d);
    d = m3dIntegrationStep(mri_in, mri_ref, mri_ref_blur, parms, m3d, dt);
    m3dComputeMetricProperties(m3d);
    sse = m3dSSE(mri_in, mri_ref, parms, m3d, &intensity_rms, &distance_rms, &area_rms);
    rms = sqrt(sse / nvox);
    if ((Gdiag & DIAG_WRITE) && (parms->write_iterations > 0) && (((n + 1) % parms->write_iterations) == 0))
      write3DSnapshot(m3d->mri_in, m3d->mri_ref, parms, m3d, n + 1);
    log3DIntegration(parms, m3d, n + 1, d, rms, intensity_rms, distance_rms, area_rms);
    if ((((n + 1) - parms->start_t) % AVG_LEN) == 0) {
      current_avg /= AVG_LEN;
      if (Gdiag & DIAG_WRITE && parms->log_fp)
        fprintf(parms->log_fp,
                "last avg=%2.2f, current avg = %2.2f, ratio=%2.3f\n\n",
                last_avg,
                current_avg,
                (last_avg - current_avg) / current_avg);
      fprintf(stdout,
              "last avg=%2.2f, current avg = %2.2f, ratio=%2.3f\n\n",
              last_avg,
              current_avg,
              (last_avg - current_avg) / current_avg);
      if ((((last_avg - current_avg) / current_avg < parms->tol) || FZERO(parms->tol)) && (last_neg <= m3d->neg)) break;
      last_neg = m3d->neg;
      last_avg = current_avg;
      current_avg = rms;
    }
    else
      current_avg += rms;
  }

  if (n < parms->niterations + parms->start_t) n++; /* broke in while, n will be one too small (never incremented) */
  parms->start_t = n;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double m3dIntegrationStep(
    MRI *mri_in, MRI *mri_ref, MRI *mri_ref_blur, MORPH_PARMS *parms, MORPH_3D *m3d, double dt)
{
  double dx, dy, dz, xgrad, ygrad, zgrad, dot, nx, ny, nz, ndx, ndy, ndz, thick, pdx, pdy, pdz;
  int width, height, depth, x, y, z, xsi, ysi, zsi, *pxi, *pyi, *pzi;
  MORPH_NODE *mn;
  MNP *mnp;
  double xin, yin, zin, scale, xref, yref, zref;

  nx = parms->in_np.neck_x0;
  ny = parms->in_np.neck_y0;
  nz = parms->in_np.neck_z0;
  ndx = parms->in_np.neck_dx;
  ndy = parms->in_np.neck_dy;
  ndz = parms->in_np.neck_dz;
  m3dclearGradient(m3d);
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  thick = mri_in->thick;
  scale = m3d->node_spacing / thick;
  pxi = mri_in->xi;
  pyi = mri_in->yi;
  pzi = mri_in->zi;
  for (x = 0; x < m3d->width; x++) {
    xin = x * scale; /* convert to voxel space of mri_in */
    xsi = pxi[nint(xin)];
    for (y = 0; y < m3d->height; y++) {
      yin = y * scale;
      ysi = pyi[nint(yin)]; /* voxel coords of mri_in */
      for (z = 0; z < m3d->depth; z++) {
        if (x == debug_x && y == debug_y && z == debug_z) DiagBreak();
        xgrad = ygrad = zgrad = 0.0;
        zin = z * scale;
        zsi = pzi[nint(zin)];        /* voxel coords of mri_in */
        mn = &m3d->nodes[x][y][z];   /* find out where this voxel went */
        mnp = &m3d->pnodes[x][y][z]; /* find out where this voxel went */
        xref = mn->x / mri_ref->thick;
        yref = mn->y / mri_ref->thick;
        zref = mn->z / mri_ref->thick;
        m3dCorrelationTerm(
            mri_in, mri_ref, mri_ref_blur, parms->l_intensity, xsi, ysi, zsi, xref, yref, zref, &dx, &dy, &dz);
        if (x == 3 * width / 4 && y == 4 * height / 4 && z == 3 * depth / 4) DiagBreak();
        pdx = (x * thick - nx);
        pdy = (y * thick - ny);
        pdz = (z * thick - nz);
        dot = pdx * ndx + pdy * ndy + pdz * ndz;
        if (dot >= 0) DiagBreak();
        dot = 1 / (1 + exp(.25 * (dot - 10))); /* 'don't care' weighting */
        dx *= dot;
        dy *= dot;
        dz *= dot;
        if (!finitep(dz) || !finitep(dy) || !finitep(dx)) DiagBreak();
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;
        if (!finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad)) DiagBreak();
        m3dDistanceTerm(m3d, parms->l_dist, x, y, z, &dx, &dy, &dz);
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;
        if (!finitep(dz) || !finitep(dy) || !finitep(dx) || !finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad))
          DiagBreak();

        m3dCompressionTerm(m3d, parms->l_compression, x, y, z, &dx, &dy, &dz);
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;

        m3dNonlinearDistanceTerm(m3d, parms->l_nldist, x, y, z, &dx, &dy, &dz);
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;

        m3dAreaTerm(m3d, parms->l_area, x, y, z, &dx, &dy, &dz);
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;
        m3dNonlinearAreaTerm(m3d, parms->l_nlarea, x, y, z, &dx, &dy, &dz, parms);
        xgrad += dx;
        ygrad += dy;
        zgrad += dz;
        if (!finitep(dz) || !finitep(dy) || !finitep(dx) || !finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad))
          DiagBreak();

        mnp->dx = xgrad;
        mnp->dy = ygrad;
        mnp->dz = zgrad;
        if (!finitep(mnp->dx) || !finitep(mnp->dy) || !finitep(mnp->dz)) DiagBreak();
      }
    }
  }

#define MAX_GRAD (m3d->node_spacing / (8.0f))
  m3dcheck(m3d);
#if USE_ITERATIVE_AVERAGING
  m3daverageGradient(m3d, parms->navgs);
#else
  m3dblurGradient(m3d, parms->sigma);
#endif
  m3dcheck(m3d);
  dt = m3dscaleDeltaT(m3d, dt, MAX_GRAD, parms);
  m3dapplyGradient(m3d, dt);
  return (dt);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dInitDistances(MORPH_3D *m3d)
{
  int i, j, k, width, height, depth, n, ndists;
  MORPH_NODE *mn;
  MNP *mnp, *mnp_nbr;
  double dist, xd, yd, zd, avg_dist;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  for (ndists = avg_dist = 0.0, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        for (n = 0; n < NEIGHBORS; n++) {
          switch (n) {
            default:
            case 0: /* i-1 */
              if (i > 0)
                mnp_nbr = &m3d->pnodes[i - 1][j][k];
              else
                mnp_nbr = NULL;
              break;
            case 1: /* i+1 */
              if (i < width - 1)
                mnp_nbr = &m3d->pnodes[i + 1][j][k];
              else
                mnp_nbr = NULL;
              break;
            case 2: /* j-1 */
              if (j > 0)
                mnp_nbr = &m3d->pnodes[i][j - 1][k];
              else
                mnp_nbr = NULL;
              break;
            case 3: /* j+1 */
              if (j < height - 1)
                mnp_nbr = &m3d->pnodes[i][j + 1][k];
              else
                mnp_nbr = NULL;
              break;
            case 4: /* k-1 */
              if (k > 0)
                mnp_nbr = &m3d->pnodes[i][j][k - 1];
              else
                mnp_nbr = NULL;
              break;
            case 5: /* k+1 */
              if (k < depth - 1)
                mnp_nbr = &m3d->pnodes[i][j][k + 1];
              else
                mnp_nbr = NULL;
              break;
          }
          if (mnp_nbr) {
            xd = mnp->ox - mnp_nbr->ox;
            yd = mnp->oy - mnp_nbr->oy;
            zd = mnp->oz - mnp_nbr->oz;
            dist = sqrt(xd * xd + yd * yd + zd * zd);
            if (!finitep(dist)) DiagBreak();
            mnp->orig_dist[n] = dist;
            avg_dist += mnp->orig_dist[n];
            ndists++;
          }
        }
      }
    }
  }
  if (Gdiag & DIAG_SHOW) {
    avg_dist /= (double)ndists;
    fprintf(stdout, "average node spacing = %2.3f\n", avg_dist);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dInitAreas(MORPH_3D *m3d)
{
  int i, j, k, width, height, depth, num, num_zero, total_num;
  MNP *mnp, *mnpi, *mnpj, *mnpk;
  VECTOR *v_i, *v_j, *v_k;
  double avg_area, min_area, max_area, area;

  v_i = VectorAlloc(3, MATRIX_REAL);
  v_j = VectorAlloc(3, MATRIX_REAL);
  v_k = VectorAlloc(3, MATRIX_REAL);
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  min_area = 100000.0;
  max_area = -1.0;
  for (avg_area = 0.0, num_zero = total_num = i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mnp = &m3d->pnodes[i][j][k];

        num = 0;
        if ((i < width - 1) && (j < height - 1) && (k < depth - 1)) {
          mnpi = &m3d->pnodes[i + 1][j][k];
          mnpj = &m3d->pnodes[i][j + 1][k];
          mnpk = &m3d->pnodes[i][j][k + 1];
          MNP_SUB(mnpi, mnp, v_i);
          MNP_SUB(mnpj, mnp, v_j);
          MNP_SUB(mnpk, mnp, v_k);
          area = VectorTripleProduct(v_j, v_k, v_i);
          num++;
        }
        else
          area = 0.0;
        if ((i > 0) && (j > 0) && (k > 0)) /* left hand coordinate system */
        {
          double a;
          mnpi = &m3d->pnodes[i - 1][j][k];
          mnpj = &m3d->pnodes[i][j - 1][k];
          mnpk = &m3d->pnodes[i][j][k - 1];

          /* invert v_i so that coordinate system is still right-handed */
          MNP_SUB(mnp, mnpi, v_i);
          MNP_SUB(mnpj, mnp, v_j);
          MNP_SUB(mnpk, mnp, v_k);
          a = VectorTripleProduct(v_j, v_k, v_i);
          if (a < 0) DiagBreak();
          area += a;
          num++;
        }
        if (!num)
          mnp->orig_area = 0;
        else {
          total_num++;
          mnp->orig_area = area / (float)num;
          avg_area += mnp->orig_area;
          if (mnp->orig_area > max_area) max_area = mnp->orig_area;
          if (mnp->orig_area < min_area) min_area = mnp->orig_area;
        }
        if (FZERO(mnp->orig_area)) num_zero++;
        if (!finitep(mnp->orig_area)) DiagBreak();
      }
    }
  }

  if (Gdiag & DIAG_SHOW) {
    avg_area /= (double)total_num;
    fprintf(
        stdout, "average node area    = %2.3f [%2.3f --> %2.3f] (%d zero)\n", avg_area, min_area, max_area, num_zero);
  }

  VectorFree(&v_i);
  VectorFree(&v_j);
  VectorFree(&v_k);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int log3DIntegration(MORPH_PARMS *parms,
                            MORPH_3D *m3d,
                            int n,
                            double dt,
                            double rms,
                            double intensity_rms,
                            double distance_rms,
                            double area_rms)
{
  if (Gdiag & DIAG_SHOW)
#if USE_ITERATIVE_AVERAGING
    fprintf(stdout,
            "%3.3d: dt = %2.4f, "
            "rms = %2.3f (%2.2f, %2.2f, %2.2f), neg=%d, avgs=%d\n",
            n,
            dt,
            rms,
            intensity_rms,
            distance_rms,
            area_rms,
            m3d->neg,
            parms->navgs);
#else
    fprintf(stdout,
            "%3.3d: dt = %2.4f, "
            "rms = %2.3f (%2.2f, %2.2f, %2.2f), neg=%d, sigma=%2.1f\n",
            n,
            dt,
            rms,
            intensity_rms,
            distance_rms,
            area_rms,
            m3d->neg,
            parms->sigma);
#endif
  if (Gdiag & DIAG_WRITE) {
#if USE_ITERATIVE_AVERAGING
    fprintf(parms->log_fp,
            "%3.3d: dt = %2.4f, "
            "rms = %2.3f (%2.2f, %2.2f, %2.2f), neg=%d, avgs=%d\n",
            n,
            dt,
            rms,
            intensity_rms,
            distance_rms,
            area_rms,
            m3d->neg,
            parms->navgs);
#else
    fprintf(parms->log_fp,
            "%3.3d: dt = %2.4f, "
            "rms = %2.3f (%2.2f, %2.2f, %2.2f), neg=%d, sigma=%2.1f\n",
            n,
            dt,
            rms,
            intensity_rms,
            distance_rms,
            area_rms,
            m3d->neg,
            parms->sigma);
#endif
    fflush(parms->log_fp);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#include "image.h"
static int write3DSnapshot(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MORPH_3D *m3d, int n)
{
  MRI *mri_tmp;
  char fname[STRLEN];

  if (!(Gdiag & DIAG_WRITE)) return (NO_ERROR);

  if (!n) {
    sprintf(fname, "in_%s", parms->base_name);
    MRIwriteImageViews(mri_in, fname, IMAGE_SIZE);
  }

  mri_tmp = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
  sprintf(fname, "%s%3.3d", parms->base_name, n);
  MRIwriteImageViews(mri_tmp, fname, IMAGE_SIZE);
  MRIfree(&mri_tmp);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRIwriteImageViews(MRI *mri, const char *base_name, int target_size)
{
  int slice_direction = getSliceDirection(mri), x, y, z;

  switch (slice_direction) {
    default:
    case MRI_CORONAL:
      x = Gsx;
      y = Gsy;
      z = Gsz;
      break;
    case MRI_SAGITTAL:
      x = Gsz;
      y = Gsy;
      z = Gsx;
      break;
    case MRI_HORIZONTAL:
      x = Gsz;
      y = Gsy;
      z = Gsx;
      break;  // not sure this works???
  }

  if (Gdiag & DIAG_SHOW) fprintf(stdout, "writing image views to ???_%s.tif...\n", base_name);
  if (z >= 0)
    mriWriteImageView(mri, base_name, target_size, MRI_CORONAL, z);
  else if (Gvz >= 0)
    mriWriteImageView(mri, base_name, target_size, MRI_CORONAL, Gvz);
  else
    mriWriteImageView(mri, base_name, target_size, MRI_CORONAL, -1);

  if (x >= 0)
    mriWriteImageView(mri, base_name, target_size, MRI_SAGITTAL, x);
  else if (Gvx >= 0)
    mriWriteImageView(mri, base_name, target_size, MRI_SAGITTAL, Gvx);
  else
    mriWriteImageView(mri, base_name, target_size, MRI_SAGITTAL, -1);

  if (y >= 0)
    mriWriteImageView(mri, base_name, target_size, MRI_HORIZONTAL, y);
  else if (Gvy >= 0 && Gvy < mri->height)
    mriWriteImageView(mri, base_name, target_size, MRI_HORIZONTAL, Gvy);
  else
    mriWriteImageView(mri, base_name, target_size, MRI_HORIZONTAL, -1);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int mriWriteImageView(MRI *mri, const char *base_name, int target_size, int view, int slice)
{
  char fname[STRLEN];
  const char *prefix;
  IMAGE *I;
  float scale;
  //  int   slice_direction ;

  switch (view) {
    default:
    case MRI_CORONAL:
      prefix = "cor";
      break;
    case MRI_SAGITTAL:
      prefix = "sag";
      break;
    case MRI_HORIZONTAL:
      prefix = "hor";
      break;
  }

  I = MRItoImageView(mri, NULL, slice, view, 0);
  if (!I) ErrorReturn(Gerror, (Gerror, "MRItoImageView failed"));

  scale = (float)target_size / (float)I->rows;
  if (!FEQUAL(scale, 1.0f)) {
    IMAGE *Itmp;

    Itmp = ImageRescale(I, NULL, scale);
    ImageFree(&I);
    I = Itmp;
  }
  sprintf(fname, "%s_%s.tif", prefix, base_name);
  ImageWrite(I, fname);
  ImageFree(&I);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dComputeMetricProperties(MORPH_3D *m3d)
{
  double area;
  int i, j, k, width, height, depth, num;
  MORPH_NODE *mn, *mni, *mnj, *mnk;
  MNP *mnp;
  VECTOR *v_i, *v_j, *v_k;

  v_i = VectorAlloc(3, MATRIX_REAL);
  v_j = VectorAlloc(3, MATRIX_REAL);
  v_k = VectorAlloc(3, MATRIX_REAL);
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  m3d->neg = 0;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        num = 0;
        if ((i < width - 1) && (j < height - 1) && (k < depth - 1)) {
          num++;
          mni = &m3d->nodes[i + 1][j][k];
          mnj = &m3d->nodes[i][j + 1][k];
          mnk = &m3d->nodes[i][j][k + 1];
          MN_SUB(mni, mn, v_i);
          MN_SUB(mnj, mn, v_j);
          MN_SUB(mnk, mn, v_k);
          area = VectorTripleProduct(v_j, v_k, v_i);
        }
        else
          area = 0;
        if ((i > 0) && (j > 0) && (k > 0)) /* left-hand coordinate system */
        {
          num++;
          mni = &m3d->nodes[i - 1][j][k];
          mnj = &m3d->nodes[i][j - 1][k];
          mnk = &m3d->nodes[i][j][k - 1];

          /* invert v_i so that coordinate system is right-handed */
          MN_SUB(mn, mni, v_i);
          MN_SUB(mnj, mn, v_j);
          MN_SUB(mnk, mn, v_k);
          area += VectorTripleProduct(v_j, v_k, v_i);
        }
        if (num > 0)
          mnp->area = area / (float)num;
        else
          mnp->area = 0;
        if ((area <= 0) && !FZERO(mnp->orig_area)) m3d->neg++;
      }
    }
  }

  VectorFree(&v_i);
  VectorFree(&v_j);
  VectorFree(&v_k);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dclearGradient(MORPH_3D *m3d)
{
  int i, j, k, width, height, depth;
  MNP *mnp;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mnp = &m3d->pnodes[i][j][k];
        mnp->dx = mnp->dy = mnp->dz = 0.0f;
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
#define DI 42
#define DJ 90
#define DK 104

static int m3dapplyGradient(MORPH_3D *m3d, double dt)
{
  int i, j, k, width, height, depth;
  MORPH_NODE *mn;
  MNP *mnp;
  double dx, dy, dz;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        dx = mnp->dx;
        dy = mnp->dy;
        dz = mnp->dz;
        mn->x += dx * dt;
        mn->y += dy * dt;
        mn->z += dz * dt;
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
static int finitep(float f)
{
  return (1);
  if (!std::isfinite(f)) return (0);
  if (fabs(f) > 1e5) return (0);
  return (1);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double m3dscaleDeltaT(MORPH_3D *m3d, double dt, double max_len, MORPH_PARMS *parms)
{
  int i, j, k, width, height, depth, maxi, maxj, maxk, mini, minj, mink;
  MORPH_NODE *mn;
  MNP *mnp;
  double len, dx, dy, dz, max_delta, new_dt, max_orig, max_current, min_ratio;

  mini = minj = mink = maxi = maxj = maxk = -1;
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  max_current = max_orig = 1.0;
  min_ratio = 10000.0;
  for (max_delta = -1.0f, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        dx = mnp->dx;
        dy = mnp->dy;
        dz = mnp->dz;
        len = sqrt(dx * dx + dy * dy + dz * dz);
        if (!FZERO(mnp->orig_area) && mnp->area / mnp->orig_area < min_ratio) {
          min_ratio = mnp->area / mnp->orig_area;
          mini = i;
          minj = j;
          mink = k;
        }
        if (len > max_delta) {
          maxi = i;
          maxj = j;
          maxk = k;
          max_delta = len;
          max_orig = mnp->orig_area;
          max_current = mnp->area;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout,
            "max delta %2.3f at (%d,%d,%d), ratio=%2.1f/%2.1f=%2.2f\n",
            max_delta,
            maxi,
            maxj,
            maxk,
            max_current,
            max_orig,
            FZERO(max_orig) ? 0.0 : max_current / max_orig);
  if (Gdiag & DIAG_WRITE && parms->log_fp)
    fprintf(parms->log_fp,
            "max delta %2.3f at (%d,%d,%d), ratio=%2.1f/%2.1f=%2.2f\n",
            max_delta,
            maxi,
            maxj,
            maxk,
            max_current,
            max_orig,
            FZERO(max_orig) ? 0.0 : max_current / max_orig);

  if ((min_ratio < .25 || max_delta > 100) && (mini != maxi || minj != maxj || mink != maxk)) {
    debug_x = mini;
    debug_y = minj;
    debug_z = mink;
    DiagBreak();
  }
  if (max_delta > 100) {
    debug_x = maxi;
    debug_y = maxj;
    debug_z = maxk;
    mn = &m3d->nodes[maxi][maxj][maxk];
    mnp = &m3d->pnodes[maxi][maxj][maxk];
  }

  if (max_delta * dt > max_len) {
    new_dt = max_len / max_delta;
    if (new_dt < dt) dt = new_dt;
  }
  return (dt);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define UNUSED_SPACE_SIZE 256

int MRI3Dwrite(MORPH_3D *m3d, char *fname)
{
  FILE *fp;
  int i, j, k;
  MORPH_NODE *mn;
  char buf[UNUSED_SPACE_SIZE + 1];

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "mris3Dwrite: could not open file %s", fname));

  if (fwriteInt(M3D_MAGIC, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
  if (fwriteInt(M3D_VERSION, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
  if (fwriteInt(m3d->width, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
  if (fwriteInt(m3d->height, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
  if (fwriteInt(m3d->depth, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
  if (fwriteFloat(m3d->node_spacing, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));

  /* so stuff can be added to the header in the future */
  memset(buf, 0, UNUSED_SPACE_SIZE * sizeof(char));
  fwrite(buf, sizeof(char), UNUSED_SPACE_SIZE, fp);
  for (i = 0; i < m3d->width; i++) {
    for (j = 0; j < m3d->height; j++) {
      for (k = 0; k < m3d->depth; k++) {
        mn = &m3d->nodes[i][j][k];
        if (fwriteFloat(mn->x, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
        if (fwriteFloat(mn->y, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
        if (fwriteFloat(mn->z, fp) != 1) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "m3dwrite: fwrite failed"));
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Read a 3d morph from a file, and initialize it,
           not allocating the properties field to conserve memory.
------------------------------------------------------*/
MORPH_3D *MRI3DreadSmall(char *fname)
{
  MORPH_3D *m3d;
  FILE *fp;
  int i, j, k, width, height, depth, magic, version;
  float node_spacing;
  MORPH_NODE *mn;
  char buf[UNUSED_SPACE_SIZE + 1];

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "MRI3DreadSmall: could not open file %s", fname));

  magic = freadInt(fp);
  if ((unsigned)magic != M3D_MAGIC) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "file %s not an old 3d morph file.\nTry a new 3d morph read routine.\n", fname));
  }
  version = freadInt(fp);
  width = freadInt(fp);
  height = freadInt(fp);
  depth = freadInt(fp);
  node_spacing = freadFloat(fp);

  /* so stuff can be added to the header in the future */
  if (fread(buf, sizeof(char), UNUSED_SPACE_SIZE, fp) != UNUSED_SPACE_SIZE && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "Could not read unused space");
  }

  m3d = m3dalloc(width, height, depth, node_spacing);
  for (i = 0; i < m3d->width; i++) {
    for (j = 0; j < m3d->height; j++) {
      for (k = 0; k < m3d->depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mn->x = freadFloat(fp);
        mn->y = freadFloat(fp);
        mn->z = freadFloat(fp);
      }
    }
  }

  fclose(fp);
  return (m3d);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Read a 3d morph from a file, and initialize it.
------------------------------------------------------*/
MORPH_3D *MRI3Dread(char *fname)
{
  MORPH_3D *m3d;

  m3d = MRI3DreadSmall(fname);
  if (!m3d) return (NULL);

  m3dAllocProperties(m3d);
  m3dComputeMetricProperties(m3d);
  return (m3d);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Read a 3d morph from a file, and initialize it.
------------------------------------------------------*/
static MORPH_3D *m3dscaleUp2(MORPH_3D *m3d_in, MORPH_3D *m3d_out, MORPH_PARMS *parms)
{
  int width, height, depth, i, j, k, alloced = 0;
  MORPH_NODE *mn;
  MNP *mnp;
  float x, y, z, node_spacing;
  VECTOR *v_X, *v_Y;

  v_X = VectorAlloc(4, MATRIX_REAL);
  v_Y = VectorAlloc(4, MATRIX_REAL);

  width = m3d_in->width * 2 - 1;
  height = m3d_in->height * 2 - 1;
  depth = m3d_in->depth * 2 - 1;
  node_spacing = m3d_in->node_spacing / 2.0;
  if (!m3d_out) {
    alloced = 1;
    m3d_out = m3dalloc(width, height, depth, node_spacing);
  }
  m3d_out->lta = m3d_in->lta;

  for (i = 0; i < width; i++) {
    x = V3_X(v_X) = (float)i * node_spacing;
    for (j = 0; j < height; j++) {
      y = V3_Y(v_X) = (float)j * node_spacing;
      for (k = 0; k < depth; k++) {
        if (k == depth - 1) DiagBreak();
        z = V3_Z(v_X) = (float)k * node_spacing;
        mn = &m3d_out->nodes[i][j][k];
        MRIsample3Dmorph(m3d_in, x, y, z, &mn->x, &mn->y, &mn->z);
      }
    }
  }
  MRI3DmorphFree(&m3d_in);

  m3dAllocProperties(m3d_out);
  for (i = 0; i < width; i++) {
    x = V3_X(v_X) = (float)i * node_spacing;
    for (j = 0; j < height; j++) {
      y = V3_Y(v_X) = (float)j * node_spacing;
      for (k = 0; k < depth; k++) {
        z = V3_Z(v_X) = (float)k * node_spacing;
        mn = &m3d_out->nodes[i][j][k];
        mnp = &m3d_out->pnodes[i][j][k];
        if (alloced) {
#if !USE_ORIGINAL_PROPTERIES
          LTAtransformPoint(parms->lta, v_X, v_Y);
#else
          MatrixCopy(v_X, v_Y);
#endif
          mnp->ox = V3_X(v_Y);
          mnp->oy = V3_Y(v_Y);
          mnp->oz = V3_Z(v_Y);
        }

        if (!finitep(mn->x) || !finitep(mn->y) || !finitep(mn->z) || !finitep(mnp->ox) || !finitep(mnp->oy) ||
            !finitep(mnp->oz))
          DiagBreak();
      }
    }
  }

  m3dComputeMetricProperties(m3d_out);
  if (alloced) {
    m3dInitDistances(m3d_out);
    m3dInitAreas(m3d_out);
  }
  return (m3d_out);
}
int MRIeraseBorders(MRI *mri, int width)
{
  int i;

  for (i = 0; i < width; i++) {
    MRIerasePlane(mri, mri->width / 2, mri->height / 2, i, 0, 0, 1, 0);
    MRIerasePlane(mri, mri->width / 2, i, mri->depth / 2, 0, 1, 0, 0);
    MRIerasePlane(mri, i, mri->height / 2, mri->depth / 2, 1, 0, 0, 0);
    MRIerasePlane(mri, mri->width / 2, mri->height / 2, mri->depth - 1 - i, 0, 0, 1, 0);
    MRIerasePlane(mri, mri->width / 2, mri->height - 1 - i, mri->depth / 2, 0, 1, 0, 0);
    MRIerasePlane(mri, mri->width - 1 - i, mri->height / 2, mri->depth / 2, 1, 0, 0, 0);
  }
  return (NO_ERROR);
}

static int mriNormalizeStds(MRI *mri)
{
  int x, y, z, width, height, depth, std, val;
  BUFTYPE *pstd, *pval;
  double mean_std;

  if (mri->nframes < 2) return (NO_ERROR);

  if (mri->nframes > 1)
    mean_std = MRImeanFrame(mri, 1);
  else
    mean_std = 1.0;
  if (FZERO(mean_std)) {
    mri->mean = 0.0;
    return (NO_ERROR);
  }
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pval = &MRIvox(mri, 0, y, z);
      pstd = &MRIseq_vox(mri, 0, y, z, 1);
      for (x = 0; x < width; x++) {
        val = *pval++;
        if (val < BACKGROUND_VAL) {
          std = *pstd;
          if (std < mean_std) std = mean_std;
          *pstd++ = std;
        }
        else
          pstd++;
      }
    }
  }
  if (mri->nframes > 1)
    mri->mean = MRImeanFrame(mri, 1);
  else
    mri->mean = 1.0;

  return (NO_ERROR);
}
#if USE_ITERATIVE_AVERAGING
static int m3daverageGradient(MORPH_3D *m3d, int niter)
{
  int i, j, k, width, height, depth, n, a, num;
  MORPH_NODE *mn;
  MNP *mnp, *mnp_nbr;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  for (a = 0; a < niter; a++) {
    for (i = 0; i < width; i++) {
      for (j = 0; j < height; j++) {
        for (k = 0; k < depth; k++) {
          mn = &m3d->nodes[i][j][k];
          mnp = &m3d->pnodes[i][j][k];
          for (num = n = 0; n < NEIGHBORS; n++) {
            switch (n) {
              default:
              case 0: /* i-1 */
                if (i > 0)
                  mnp_nbr = &m3d->pnodes[i - 1][j][k];
                else
                  mnp_nbr = NULL;
                break;
              case 1: /* i+1 */
                if (i < width - 1)
                  mnp_nbr = &m3d->pnodes[i + 1][j][k];
                else
                  mnp_nbr = NULL;
                break;
              case 2: /* j-1 */
                if (j > 0)
                  mnp_nbr = &m3d->pnodes[i][j - 1][k];
                else
                  mnp_nbr = NULL;
                break;
              case 3: /* j+1 */
                if (j < height - 1)
                  mnp_nbr = &m3d->pnodes[i][j + 1][k];
                else
                  mnp_nbr = NULL;
                break;
              case 4: /* k-1 */
                if (k > 0)
                  mnp_nbr = &m3d->pnodes[i][j][k - 1];
                else
                  mnp_nbr = NULL;
                break;
              case 5: /* k+1 */
                if (k < depth - 1)
                  mnp_nbr = &m3d->pnodes[i][j][k + 1];
                else
                  mnp_nbr = NULL;
                break;
            }
            if (mnp_nbr) {
              num++;
              mnp->tdx += mnp_nbr->dx;
              mnp->tdy += mnp_nbr->dy;
              mnp->tdz += mnp_nbr->dz;
            }
          }
          if (num >= 0) {
            mnp->tdx /= (float)num;
            mnp->tdy /= (float)num;
            mnp->tdz /= (float)num;
          }
        }
      }
    }
    for (i = 0; i < width; i++) {
      for (j = 0; j < height; j++) {
        for (k = 0; k < depth; k++) {
          mnp = &m3d->pnodes[i][j][k];
          mnp->dx = mnp->tdx;
          mnp->dy = mnp->tdy;
          mnp->dz = mnp->tdz;
          mnp->tdx = mnp->tdy = mnp->tdz = 0.0;
        }
      }
    }
  }
  return (NO_ERROR);
}
#else
static int m3dblurGradient(MORPH_3D *m3d, float sigma)
{
  int i, j, k, i1, j1, k1, whalf, width, height, depth;
  MORPH_NODE *mn;
  MNP *mnp, *mnp_nbr;
  double dx_total, dy_total, dz_total, ktotal, two_sigma_sq, dist, kernel;

  if (FZERO(sigma)) return (NO_ERROR);

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  whalf = nint(1.5 * sigma) /*nint(2*sigma)*/;
  two_sigma_sq = 2.0 * sigma * sigma;

  /* convolve in width */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (i1 = i - whalf; i1 <= i + whalf; i1++) {
          if (i1 < 0 || i1 >= width) continue;
          mnp_nbr = &m3d->pnodes[i1][j][k];
          dist = (i1 - i);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
          dy_total += mnp_nbr->dy * kernel;
          dz_total += mnp_nbr->dz * kernel;
        }
        mnp->tdx = dx_total / ktotal;
        mnp->tdy = dy_total / ktotal;
        mnp->tdz = dz_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);

  /* convolve in height */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (j1 = j - whalf; j1 <= j + whalf; j1++) {
          if (j1 < 0 || j1 >= height) continue;
          mnp_nbr = &m3d->pnodes[i][j1][k];
          dist = (j1 - j);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
          dy_total += mnp_nbr->dy * kernel;
          dz_total += mnp_nbr->dz * kernel;
        }
        mnp->tdx = dx_total / ktotal;
        mnp->tdy = dy_total / ktotal;
        mnp->tdz = dz_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);

  /* convolve in depth */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (k1 = k - whalf; k1 <= k + whalf; k1++) {
          if (k1 < 0 || k1 >= depth) continue;
          mnp_nbr = &m3d->pnodes[i][j][k1];
          dist = (k1 - k);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
          dy_total += mnp_nbr->dy * kernel;
          dz_total += mnp_nbr->dz * kernel;
        }
        mnp->tdx = dx_total / ktotal;
        mnp->tdy = dy_total / ktotal;
        mnp->tdz = dz_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);
  return (NO_ERROR);
}
#endif

static int m3dCopyTempToGradient(MORPH_3D *m3d)
{
  int i, j, k, width, height, depth;
  MNP *mnp;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mnp = &m3d->pnodes[i][j][k];
        mnp->dx = mnp->tdx;
        mnp->dy = mnp->tdy;
        mnp->dz = mnp->tdz;
      }
    }
  }
  return (NO_ERROR);
}
static int m3dcheck(MORPH_3D *m3d)
{
  MNP *mnp;
  int i, j, k;

  for (i = 0; i < m3d->width; i++) {
    for (j = 0; j < m3d->height; j++) {
      for (k = 0; k < m3d->depth; k++) {
        mnp = &m3d->pnodes[i][j][k];
        if (FZERO(mnp->orig_area)) DiagBreak();
      }
    }
  }
  return (NO_ERROR);
}

static float mriMaxRadius(MRI *mri, int low_val, int hi_val, float *px0, float *py0, float *pz0)
{
  int x, y, z, width, depth, height, nvox;
  double xw, yw, zw, max_radius, radius, xd, yd, zd;
  float x0, y0, z0;
  BUFTYPE val;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  x0 = y0 = z0 = 0.0f;
  nvox = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        val = MRIvox(mri, x, y, z);
        if (val >= low_val && val <= hi_val) {
          nvox++;
          MRIvoxelToWorld(mri, (double)x, (double)y, (double)z, &xw, &yw, &zw);
          x0 += xw;
          y0 += yw;
          z0 += zw;
        }
      }
    }
  }

  x0 /= (float)nvox;
  y0 /= (float)nvox;
  z0 /= (float)nvox;
  max_radius = 0.0f;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        val = MRIvox(mri, x, y, z);
        if (val >= low_val && val <= hi_val) {
          MRIvoxelToWorld(mri, (double)x, (double)y, (double)z, &xw, &yw, &zw);
          xd = xw - x0;
          yd = yw - y0;
          zd = zw - z0;
          radius = sqrt(xd * xd + yd * yd + zd * zd);
          if (radius > max_radius) max_radius = radius;
        }
      }
    }
  }
  *px0 = x0;
  *py0 = y0;
  *pz0 = z0;
  return (max_radius);
}

MRI_SURFACE *MRISshrinkWrapSkull(MRI *mri, MORPH_PARMS *parms)
{
  MRI_SURFACE *mris;
  MRI *mri_smooth, *mri_kernel;
  float radius, max_radius, scale, x0, y0, z0, sx, sy, sz, thick;
  INTEGRATION_PARMS lparms;
  int diag;
  static int ncalls = 0;
  MRI_REGION bbox;

  thick = mri->thick;
  mri_kernel = MRIgaussian1d(0.5, 13);
  mri_smooth = MRIconvolveGaussian(mri, NULL, mri_kernel);
  diag = Gdiag;
  lparms.momentum = .75;
  lparms.dt = .75;
  lparms.l_spring_norm = 1.0f;
  lparms.write_iterations = 0 * parms->write_iterations;
  if (ncalls++ == 0)
    strcpy(lparms.base_name, "in_skull.geo");
  else
    strcpy(lparms.base_name, "ref_skull.geo");

  mris = ic2562_make_surface(0, 0);
  MRISsetNeighborhoodSizeAndDist(mris, 2);
  max_radius = mriMaxRadius(mri, 120, 255, &x0, &y0, &z0);
  radius = MRISaverageRadius(mris);
  scale = (max_radius + 1.0f) / radius;
  fprintf(stdout, "origin = (%2.1f, %2.1f, %2.1f), max radius = %2.1f.\n", x0, y0, z0, max_radius);
  MRISsetVals(mris, 150);
  /*  MRISscaleBrain(mris, mris, scale) ;*/
  MRIboundingBox(mri, 80, &bbox);
  fprintf(stdout, "bounding box size (%d, %d, %d)\n", bbox.dx, bbox.dy, bbox.dz);

  /* convert from MRI to surface coordinate system (y<-->z) */
  sx = bbox.dx * thick * .7 / radius;  /* ear to ear */
  sy = bbox.dy * thick * .7 / radius;  /* inferior->superior */
  sz = bbox.dz * thick * .65 / radius; /* anterior->posterior */
  fprintf(stdout, "scaling by (%2.1f, %2.1f, %2.1f)\n", sx, sy, sz);
  MRISanisotropicScale(mris, sx, sz, sy);
  MRISmoveOrigin(mris, x0, y0, z0);
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "start.geo");
  strcpy(mris->fname, "skull.geo");
  lparms.flags |= IPFLAG_NO_SELF_INT_TEST;
  /*  lparms.write_iterations = 5 ;*/
  if (!lparms.write_iterations) Gdiag = 0;
  MRISsetMarks(mris, 1);
  lparms.l_intensity = 0.025;
  lparms.niterations = 300;
  MRISpositionSurface(mris, mri, mri_smooth, &lparms); /* outer skin */
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "outer_skin1.geo");
  lparms.l_intensity = 0.015;
  lparms.niterations = 200;
  MRISpositionSurface(mris, mri, mri_smooth, &lparms); /* outer skin */
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "outer_skin2.geo");
  lparms.l_intensity = 0.01;
  lparms.niterations = 200;
  MRISpositionSurface(mris, mri, mri_smooth, &lparms); /* outer skin */
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "outer_skin3.geo");
  MRISexpandSurface(mris, -5, NULL, 0, 1);
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "outer_skin_expanded.geo");
  MRISsetVals(mris, 0);
  lparms.niterations = 200;
  lparms.l_intensity = 0.1;
  lparms.l_spring_norm = 1.0;
  MRISpositionSurface(mris, mri, mri_smooth, &lparms);

  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "inner_skin_expanded.geo");
  /* smooth final surface */
  lparms.l_intensity = 0.005;
  lparms.l_spring_norm = 1.0;
  lparms.niterations = 50;
  MRISpositionSurface(mris, mri, mri_smooth, &lparms);
  if (ncalls == 1 && (Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) MRISwrite(mris, "inner_skin_smoothed.geo");
  Gdiag = diag;

  MRIfree(&mri_smooth);
  MRIfree(&mri_kernel);
  return (mris);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the directional derivative of the distance term
          w.r.t the i,j,kth node.
------------------------------------------------------*/
static int m3dCompressionTerm(
    MORPH_3D *m3d, double l_compression, int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  int n, num;
  double dx, dy, dz, dist, delta, xgrad, ygrad, zgrad, expansion;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;

  if (FZERO(l_compression)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }

  l_compression /= (float)NEIGHBORS;
  mn = &m3d->nodes[i][j][k];
  mnp = &m3d->pnodes[i][j][k];
  for (xgrad = ygrad = zgrad = 0.0, num = n = 0; n < NEIGHBORS; n++) {
    switch (n) {
      default:
      case 0: /* i-1 */
        if (i > 0)
          mn_nbr = &m3d->nodes[i - 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 1: /* i+1 */
        if (i < m3d->width - 1)
          mn_nbr = &m3d->nodes[i + 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 2: /* j-1 */
        if (j > 0)
          mn_nbr = &m3d->nodes[i][j - 1][k];
        else
          mn_nbr = NULL;
        break;
      case 3: /* j+1 */
        if (j < m3d->height - 1)
          mn_nbr = &m3d->nodes[i][j + 1][k];
        else
          mn_nbr = NULL;
        break;
      case 4: /* k-1 */
        if (k > 0)
          mn_nbr = &m3d->nodes[i][j][k - 1];
        else
          mn_nbr = NULL;
        break;
      case 5: /* k+1 */
        if (k < m3d->depth - 1)
          mn_nbr = &m3d->nodes[i][j][k + 1];
        else
          mn_nbr = NULL;
        break;
    }
    if (mn_nbr) {
      num++;
      expansion = m3dNodeAverageExpansion(m3d, i, j, k);
      dx = mn_nbr->x - mn->x;
      dy = mn_nbr->y - mn->y;
      dz = mn_nbr->z - mn->z;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      delta = dist - expansion * mnp->orig_dist[n];
      if (!FZERO(dist)) /* make it a unit vector */
      {
        dx /= dist;
        dx /= dist;
        dz /= dist;
      }
#if SCALE_INVARIANT
      delta /= m3d->node_spacing;
#endif
      xgrad += delta * dx;
      ygrad += delta * dy;
      zgrad += delta * dz;
      if (!finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad) || (fabs(dx) > 1e5) || (fabs(dy) > 1e5) ||
          (fabs(dz) > 1e5))
        DiagBreak();
    }
  }

  l_compression /= (float)num;
  *pdx = l_compression * xgrad;
  *pdy = l_compression * ygrad;
  *pdz = l_compression * zgrad;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Compute the directional derivative of the distance term
          w.r.t the i,j,kth node.
------------------------------------------------------*/
static int m3dNonlinearDistanceTerm(
    MORPH_3D *m3d, double l_nldist, int i, int j, int k, double *pdx, double *pdy, double *pdz)
{
  int n, num;
  double dx, dy, dz, dist, delta, xgrad, ygrad, zgrad, expansion;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;

  if (FZERO(l_nldist)) {
    *pdx = *pdy = *pdz = 0.0;
    return (NO_ERROR);
  }

  l_nldist /= (float)NEIGHBORS;
  mn = &m3d->nodes[i][j][k];
  mnp = &m3d->pnodes[i][j][k];
  for (xgrad = ygrad = zgrad = 0.0, num = n = 0; n < NEIGHBORS; n++) {
    switch (n) {
      default:
      case 0: /* i-1 */
        if (i > 0)
          mn_nbr = &m3d->nodes[i - 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 1: /* i+1 */
        if (i < m3d->width - 1)
          mn_nbr = &m3d->nodes[i + 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 2: /* j-1 */
        if (j > 0)
          mn_nbr = &m3d->nodes[i][j - 1][k];
        else
          mn_nbr = NULL;
        break;
      case 3: /* j+1 */
        if (j < m3d->height - 1)
          mn_nbr = &m3d->nodes[i][j + 1][k];
        else
          mn_nbr = NULL;
        break;
      case 4: /* k-1 */
        if (k > 0)
          mn_nbr = &m3d->nodes[i][j][k - 1];
        else
          mn_nbr = NULL;
        break;
      case 5: /* k+1 */
        if (k < m3d->depth - 1)
          mn_nbr = &m3d->nodes[i][j][k + 1];
        else
          mn_nbr = NULL;
        break;
    }
    if (mn_nbr) {
      num++;
      expansion = m3dNodeAverageExpansion(m3d, i, j, k);
      dx = mn_nbr->x - mn->x;
      dy = mn_nbr->y - mn->y;
      dz = mn_nbr->z - mn->z;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      delta = dist - expansion * mnp->orig_dist[n];
      if (!FZERO(dist)) /* make it a unit vector */
      {
        dx /= dist;
        dx /= dist;
        dz /= dist;
      }
#if SCALE_INVARIANT
      delta /= m3d->node_spacing;
#endif
      xgrad += delta * dx;
      ygrad += delta * dy;
      zgrad += delta * dz;
      if (!finitep(xgrad) || !finitep(ygrad) || !finitep(zgrad) || (fabs(dx) > 1e5) || (fabs(dy) > 1e5) ||
          (fabs(dz) > 1e5))
        DiagBreak();
    }
  }

  l_nldist /= (float)num;
  *pdx = l_nldist * xgrad;
  *pdy = l_nldist * ygrad;
  *pdz = l_nldist * zgrad;
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static float m3dNodeAverageExpansion(MORPH_3D *m3d, int i, int j, int k)
{
  int n, num;
  double dx, dy, dz, dist, xgrad, ygrad, zgrad, avg_expansion;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;

  mn = &m3d->nodes[i][j][k];
  mnp = &m3d->pnodes[i][j][k];
  avg_expansion = 0;
  for (xgrad = ygrad = zgrad = 0.0, num = n = 0; n < NEIGHBORS; n++) {
    switch (n) {
      default:
      case 0: /* i-1 */
        if (i > 0)
          mn_nbr = &m3d->nodes[i - 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 1: /* i+1 */
        if (i < m3d->width - 1)
          mn_nbr = &m3d->nodes[i + 1][j][k];
        else
          mn_nbr = NULL;
        break;
      case 2: /* j-1 */
        if (j > 0)
          mn_nbr = &m3d->nodes[i][j - 1][k];
        else
          mn_nbr = NULL;
        break;
      case 3: /* j+1 */
        if (j < m3d->height - 1)
          mn_nbr = &m3d->nodes[i][j + 1][k];
        else
          mn_nbr = NULL;
        break;
      case 4: /* k-1 */
        if (k > 0)
          mn_nbr = &m3d->nodes[i][j][k - 1];
        else
          mn_nbr = NULL;
        break;
      case 5: /* k+1 */
        if (k < m3d->depth - 1)
          mn_nbr = &m3d->nodes[i][j][k + 1];
        else
          mn_nbr = NULL;
        break;
    }
    if (mn_nbr) {
      dx = mn_nbr->x - mn->x;
      dy = mn_nbr->y - mn->y;
      dz = mn_nbr->z - mn->z;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      num++;
      if (!FZERO(mnp->orig_dist[n])) avg_expansion += dist / mnp->orig_dist[n];
    }
  }
  avg_expansion /= (float)num;
  return (avg_expansion);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static double m3dCompressionSSE(MRI *mri_in, MRI *mri_ref, MORPH_3D *m3d)
{
  int i, j, k, width, height, depth, n;
  MORPH_NODE *mn, *mn_nbr;
  MNP *mnp;
  double dist, xd, yd, zd, delta, sse, node_spacing, expansion;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  node_spacing = m3d->node_spacing;
  for (sse = 0.0, i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        expansion = m3dNodeAverageExpansion(m3d, i, j, k);
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        for (n = 0; n < NEIGHBORS; n++) {
          switch (n) {
            default:
            case 0: /* i-1 */
              if (i > 0)
                mn_nbr = &m3d->nodes[i - 1][j][k];
              else
                mn_nbr = NULL;
              break;
            case 1: /* i+1 */
              if (i < width - 1)
                mn_nbr = &m3d->nodes[i + 1][j][k];
              else
                mn_nbr = NULL;
              break;
            case 2: /* j-1 */
              if (j > 0)
                mn_nbr = &m3d->nodes[i][j - 1][k];
              else
                mn_nbr = NULL;
              break;
            case 3: /* j+1 */
              if (j < height - 1)
                mn_nbr = &m3d->nodes[i][j + 1][k];
              else
                mn_nbr = NULL;
              break;
            case 4: /* k-1 */
              if (k > 0)
                mn_nbr = &m3d->nodes[i][j][k - 1];
              else
                mn_nbr = NULL;
              break;
            case 5: /* k+1 */
              if (k < depth - 1)
                mn_nbr = &m3d->nodes[i][j][k + 1];
              else
                mn_nbr = NULL;
              break;
          }
          if (mn_nbr) {
            xd = mn->x - mn_nbr->x;
            yd = mn->y - mn_nbr->y;
            zd = mn->z - mn_nbr->z;
            dist = sqrt(xd * xd + yd * yd + zd * zd);
            delta = (dist - expansion * mnp->orig_dist[n]);
#if SCALE_INVARIANT
            delta /= node_spacing;
#endif
            sse += delta * delta * mri_in->thick;
          }
        }
      }
    }
  }

  return (sse);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIeraseNeck(MRI *mri, NECK_PARMS *np)
{
  float dx, dy, dz, y, thick, x0, y0, z0;

  thick = mri->thick;
  x0 = np->neck_x0 / thick;
  y0 = np->neck_y0 / thick;
  z0 = np->neck_z0 / thick;
  dx = np->neck_dx;
  dy = np->neck_dy;
  dz = np->neck_dz;
  for (y = y0; y < mri->height + 25; y++) MRIerasePlane(mri, x0, y, z0, np->neck_dx, np->neck_dy, np->neck_dz, 0);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dPositionBorderNodes(MORPH_3D *m3d)
{
  return (NO_ERROR);
}

static int m3dblurDx(MORPH_3D *m3d, float sigma)
{
  int i, j, k, i1, j1, k1, whalf, width, height, depth;
  MORPH_NODE *mn;
  MNP *mnp, *mnp_nbr;
  double dx_total, dy_total, dz_total, ktotal, two_sigma_sq, dist, kernel;

  if (FZERO(sigma)) return (NO_ERROR);

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  whalf = nint(1.5 * sigma) /*nint(2*sigma)*/;
  two_sigma_sq = 2.0 * sigma * sigma;

  /* convolve in width */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (i1 = i - whalf; i1 <= i + whalf; i1++) {
          if (i1 < 0 || i1 >= width) continue;
          mnp_nbr = &m3d->pnodes[i1][j][k];
          dist = (i1 - i);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
        }
        mnp->tdx = dx_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);

  /* convolve in height */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (j1 = j - whalf; j1 <= j + whalf; j1++) {
          if (j1 < 0 || j1 >= height) continue;
          mnp_nbr = &m3d->pnodes[i][j1][k];
          dist = (j1 - j);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
        }
        mnp->tdx = dx_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);

  /* convolve in depth */
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        dx_total = dy_total = dz_total = ktotal = 0.0;
        for (k1 = k - whalf; k1 <= k + whalf; k1++) {
          if (k1 < 0 || k1 >= depth) continue;
          mnp_nbr = &m3d->pnodes[i][j][k1];
          dist = (k1 - k);
          dist *= dist;
          kernel = exp(-dist / two_sigma_sq);
          ktotal += kernel;
          dx_total += mnp_nbr->dx * kernel;
        }
        mnp->tdx = dx_total / ktotal;
      }
    }
  }
  m3dCopyTempToGradient(m3d);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int m3dTranslate(MORPH_3D *m3d, float dx, float dy, float dz)
{
  MORPH_NODE *mn;
  int i, j, k, width, height, depth;

  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        if (i == width / 2 && j == height / 2 && k == depth / 2) {
          DiagBreak();
          i = width / 2;
          j = height / 2;
          k = depth / 2;
        }

        mn = &m3d->nodes[i][j][k];
        mn->x += dx;
        mn->y += dy;
        mn->z += dz;
      }
    }
  }
  return (NO_ERROR);
}
/*
   note that mri_in and mri_ref will be copied so that the neck can
   be erased from both without changing the caller's volumes.
   */
static int m3RecomputeTranslation(MORPH_3D *m3d, MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  float in_means[4], ref_means[4];
  double dx, dy, dz;
  MRI *mri_in_xformed;

  mri_in = MRIcopy(mri_in, NULL);
  mri_ref = MRIcopy(mri_ref, NULL);
  MRIeraseNeck(mri_ref, &parms->ref_np);
  MRIeraseNeck(mri_in, &parms->ref_np);

  if (Gdiag & DIAG_SHOW) fprintf(stdout, "recomputing translation values based on skull morph...\n");
  mri_in_xformed = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
  if (Gdiag & DIAG_WRITE) MRIwriteImageViews(mri_in_xformed, "postxform", IMAGE_SIZE);
  MRIfindCenterOfBrain(mri_in, in_means, in_means + 1, in_means + 2);
  MRIfindCenterOfBrain(mri_in_xformed, ref_means, ref_means + 1, ref_means + 2);
  dx = (double)(ref_means[0] - in_means[0]) * mri_in->thick;
  dy = (double)(ref_means[1] - in_means[1]) * mri_in->thick;
  dz = (double)(ref_means[2] - in_means[2]) * mri_in->thick;
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "translating morph by (%2.1f,%2.1f,%2.1f)\n", dx, dy, dz);
  m3dTranslate(m3d, dx, dy, dz);
  MRIfree(&mri_in_xformed);
  if (Gdiag & DIAG_WRITE) {
    MRI *mri_tmp;
    MRIwriteImageViews(mri_in, "input", IMAGE_SIZE);
    mri_tmp = MRIapplyInverse3DMorph(mri_ref, m3d, NULL);
    MRIwriteImageViews(mri_tmp, "posttrans", IMAGE_SIZE);
    MRIfree(&mri_tmp);
  }
  MRIfree(&mri_in);
  MRIfree(&mri_ref);
  return (NO_ERROR);
}
/*
  transform in volume using mapping specified by mri_in_skull -> mri_ref_skull
*/
static int m3dMorphSkull(MORPH_3D *m3d, MRI_SURFACE *mris_in_skull, MRI_SURFACE *mris_ref_skull, MRI *mri)
{
  MHT *mht_in, *mht_ref;
  MORPH_NODE *mn;
  MNP *mnp;
  int i, j, k, width, height, depth, n;
  float in_x0, in_y0, in_z0, ref_x0, ref_y0, ref_z0, scale, nx, ny, nz, ref_v_x0, ref_v_y0, ref_v_z0, thick, dx, dy, dz,
      ref_dist, in_dist;
  double xw, yw, zw, xv, yv, zv, mean_scale;

  thick = mri->thick;

  /* build hash tables with in and ref skull position */

  in_x0 = mris_in_skull->xctr;
  in_y0 = mris_in_skull->yctr;
  in_z0 = mris_in_skull->zctr;
  ref_x0 = mris_ref_skull->xctr;
  ref_y0 = mris_ref_skull->yctr;
  ref_z0 = mris_ref_skull->zctr;

  dx = ref_x0 - in_x0;
  dy = ref_y0 - in_y0;
  dz = ref_z0 - in_z0;
/*  MRIStranslate(mris_in_skull, -dx, -dy, -dz) ;*/
  mht_in  = MHTcreateFaceTable(mris_in_skull);
  mht_ref = MHTcreateFaceTable(mris_ref_skull);

  /* find morph coordinates of ref origin */
  // MRIworldToVoxel(mri, ref_x0, ref_y0, ref_z0, &xv, &yv, &zv) ;
  MRIsurfaceRASToVoxel(mri, ref_x0, ref_y0, ref_z0, &xv, &yv, &zv);
  ref_v_x0 = xv * thick;
  ref_v_y0 = yv * thick;
  ref_v_z0 = zv * thick;
  width = m3d->width;
  height = m3d->height;
  depth = m3d->depth;
  mean_scale = 0.0;
  n = 0;
  for (i = 0; i < width; i++) {
    DiagHeartbeat((float)i / (float)(width - 1));
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        if (i == width / 2 && j == height / 2 && k == depth / 2) {
          DiagBreak();
          i = width / 2;
          j = height / 2;
          k = depth / 2;
        }

        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];

        /*
           find length of the ray from the center of ref surface to the
           ref skull passing through the morphed point
           */
        xv = mn->x / thick;
        yv = mn->y / thick;
        zv = mn->z / thick;
        MRIvoxelToWorld(mri, (double)xv, (double)yv, (double)zv, &xw, &yw, &zw);
        nx = xw - ref_x0;
        ny = yw - ref_y0;
        nz = zw - ref_z0;
        ref_dist = MRISdistanceToSurface(mris_ref_skull, mht_ref, ref_x0, ref_y0, ref_z0, nx, ny, nz);

        /*
           find length of the ray from the center of in surface to the
           in skull passing through the original point
           */
        in_dist = MRISdistanceToSurface(mris_in_skull, mht_in, ref_x0, ref_y0, ref_z0, nx, ny, nz);

        /* scale point by ratio */
        if (!FZERO(in_dist) && !FZERO(ref_dist)) {
          scale = ref_dist / in_dist;
          mean_scale += scale;
          n++;
          if (scale > MAX_DX)
            scale = MAX_DX;
          else if (scale < MIN_DX)
            scale = MIN_DX;
        }
        else {
          scale = 1.0;
          fprintf(stdout, "warning: node (%d,%d,%d) has zero in dist\n", i, j, k);
        }
        mnp->dx = scale;
      }
    }
  }

  /* set unassigned values to mean of assigned ones */
  mean_scale /= (double)n;
  fprintf(stdout, "mean brain scaling = %2.1f\n", mean_scale);

  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        if (i == width / 2 && j == height / 2 && k == depth / 2) {
          DiagBreak();
          i = width / 2;
          j = height / 2;
          k = depth / 2;
        }

        mnp = &m3d->pnodes[i][j][k];
        if (FZERO(mnp->dx)) mnp->dx = mean_scale;
      }
    }
  }
  {
    const char *cp;
    float sigma;
    cp = getenv("SIGMA");
    if (!cp) cp = "1.0";
    sigma = atof(cp);
    m3dblurDx(m3d, sigma);
  }
  for (i = 0; i < width; i++) {
    for (j = 0; j < height; j++) {
      for (k = 0; k < depth; k++) {
        if (i == width / 2 && j == height / 2 && k == depth / 2) {
          DiagBreak();
          i = width / 2;
          j = height / 2;
          k = depth / 2;
        }

        mn = &m3d->nodes[i][j][k];
        mnp = &m3d->pnodes[i][j][k];
        scale = mnp->dx;
        mn->x = (mn->x - ref_v_x0) * scale + ref_v_x0;
        mn->y = (mn->y - ref_v_y0) * scale + ref_v_y0;
        mn->z = (mn->z - ref_v_z0) * scale + ref_v_z0;
      }
    }
  }

  MHTfree(&mht_in);
  MHTfree(&mht_ref);
  m3dComputeMetricProperties(m3d);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIrigidAlign(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MATRIX *m_L)
{
  int nlevels, i;
  MRI *mri_in_pyramid[MAX_LEVELS], *mri_ref_pyramid[MAX_LEVELS];
  char base_name[STRLEN];
  /*  double tmul ;*/

  parms->rigid = 1;
  if (!parms->lta) parms->lta = LTAalloc(1, NULL);

  if (m_L) parms->lta->xforms[0].m_L = MatrixCopy(m_L, parms->lta->xforms[0].m_L);
  mriNormalizeStds(mri_ref);
  if (DZERO(parms->dt)) parms->dt = 1e-6;

  if (mri_ref->nframes > 1)
    mri_ref->mean = MRImeanFrame(mri_ref, 1);
  else
    mri_ref->mean = 1.0;

  strcpy(base_name, parms->base_name);
  openLogFile(parms);


  /* disable all the neck "don't care" stuff */
  parms->ref_np.neck_x0 = parms->ref_np.neck_y0 = parms->ref_np.neck_z0 = 1000;
  parms->ref_np.neck_dx = parms->ref_np.neck_dy = parms->ref_np.neck_dz = 1;
  parms->in_np.neck_x0 = parms->in_np.neck_y0 = parms->in_np.neck_z0 = 1000;
  parms->in_np.neck_dx = parms->in_np.neck_dy = parms->in_np.neck_dz = 1;

  /* build Gaussian pyramid */
  mri_in_pyramid[0] = mri_in;
  mri_ref_pyramid[0] = mri_ref;
  for (nlevels = 1; nlevels <= parms->max_levels; nlevels++) {
    if (mri_in_pyramid[nlevels - 1]->width <= MIN_PYR_WIDTH) break;
    mri_in_pyramid[nlevels] = MRIreduceByte(mri_in_pyramid[nlevels - 1], NULL);
    mri_ref_pyramid[nlevels] = MRIreduceMeanAndStdByte(mri_ref_pyramid[nlevels - 1], NULL);
  }

  for (i = nlevels - 1; i >= 0; i--) {
    /* convert transform to voxel coordinates for this level */
    MRIrasXformToVoxelXform(
        mri_in_pyramid[i], mri_ref_pyramid[i], parms->lta->xforms[0].m_L, parms->lta->xforms[0].m_L);
    if (parms->m_xform_mean) {
      MRIrasXformToVoxelXform(mri_in_pyramid[i], mri_ref_pyramid[i], parms->m_xform_mean, parms->m_xform_mean);
    }
    if (parms->m_xforms) {
      int sno;
      VECTOR *v = NULL, *vT = NULL;
      MATRIX *m_vvT = NULL, *m_voxel;

      if (parms->m_xform_covariance) MatrixFree(&parms->m_xform_covariance);
      if (parms->m_inv_cov) MatrixFree(&parms->m_inv_cov);

      /* compute covariance matrix in specified coordinate system */
      for (sno = 0; sno < parms->nxforms; sno++) {
        m_voxel = MRIrasXformToVoxelXform(mri_in_pyramid[i], mri_ref_pyramid[i], parms->m_xforms[sno], NULL);
        v = MatrixReshape(m_voxel, v, (m_voxel->rows - 1) * m_voxel->cols, 1);
        vT = MatrixTranspose(v, vT);
        m_vvT = MatrixMultiply(v, vT, m_vvT);
        if (!parms->m_xform_covariance)
          parms->m_xform_covariance = MatrixCopy(m_vvT, NULL);
        else
          MatrixAdd(m_vvT, parms->m_xform_covariance, parms->m_xform_covariance);
      }
      MatrixScalarMul(parms->m_xform_covariance, 1 / (double)parms->nxforms, parms->m_xform_covariance);
      parms->m_inv_cov = MatrixInverse(parms->m_xform_covariance, NULL);
      if (!parms->m_inv_cov) ErrorExit(ERROR_BADPARM, "%s: could not invert covariance matrix", Progname);
    }

    if (Gdiag & DIAG_SHOW) {
      printf("initial voxel transform:\n");
      MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);
    }

    fprintf(stdout, "aligning pyramid level %d.\n", i);
    if ((Gdiag & DIAG_WRITE) && parms->log_fp) fprintf(parms->log_fp, "aligning pyramid level %d.\n", i);
    mriQuasiNewtonLinearAlignPyramidLevel(mri_in_pyramid[i], mri_ref_pyramid[i], parms);

    /* convert transform to RAS coordinate representation */
    MRIvoxelXformToRasXform(
        mri_in_pyramid[i], mri_ref_pyramid[i], parms->lta->xforms[0].m_L, parms->lta->xforms[0].m_L);
    if (parms->m_xform_mean) {
      MRIvoxelXformToRasXform(mri_in_pyramid[i], mri_ref_pyramid[i], parms->m_xform_mean, parms->m_xform_mean);
    }
  }

  /* free Gaussian pyramid */
  for (i = 1; i < nlevels; i++) {
    MRIfree(&mri_in_pyramid[i]);
    MRIfree(&mri_ref_pyramid[i]);
  }
  strcpy(parms->base_name, base_name);
  if (parms->log_fp) {
    fclose(parms->log_fp);
    parms->log_fp = NULL;
  }

  /*  mriOrthonormalizeTransform(parms->lta->xforms[0].m_L) ;*/
  return (NO_ERROR);
}

/* for 3d vector macros */
#include "tritri.h"
static int mriOrthonormalizeTransform(MATRIX *m_L)
{
  double dot, c1[3], c2[3], c3[3], len;
  int i;

  for (i = 0; i < 3; i++) {
    c1[i] = *MATRIX_RELT(m_L, i + 1, 1);
    c2[i] = *MATRIX_RELT(m_L, i + 1, 2);
    c3[i] = *MATRIX_RELT(m_L, i + 1, 3);
  }

  /* make 1st column vector unit length */
  len = VLEN(c1);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c1, 1.0 / len, c1);

  /* project out component of 2nd vector in direction of 1st column vector */
  dot = DOT(c1, c2);
  for (i = 0; i < 3; i++) c2[i] -= dot * c1[i];

  /* make 2nd column vector unit length */
  len = VLEN(c2);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c2, 1.0 / len, c2);

  /* project out component of 3rd vector in direction of 1st column vector */
  dot = DOT(c1, c3);
  for (i = 0; i < 3; i++) c3[i] -= dot * c1[i];

  /* project out component of 3rd vector in direction of 2nd column vector */
  dot = DOT(c2, c3);
  for (i = 0; i < 3; i++) c3[i] -= dot * c2[i];

  /* make 3rd column vector unit length */
  len = VLEN(c3);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c3, 1.0 / len, c3);

  for (i = 0; i < 3; i++) {
    *MATRIX_RELT(m_L, i + 1, 1) = c1[i];
    *MATRIX_RELT(m_L, i + 1, 2) = c2[i];
    *MATRIX_RELT(m_L, i + 1, 3) = c3[i];
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#include "gca.h"

static MRI *g_mri_in, *g_mri_ref;
static MP *g_parms;
static GCA *g_gca;
static double g_clamp;
extern void (*user_call_func)(float[]);


static void dfp_step_func(int itno, float rms, void *vparms, float *p);
static void dfp_em_step_func(int itno, float rms, void *vparms, float *p);

static void dfp_step_func(int itno, float sse, void *vparms, float *p)
{
  MP *parms = (MP *)vparms;
  int i, row, col;
  float rms;

  printf("dfp_step_func: %03d: -log(p) = %2.1f\n", itno, sse);
  printf("transform: ( %.2f, %.2f, %.2f, %.2f)\n", p[1], p[2], p[3], p[4]);
  printf("transform: ( %.2f, %.2f, %.2f, %.2f)\n", p[5], p[6], p[7], p[8]);
  printf("transform: ( %.2f, %.2f, %.2f, %.2f)\n", p[9], p[10], p[11], p[12]);
  fflush(stdout);
  rms = sqrt(sse);
  if (parms->l_priors > 0) {
    MATRIX *m_save;
    float intensity_rms;

    m_save = parms->m_xform_mean;
    parms->m_xform_mean = NULL;
    intensity_rms = computeRigidAlignmentErrorFunctional(p);
    parms->m_xform_mean = m_save;
    intensity_rms = sqrt(intensity_rms);
    fprintf(stdout, "%03d: %2.3f (intensity=%2.3f)\n", parms->start_t + itno, rms, intensity_rms);
  }
  else {
    fprintf(stdout, "%03d: %2.3f\n", parms->start_t + itno, rms);
  }

  /* read out current transform */
  for (i = row = 1; row <= 4; row++) {
    for (col = 1; col <= 4; col++) {
      parms->lta->xforms[0].m_L->rptr[row][col] = p[i++];
    }
  }
  if ((parms->write_iterations > 0) && Gdiag & DIAG_WRITE && (((itno + 1) % parms->write_iterations) == 0)) {
    MATRIX *m_tmp, *m_save, *m_voxel;

    m_save = parms->lta->xforms[0].m_L;
    m_tmp = MRIvoxelXformToRasXform(g_mri_in, g_mri_ref, m_save, NULL);
    m_voxel = MRIrasXformToVoxelXform(parms->mri_in, parms->mri_ref, m_tmp, NULL);
    parms->lta->xforms[0].m_L = m_voxel;
    writeSnapshot(parms->mri_in, parms, parms->start_t + itno);
    MatrixFree(&m_voxel);
    MatrixFree(&m_tmp);
    parms->lta->xforms[0].m_L = m_save;
  }
  logIntegration(parms, parms->start_t + itno, (double)rms);
  fflush(stdout);
}

static void dfp_em_step_func(int itno, float sse, void *vparms, float *p)
{
  MP *parms = (MP *)vparms;
  int i, row, col;

  printf("dfp_em_step_func: %03d: -log(p) = %2.1f\n", parms->start_t + itno, sse);
  fflush(stdout);
  /* read out current transform */
  for (i = row = 1; row <= 3; row++)  // used to be 4 here ... tosa
  {
    for (col = 1; col <= 4; col++) {
      parms->lta->xforms[0].m_L->rptr[row][col] = p[i++];
    }
  }
//////////////////////////// new //////////////////////////
  ///////////////////////////////////////////////////////////
  if ((parms->write_iterations > 0) && Gdiag & DIAG_WRITE && (((itno + 1) % parms->write_iterations) == 0)) {
    MATRIX *m_tmp, *m_save, *m_voxel;

    m_save = parms->lta->xforms[0].m_L;
    m_tmp = MRIvoxelXformToRasXform(g_mri_in, g_mri_ref, m_save, NULL);
    m_voxel = MRIrasXformToVoxelXform(parms->mri_in, parms->mri_ref, m_tmp, NULL);
    parms->lta->xforms[0].m_L = m_voxel;
    writeSnapshot(parms->mri_in, parms, parms->start_t + itno);
    MatrixFree(&m_voxel);
    MatrixFree(&m_tmp);
    parms->lta->xforms[0].m_L = m_save;
  }
  logIntegration(parms, parms->start_t + itno, (double)sse);
  fflush(stdout);
}

static int mriQuasiNewtonLinearAlignPyramidLevel(MRI *mri_in, MRI *mri_ref, MP *parms)
{
  float p[5 * 4], fold, fnew;
  int row, col, i, iter, steps;
  MATRIX *m_L;

  /*  user_call_func = integration_step ;*/
  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      p[i++] = parms->lta->xforms[0].m_L->rptr[row][col];
    }
  }
  g_mri_in = mri_in;
  g_mri_ref = mri_ref;
  g_parms = parms;
  parms->mri_red_in = mri_in;
  parms->mri_red_ref = mri_ref;
  fnew = computeRigidAlignmentErrorFunctional(p);
  steps = 0;
  do {
    if (steps++ > 5) break;
    if (steps > 1) fprintf(stdout, "pass %d through quasi-newton minimization...\n", steps);
    fold = fnew;
    dfp_step_func(0, fnew, parms, p);

    OpenDFPMin(p,
               12,
               parms->tol,
               &iter,
               &fnew,
               computeRigidAlignmentErrorFunctional,
               computeRigidAlignmentGradient,
               dfp_step_func,
               parms,
               NULL);

    parms->start_t += iter;

    /* read out current transform */
    m_L = parms->lta->xforms[0].m_L;
    for (i = row = 1; row <= 3; row++) {
      for (col = 1; col <= 4; col++) {
        m_L->rptr[row][col] = p[i++];
      }
    }
    if (DIAG_VERBOSE_ON) {
      printf("after pass:transform: ( %.2f, %.2f, %.2f, %.2f)\n", p[1], p[2], p[3], p[4]);
      printf("                      ( %.2f, %.2f, %.2f, %.2f)\n", p[5], p[6], p[7], p[8]);
      printf("                      ( %.2f, %.2f, %.2f, %.2f)\n", p[9], p[10], p[11], p[12]);
    }

  } while ((fold - fnew) / fold > parms->tol);

  fflush(stdout);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void computeRigidAlignmentGradient(float *p, float *g)
{
  int width, height, depth, row, col, i, y1, y2, y3;
  VECTOR *v_X, *v_Y, *v_Yk; /* original and transformed coordinate systems */
  VECTOR *v_crop;
  MATRIX *m_tmp, *m_L, *m_dL, *m_dT_X_T, *m_L_inv, *m_dw_X_T, *m_crop;
  VECTOR *v_dT, *v_X_T, *v_dw; /* gradient of mri_ref */
  double ref_val, in_val, x1, x2, x3 /*, len*/;
  double delta = 0.0, dx, dy, dz, std, mean_std, ref_wt;
  MRI *mri_ref, *mri_in;
  MP *parms;

  mri_ref = g_mri_ref;
  parms = g_parms;
  mri_in = g_mri_in;
  mean_std = mri_ref->mean;
  if (FZERO(mean_std)) mean_std = 1.0;

  /* copy current matrix out of p into matrix format */
  m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  m_dL = MatrixAlloc(4, 4, MATRIX_REAL);
  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      m_L->rptr[row][col] = p[i++];
    }
  }
  m_L->rptr[4][1] = m_L->rptr[4][2] = m_L->rptr[4][3] = 0;
  m_L->rptr[4][4] = 1.0;
  m_L_inv = MatrixInverse(m_L, NULL);

  if (parms->mri_crop) {
    MATRIX *m_tmp;

    m_tmp = MRIvoxelXformToRasXform(mri_in, mri_ref, m_L_inv, NULL);
    m_crop = MRIvoxelXformToRasXform(mri_ref, parms->mri_crop, m_tmp, NULL);
    MatrixFree(&m_tmp);
    v_crop = VectorAlloc(4, MATRIX_REAL);
  }
  else
    v_crop = m_crop = NULL;

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  m_tmp = MatrixAlloc(4, 4, MATRIX_REAL);
  m_dT_X_T = MatrixAlloc(4, 4, MATRIX_REAL);
  m_dw_X_T = MatrixAlloc(4, 4, MATRIX_REAL);
  v_X = VectorAlloc(4, MATRIX_REAL);    /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL);    /* transformed (dst) coordinates */
  v_Yk = VectorAlloc(4, MATRIX_REAL);   /* transformed (dst) coordinates */
  v_X_T = RVectorAlloc(4, MATRIX_REAL); /* transpose of input coords */
  v_dT = VectorAlloc(4, MATRIX_REAL);   /* gradient of target image */
  v_dw = VectorAlloc(4, MATRIX_REAL);   /* gradient of weighting */

  v_Y->rptr[4][1] = 1.0f;

  for (y3 = 0; y3 < depth; y3++) {
    V3_Z(v_Y) = (double)y3;
    for (y2 = 0; y2 < height; y2++) {
      V3_Y(v_Y) = (double)y2;
      for (y1 = 0; y1 < width; y1++) {
        ref_val = (double)MRIgetVoxVal(mri_ref, y1, y2, y3, 0);
        ref_wt = MRIsampleReferenceWeighting(mri_ref, y1, y2, y3);
        if (FZERO(ref_wt)) /* "don't care" weighting */
          continue;        /* gradient will be 0 */

        V3_X(v_Y) = (double)y1;
        MatrixMultiply(m_L_inv, v_Y, v_X);
        MatrixTranspose(v_X, v_X_T);

        x1 = V3_X(v_X);
        x2 = V3_Y(v_X);
        x3 = V3_Z(v_X);
        if (parms->scout_flag && ((nint(x1) != width / 2) && (nint(x2) != height / 2) && (nint(x3) != depth / 2)))
          continue;

        if (mri_ref->nframes > 1)
          std = MRIgetVoxVal(mri_ref, y1, y2, y3, 1);
        else
          std = mean_std;
        std /= mean_std;
        if (DZERO(std)) std = 1;
        MRIsampleVolumeGradient(mri_ref, y1, y2, y3, &dx, &dy, &dz);
        if (x1 > -1 && x1 < width && x2 > -1 && x2 < height && x3 > -1 && x3 < depth) {
          if (parms->mri_crop) {
            MatrixMultiply(m_crop, v_X, v_crop);
            MRIsampleVolume(parms->mri_crop, V3_X(v_crop), V3_Y(v_crop), V3_Z(v_crop), &in_val);
            if (in_val > 0) continue;
          }
          MRIsampleVolume(mri_in, x1, x2, x3, &in_val);
        }
        else
          in_val = 0.0;
        delta = ref_wt * (ref_val - in_val) / std;

        V3_X(v_dT) = dx;
        V3_Y(v_dT) = dy;
        V3_Z(v_dT) = dz;
        MatrixMultiply(v_dT, v_X_T, m_dT_X_T);

        MatrixScalarMul(m_dT_X_T, ref_wt * delta, m_tmp);
        MatrixCheck(m_tmp);
        MatrixAdd(m_tmp, m_dL, m_dL);

        /* second term in product rule */
        MRIsampleReferenceWeightingGradient(mri_ref, y1, y2, y3, &dx, &dy, &dz);
        if (!FZERO(dx) || !FZERO(dy) || !FZERO(dx)) DiagBreak();
        V3_X(v_dw) = dx;
        V3_Y(v_dw) = dy;
        V3_Z(v_dw) = dz;
        MatrixMultiply(v_dw, v_X_T, m_dw_X_T);
        MatrixScalarMul(m_dw_X_T, delta * delta, m_tmp);
        MatrixAdd(m_tmp, m_dL, m_dL);

        MatrixCheck(m_dL);
      }
    }
  }

  if (parms->mri_crop) {
    MatrixFree(&m_crop);
    MatrixFree(&v_crop);
  }
  MatrixScalarMul(m_dL, parms->l_intensity / (double)(width * height * depth), m_dL);
  MatrixFree(&m_dT_X_T);
  MatrixFree(&v_X);
  MatrixFree(&v_Y);
  MatrixFree(&v_X_T);
  MatrixFree(&v_dT);
  MatrixFree(&m_tmp);
  MatrixFree(&m_L_inv);
  MatrixFree(&v_dw);
  MatrixFree(&m_dw_X_T);
  if (!FZERO(parms->factor) && !FEQUAL(parms->factor, 1.0f)) {
    MatrixScalarMul(m_dL, parms->factor, m_dL);
  }


  if (parms->m_xform_mean) {
    MATRIX *m_diff;
    VECTOR *v_diff, *v_C_inv_diff;

    m_diff = MatrixSubtract(m_L, parms->m_xform_mean, NULL);
    *MATRIX_RELT(m_diff, 1, 4) = *MATRIX_RELT(m_diff, 2, 4) = *MATRIX_RELT(m_diff, 3, 4) =
        0.0; /* ignore translation priors */
    v_diff = MatrixReshape(m_diff, NULL, 12, 1);
    v_C_inv_diff = MatrixMultiply(parms->m_inv_cov, v_diff, NULL);
    MatrixReshape(v_C_inv_diff, m_diff, 0, 0);
    MatrixScalarMul(m_diff, parms->l_priors, m_diff);
    MatrixAdd(m_diff, m_dL, m_dL);
    MatrixFree(&m_diff);
    VectorFree(&v_diff);
    VectorFree(&v_C_inv_diff);
  }

  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      /*      p[i] = m_L->rptr[row][col] ;*/
      g[i++] = m_dL->rptr[row][col];
    }
  }
  MatrixFree(&m_L);
  MatrixFree(&m_dL);
}
static float computeRigidAlignmentErrorFunctional(float *p)
{
  int y1, y2, y3, width, height, depth, row, col, i;
  VECTOR *v_X, *v_Y, *v_crop; /* original and transformed coordinate systems */
  double ref_val, in_val, x1, x2, x3;
  double sse = 0.0f, delta = 0.0, std, mean_std, wt;
  MATRIX *m_L, *m_L_inv, *m_crop;
  MRI *mri_ref, *mri_in;
  MP *parms;

  mri_ref = g_mri_ref;
  parms = g_parms;
  mri_in = g_mri_in;

  /* copy current matrix out of p into matrix format */
  m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      m_L->rptr[row][col] = p[i++];
    }
  }
  m_L->rptr[4][1] = m_L->rptr[4][2] = m_L->rptr[4][3] = 0;
  m_L->rptr[4][4] = 1.0;
  m_L_inv = MatrixInverse(m_L, NULL);

  if (parms->mri_crop) {
    MATRIX *m_tmp;

    m_tmp = MRIvoxelXformToRasXform(mri_in, mri_ref, m_L_inv, NULL);
    m_crop = MRIvoxelXformToRasXform(mri_ref, parms->mri_crop, m_tmp, NULL);
    MatrixFree(&m_tmp);
    v_crop = VectorAlloc(4, MATRIX_REAL);
  }
  else
    v_crop = m_crop = NULL;

  mean_std = mri_ref->mean;
  if (FZERO(mean_std)) mean_std = 1.0;
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;

  v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f;
  for (y3 = 0; y3 < depth; y3++) {
    V3_Z(v_Y) = y3;
    for (y2 = 0; y2 < height; y2++) {
      V3_Y(v_Y) = y2;
      for (y1 = 0; y1 < width; y1++) {
        V3_X(v_Y) = y1;

        v_X = MatrixMultiply(m_L_inv, v_Y, v_X);

        x1 = V3_X(v_X);
        x2 = V3_Y(v_X);
        x3 = V3_Z(v_X);
        if (x1 == 3 * width / 4 && x2 == 3 * height / 4 && x3 == 3 * depth / 4) DiagBreak();

        if (parms->scout_flag && ((nint(x1) != width / 2) && (nint(x2) != height / 2) && (nint(x3) != depth / 2)))
          continue;
        ref_val = (double)MRIgetVoxVal(mri_ref, y1, y2, y3, 0);
        if (mri_ref->nframes > 1)
          std = MRIseq_vox(mri_ref, y1, y2, y3, 1);
        else
          std = mean_std;
        std /= mean_std;
        if (DZERO(std)) std = 1.0;
        if (x1 > -1 && x1 < width && x2 > -1 && x2 < height && x3 > -1 && x3 < depth) {
          if (parms->mri_crop) {
            MatrixMultiply(m_crop, v_X, v_crop);
            MRIsampleVolume(parms->mri_crop, V3_X(v_crop), V3_Y(v_crop), V3_Z(v_crop), &in_val);
            if (in_val > 0) continue;
          }
          MRIsampleVolume(mri_in, x1, x2, x3, &in_val);
        }
        else /* out of bounds, assume in val is 0 */
          in_val = 0;

        wt = MRIsampleReferenceWeighting(mri_ref, y1, y2, y3);
        delta = wt * (ref_val - in_val) / std;
        if (fabs(delta) > 20) DiagBreak();
        sse += delta * delta;
        if (sse > 2000) DiagBreak();
        if (!std::isfinite(sse)) DiagBreak();
      }
    }
  }

  sse = sse * parms->l_intensity / (double)(width * height * depth);

  if (parms->mri_crop) {
    MatrixFree(&m_crop);
    MatrixFree(&v_crop);
  }
  if (parms->m_xform_mean) {
    MATRIX *m_diff, *m_tmp, *m_scalar;
    VECTOR *v_diff, *v_diff_T;

    m_diff = MatrixSubtract(m_L, parms->m_xform_mean, NULL);
    *MATRIX_RELT(m_diff, 1, 4) = *MATRIX_RELT(m_diff, 2, 4) = *MATRIX_RELT(m_diff, 3, 4) =
        0.0; /* ignore translation priors */
    v_diff = MatrixReshape(m_diff, NULL, 12, 1);
    v_diff_T = MatrixTranspose(v_diff, NULL);
    m_tmp = MatrixMultiply(parms->m_inv_cov, v_diff, NULL);
    m_scalar = MatrixMultiply(v_diff_T, m_tmp, NULL);
    sse += *MATRIX_RELT(m_scalar, 1, 1);
    MatrixFree(&m_diff);
    MatrixFree(&m_scalar);
    MatrixFree(&m_tmp);
    VectorFree(&v_diff);
    VectorFree(&v_diff_T);
  }

  MatrixFree(&v_X);
  MatrixFree(&v_Y);
  MatrixFree(&m_L);
  MatrixFree(&m_L_inv);

  return ((float)sse); /* rms */
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void  // in         out
computeEMAlignmentGradient(float *p, float *g)
{
  int width, height, depth, row, col, i, xn, yn, zn, xv, yv, zv, n, ndets = 0;
  VECTOR *v_X, *v_Y, *v_means, *v_X_T; /* original and transformed coordinate systems */
  MATRIX *m_L, *m_dL, *m_dI_X_T, *m_delI, *m_inv_cov, *m_tmp1, *m_tmp2;
  double dx, dy, dz, det, det_total;
  MRI *mri_in;
  MP *parms;
  GCA *gca;
  GCA_SAMPLE *gcas;
  TRANSFORM *transform;
  float vals[MAX_GCA_INPUTS];
  GC1D *gc;

  gca = g_gca;
  parms = g_parms;
  mri_in = g_mri_in;
  gcas = parms->gcas;

  /* copy current matrix out of p into matrix format */
  m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  m_dL = MatrixAlloc(4, 4, MATRIX_REAL);
  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      m_L->rptr[row][col] = p[i++];
    }
  }
  m_L->rptr[4][1] = m_L->rptr[4][2] = m_L->rptr[4][3] = 0;
  m_L->rptr[4][4] = 1.0;
#ifndef __OPTIMIZE__
#endif
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  m_dI_X_T = MatrixAlloc(4, 4, MATRIX_REAL);
  v_X = VectorAlloc(4, MATRIX_REAL);                 /* input (src) coordinates */
  v_X_T = RVectorAlloc(4, MATRIX_REAL);              /* transpose of input coords */
  v_Y = VectorAlloc(4, MATRIX_REAL);                 /* transformed (dst) coordinates */
  v_means = VectorAlloc(gca->ninputs, MATRIX_REAL);  // ninputs (x) 1
  m_tmp1 = m_tmp2 = NULL;

  m_delI = MatrixAlloc(4, gca->ninputs, MATRIX_REAL);                // 4 (x) 4
  m_inv_cov = MatrixAlloc(gca->ninputs, gca->ninputs, MATRIX_REAL);  // ninputs (x) ninputs

  transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL);
  MatrixCopy(m_L, ((LTA *)(transform->xform))->xforms[0].m_L);
  TransformInvert(transform, NULL);
  det_total = 0.0;
  ndets = 0;
  for (i = 0; i < parms->nsamples; i++) {
    xn = gcas[i].xp;
    yn = gcas[i].yp;
    zn = gcas[i].zp;
    // if xv is inside the volume ///////////////////////////////////////////////
    if (!GCApriorToSourceVoxel(gca, mri_in, transform, xn, yn, zn, &xv, &yv, &zv)) {
      gc = GCAfindSourceGC(gca, mri_in, transform, xv, yv, zv, gcas[i].label);

      /////////////////////// here is the one makes translation ///////////////
      /*      VECTOR_ELT(v_X, 4) = 1.;*/  ////////////////////////////////////////////////
      V3_X(v_X) = xv;
      V3_Y(v_X) = yv;
      V3_Z(v_X) = zv;
      MatrixTranspose(v_X, v_X_T);

      if (gc) {
        load_mean_vector(gc, v_means, gca->ninputs);
        load_inverse_covariance_matrix(gc, m_inv_cov, gca->ninputs);
        det = MatrixDeterminant(m_inv_cov);
        det_total += det;
        ndets++;
      }
      else {
        MatrixClear(v_means);
        MatrixIdentity(gca->ninputs, m_inv_cov);
      }

      load_vals(mri_in, xv, yv, zv, vals, gca->ninputs);

      /* construct gradient tensor */
      for (n = 0; n < gca->ninputs; n++) {
        // taking x+1, x-1 -> dx, y+1, y-1 -> dy, z+1, z-1 -> dz
        MRIsampleVolumeGradientFrame(mri_in, xv, yv, zv, &dx, &dy, &dz, n);
        *MATRIX_RELT(m_delI, 1, n + 1) = dx;
        *MATRIX_RELT(m_delI, 2, n + 1) = dy;
        *MATRIX_RELT(m_delI, 3, n + 1) = dz;
        *MATRIX_RELT(m_delI, 4, n + 1) = 1.;
        VECTOR_ELT(v_means, n + 1) -= vals[n];
      }
      // m_tmp1 = ninputs (x) 4
      m_tmp1 = MatrixMultiply(v_means, v_X_T, m_tmp1); /* (I(Lr) - A(r)) r' */
      MatrixMultiply(m_inv_cov, m_tmp1, m_tmp1);       /* inv(C) * (I(Lr) - A(r)) r' */
      m_tmp2 = MatrixMultiply(m_delI, m_tmp1, m_tmp2);

      MatrixCheck(m_tmp2);
      MatrixAdd(m_tmp2, m_dL, m_dL);

      MatrixCheck(m_dL);
    }
  }

  if (ndets > 0)
    det_total /= (double)ndets;
  else
    det_total = 0;
  MatrixScalarMul(m_dL, parms->dt * parms->l_intensity / (double)(det_total * gca->ninputs * parms->nsamples), m_dL);
  MatrixFree(&m_dI_X_T);
  MatrixFree(&v_X);
  MatrixFree(&m_inv_cov);
  if (m_tmp1) MatrixFree(&m_tmp1);
  if (m_tmp2) MatrixFree(&m_tmp2);
  MatrixFree(&v_X_T);

  for (i = row = 1; row <= 3; row++) {
    for (col = 1; col <= 4; col++) {
      /*      p[i] = m_L->rptr[row][col] ;*/
      g[i++] = m_dL->rptr[row][col];
    }
  }
  // copied only up to 12
  MatrixFree(&m_L);
  MatrixFree(&m_dL);
  MatrixFree(&m_delI);
  VectorFree(&v_means);
}
static float computeEMAlignmentErrorFunctional(float *p)
{
  float log_p;
  LTA *lta, *old_lta;
  MATRIX *m_L;
  int row, col, i;
  GCA *gca;
  MRI *mri_inputs;
  MP *parms;
  double clamp;

  parms = g_parms;
  gca = g_gca;
  clamp = g_clamp;
  mri_inputs = g_mri_in;

  /* copy current matrix out of p into matrix format */
  lta = LTAalloc(1, NULL);
  m_L = lta->xforms[0].m_L;
  for (i = row = 1; row <= 3; row++)  // row up to 3
  {
    for (col = 1; col <= 4; col++)  // column up to 4
    {
      m_L->rptr[row][col] = p[i++];
    }
  }
  m_L->rptr[4][1] = m_L->rptr[4][2] = m_L->rptr[4][3] = 0;
  m_L->rptr[4][4] = 1.0;
#ifndef __OPTIMIZE__
#endif
  old_lta = (LTA *)parms->transform->xform;
  parms->transform->xform = (void *)lta;
  log_p = GCAcomputeLogSampleProbability(gca, parms->gcas, g_mri_in, parms->transform, parms->nsamples, clamp);
  parms->transform->xform = (void *)old_lta;
  LTAfree(&lta);
  return (-log_p);
}


#define HANNING_RADIUS 100
int MRIsampleReferenceWeightingGradient(MRI *mri, int x, int y, int z, double *pdx, double *pdy, double *pdz)
{
  double dx, dy, dz, vp1, vm1;

  vp1 = MRIsampleReferenceWeighting(mri, x + 1, y, z);
  vm1 = MRIsampleReferenceWeighting(mri, x - 1, y, z);
  dx = (vp1 - vm1) / (2 * mri->xsize);

  vp1 = MRIsampleReferenceWeighting(mri, x, y + 1, z);
  vm1 = MRIsampleReferenceWeighting(mri, x, y - 1, z);
  dy = (vp1 - vm1) / (2 * mri->ysize);

  vp1 = MRIsampleReferenceWeighting(mri, x, y, z + 1);
  vm1 = MRIsampleReferenceWeighting(mri, x, y, z - 1);
  dz = (vp1 - vm1) / (2 * mri->zsize);

  *pdx = dx;
  *pdy = dy;
  *pdz = dz;
  return (NO_ERROR);
}
#define T0_X -1
#define T0_Y -10
#define T0_Z -17
double MRIsampleReferenceWeighting(MRI *mri, int x, int y, int z)
{
  return (1.0);
}

#define MAX_CLASSES 4
int MRIemAlign(MRI *mri_in, GCA *gca, MORPH_PARMS *parms, MATRIX *m_L)
{
  int i;
  char base_name[STRLEN];
  float pcurrent, pold;

  // if no transform
  if (!parms->transform) {
    parms->transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL);
    parms->lta = (LTA *)parms->transform->xform;
  }

  // if matrix is given
  if (m_L) {
    // copy to lta matrix                  m_L copied to parms->lta->xforms[0].m_L
    // this process is redundant, since mri_em_align calls as
    //             ...                            &parms, parms.lta->xforms[0].m_L)
    parms->lta->xforms[0].m_L = MatrixCopy(m_L, parms->lta->xforms[0].m_L);
    // assign transform->xform to the same one
    /* make sure transform and lta are the same (sorry - retrofitting!).... mri_em_register comment */
    parms->transform->xform = (void *)parms->lta;
  }
  // check transform
  if (parms->transform->type != LINEAR_VOX_TO_VOX) {
    fprintf(stdout, "parms->transform->type is %d\n", parms->transform->type);
    ErrorExit(-1, "MRIemAlign: transform type must be LINEAR_VOX_TO_VOX");
  }

  if (DZERO(parms->dt)) parms->dt = 1e-6;

  strcpy(base_name, parms->base_name);
  openLogFile(parms);

  /* disable all the neck "don't care" stuff */
  parms->ref_np.neck_x0 = parms->ref_np.neck_y0 = parms->ref_np.neck_z0 = 1000;
  parms->ref_np.neck_dx = parms->ref_np.neck_dy = parms->ref_np.neck_dz = 1;
  parms->in_np.neck_x0 = parms->in_np.neck_y0 = parms->in_np.neck_z0 = 1000;
  parms->in_np.neck_dx = parms->in_np.neck_dy = parms->in_np.neck_dz = 1;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("initial voxel transform:\n");
    MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);
  }

  fprintf(stdout, "Aligning input volume to GCA...\n");
  if ((Gdiag & DIAG_WRITE) && parms->log_fp) fprintf(parms->log_fp, "aligning input volume to GCA.\n");

  fprintf(stdout, "Transform matrix\n");
  MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);
  fprintf(stdout, "nsamples %d\n", parms->nsamples);
  /* E step */
  pcurrent = -GCAcomputeLogSampleProbability(gca, parms->gcas, mri_in, parms->transform, parms->nsamples, parms->clamp);

  i = 0;
  do {
    pold = pcurrent;
    /* M step */
    mriQuasiNewtonEMAlignPyramidLevel(mri_in, gca, parms);

    pcurrent =
        -GCAcomputeLogSampleProbability(gca, parms->gcas, mri_in, parms->transform, parms->nsamples, parms->clamp);
    i++;
    printf("outof QuasiNewtonEMA: %03d: -log(p) = %6.1f  tol %f\n", parms->start_t + i, pcurrent, parms->tol);
  } while (((pcurrent - pold) / (pold)) > parms->tol);

  strcpy(parms->base_name, base_name);
  if (parms->log_fp) {
    fclose(parms->log_fp);
    parms->log_fp = NULL;
  }

  /*  mriOrthonormalizeTransform(parms->lta->xforms[0].m_L) ;*/
  return (NO_ERROR);
}
int MRIquasiNewtonAlignVolumes(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms, MATRIX *m_L)
{
  int i;
  char base_name[STRLEN];
  float pcurrent, pold;

  if (parms->mri_in == NULL) parms->mri_in = mri_in;
  if (parms->mri_ref == NULL) parms->mri_ref = mri_ref;
  parms->l_intensity = 1.0;

  // if no transform
  if (!parms->transform) {
    parms->transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL);
    parms->lta = (LTA *)parms->transform->xform;
  }

  // if matrix is given
  if (m_L) {
    // copy to lta matrix                  m_L copied to parms->lta->xforms[0].m_L
    // this process is redundant, since mri_em_align calls as
    //             ...                            &parms, parms.lta->xforms[0].m_L)
    parms->lta->xforms[0].m_L = MatrixCopy(m_L, parms->lta->xforms[0].m_L);
    // assign transform->xform to the same one
    /* make sure transform and lta are the same (sorry - retrofitting!).... mri_em_register comment */
    parms->transform->xform = (void *)parms->lta;
  }
  // check transform
  if (parms->transform->type != LINEAR_VOX_TO_VOX) {
    fprintf(stdout, "parms->transform->type is %d\n", parms->transform->type);
    ErrorExit(-1, "MRIemAlign: transform type must be LINEAR_VOX_TO_VOX");
  }

  if (DZERO(parms->dt)) parms->dt = 1e-6;

  strcpy(base_name, parms->base_name);
  openLogFile(parms);

  /* disable all the neck "don't care" stuff */
  parms->ref_np.neck_x0 = parms->ref_np.neck_y0 = parms->ref_np.neck_z0 = 1000;
  parms->ref_np.neck_dx = parms->ref_np.neck_dy = parms->ref_np.neck_dz = 1;
  parms->in_np.neck_x0 = parms->in_np.neck_y0 = parms->in_np.neck_z0 = 1000;
  parms->in_np.neck_dx = parms->in_np.neck_dy = parms->in_np.neck_dz = 1;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    printf("initial voxel transform:\n");
    MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);
  }

  fprintf(stdout, "Aligning input volume to reference...\n");
  if ((Gdiag & DIAG_WRITE) && parms->log_fp) fprintf(parms->log_fp, "aligning input volume to reference.\n");

  fprintf(stdout, "Transform matrix\n");
  MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);

  /* E step */
  pcurrent = mriIntensitySSE(mri_in, mri_ref, parms->lta->xforms[0].m_L);

  i = 0;
  do {
    pold = pcurrent;
    /* M step */
    mriQuasiNewtonLinearAlignPyramidLevel(mri_in, mri_ref, parms);

    pcurrent = mriIntensitySSE(mri_in, mri_ref, parms->lta->xforms[0].m_L);

    i++;
    printf("outof QuasiNewtonEMA: %03d: -log(p) = %6.1f  tol %f\n", parms->start_t + i, pcurrent, parms->tol);
  } while (((pcurrent - pold) / (pold)) > parms->tol);

  strcpy(parms->base_name, base_name);
  if (parms->log_fp) {
    fclose(parms->log_fp);
    parms->log_fp = NULL;
  }

  /*  mriOrthonormalizeTransform(parms->lta->xforms[0].m_L) ;*/
  return (NO_ERROR);
}
#include "voxlist.h"
static int powell_minimize(VOXEL_LIST *vl_source,
                           VOXEL_LIST *vl_target,
                           MATRIX *mat,
                           float *pscale_factor,
                           MATRIX *m_constraint,
                           int map_both_ways);
static float compute_powell_sse(float *p);
static double compute_likelihood(
    VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L, float scale_factor, int map_both_ways);
/* compute mean squared error of two images with a transform */
static double compute_likelihood(
    VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L, float scale_factor, int map_both_ways)
{
  int x, y, z, width, height, depth, hwidth, hheight, hdepth, i;
  VECTOR *v1, *v2;
  MRI *mri_target, *mri_source;
  double sse, error;
  double d1, d2, xd, yd, zd;
  MATRIX *m_L_inv;

  m_L_inv = MatrixInverse(m_L, NULL);
  if (m_L_inv == NULL) {
    //    ErrorPrintf(ERROR_BADPARM, "compute_likelihood: singular matrix.") ;
    return (-(vl_source->nvox + vl_target->nvox) * 1e4);
  }

  mri_target = vl_target->mri2;
  mri_source = vl_source->mri2;

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  *MATRIX_RELT(v2, 4, 1) = 1.0;

  width = mri_target->width;
  height = mri_target->height;
  depth = mri_target->depth;
  hwidth = mri_source->width;
  hheight = mri_source->height;
  hdepth = mri_source->depth;

  /* go through both voxel lists and compute the sse
     map it to the source, and if the source hasn't been counted yet, count it.
  */

  sse = 0.0;
  for (i = 0; i < vl_source->nvox; i++) {
    x = vl_source->xi[i];
    y = vl_source->yi[i];
    z = vl_source->zi[i];

    V3_X(v1) = x;
    V3_Y(v1) = y;
    V3_Z(v1) = z;
    MatrixMultiply(m_L, v1, v2);
    d1 = MRIgetVoxVal(vl_source->mri2, x, y, z, 0);
    xd = V3_X(v2);
    yd = V3_Y(v2);
    zd = V3_Z(v2);
    if (xd < 0)
      xd = 0;
    else if (xd >= width - 1)
      xd = width - 1;
    if (yd < 0)
      yd = 0;
    else if (yd >= height - 1)
      yd = height - 1;
    if (zd < 0)
      zd = 0;
    else if (zd >= depth - 1)
      zd = depth - 1;
    MRIsampleVolume(vl_target->mri2, xd, yd, zd, &d2);
    error = scale_factor * d1 - d2;
    sse += error * error;
  }

  /* now count target voxels that weren't mapped to in union */
  if (map_both_ways)
    for (i = 0; i < vl_target->nvox; i++) {
      x = vl_target->xi[i];
      y = vl_target->yi[i];
      z = vl_target->zi[i];
      V3_X(v1) = x;
      V3_Y(v1) = y;
      V3_Z(v1) = z;
      MatrixMultiply(m_L_inv, v1, v2);
      d1 = MRIgetVoxVal(vl_target->mri2, x, y, z, 0);

      xd = V3_X(v2);
      yd = V3_Y(v2);
      zd = V3_Z(v2);
      if (xd < 0)
        xd = 0;
      else if (xd >= hwidth - 1)
        xd = hwidth - 1;
      if (yd < 0)
        yd = 0;
      else if (yd >= hheight - 1)
        yd = hheight - 1;
      if (zd < 0)
        zd = 0;
      else if (zd >= hdepth - 1)
        zd = hdepth - 1;
      MRIsampleVolume(vl_source->mri2, xd, yd, zd, &d2);
      error = d1 - d2;
      sse += error * error;
    }

  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&m_L_inv);
  return (-sqrt(sse / (double)(vl_target->nvox + vl_source->nvox)));
}
MATRIX *MRIpowellAlignImages(MRI *mri_source,
                             MRI *mri_target,
                             MATRIX *m_L,
                             float *pscale_factor,
                             MATRIX *m_constraint,
                             MRI *mri_source_mask,
                             MRI *mri_target_mask,
                             int map_both_ways)
{
  int i;
  float fmin, fmax;
  VOXEL_LIST *vl_target, *vl_source;

  if (m_L == NULL) m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target);

  if (mri_target_mask)
    vl_target = VLSTcreate(mri_target_mask, 1, 2, NULL, 0, 0);
  else {
    MRIvalRange(mri_target, &fmin, &fmax);
    vl_target = VLSTcreate(mri_target, 1, fmax + 1, NULL, 0, 0);
  }
  vl_target->mri2 = mri_target;
  for (i = 0; i < 1; i++) {
    if (mri_source_mask)
      vl_source = VLSTcreate(mri_source_mask, 1, 2, NULL, 0, 0);
    else {
      MRIvalRange(mri_source, &fmin, &fmax);
      vl_source = VLSTcreate(mri_source, 1, fmax + 1, NULL, 0, 0);
    }
    vl_source->mri2 = mri_source;
    powell_minimize(vl_source, vl_target, m_L, pscale_factor, m_constraint, map_both_ways);
    VLSTfree(&vl_source);
  }
  if (pscale_factor) printf("best matching intensity scaling = %2.4f\n", *pscale_factor);
  VLSTfree(&vl_target);
  //  write_snapshot(mri_source, mri_target, m_L, "after_powell", 0) ;
  return (m_L);
}
#define NPARMS (3 * 4)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-5

static VOXEL_LIST *Gvl_target, *Gvl_source;
static int Gmap_both_ways = 0;
static int powell_minimize(VOXEL_LIST *vl_source,
                           VOXEL_LIST *vl_target,
                           MATRIX *mat,
                           float *pscale_factor,
                           MATRIX *m_constraint,
                           int map_both_ways)
{
  float *p, **xi, fret, fstart, scale_factor, min_sse;
  int i, r, c, iter;
  p = vector(1, NPARMS + 1);
  xi = matrix(1, NPARMS + 1, 1, NPARMS + 1);
  i = 1;
  Gmap_both_ways = map_both_ways;
  p[i++] = *MATRIX_RELT(mat, 1, 4);
  p[i++] = *MATRIX_RELT(mat, 2, 4);
  p[i++] = *MATRIX_RELT(mat, 3, 4);
  for (r = 1; r <= 3; r++) {
    for (c = 1; c <= 3; c++) {
      p[i++] = *MATRIX_RELT(mat, r, c);
    }
  }
  p[i] = 1.0;  // scaling parameter

  Gvl_target = vl_target;
  Gvl_source = vl_source;
  for (r = 1; r <= NPARMS + 1; r++) {
    for (c = 1; c <= NPARMS + 1; c++) {
      xi[r][c] = r == c ? 1 : 0;
    }
  }

  if (m_constraint) {
    for (i = r = 1; r <= 4; r++) {
      for (c = 1; c <= 4; c++, i++) {
        if (*MATRIX_RELT(m_constraint, r, c) > 0) xi[i][i] = 0;  // remove this direction from powell set
      }
    }
  }


  min_sse = compute_powell_sse(p);
  if (pscale_factor)
    OpenPowell(p, xi, NPARMS + 1, TOL, &iter, &fret, compute_powell_sse);
  else
    OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
  scale_factor = p[NPARMS + 1];
  do {
    // reinitialize powell directions
    for (r = 1; r <= NPARMS + 1; r++) {
      for (c = 1; c <= NPARMS + 1; c++) {
        xi[r][c] = r == c ? 1 : 0;
      }
    }

    if (m_constraint) {
      for (i = r = 1; r <= 4; r++) {
        for (c = 1; c <= 4; c++, i++) {
          if (*MATRIX_RELT(m_constraint, r, c) > 0) xi[i][i] = 0;  // remove this direction from powell set
        }
      }
    }
    fstart = fret;
    if (pscale_factor)
      OpenPowell(p, xi, NPARMS + 1, TOL, &iter, &fret, compute_powell_sse);
    else
      OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);

    i = 1;
    *MATRIX_RELT(mat, 1, 4) = p[i++];
    *MATRIX_RELT(mat, 2, 4) = p[i++];
    *MATRIX_RELT(mat, 3, 4) = p[i++];
    for (r = 1; r <= 3; r++) {
      for (c = 1; c <= 3; c++) {
        *MATRIX_RELT(mat, r, c) = p[i++];
      }
    }
    scale_factor = p[i];
    *MATRIX_RELT(mat, 4, 1) = 0.0;
    *MATRIX_RELT(mat, 4, 2) = 0.0;
    *MATRIX_RELT(mat, 4, 3) = 0.0;
    *MATRIX_RELT(mat, 4, 4) = 1.0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("best alignment after powell: %2.3f (%d steps)\n", fret, iter);
    if ((fstart - fret) / fstart < TOL) break;
  } while (fret < fstart);

  free_matrix(xi, 1, NPARMS, 1, NPARMS);
  free_vector(p, 1, NPARMS);
  if (pscale_factor) *pscale_factor = scale_factor;
  return (NO_ERROR);
}
static float compute_powell_sse(float *p)
{
  static MATRIX *mat = NULL;
  float error, scale_factor;
  int i, r, c;

  if (mat == NULL) mat = MatrixAlloc(4, 4, MATRIX_REAL);
  i = 1;
  *MATRIX_RELT(mat, 1, 4) = p[i++];
  *MATRIX_RELT(mat, 2, 4) = p[i++];
  *MATRIX_RELT(mat, 3, 4) = p[i++];
  for (r = 1; r <= 3; r++) {
    for (c = 1; c <= 3; c++) {
      *MATRIX_RELT(mat, r, c) = p[i++];
    }
  }
  scale_factor = p[i];
  *MATRIX_RELT(mat, 4, 1) = 0.0;
  *MATRIX_RELT(mat, 4, 2) = 0.0;
  *MATRIX_RELT(mat, 4, 3) = 0.0;
  *MATRIX_RELT(mat, 4, 4) = 1.0;
  error = -compute_likelihood(Gvl_source, Gvl_target, mat, scale_factor, Gmap_both_ways);
  return (error);
}

extern void (*user_call_func)(float[]);

static int mriQuasiNewtonEMAlignPyramidLevel(MRI *mri_in, GCA *gca, MP *parms)
{
  float p[5 * 5], fold, fnew;
  int row, col, i, iter, steps, retval;
  MATRIX *m_L;

  user_call_func = showTran;

  fprintf(stdout, "Quasinewton: input matrix\n");
  MatrixPrintHires(stdout, parms->lta->xforms[0].m_L);

  /*  user_call_func = integration_step ;*/
  // copy transform into an array
  for (i = row = 1; row <= 4; row++) {
    for (col = 1; col <= 4; col++) {
      p[i++] = parms->lta->xforms[0].m_L->rptr[row][col];
    }
  }
  //////
  g_mri_in = mri_in;
  g_clamp = parms->clamp;
  g_gca = gca;
  g_parms = parms;
  //////
  parms->mri_red_in = mri_in;
  parms->gca_red = gca;
  //////
  fnew = computeEMAlignmentErrorFunctional(p);
  steps = 0;
  do {
    if (steps++ > 5) break;
    if (steps > 1) fprintf(stdout, "pass %d through quasi-newton minimization...\n", steps);
    fold = fnew;

    // retval is not to be trusted, not sure why
    retval = OpenDFPMin(p,
                        12,
                        parms->tol,
                        &iter,
                        &fnew,
                        //    void (*func)(float [])
                        computeEMAlignmentErrorFunctional,
                        //  void (*dfunc)(float [], float []), void (*stepfunc), parms
                        computeEMAlignmentGradient,
                        dfp_em_step_func,
                        parms,
                        user_call_func);

    if (iter == 0)  // parms can be corrupted if no step taken
      break;
    if (parms->rigid) mriOrthonormalizeTransform(parms->lta->xforms[0].m_L);
    parms->start_t += iter;
    /* read out current transform */
    m_L = parms->lta->xforms[0].m_L;

    for (i = row = 1; row <= 3; row++) {
      for (col = 1; col <= 4; col++) {
        m_L->rptr[row][col] = p[i++];
      }
    }
    ///////////// missing
    m_L->rptr[4][1] = m_L->rptr[4][2] = m_L->rptr[4][3] = 0.;
    m_L->rptr[4][4] = 1.;
  } while ((fold - fnew) / fold > parms->tol);

  if (parms->rigid) mriOrthonormalizeTransform(parms->lta->xforms[0].m_L);
  return (NO_ERROR);
}

MATRIX *MRIfaridAlignImages(MRI *mri_source, MRI *mri_target, MATRIX *m_L)
{
  // compute linear registration using Hanry Farid's method
  float fmin, fmax;
  VOXEL_LIST *vl_target, *vl_source;
  MRI *mri_tmp;

  if (m_L == NULL) m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target);

  MRIvalRange(mri_target, &fmin, &fmax);
  vl_target = VLSTcreate(mri_target, 1, fmax + 1, NULL, 0, 0);
  mri_tmp = MRIcopy(mri_target, NULL);  // store the transformed target
  if (mri_tmp == NULL) {
    fprintf(stdout, "Unable to allocate memory for linear alignment. Exit.\n");
    exit(0);
  }
  vl_target->mri = mri_target;
  vl_target->mri2 = mri_tmp;

  MRIvalRange(mri_target, &fmin, &fmax);
  vl_source = VLSTcreate(mri_source, 1, fmax + 1, NULL, 0, 0);
  vl_source->mri2 = mri_source;
  VLSTcomputeStats(vl_source);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("source mean =%g, std = %g\n", vl_source->mean, vl_source->std);
  fflush(stdout);
  farid_align(vl_source, vl_target, m_L);

  VLSTfree(&vl_source);

  VLSTfree(&vl_target);

  MRIfree(&mri_tmp);

  return (m_L);
}

static int ComputeStepTransform(
    VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, float cx, float cy, float cz, MATRIX *Minc)
{
  //(cx, cy, cz) is used to translate the coordinates system to be centered at it; may help reduce numerical errors
  int i, j, k, width, height, depth;
  int x, y, z;
  double d1, d2;
  float x1, y1, z1;

  double fx, fy, fz, ft, ck;
  MATRIX *sumCC, *invCC;
  VECTOR *sumCk;
#ifdef NPARMS
#undef NPARMS
#endif
#define NPARMS 13  // remove offset at 14
  float Mvec[NPARMS];

  //  float  evalues[4] ;
  // MATRIX *m_evectors ;

  double c[NPARMS];

  sumCC = MatrixAlloc(NPARMS, NPARMS, MATRIX_REAL);

  sumCk = VectorAlloc(NPARMS, MATRIX_REAL);

  for (j = 1; j <= NPARMS; j++) {
    VECTOR_ELT(sumCk, j) = 0;
    for (k = 1; k <= NPARMS; k++) {
      sumCC->rptr[j][k] = 0.0;
    }
  }
  for (i = 0; i < vl_source->nvox; i++) {
    x = vl_source->xi[i];
    y = vl_source->yi[i];
    z = vl_source->zi[i];
    d1 = MRIgetVoxVal(vl_source->mri2, x, y, z, 0);
    d2 = MRIgetVoxVal(vl_target->mri2, x, y, z, 0);

    ft = d1 - d2;

    /*    MRIsampleVolumeDerivativeScale(vl_source->mri2, x, y, z, 1, 0, 0, &fx, 0.5);
     MRIsampleVolumeDerivativeScale(vl_source->mri2, x, y, z, 0, 1, 0, &fy, 0.5);
     MRIsampleVolumeDerivativeScale(vl_source->mri2, x, y, z, 0, 0, 1, &fz, 0.5);
    */
    fx = (MRIgetVoxVal(vl_source->mri2, x + 1, y, z, 0) - MRIgetVoxVal(vl_source->mri2, x - 1, y, z, 0)) * 0.5;
    fy = (MRIgetVoxVal(vl_source->mri2, x, y + 1, z, 0) - MRIgetVoxVal(vl_source->mri2, x, y - 1, z, 0)) * 0.5;
    fz = (MRIgetVoxVal(vl_source->mri2, x, y, z + 1, 0) - MRIgetVoxVal(vl_source->mri2, x, y, z - 1, 0)) * 0.5;

    x1 = x;
    y1 = y;
    z1 = z;
    c[0] = x1 * fx;
    c[1] = y1 * fx;
    c[2] = z1 * fx;
    c[3] = x1 * fy;
    c[4] = y1 * fy;
    c[5] = z1 * fy;
    c[6] = x1 * fz;
    c[7] = y1 * fz;
    c[8] = z1 * fz;
    c[9] = fx;
    c[10] = fy;
    c[11] = fz;
    c[12] = -d1;  // c[13] = -1.0;
    ck = ft - d1 + x1 * fx + y1 * fy + z1 * fz;

    for (j = 1; j <= NPARMS; j++) {
      VECTOR_ELT(sumCk, j) += c[j - 1] * ck;
      for (k = 1; k <= NPARMS; k++) {
        sumCC->rptr[j][k] += c[j - 1] * c[k - 1];
      }
    }
  }

  for (j = 1; j <= NPARMS; j++) {
    VECTOR_ELT(sumCk, j) /= (float)vl_source->nvox;
    for (k = 1; k <= NPARMS; k++) {
      sumCC->rptr[j][k] /= (float)vl_source->nvox;
    }
  }

  invCC = MatrixInverse(sumCC, NULL);

  for (i = 0; i < NPARMS; i++) Mvec[i] = 0;

  if (!invCC) {
    Mvec[0] = 1.0;
    Mvec[4] = 1.0;
    Mvec[8] = 1.0;
    Mvec[12] = 1.0;

    Mvec[0] = 1.0;
    Mvec[4] = 1.0;
    Mvec[8] = 1.0;
    Mvec[12] = 1.0;
    if (NPARMS == 14) Mvec[13] = 0.0;  // offset
  }
  else {
    for (i = 1; i <= NPARMS; i++)
      for (j = 1; j <= NPARMS; j++) Mvec[i - 1] += invCC->rptr[i][j] * VECTOR_ELT(sumCk, j);
  }

  if (!Minc) {
    printf("This shouldn't happen, be sure to alloc Minc before calling this function\n");
  }

  Minc->rptr[1][1] = Mvec[0];
  Minc->rptr[1][2] = Mvec[1];
  Minc->rptr[1][3] = Mvec[2];
  Minc->rptr[2][1] = Mvec[3];
  Minc->rptr[2][2] = Mvec[4];
  Minc->rptr[2][3] = Mvec[5];
  Minc->rptr[3][1] = Mvec[6];
  Minc->rptr[3][2] = Mvec[7];
  Minc->rptr[3][3] = Mvec[8];
  Minc->rptr[4][1] = 0;
  Minc->rptr[4][2] = 0;
  Minc->rptr[4][3] = 0;

  Minc->rptr[1][4] = Mvec[9];
  Minc->rptr[2][4] = Mvec[10];
  Minc->rptr[3][4] = Mvec[11];
  Minc->rptr[4][4] = 1.0;

  //  printf("intensity scaling = %g, offset = %g\n", Mvec[12], Mvec[13]);

  MatrixFree(&sumCC);
  MatrixFree(&invCC);
  VectorFree(&sumCk);

  width = 128;
  height = 128;
  depth = 128;
  x1 = Mvec[0] * width + Mvec[1] * height + Mvec[2] * depth + Mvec[9];
  y1 = Mvec[3] * width + Mvec[4] * height + Mvec[5] * depth + Mvec[10];
  z1 = Mvec[6] * width + Mvec[7] * height + Mvec[8] * depth + Mvec[11];

  d2 = (x1 - width) * (x1 - width) + (y1 - height) * (y1 - height) + (z1 - depth) * (z1 - depth);

  if (d2 < 0.01) {
    if (NPARMS == 14)
      printf("final contrast scaling = %2.2f + %2.2f\n", Mvec[12], Mvec[13]);
    else
      printf("final contrast scaling = %2.2f\n", Mvec[12]);
    return 1;  // converged
  }
  else
    return 0;
}

static int farid_align(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L)
{
  int iterations;
  int i;
  int width, height, depth, x, y, z;
  int minX, maxX, minY, maxY, minZ, maxZ;
  float cx, cy, cz;
  VECTOR *v1, *v2;
  double d2, xd, yd, zd;
  int flag;
  MATRIX *Minc, *Mtmp;
  //  float  evalues[4] ;
  // MATRIX *m_evectors ;
  //  int  out_of_range ;

  Minc = MatrixAlloc(4, 4, MATRIX_REAL);

  width = vl_source->mri2->width;
  height = vl_source->mri2->height;
  depth = vl_source->mri2->depth;
  minX = width;
  maxX = 0;
  minY = height;
  maxY = 0;
  minZ = depth;
  maxZ = 0;

  // to save time, first compute the range of voxels that need to be filled
  // for transforming target volume
  cx = 0;
  cy = 0;
  cz = 0;
  for (i = 0; i < vl_source->nvox; i++) {
    x = vl_source->xi[i];
    y = vl_source->yi[i];
    z = vl_source->zi[i];
    cx += x;
    cy += y;
    cz += z;
    if (minX > x) minX = x;
    if (maxX < x) maxX = x;
    if (minY > y) minY = y;
    if (maxY < y) maxY = y;
    if (minZ > z) minZ = z;
    if (maxZ < z) maxZ = z;
  }

  // centroid of structure
  cx /= (float)vl_source->nvox;
  cy /= (float)vl_source->nvox;
  cz /= (float)vl_source->nvox;

  if (minX > 0) minX -= 1;
  if (maxX < width - 1) maxX += 1;
  if (minY > 0) minY -= 1;
  if (maxY < height - 1) maxY += 1;
  if (minZ > 0) minZ -= 1;
  if (maxZ < depth - 1) maxZ += 1;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("X:%d-%d; Y:%d-%d; Z:%d-%d\n", minX, maxX, minY, maxY, minZ, maxZ);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  *MATRIX_RELT(v2, 4, 1) = 1.0;
  flag = 0;
  for (iterations = 1; iterations <= 30; iterations++) {
    // apply current registration to target volume
    for (z = minZ; z <= maxZ; z++)
      for (y = minY; y <= maxY; y++)
        for (x = minX; x <= maxX; x++) {
          V3_X(v1) = x;
          V3_Y(v1) = y;
          V3_Z(v1) = z;
          MatrixMultiply(m_L, v1, v2);
          xd = V3_X(v2);
          yd = V3_Y(v2);
          zd = V3_Z(v2);
          if (xd < 0 || xd >= width - 1 || yd < 0 || yd >= height - 1 || zd < 0 || zd >= depth - 1)
            d2 = 0;
          else
            MRIsampleVolume(vl_target->mri, xd, yd, zd, &d2);

          MRIsetVoxVal(vl_target->mri2, x, y, z, 0, d2);
        }

    // find the alignment between src and transformed target
    flag = ComputeStepTransform(vl_source, vl_target, cx, cy, cz, Minc);
    // printf("Minc\n");
    // MatrixPrint(stdout, Minc);
    Mtmp = MatrixCopy(m_L, NULL);
    m_L = MatrixMultiply(Mtmp, Minc, m_L);
    MatrixFree(&Mtmp);

    //    printf("M_L\n");
    // MatrixPrint(stdout, m_L);

    if (flag) break;
  }

  MatrixFree(&Minc);
  VectorFree(&v1);
  VectorFree(&v2);
  if (iterations == 31) {
    //    printf("Not converged yet after 30 iterations \n");
    return 0;
  }
  else
    return 1;
}

static int powell_minimize_label(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *mat);
static float compute_powell_label_sse(float *p);
static double compute_label_likelihood(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L);

/* compute mean squared error of two label images with a transform */
static double compute_label_likelihood(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *m_L)
{
  int x, y, z, width, height, depth, hwidth, hheight, hdepth, i;
  VECTOR *v1, *v2;
  MRI *mri_target, *mri_source;
  double sse;
  double xd, yd, zd;
  int d1, d2;
  MATRIX *m_L_inv;

  m_L_inv = MatrixInverse(m_L, NULL);
  if (m_L_inv == NULL) ErrorExit(ERROR_BADPARM, "compute_label_likelihood: singular matrix.");

  mri_target = vl_target->mri2;
  mri_source = vl_source->mri2;

  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v1, 4, 1) = 1.0;
  *MATRIX_RELT(v2, 4, 1) = 1.0;

  width = mri_target->width;
  height = mri_target->height;
  depth = mri_target->depth;
  hwidth = mri_source->width;
  hheight = mri_source->height;
  hdepth = mri_source->depth;

  /* go through both voxel lists and compute the sse
     map it to the source, and if the source hasn't been counted yet, count it.
  */

  sse = 0.0;
  for (i = 0; i < vl_source->nvox; i++) {
    x = vl_source->xi[i];
    y = vl_source->yi[i];
    z = vl_source->zi[i];

    V3_X(v1) = x;
    V3_Y(v1) = y;
    V3_Z(v1) = z;
    MatrixMultiply(m_L, v1, v2);
    d1 = (int)MRIgetVoxVal(vl_source->mri2, x, y, z, 0);
    xd = (int)(V3_X(v2) + 0.5);
    yd = (int)(V3_Y(v2) + 0.5);
    zd = (int)(V3_Z(v2) + 0.5);
    if (xd < 0 || xd >= width - 1 || yd < 0 || yd >= height - 1 || zd < 0 || zd >= depth - 1)
      d2 = 0;
    else
      d2 = (int)MRIgetVoxVal(vl_target->mri2, xd, yd, zd, 0);

    if (d1 != d2) sse += 1;  // may need to weighted with the probability
  }

  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&m_L_inv);
  return (-sqrt(sse / (double)(vl_target->nvox + vl_source->nvox)));
}

MATRIX *MRIpowellAlignLabels(MRI *mri_source, MRI *mri_target, MATRIX *m_L)
{
  // align two label volumes
  VOXEL_LIST *vl_target, *vl_source;

  if (m_L == NULL) m_L = MRIgetVoxelToVoxelXform(mri_source, mri_target);

  //  MRIvalRange(mri_target, &fmin, &fmax) ;
  vl_target = VLSTcreate(mri_target, 1, 255, NULL, 0, 0);
  vl_target->mri2 = mri_target;

  //  MRIvalRange(mri_target, &fmin, &fmax) ;
  vl_source = VLSTcreate(mri_source, 1, 255, NULL, 0, 0);
  vl_source->mri2 = mri_source;
  powell_minimize_label(vl_source, vl_target, m_L);
  VLSTfree(&vl_source);

  VLSTfree(&vl_target);
  // write_snapshot(mri_source, mri_target, m_L, "after_powell", 0) ;
  return (m_L);
}

#ifdef NPARMS
#undef NPARMS
#endif
#define NPARMS (3 * 4)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-5

static int powell_minimize_label(VOXEL_LIST *vl_source, VOXEL_LIST *vl_target, MATRIX *mat)
{
  // used for registration of two label volumes
  float *p, **xi, fret, fstart;
  int i, r, c, iter;
  p = vector(1, NPARMS);
  xi = matrix(1, NPARMS, 1, NPARMS);
  i = 1;
  p[i++] = *MATRIX_RELT(mat, 1, 4);
  p[i++] = *MATRIX_RELT(mat, 2, 4);
  p[i++] = *MATRIX_RELT(mat, 3, 4);
  for (r = 1; r <= 3; r++) {
    for (c = 1; c <= 3; c++) {
      p[i++] = *MATRIX_RELT(mat, r, c);
    }
  }

  Gvl_target = vl_target;
  Gvl_source = vl_source;
  for (r = 1; r <= NPARMS; r++) {
    for (c = 1; c <= NPARMS; c++) {
      xi[r][c] = r == c ? 1 : 0;
    }
  }

  OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_label_sse);

  do {
    // reinitialize powell directions
    for (r = 1; r <= NPARMS; r++) {
      for (c = 1; c <= NPARMS; c++) {
        xi[r][c] = r == c ? 1 : 0;
      }
    }

    fstart = fret;
    OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_label_sse);

    i = 1;
    *MATRIX_RELT(mat, 1, 4) = p[i++];
    *MATRIX_RELT(mat, 2, 4) = p[i++];
    *MATRIX_RELT(mat, 3, 4) = p[i++];
    for (r = 1; r <= 3; r++) {
      for (c = 1; c <= 3; c++) {
        *MATRIX_RELT(mat, r, c) = p[i++];
      }
    }
    *MATRIX_RELT(mat, 4, 1) = 0.0;
    *MATRIX_RELT(mat, 4, 2) = 0.0;
    *MATRIX_RELT(mat, 4, 3) = 0.0;
    *MATRIX_RELT(mat, 4, 4) = 1.0;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("best alignment after powell: %2.3f (%d steps)\n", fret, iter);
    if ((fstart - fret) / fstart < TOL) break;
  } while (fret < fstart);

  free_matrix(xi, 1, NPARMS, 1, NPARMS);
  free_vector(p, 1, NPARMS);

  return (NO_ERROR);
}
static float compute_powell_label_sse(float *p)
{
  static MATRIX *mat = NULL;
  float error;
  int i, r, c;

  if (mat == NULL) mat = MatrixAlloc(4, 4, MATRIX_REAL);
  i = 1;
  *MATRIX_RELT(mat, 1, 4) = p[i++];
  *MATRIX_RELT(mat, 2, 4) = p[i++];
  *MATRIX_RELT(mat, 3, 4) = p[i++];
  for (r = 1; r <= 3; r++) {
    for (c = 1; c <= 3; c++) {
      *MATRIX_RELT(mat, r, c) = p[i++];
    }
  }

  *MATRIX_RELT(mat, 4, 1) = 0.0;
  *MATRIX_RELT(mat, 4, 2) = 0.0;
  *MATRIX_RELT(mat, 4, 3) = 0.0;
  *MATRIX_RELT(mat, 4, 4) = 1.0;
  error = -compute_label_likelihood(Gvl_source, Gvl_target, mat);
  return (error);
}
double MRIcomputeOptimalLinearXform(MRI *mri_source,
                                    MRI *mri_target,
                                    MATRIX *m_L,
                                    float min_angle,
                                    float max_angle,
                                    float min_scale,
                                    float max_scale,
                                    float min_trans,
                                    float max_trans,
                                    float angle_steps,
                                    float scale_steps,
                                    float trans_steps,
                                    int nreductions,
                                    char *base_name,
                                    int map_both_ways)
{
  MATRIX *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp, *m_L_tmp, *m_origin_inv, *m_tmp2, *m_scale, *m_trans,
      *m_tmp3 = NULL, *m_origin;
  double x_angle, y_angle, z_angle, x_max_rot;
  double y_max_rot, z_max_rot, delta_rot, src_means[3];
  double x_max_scale, y_max_scale, z_max_scale;
  double delta_scale, x_trans, delta_trans, y_trans, z_trans;
  double log_p, max_log_p, mean_angle, x_scale, y_scale, z_scale;
  double mean_scale, x_max_trans, y_max_trans, z_max_trans, mean_trans;
  int i, found, nreduced;
  VOXEL_LIST *vl_source, *vl_target;
  float fmin, fmax;
  MRI *mri_tmp;

  Gmap_both_ways = map_both_ways;
  MRIcenterOfMass(mri_source, src_means, 0);
  // unit matrix
  m_origin = MatrixIdentity(4, NULL);
  // set the translation  (center = 1)
  *MATRIX_RELT(m_origin, 1, 4) = src_means[0];
  *MATRIX_RELT(m_origin, 2, 4) = src_means[1];
  *MATRIX_RELT(m_origin, 3, 4) = src_means[2];
  *MATRIX_RELT(m_origin, 4, 4) = 1;

  m_trans = MatrixIdentity(4, NULL);
  m_origin_inv = MatrixCopy(m_origin, NULL);
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL;
  x_max_trans = y_max_trans = z_max_trans = x_max_rot = y_max_rot = z_max_rot = 0.0;
  x_max_scale = y_max_scale = z_max_scale = 1.0f;
  m_scale = MatrixIdentity(4, NULL);
  MRIvalRange(mri_source, &fmin, &fmax);
  vl_source = VLSTcreate(mri_source, 1, fmax + 1, NULL, 0, 0);
  vl_source->mri2 = mri_source;
  MRIvalRange(mri_target, &fmin, &fmax);
  vl_target = VLSTcreate(mri_target, 1, fmax + 1, NULL, 0, 0);
  mri_tmp = MRIcopy(mri_target, NULL);  // store the transformed target
  vl_target->mri2 = mri_tmp;
  VLSTcomputeStats(vl_source);
  if (Gdiag & DIAG_SHOW) printf("source mean =%g, std = %g\n", vl_source->mean, vl_source->std);
  max_log_p = compute_likelihood(vl_source, vl_target, m_L, 1.0, map_both_ways);
  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    MRI *mri_aligned = MRIlinearTransform(mri_source, NULL, m_L);
    sprintf(fname, "%s_init.mgz", base_name);
    printf("writing snapshot to %s\n", fname);
    MRIwrite(mri_aligned, fname);
    MRIfree(&mri_aligned);
  }

  found = nreduced = i = 0;
  do {
    found = 0;
    delta_trans = (max_trans - min_trans) / (trans_steps - 1);
    delta_scale = (max_scale - min_scale) / (scale_steps - 1);
    if (FZERO(delta_scale)) delta_scale = max_scale;

    delta_rot = (max_angle - min_angle) / (angle_steps - 1);
    if (Gdiag & DIAG_SHOW) {
      printf(
          "  scanning %2.2f degree nbhd (%2.1f)\n"
          "  scale %2.3f->%2.3f (step %2.3f), "
          "trans %2.2f->%2.2f (step %2.2f)\n",
          (float)DEGREES(max_angle),
          (float)DEGREES(delta_rot),
          min_scale,
          max_scale,
          delta_scale,
          min_trans,
          max_trans,
          delta_trans);
      fflush(stdout);
    }

    // scale /////////////////////////////////////////////////////////////
    for (x_scale = min_scale; x_scale <= max_scale; x_scale += delta_scale) {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale;
      for (y_scale = min_scale; y_scale <= max_scale; y_scale += delta_scale) {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale;
        for (z_scale = min_scale; z_scale <= max_scale; z_scale += delta_scale) {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) = *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp);
          MatrixMultiply(m_origin, m_tmp, m_scale);

          // angle //////////////////////////////
          for (x_angle = min_angle; x_angle <= max_angle; x_angle += delta_rot) {
            m_x_rot = MatrixReallocRotation(4, x_angle, X_ROTATION, m_x_rot);
            for (y_angle = min_angle; y_angle <= max_angle; y_angle += delta_rot) {
              m_y_rot = MatrixReallocRotation(4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp);
              for (z_angle = min_angle; z_angle <= max_angle; z_angle += delta_rot) {
                m_z_rot = MatrixReallocRotation(4, z_angle, Z_ROTATION, m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot);
                m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2);
                MatrixMultiply(m_origin, m_tmp2, m_rot);

                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2);
                m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3);

                // translation //////////
                for (x_trans = min_trans; x_trans <= max_trans; x_trans += delta_trans) {
                  *MATRIX_RELT(m_trans, 1, 4) = x_trans;
                  for (y_trans = min_trans; y_trans <= max_trans; y_trans += delta_trans) {
                    *MATRIX_RELT(m_trans, 2, 4) = y_trans;
                    for (z_trans = min_trans; z_trans <= max_trans; z_trans += delta_trans) {
                      *MATRIX_RELT(m_trans, 3, 4) = z_trans;

                      m_L_tmp = MatrixMultiply(m_trans, m_tmp3, m_L_tmp);
                      log_p = compute_likelihood(vl_source, vl_target, m_L_tmp, 1.0, Gmap_both_ways);

                      if (log_p > max_log_p) {
                        found = 1;
                        max_log_p = log_p;
                        x_max_scale = x_scale;
                        y_max_scale = y_scale;
                        z_max_scale = z_scale;
                        x_max_rot = x_angle;
                        y_max_rot = y_angle;
                        z_max_rot = z_angle;
                        x_max_trans = x_trans;
                        y_max_trans = y_trans;
                        z_max_trans = z_trans;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (Gdiag & DIAG_SHOW) {
      printf(
          "  max log p = %2.1f @ R=(%2.3f,%2.3f,%2.3f),"
          "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f,%2.1f,%2.1f)\n",
          max_log_p,
          DEGREES(x_max_rot),
          DEGREES(y_max_rot),
          DEGREES(z_max_rot),
          x_max_scale,
          y_max_scale,
          z_max_scale,
          x_max_trans,
          y_max_trans,
          z_max_trans);
    }

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) = *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f;
    *MATRIX_RELT(m_scale, 1, 1) = x_max_scale;
    *MATRIX_RELT(m_scale, 2, 2) = y_max_scale;
    *MATRIX_RELT(m_scale, 3, 3) = z_max_scale;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp);
    MatrixMultiply(m_origin, m_tmp, m_scale);

    x_max_scale = y_max_scale = z_max_scale = 1.0;

    /* update L to reflect new maximum and search around it */
    MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot);
    MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot);
    MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot);
    MatrixMultiply(m_y_rot, m_x_rot, m_tmp);
    MatrixMultiply(m_z_rot, m_tmp, m_rot);
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2);
    MatrixMultiply(m_origin, m_tmp2, m_rot);

    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2);
    m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3);

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max_trans;
    *MATRIX_RELT(m_trans, 2, 4) = y_max_trans;
    *MATRIX_RELT(m_trans, 3, 4) = z_max_trans;
    MatrixMultiply(m_trans, m_tmp3, m_L_tmp);

    MatrixCopy(m_L_tmp, m_L);

    x_max_trans = y_max_trans = z_max_trans = 0.0;

    /* we've translated transform by old maxs */
    x_max_rot = y_max_rot = z_max_rot = 0.0;
    if (Gdiag & DIAG_WRITE /* && found*/) {
      char fname[STRLEN];
      MRI *mri_aligned = MRIlinearTransform(mri_source, NULL, m_L);

      sprintf(fname, "%s_iter%d.mgz", base_name, i);
      printf("writing snapshot to %s\n", fname);
      MRIwrite(mri_aligned, fname);
      MRIfree(&mri_aligned);
    }
    if (found == 0)  // reduce scale
    {
      nreduced++;
      if (nreduced <= nreductions)
        printf("no new mimimum found, reduction %d of %d\n", nreduced, nreductions);
      else
        printf("terminating search\n");
      mean_scale = (max_scale + min_scale) / 2;
      delta_scale = (max_scale - min_scale) / 4;
      if (mean_scale - delta_scale < 1) min_scale = mean_scale - delta_scale;
      if (mean_scale + delta_scale > 1) max_scale = mean_scale + delta_scale;

      mean_trans = (max_trans + min_trans) / 2;
      delta_trans = (max_trans - min_trans) / 4;
      min_trans = mean_trans - delta_trans;
      max_trans = mean_trans + delta_trans;

      /* we've rotated transform to old max */

      mean_angle = (max_angle + min_angle) / 2;
      delta_rot = (max_angle - min_angle) / 4;
      min_angle = mean_angle - delta_rot;
      max_angle = mean_angle + delta_rot;
    }
    i++;
  } while (nreduced <= nreductions);

  MatrixFree(&m_x_rot);
  MatrixFree(&m_y_rot);
  MatrixFree(&m_z_rot);
  MatrixFree(&m_rot);
  MatrixFree(&m_tmp);
  MatrixFree(&m_origin_inv);
  MatrixFree(&m_tmp2);
  MatrixFree(&m_trans);
  MatrixFree(&m_tmp3);

  return (max_log_p);
}

static double mriIntensitySSE(MRI *mri_in, MRI *mri_ref, MATRIX *m_L)
{
  int x, y, z;
  VECTOR *v1, *v2;
  double xd, yd, zd, sse, width, height, depth;
  MATRIX *m_L_inv;
  double d1, d2, error;

  m_L_inv = MatrixInverse(m_L, NULL);
  if (m_L_inv == NULL) {
    int mx;
    mx = MAX(MAX(mri_in->width, mri_in->height), mri_in->depth);
    mx = MAX(mx, MAX(MAX(mri_ref->width, mri_ref->height), mri_ref->depth));
    return (mx * mx * mx * 10);
  }
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0;

  width = mri_ref->width;
  height = mri_ref->height;
  depth = mri_ref->depth;
  for (sse = 0.0, x = 0; x < mri_in->width; x++) {
    V3_X(v1) = x;
    for (y = 0; y < mri_in->height; y++) {
      V3_Y(v1) = y;
      for (z = 0; z < mri_in->depth; z++) {
        V3_Z(v1) = z;
        d1 = MRIgetVoxVal(mri_in, x, y, z, 0);
        MatrixMultiply(m_L, v1, v2);
        xd = V3_X(v2);
        yd = V3_Y(v2);
        zd = V3_Z(v2);
        if (xd < 0)
          xd = 0;
        else if (xd >= width - 1)
          xd = width - 1;
        if (yd < 0)
          yd = 0;
        else if (yd >= height - 1)
          yd = height - 1;
        if (zd < 0)
          zd = 0;
        else if (zd >= depth - 1)
          zd = depth - 1;
        MRIsampleVolume(mri_ref, xd, yd, zd, &d2);
        error = d1 - d2;
        sse += error * error;
      }
    }
  }

  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&m_L_inv);
  return (sse);
}
