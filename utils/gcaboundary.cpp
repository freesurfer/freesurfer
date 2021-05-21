/**
 * @brief Routines for boundary deformations based on gca.
 *
 * Reference:
 * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
 * Structures in the Human Brain", Fischl et al.
 * (2002) Neuron, 33:341-355.
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
#include <stdlib.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "gcaboundary.h"
#include "macros.h"
#include "matrix.h"
#include "talairachex.h"
#include "transform.h"
#include "tritri.h"

#define GCAB_IIN 0
#define GCAB_IOUT 1
#define GCAB_GRAD 2

#define GRAD_SIGMA 1.0

static double gcabSamplePDF(GCAB *gcab, GCABS *gcabs, int which, int vno, double val);

static int gcabUpdateDistributions(GCAB *gcab,
                                   MRI *mri_int,
                                   MRI *mri_seg,
                                   MRI *mri_dist,
                                   TRANSFORM *transform,
                                   int x,
                                   int y,
                                   int z,
                                   float xn,
                                   float yn,
                                   float zn,
                                   int target_label);

static int gcabUpdateNodeIin(GCAB *gcab, int xn, int yn, int zn, double val1, int vno, float wt);
static int gcabUpdateNodeIout(GCAB *gcab, int xn, int yn, int zn, double val2, int vno, float wt);
static MRI *gcabWritePDFToMRI(GCAB *gcab, MRI *mri, int x, int y, int z);
static int gcabSmoothPDF(GCAB *gcab, int x, int y, int z, MRI *mri_kernel);
static MRI *gcabWriteMRIToPDF(GCAB *gcab, MRI *mri, int x, int y, int z);
static int gcabComputeBorderNormal(GCAB *gcab,
                                   MRI *mri_seg,
                                   TRANSFORM *transform,
                                   int x,
                                   int y,
                                   int z,
                                   float *pnx,
                                   float *pny,
                                   float *pnz,
                                   int target_label);

static int gcabUpdateNode(GCAB *gcab, int xn, int yn, int zn, double val1, double val2, float grad, int vno, float wt);
static int gcabFindClosestIcoVertex(
    GCAB *gcab, TRANSFORM *transform, float x0, float y0, float z0, float nx, float ny, float nz);
static int gcabValsToBins(GCAB *gcab, double val1, double val2, float grad, float *pi1f, float *pi2f, float *pgf);
GCAB *GCABalloc(GCA *gca, int spacing, int ico_order, int nintensity_bins, int ngrad_bins, int target_label)
{
  GCAB *gcab;
  GCABS *gcabs;
  int width, height, depth, x, y, z, vno, i1, x1, y1, z1;

  gcab = (GCAB *)calloc(1, sizeof(GCAB));
  if (gcab == NULL) ErrorExit(ERROR_NOMEMORY, "GCABalloc: could not alloc GCAB");

  gcab->target_label = target_label;
#define PAD 5
  // find bounding box for label, then expand it and crop it to be in fov
  GCAstructureBoundingBox(gca, target_label, &gcab->bounding_box);
  gcab->bounding_box.x -= PAD;
  gcab->bounding_box.y -= PAD;
  gcab->bounding_box.z -= PAD;
  gcab->bounding_box.dx += 2 * PAD;
  gcab->bounding_box.dy += 2 * PAD;
  gcab->bounding_box.dz += 2 * PAD;
  x1 = gcab->bounding_box.x + gcab->bounding_box.dx;
  y1 = gcab->bounding_box.y + gcab->bounding_box.dy;
  z1 = gcab->bounding_box.z + gcab->bounding_box.dz;
  gcab->bounding_box.x = MAX(0, gcab->bounding_box.x);
  gcab->bounding_box.y = MAX(0, gcab->bounding_box.y);
  gcab->bounding_box.z = MAX(0, gcab->bounding_box.z);
  x1 = MIN(x1, gca->width);
  y1 = MIN(y1, gca->height);
  z1 = MIN(z1, gca->depth);
  gcab->bounding_box.dx = x1 - gcab->bounding_box.x;
  gcab->bounding_box.dy = y1 - gcab->bounding_box.y;
  gcab->bounding_box.dz = z1 - gcab->bounding_box.z;

  gcab->min_intensity = 0.0f;
  gcab->max_intensity = 255.0f;
  gcab->max_grad = 15;
  gcab->min_grad = -15;
  gcab->spacing = spacing;
  gcab->nintensity_bins = nintensity_bins;
  gcab->ngrad_bins = ngrad_bins;
  gcab->gca = gca;
  gcab->ico_order = ico_order;
  gcab->ico = read_icosahedron_by_order(ico_order);

  gcab->iscale = (gcab->nintensity_bins - 1) / (gcab->max_intensity - gcab->min_intensity);
  gcab->gscale = (gcab->ngrad_bins - 1) / (gcab->max_grad - gcab->min_grad);

  if (gcab->ico == NULL) ErrorExit(ERROR_NOMEMORY, "GCABalloc: could not read ico");
  gcab->nvertices = gcab->ico->nvertices;

  gcab->width = width = (int)(((float)gcab->bounding_box.dx / spacing) + .99);
  gcab->height = height = (int)((float)gcab->bounding_box.dy / spacing + .99);
  gcab->depth = depth = (int)(((float)gcab->bounding_box.dz / spacing) + .99);
  gcab->bs = (GCABS ***)calloc(width, sizeof(GCABS **));
  if (!gcab->bs) ErrorExit(ERROR_NOMEMORY, "GCABalloc: could not allocate BS");

  // setting vlaues gcab->bs volume
  for (x = 0; x < width; x++) {
    gcab->bs[x] = (GCABS **)calloc(height, sizeof(GCABS *));
    if (!gcab->bs[x]) ErrorExit(ERROR_NOMEMORY, "GCABalloc: could not allocate %dth **", x);

    for (y = 0; y < height; y++) {
      gcab->bs[x][y] = (GCABS *)calloc(depth, sizeof(GCABS));
      if (!gcab->bs[x][y]) ErrorExit(ERROR_NOMEMORY, "GCABalloc: could not allocate %d,%dth *", x, y);
      for (z = 0; z < depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        gcabs = &gcab->bs[x][y][z];
        gcabs->pdfs = (double ***)calloc(gcab->nvertices, sizeof(double **));
        gcabs->h_grad_pdfs = (HISTOGRAM **)calloc(gcab->nvertices, sizeof(HISTOGRAM *));
        gcabs->h_Iin_pdfs = (HISTOGRAM **)calloc(gcab->nvertices, sizeof(HISTOGRAM *));
        gcabs->h_Iout_pdfs = (HISTOGRAM **)calloc(gcab->nvertices, sizeof(HISTOGRAM *));
        gcabs->ntraining = (float *)calloc(gcab->nvertices, sizeof(float));
        for (vno = 0; vno < gcab->nvertices; vno++) {
          gcabs->pdfs[vno] = (double **)calloc(nintensity_bins, sizeof(double *));
          gcabs->h_grad_pdfs[vno] = HISTOinit(NULL, gcab->ngrad_bins, gcab->min_grad, gcab->max_grad);
          gcabs->h_Iin_pdfs[vno] = HISTOinit(NULL, gcab->nintensity_bins, gcab->min_intensity, gcab->max_intensity);
          gcabs->h_Iout_pdfs[vno] = HISTOinit(NULL, gcab->nintensity_bins, gcab->min_intensity, gcab->max_intensity);
          for (i1 = 0; i1 < nintensity_bins; i1++) {
            gcabs->pdfs[vno][i1] = (double *)calloc(nintensity_bins, sizeof(double));
          }
        }
      }
    }
  }
  return (gcab);
}

int GCABfree(GCAB **pgcab)
{
  GCAB *gcab;

  gcab = *pgcab;
  *pgcab = NULL;

  free(gcab);
  return (NO_ERROR);
}

int GCABcompleteTraining(GCAB *gcab)
{
  int x, y, z, i1, i2, vno;
  GCABS *gcabs;

  for (x = 0; x < gcab->width; x++)
    for (y = 0; y < gcab->height; y++)
      for (z = 0; z < gcab->depth; z++) {
        gcabs = &gcab->bs[x][y][z];
        for (vno = 0; vno < gcab->ico->nvertices; vno++) {
          if (gcabs->ntraining[vno] <= 0.0) continue;

          //          HISTOsmooth(gcabs->h_Iin_pdfs[vno], gcabs->h_Iin_pdfs[vno], .5) ;
          //          HISTOsmooth(gcabs->h_Iout_pdfs[vno], gcabs->h_Iout_pdfs[vno], .5) ;
          //          HISTOsmooth(gcabs->h_grad_pdfs[vno], gcabs->h_grad_pdfs[vno], .5) ;
          HISTOmakePDF(gcabs->h_grad_pdfs[vno], gcabs->h_grad_pdfs[vno]);
          HISTOmakePDF(gcabs->h_Iin_pdfs[vno], gcabs->h_Iin_pdfs[vno]);
          HISTOmakePDF(gcabs->h_Iout_pdfs[vno], gcabs->h_Iout_pdfs[vno]);
          for (i1 = 0; i1 < gcab->nintensity_bins; i1++)
            for (i2 = 0; i2 < gcab->nintensity_bins; i2++) {
              gcabs->pdfs[vno][i1][i2] /= gcabs->ntraining[vno];
            }
        }
      }

  GCABsmoothPDFs(gcab, 1);
  return (NO_ERROR);
}

#define SDIST .25

double GCABgetProbability(GCAB *gcab,
                          MRI *mri_int,
                          MRI *mri_dist,
                          TRANSFORM *transform,
                          float x0,
                          float y0,
                          float z0,
                          float nx,
                          float ny,
                          float nz)
{

  // Note: warnings fix: commented out variables should be refactored into the #else below if ever used again - clarsen

  float xn, yn, zn, xi, yi, zi, xo, yo, zo, xi0, yi0, zi0, d, dx, dy, dz, xo0, yo0, zo0;
  // float  xd, yd, zd;
  int vno, oob,  x, y, z, num;
  // int xu, yu, zu, xb, yb, zb;
  double grad, p, pdelI, mag, val1, val2, pI_in, pI_out, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2;
  // double pI_total;
  GCABS *gcabs;

  oob = GCABsourceVoxelToNode(gcab, mri_int, transform, x0, y0, z0, &xn, &yn, &zn);
  if (oob != NO_ERROR) return (0.0);
  xo = x0 + SDIST * nx;
  yo = y0 + SDIST * ny;
  zo = z0 + SDIST * nz;
  xi = x0 - SDIST * nx;
  yi = y0 - SDIST * ny;
  zi = z0 - SDIST * nz;
  dx = xo - xi;
  dy = yo - yi;
  dz = zo - zi;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(d)) DiagBreak();
  dx /= d;
  dy /= d;
  dz /= d;
  MRIsampleVolumeDerivativeScale(mri_int, x0, y0, z0, nx, ny, nz, &grad, GRAD_SIGMA);
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  if (grad < gcab->min_grad) grad = gcab->min_grad;
  MRIsampleVolume(mri_int, xo, yo, zo, &val2);
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
  MRIsampleVolume(mri_int, xi, yi, zi, &val1);
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  vno = gcabFindClosestIcoVertex(gcab, transform, x0, y0, z0, nx, ny, nz);
  if (vno < 0) return (0);

  // xb = (int)xn;
  // yb = (int)yn;
  // zb = (int)zn;
  // xd = xn - xb;
  // yd = yn - yb;
  // zd = zn - zb;
  // xu = xb + 1;
  // yu = yb + 1;
  // zu = zb + 1;

  // pI_total = pdelI = 0;
#if 1
  x = nint(xn);
  y = nint(yn);
  z = nint(zn);
  gcabs = &gcab->bs[x][y][z];
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri;

    char fname[STRLEN];

    mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
    sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
    MRIwrite(mri, fname);
    MRIfree(&mri);
    sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_grad_pdfs[vno], fname);
    sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
    sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
  }
  pdelI = gcabSamplePDF(gcab, gcabs, GCAB_GRAD, vno, grad);  // p(gradient)
  pdelI *= exp(-((gcab->max_grad + 1) - fabs(grad)));
  pI_out = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iout)

  // compute probability of interior
  xi0 = xi;
  yi0 = yi;
  zi0 = zi;
  pI_in = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
  //  pI_in = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
  for (d = 0.0, num = 0; d > -10; d -= 0.5)  // start 1/2mm inside border
  {
    xi = xi0 + d * dx;
    yi = yi0 + d * dy;
    zi = zi0 + d * dz;
    MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
    if (mag < 0)  // closer to the opposite side of the surface
      break;
    num++;
    MRIsampleVolume(mri_int, xi, yi, zi, &val1);
    if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
    if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
    p = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
    //    p = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
    pI_in = MIN(pI_in, p);
  }

  // build two tangent vectors for sampling exterior values
  e1x = -ny;
  e1y = -nz;
  e1z = nx;  // non-parallel vector
  CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z);
  CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z);
  mag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  if (FZERO(mag)) mag = 1e-5;
  e1x /= mag;
  e1y /= mag;
  e1z /= mag;
  mag = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
  if (FZERO(mag)) mag = 1e-5;
  e2x /= mag;
  e2y /= mag;
  e2z /= mag;

  xo0 = xo;
  yo0 = yo;
  zo0 = zo;
  for (d1 = -1; d1 <= 1.0; d1++)
    for (d2 = -1.0; d2 <= 1.0; d2++) {
      xo = xo0 + d1 * e1x + d2 * e2x;
      yo = yo0 + d1 * e1y + d2 * e2y;
      zo = zo0 + d1 * e1z + d2 * e2z;
      MRIsampleVolume(mri_int, xo, yo, zo, &val2);
      if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
      if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
      p = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iin)
      pI_out = MIN(pI_out, p);
    }

  p = pow(pI_out * pI_in * pdelI, 0.3333);

#else
  for (x = xb; x <= xu; x++) {
    for (y = yb; y <= yu; y++) {
      for (z = zb; z <= zu; z++) {
        wt = (x == xb ? (1 - xd) : xd);
        wt *= (y == yb ? (1 - yd) : yd);
        wt *= (z == zb ? (1 - zd) : zd);

        gcabs = &gcab->bs[x][y][z];
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
          MRI *mri;

          char fname[STRLEN];

          mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
          sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
          MRIwrite(mri, fname);
          MRIfree(&mri);
          sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_grad_pdfs[vno], fname);
          sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
          sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
        }
        xi0 = xi;
        yi0 = yi;
        zi0 = zi;
        pI_in = 1;
        pI_out = pI = 0.0;
        // compute joint probability of all interior voxels with exterior one
        for (d = 0.0, num = 0; d > -10; d--)  // start 1/2mm inside border
        {
          xi = xi0 + d * dx;
          yi = yi0 + d * dy;
          zi = zi0 + d * dz;
          MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
          if (mag < 0)  // closer to the opposite side of the surface
            break;
          num++;
          MRIsampleVolume(mri_int, xi, yi, zi, &val1);
          if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
          if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
          gcabValsToBins(gcab, val1, val2, grad, &i1f, &i2f, &gbinf);
          i1b = (int)i1f;
          i2b = (int)i2f;
          gbinb = (int)gbinf;
          i1u = i1b + 1;
          i2u = i2b + 1;
          gbinu = gbinb + 1;
          i1d = i1f - i1b;
          i2d = i2f - i2b;
          gbind = gbinf - gbinb;

#if 0
          pI += gcabs->pdfs[vno][i1b][i2b] * wt * (1.0-i1d) * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1b][i2u] * wt * (1.0-i1d) * i2d;
          pI += gcabs->pdfs[vno][i1u][i2b] * wt * i1d       * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1u][i2u] * wt * i1d       * i2d;
#else
          pI += gcabs->pdfs[vno][nint(i1f)][nint(i2f)];
          pI_in = MIN(pI_in, gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]);
#endif
        }

        if (num > 0) pI_total += pI / num;
#if 0
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinb] * (1.0-gbind);
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinu] * gbind;
#else
        break;
#endif
      }
    }
  }
  p = pI_total * pdelI;
  p = pI_out * pI_in;
#endif

  return (p);
}

double GCABgetPin(GCAB *gcab,
                  MRI *mri_int,
                  MRI *mri_dist,
                  TRANSFORM *transform,
                  float x0,
                  float y0,
                  float z0,
                  float nx,
                  float ny,
                  float nz)
{

  float xn, yn, zn, xi, yi, zi, xo, yo, zo, xi0, yi0, zi0, d, dx, dy, dz, xo0, yo0, zo0;
  // float xd, yd, zd;
  int vno, oob,  x, y, z, num;
  // int xb, yb, zb, xu, yu, zu;
  double grad, p, pdelI, mag, val1, val2, pI_in, pI_out, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2;
  // double pI_total;
  GCABS *gcabs;

  oob = GCABsourceVoxelToNode(gcab, mri_int, transform, x0, y0, z0, &xn, &yn, &zn);
  if (oob != NO_ERROR) return (0.0);
  xo = x0 + SDIST * nx;
  yo = y0 + SDIST * ny;
  zo = z0 + SDIST * nz;
  xi = x0 - SDIST * nx;
  yi = y0 - SDIST * ny;
  zi = z0 - SDIST * nz;
  dx = xo - xi;
  dy = yo - yi;
  dz = zo - zi;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(d)) DiagBreak();
  dx /= d;
  dy /= d;
  dz /= d;
  MRIsampleVolumeDerivativeScale(mri_int, x0, y0, z0, nx, ny, nz, &grad, GRAD_SIGMA);
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  if (grad < gcab->min_grad) grad = gcab->min_grad;
  MRIsampleVolume(mri_int, xo, yo, zo, &val2);
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
  MRIsampleVolume(mri_int, xi, yi, zi, &val1);
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  vno = gcabFindClosestIcoVertex(gcab, transform, x0, y0, z0, nx, ny, nz);
  if (vno < 0) return (0);

  // xb = (int)xn;
  // yb = (int)yn;
  // zb = (int)zn;
  // xd = xn - xb;
  // yd = yn - yb;
  // zd = zn - zb;
  // xu = xb + 1;
  // yu = yb + 1;
  // zu = zb + 1;

  // pI_total = pdelI = 0;
#if 1
  x = nint(xn);
  y = nint(yn);
  z = nint(zn);
  gcabs = &gcab->bs[x][y][z];
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri;

    char fname[STRLEN];

    mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
    sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
    MRIwrite(mri, fname);
    MRIfree(&mri);
    sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_grad_pdfs[vno], fname);
    sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
    sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
  }
  pdelI = gcabSamplePDF(gcab, gcabs, GCAB_GRAD, vno, grad);  // p(gradient)
  pdelI *= exp(-((gcab->max_grad + 1) - fabs(grad)));
  pI_out = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iout)

  // compute probability of interior
  xi0 = xi;
  yi0 = yi;
  zi0 = zi;
  pI_in = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
  //  pI_in = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
  for (d = 0.0, num = 0; d > -10; d -= 0.1)  // start 1/2mm inside border
  {
    xi = xi0 + d * dx;
    yi = yi0 + d * dy;
    zi = zi0 + d * dz;
    MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
    if (mag < 0)  // closer to the opposite side of the surface
      break;
    num++;
    MRIsampleVolume(mri_int, xi, yi, zi, &val1);
    if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
    if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
    p = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
    p = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1);  // p(Iin)
    pI_in = MIN(pI_in, p);
  }

  // build two tangent vectors for sampling exterior values
  e1x = -ny;
  e1y = -nz;
  e1z = nx;  // non-parallel vector
  CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z);
  CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z);
  mag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  if (FZERO(mag)) mag = 1e-5;
  e1x /= mag;
  e1y /= mag;
  e1z /= mag;
  mag = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
  if (FZERO(mag)) mag = 1e-5;
  e2x /= mag;
  e2y /= mag;
  e2z /= mag;

  xo0 = xo;
  yo0 = yo;
  zo0 = zo;
  for (d1 = -1; d1 <= 1.0; d1++)
    for (d2 = -1.0; d2 <= 1.0; d2++) {
      xo = xo0 + d1 * e1x + d2 * e2x;
      yo = yo0 + d1 * e1y + d2 * e2y;
      zo = zo0 + d1 * e1z + d2 * e2z;
      MRIsampleVolume(mri_int, xo, yo, zo, &val2);
      if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
      if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
      p = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iin)
      pI_out = MIN(pI_out, p);
    }

  p = pI_in;

#else
  for (x = xb; x <= xu; x++) {
    for (y = yb; y <= yu; y++) {
      for (z = zb; z <= zu; z++) {
        wt = (x == xb ? (1 - xd) : xd);
        wt *= (y == yb ? (1 - yd) : yd);
        wt *= (z == zb ? (1 - zd) : zd);

        gcabs = &gcab->bs[x][y][z];
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
          MRI *mri;

          char fname[STRLEN];

          mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
          sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
          MRIwrite(mri, fname);
          MRIfree(&mri);
          sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_grad_pdfs[vno], fname);
          sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
          sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
        }
        xi0 = xi;
        yi0 = yi;
        zi0 = zi;
        pI_in = 1;
        pI_out = pI = 0.0;
        // compute joint probability of all interior voxels with exterior one
        for (d = 0.0, num = 0; d > -10; d--)  // start 1/2mm inside border
        {
          xi = xi0 + d * dx;
          yi = yi0 + d * dy;
          zi = zi0 + d * dz;
          MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
          if (mag < 0)  // closer to the opposite side of the surface
            break;
          num++;
          MRIsampleVolume(mri_int, xi, yi, zi, &val1);
          if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
          if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
          gcabValsToBins(gcab, val1, val2, grad, &i1f, &i2f, &gbinf);
          i1b = (int)i1f;
          i2b = (int)i2f;
          gbinb = (int)gbinf;
          i1u = i1b + 1;
          i2u = i2b + 1;
          gbinu = gbinb + 1;
          i1d = i1f - i1b;
          i2d = i2f - i2b;
          gbind = gbinf - gbinb;

#if 0
          pI += gcabs->pdfs[vno][i1b][i2b] * wt * (1.0-i1d) * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1b][i2u] * wt * (1.0-i1d) * i2d;
          pI += gcabs->pdfs[vno][i1u][i2b] * wt * i1d       * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1u][i2u] * wt * i1d       * i2d;
#else
          pI += gcabs->pdfs[vno][nint(i1f)][nint(i2f)];
          pI_in = MIN(pI_in, gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]);
#endif
        }

        if (num > 0) pI_total += pI / num;
#if 0
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinb] * (1.0-gbind);
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinu] * gbind;
#else
        break;
#endif
      }
    }
  }
  p = pI_total * pdelI;
  p = pI_out * pI_in;
#endif

  return (p);
}

double GCABgetPout(GCAB *gcab,
                   MRI *mri_int,
                   MRI *mri_dist,
                   TRANSFORM *transform,
                   float x0,
                   float y0,
                   float z0,
                   float nx,
                   float ny,
                   float nz)
{
  float xn, yn, zn, xi, yi, zi, xo, yo, zo, xi0, yi0, zi0, d, dx, dy, dz, xo0, yo0, zo0;
  // float xd, yd, zd;
  int vno, oob, x, y, z, num;
  // int xb, yb, zb, xu, yu, zu;
  double grad, p, pdelI, mag, val1, val2, pI_in, pI_out, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2;
  // double pI_total;
  GCABS *gcabs;

  oob = GCABsourceVoxelToNode(gcab, mri_int, transform, x0, y0, z0, &xn, &yn, &zn);
  if (oob != NO_ERROR) return (0.0);
  xo = x0 + SDIST * nx;
  yo = y0 + SDIST * ny;
  zo = z0 + SDIST * nz;
  xi = x0 - SDIST * nx;
  yi = y0 - SDIST * ny;
  zi = z0 - SDIST * nz;
  dx = xo - xi;
  dy = yo - yi;
  dz = zo - zi;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(d)) DiagBreak();
  dx /= d;
  dy /= d;
  dz /= d;
  MRIsampleVolumeDerivativeScale(mri_int, x0, y0, z0, nx, ny, nz, &grad, GRAD_SIGMA);
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  if (grad < gcab->min_grad) grad = gcab->min_grad;
  MRIsampleVolume(mri_int, xo, yo, zo, &val2);
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
  MRIsampleVolume(mri_int, xi, yi, zi, &val1);
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  vno = gcabFindClosestIcoVertex(gcab, transform, x0, y0, z0, nx, ny, nz);
  if (vno < 0) return (0);

  // xb = (int)xn;
  // yb = (int)yn;
  // zb = (int)zn;
  // xd = xn - xb;
  // yd = yn - yb;
  // zd = zn - zb;
  // xu = xb + 1;
  // yu = yb + 1;
  // zu = zb + 1;

  // pI_total = pdelI = 0;
#if 1
  x = nint(xn);
  y = nint(yn);
  z = nint(zn);
  gcabs = &gcab->bs[x][y][z];
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri;

    char fname[STRLEN];

    mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
    sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
    MRIwrite(mri, fname);
    MRIfree(&mri);
    sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_grad_pdfs[vno], fname);
    sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
    sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
  }
  pdelI = gcabSamplePDF(gcab, gcabs, GCAB_GRAD, vno, grad);  // p(gradient)
  pdelI *= exp(-((gcab->max_grad + 1) - fabs(grad)));
  pI_out = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iout)

  // compute probability of interior
  xi0 = xi;
  yi0 = yi;
  zi0 = zi;
  pI_in = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
  //  pI_in = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
  for (d = 0.0, num = 0; d > -10; d -= 0.5)  // start 1/2mm inside border
  {
    xi = xi0 + d * dx;
    yi = yi0 + d * dy;
    zi = zi0 + d * dz;
    MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
    if (mag < 0)  // closer to the opposite side of the surface
      break;
    num++;
    MRIsampleVolume(mri_int, xi, yi, zi, &val1);
    if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
    if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
    p = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
    //    p = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
    pI_in = MIN(pI_in, p);
  }

  // build two tangent vectors for sampling exterior values
  e1x = -ny;
  e1y = -nz;
  e1z = nx;  // non-parallel vector
  CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z);
  CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z);
  mag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  if (FZERO(mag)) mag = 1e-5;
  e1x /= mag;
  e1y /= mag;
  e1z /= mag;
  mag = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
  if (FZERO(mag)) mag = 1e-5;
  e2x /= mag;
  e2y /= mag;
  e2z /= mag;

  xo0 = xo;
  yo0 = yo;
  zo0 = zo;
  for (d1 = -1; d1 <= 1.0; d1++)
    for (d2 = -1.0; d2 <= 1.0; d2++) {
      xo = xo0 + d1 * e1x + d2 * e2x;
      yo = yo0 + d1 * e1y + d2 * e2y;
      zo = zo0 + d1 * e1z + d2 * e2z;
      MRIsampleVolume(mri_int, xo, yo, zo, &val2);
      if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
      if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
      p = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iin)
      pI_out = MIN(pI_out, p);
    }

  p = pI_out;

#else
  for (x = xb; x <= xu; x++) {
    for (y = yb; y <= yu; y++) {
      for (z = zb; z <= zu; z++) {
        wt = (x == xb ? (1 - xd) : xd);
        wt *= (y == yb ? (1 - yd) : yd);
        wt *= (z == zb ? (1 - zd) : zd);

        gcabs = &gcab->bs[x][y][z];
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
          MRI *mri;

          char fname[STRLEN];

          mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
          sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
          MRIwrite(mri, fname);
          MRIfree(&mri);
          sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_grad_pdfs[vno], fname);
          sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
          sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
        }
        xi0 = xi;
        yi0 = yi;
        zi0 = zi;
        pI_in = 1;
        pI_out = pI = 0.0;
        // compute joint probability of all interior voxels with exterior one
        for (d = 0.0, num = 0; d > -10; d--)  // start 1/2mm inside border
        {
          xi = xi0 + d * dx;
          yi = yi0 + d * dy;
          zi = zi0 + d * dz;
          MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
          if (mag < 0)  // closer to the opposite side of the surface
            break;
          num++;
          MRIsampleVolume(mri_int, xi, yi, zi, &val1);
          if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
          if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
          gcabValsToBins(gcab, val1, val2, grad, &i1f, &i2f, &gbinf);
          i1b = (int)i1f;
          i2b = (int)i2f;
          gbinb = (int)gbinf;
          i1u = i1b + 1;
          i2u = i2b + 1;
          gbinu = gbinb + 1;
          i1d = i1f - i1b;
          i2d = i2f - i2b;
          gbind = gbinf - gbinb;

#if 0
          pI += gcabs->pdfs[vno][i1b][i2b] * wt * (1.0-i1d) * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1b][i2u] * wt * (1.0-i1d) * i2d;
          pI += gcabs->pdfs[vno][i1u][i2b] * wt * i1d       * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1u][i2u] * wt * i1d       * i2d;
#else
          pI += gcabs->pdfs[vno][nint(i1f)][nint(i2f)];
          pI_in = MIN(pI_in, gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]);
#endif
        }

        if (num > 0) pI_total += pI / num;
#if 0
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinb] * (1.0-gbind);
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinu] * gbind;
#else
        break;
#endif
      }
    }
  }
  p = pI_total * pdelI;
  p = pI_out * pI_in;
#endif

  return (p);
}

double GCABgetPgrad(GCAB *gcab,
                    MRI *mri_int,
                    MRI *mri_dist,
                    TRANSFORM *transform,
                    float x0,
                    float y0,
                    float z0,
                    float nx,
                    float ny,
                    float nz)
{
  float xn, yn, zn, xi, yi, zi, xo, yo, zo, xi0, yi0, zi0, d, dx, dy, dz, xo0, yo0, zo0;
  // float xd, yd, zd;
  int vno, oob, x, y, z, num;
  // int xb, yb, zb, xu, yu, zu;
  double grad, p, pdelI, mag, val1, val2, pI_in, pI_out, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2;
  // double pI_total;
  GCABS *gcabs;

  oob = GCABsourceVoxelToNode(gcab, mri_int, transform, x0, y0, z0, &xn, &yn, &zn);
  if (oob != NO_ERROR) return (0.0);
  xo = x0 + SDIST * nx;
  yo = y0 + SDIST * ny;
  zo = z0 + SDIST * nz;
  xi = x0 - SDIST * nx;
  yi = y0 - SDIST * ny;
  zi = z0 - SDIST * nz;
  dx = xo - xi;
  dy = yo - yi;
  dz = zo - zi;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(d)) DiagBreak();
  dx /= d;
  dy /= d;
  dz /= d;
  MRIsampleVolumeDerivativeScale(mri_int, x0, y0, z0, nx, ny, nz, &grad, GRAD_SIGMA);
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  if (grad < gcab->min_grad) grad = gcab->min_grad;
  MRIsampleVolume(mri_int, xo, yo, zo, &val2);
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
  MRIsampleVolume(mri_int, xi, yi, zi, &val1);
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  vno = gcabFindClosestIcoVertex(gcab, transform, x0, y0, z0, nx, ny, nz);
  if (vno < 0) return (0);

  // xb = (int)xn;
  // yb = (int)yn;
  // zb = (int)zn;
  // xd = xn - xb;
  // yd = yn - yb;
  // zd = zn - zb;
  // xu = xb + 1;
  // yu = yb + 1;
  // zu = zb + 1;

  // pI_total = pdelI = 0;
#if 1
  x = nint(xn);
  y = nint(yn);
  z = nint(zn);
  gcabs = &gcab->bs[x][y][z];
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRI *mri;

    char fname[STRLEN];

    mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
    sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
    MRIwrite(mri, fname);
    MRIfree(&mri);
    sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_grad_pdfs[vno], fname);
    sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
    sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
    HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
  }
  pdelI = gcabSamplePDF(gcab, gcabs, GCAB_GRAD, vno, grad);  // p(gradient)
  pdelI *= exp(-((gcab->max_grad + 1) - fabs(grad)));
  pI_out = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iout)

  // compute probability of interior
  xi0 = xi;
  yi0 = yi;
  zi0 = zi;
  pI_in = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
  //  pI_in = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
  for (d = 0.0, num = 0; d > -10; d -= 0.5)  // start 1/2mm inside border
  {
    xi = xi0 + d * dx;
    yi = yi0 + d * dy;
    zi = zi0 + d * dz;
    MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
    if (mag < 0)  // closer to the opposite side of the surface
      break;
    num++;
    MRIsampleVolume(mri_int, xi, yi, zi, &val1);
    if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
    if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
    p = GCAcomputeLabelLikelihood(gcab->gca, transform, mri_int, xi, yi, zi, gcab->target_label);
    //    p = gcabSamplePDF(gcab, gcabs, GCAB_IIN, vno, val1) ;  // p(Iin)
    pI_in = MIN(pI_in, p);
  }

  // build two tangent vectors for sampling exterior values
  e1x = -ny;
  e1y = -nz;
  e1z = nx;  // non-parallel vector
  CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z);
  CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z);
  mag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  if (FZERO(mag)) mag = 1e-5;
  e1x /= mag;
  e1y /= mag;
  e1z /= mag;
  mag = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
  if (FZERO(mag)) mag = 1e-5;
  e2x /= mag;
  e2y /= mag;
  e2z /= mag;

  xo0 = xo;
  yo0 = yo;
  zo0 = zo;
  for (d1 = -1; d1 <= 1.0; d1++)
    for (d2 = -1.0; d2 <= 1.0; d2++) {
      xo = xo0 + d1 * e1x + d2 * e2x;
      yo = yo0 + d1 * e1y + d2 * e2y;
      zo = zo0 + d1 * e1z + d2 * e2z;
      MRIsampleVolume(mri_int, xo, yo, zo, &val2);
      if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
      if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
      p = gcabSamplePDF(gcab, gcabs, GCAB_IOUT, vno, val2);  // p(Iin)
      pI_out = MIN(pI_out, p);
    }

  p = pdelI;

#else
  for (x = xb; x <= xu; x++) {
    for (y = yb; y <= yu; y++) {
      for (z = zb; z <= zu; z++) {
        wt = (x == xb ? (1 - xd) : xd);
        wt *= (y == yb ? (1 - yd) : yd);
        wt *= (z == zb ? (1 - zd) : zd);

        gcabs = &gcab->bs[x][y][z];
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
          MRI *mri;

          char fname[STRLEN];

          mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
          sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
          MRIwrite(mri, fname);
          MRIfree(&mri);
          sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_grad_pdfs[vno], fname);
          sprintf(fname, "ihisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
          sprintf(fname, "ohisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
          HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
        }
        xi0 = xi;
        yi0 = yi;
        zi0 = zi;
        pI_in = 1;
        pI_out = pI = 0.0;
        // compute joint probability of all interior voxels with exterior one
        for (d = 0.0, num = 0; d > -10; d--)  // start 1/2mm inside border
        {
          xi = xi0 + d * dx;
          yi = yi0 + d * dy;
          zi = zi0 + d * dz;
          MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
          if (mag < 0)  // closer to the opposite side of the surface
            break;
          num++;
          MRIsampleVolume(mri_int, xi, yi, zi, &val1);
          if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
          if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
          gcabValsToBins(gcab, val1, val2, grad, &i1f, &i2f, &gbinf);
          i1b = (int)i1f;
          i2b = (int)i2f;
          gbinb = (int)gbinf;
          i1u = i1b + 1;
          i2u = i2b + 1;
          gbinu = gbinb + 1;
          i1d = i1f - i1b;
          i2d = i2f - i2b;
          gbind = gbinf - gbinb;

#if 0
          pI += gcabs->pdfs[vno][i1b][i2b] * wt * (1.0-i1d) * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1b][i2u] * wt * (1.0-i1d) * i2d;
          pI += gcabs->pdfs[vno][i1u][i2b] * wt * i1d       * (1.0-i2d);
          pI += gcabs->pdfs[vno][i1u][i2u] * wt * i1d       * i2d;
#else
          pI += gcabs->pdfs[vno][nint(i1f)][nint(i2f)];
          pI_in = MIN(pI_in, gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]);
#endif
        }

        if (num > 0) pI_total += pI / num;
#if 0
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinb] * (1.0-gbind);
        pdelI += wt * gcabs->h_grad_pdfs[vno]->counts[gbinu] * gbind;
#else
        break;
#endif
      }
    }
  }
  p = pI_total * pdelI;
  p = pI_out * pI_in;
#endif

  return (p);
}

int GCABtrain(GCAB *gcab, MRI *mri_int, MRI *mri_seg, TRANSFORM *transform, int target_label)
{
  int x, y, z, oob, border;
  float xn, yn, zn;
  MRI *mri_border, *mri_dist, *mri_tmp;

  mri_tmp = MRIclone(mri_seg, NULL);
  MRIcopyLabel(mri_seg, mri_tmp, target_label);
  mri_dist = MRIdistanceTransform(mri_tmp, NULL, target_label, 10, DTRANS_MODE_SIGNED, NULL);
  MRIfree(&mri_tmp);

  gcab->ntraining++;  // # of volumes that went into this atlas
  mri_border = MRImarkLabelBorderVoxels(mri_seg, NULL, target_label, 1, 1);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_border, "b.mgz");
    MRIwrite(mri_dist, "d.mgz");
  }
  for (x = 0; x < mri_int->width; x++) {
    for (y = 0; y < mri_int->height; y++) {
      for (z = 0; z < mri_int->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        border = (int)MRIgetVoxVal(mri_border, x, y, z, 0);
        if (border == 0) continue;
        oob = GCABsourceVoxelToNode(gcab, mri_int, transform, x, y, z, &xn, &yn, &zn);
        if (oob != NO_ERROR)  // not in bounding box
          continue;
        DiagBreak();
        gcabUpdateDistributions(gcab, mri_int, mri_seg, mri_dist, transform, x, y, z, xn, yn, zn, target_label);
      }
    }
  }

  MRIfree(&mri_border);
  MRIfree(&mri_dist);
  return (NO_ERROR);
}

int GCABsourceVoxelToNode(
    GCAB *gcab, MRI *mri, TRANSFORM *transform, float xv, float yv, float zv, float *pxn, float *pyn, float *pzn)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt;
  LTA *lta;

  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX)  // from src to talairach volume
    {
      lta = (LTA *)transform->xform;
      // get the talairach position
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, &xrt, &yrt, &zrt);
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
      xt = xrt;
      yt = yrt;
      zt = zrt;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCAsourceVoxelToNode: needs vox-to-vox transform");
  }
  else {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
  }
  if (Ggca_x == xv && Ggca_y == yv && Ggca_z == zv && DIAG_VERBOSE_ON)
    fprintf(stdout, "source (%2.1f, %2.1f, %2.1f) to talposition (%.2f, %.2f, %.2f)\n", xv, yv, zv, xt, yt, zt);
  fflush(stdout);
  xt -= gcab->bounding_box.x;
  yt -= gcab->bounding_box.y;
  zt -= gcab->bounding_box.z;

  // get the position in node from the talairach position
  return GCABvoxelToNode(gcab, xt, yt, zt, pxn, pyn, pzn);
}
////////////////////////////////////////////////////////////////
// transform from template -> node
////////////////////////////////////////////////////////////////
int GCABvoxelToNode(GCAB *gcab, float xv, float yv, float zv, float *pxn, float *pyn, float *pzn)
{
  float xn, yn, zn;

  xn = *pxn = xv / gcab->spacing;
  yn = *pyn = yv / gcab->spacing;
  zn = *pzn = zv / gcab->spacing;
  if (xn < 0 || yn < 0 || zn < 0 || xn >= gcab->width || yn >= gcab->height || zn >= gcab->depth)
    return (ERROR_OUT_OF_BOUNDS);
  return NO_ERROR;
}
int GCABwrite(GCAB *gcab, char *fname)
{
  FILE *fp;
  int x, y, z, vno, i1, i2;
  GCABS *gcabs;
  long where;

  strcpy(gcab->fname, fname);

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "GCABwrite(%s) fopen failed", fname));

  fwriteFloat(GCAB_VERSION, fp);
  fwrite(gcab->gca_fname, sizeof(char), STRLEN - 1, fp);
  where = ftell(fp);
  if (where == -1L){
    ErrorPrintf(ERROR_BADFILE, "GCABwrite(%s) ftell returned invalid position indicator", fname);
  }
  fwriteInt(gcab->target_label, fp);
  fwriteInt(gcab->spacing, fp);
  fwriteInt(gcab->nvertices, fp);
  fwriteInt(gcab->ico_order, fp);
  fwriteInt(gcab->nintensity_bins, fp);
  fwriteInt(gcab->ngrad_bins, fp);
  fwriteInt(gcab->width, fp);
  fwriteInt(gcab->height, fp);
  fwriteInt(gcab->depth, fp);

  fwriteFloat(gcab->min_intensity, fp);
  fwriteFloat(gcab->max_intensity, fp);
  fwriteFloat(gcab->min_grad, fp);
  fwriteFloat(gcab->max_grad, fp);
  fwriteFloat(gcab->iscale, fp);
  fwriteFloat(gcab->gscale, fp);
  fwriteInt(gcab->ntraining, fp);
  fwriteInt(gcab->bounding_box.x, fp);
  fwriteInt(gcab->bounding_box.y, fp);
  fwriteInt(gcab->bounding_box.z, fp);

  fwriteInt(gcab->bounding_box.dx, fp);
  fwriteInt(gcab->bounding_box.dy, fp);
  fwriteInt(gcab->bounding_box.dz, fp);

  fwrite(gcab->fname, sizeof(char), STRLEN - 1, fp);

  for (x = 0; x < gcab->width; x++) {
    for (y = 0; y < gcab->height; y++) {
      for (z = 0; z < gcab->depth; z++) {
        gcabs = &gcab->bs[x][y][z];
        if (Gdiag & DIAG_WRITE /* && DIAG_VERBOSE_ON*/) {
          MRI *mri;

          char fname[STRLEN];

          mri = gcabWritePDFToMRI(gcab, NULL, x, y, z);
          sprintf(fname, "Ihisto_%d_%d_%d.mgz", x, y, z);
          MRIwrite(mri, fname);
          MRIfree(&mri);
        }
        for (vno = 0; vno < gcab->nvertices; vno++) {
          fwriteFloat(gcabs->ntraining[vno], fp);
          HISTOwriteInto(gcabs->h_grad_pdfs[vno], fp);
          HISTOwriteInto(gcabs->h_Iin_pdfs[vno], fp);
          HISTOwriteInto(gcabs->h_Iout_pdfs[vno], fp);
          if (Gdiag & DIAG_WRITE /*&& DIAG_VERBOSE_ON*/) {
            char fname[STRLEN];
            sprintf(fname, "ghisto_%d_%d_%d_vno%d.plt", x, y, z, vno);
            HISTOplot(gcabs->h_grad_pdfs[vno], fname);
            HISTOplot(gcabs->h_Iin_pdfs[vno], fname);
            HISTOplot(gcabs->h_Iout_pdfs[vno], fname);
          }
          for (i1 = 0; i1 < gcab->nintensity_bins; i1++)
            for (i2 = 0; i2 < gcab->nintensity_bins; i2++) fwriteFloat(gcabs->pdfs[vno][i1][i2], fp);
        }
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}

GCAB *GCABread(char *fname, GCA *gca)
{
  FILE *fp;
  int x, y, z, vno, i1, i2, spacing, ico_order, nibins, ngbins, label;
  // int nvertices;
  float version;
  GCABS *gcabs;
  GCAB *gcab;
  char gca_fname[STRLEN];
  long where;

  fp = fopen(fname, "rb");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_BADFILE, "GCABread(%s) fopen failed", fname));

  version = freadFloat(fp);
  if (version != GCAB_VERSION)
    ErrorReturn(NULL,
                (ERROR_BADFILE,
                 "GCABread(%s): version %2.1f does not match expected version %2.1f",
                 fname,
                 version,
                 GCAB_VERSION));
  if(!fread(gca_fname, sizeof(char), STRLEN - 1, fp)){
    ErrorPrintf(ERROR_BADFILE, "GCABRead(%s) failed fread", fname);
  }
  if (gca == NULL) {
    gca = GCAread(gca_fname);
    if (gca == NULL) ErrorReturn(NULL, (ERROR_BADFILE, "GCABread(%s): could not read gca from %s", fname, gca_fname));
  }
  where = ftell(fp);
  if(where == -1L){
    ErrorPrintf(ERROR_BADFILE, "GCABread(%s) ftell returned invalid position indicator", fname);
  }
  label = freadInt(fp);
  spacing = freadInt(fp);
  // nvertices = 
  freadInt(fp);
  ico_order = freadInt(fp);
  nibins = freadInt(fp);
  ngbins = freadInt(fp);
  gcab = GCABalloc(gca, spacing, ico_order, nibins, ngbins, label);

  freadInt(fp);
  freadInt(fp);
  freadInt(fp);  // width, height and depth
  gcab->min_intensity = freadFloat(fp);
  gcab->max_intensity = freadFloat(fp);
  gcab->min_grad = freadFloat(fp);
  gcab->max_grad = freadFloat(fp);
  gcab->iscale = freadFloat(fp);
  gcab->gscale = freadFloat(fp);
  gcab->ntraining = freadInt(fp);
  gcab->bounding_box.x = freadInt(fp);
  gcab->bounding_box.y = freadInt(fp);
  gcab->bounding_box.z = freadInt(fp);

  gcab->bounding_box.dx = freadInt(fp);
  gcab->bounding_box.dy = freadInt(fp);
  gcab->bounding_box.dz = freadInt(fp);

  if(!fread(gcab->fname, sizeof(char), STRLEN - 1, fp)){
    ErrorPrintf(ERROR_BADFILE, "GCABread(%s) failed fread", fname);
  }
  strcpy(gcab->fname, fname);

  for (x = 0; x < gcab->width; x++) {
    for (y = 0; y < gcab->height; y++) {
      for (z = 0; z < gcab->depth; z++) {
        gcabs = &gcab->bs[x][y][z];
        for (vno = 0; vno < gcab->nvertices; vno++) {
          gcabs->ntraining[vno] = freadFloat(fp);
          HISTOfree(&gcabs->h_grad_pdfs[vno]);  // alloced in GCABalloc
          HISTOfree(&gcabs->h_Iin_pdfs[vno]);   // alloced in GCABalloc
          HISTOfree(&gcabs->h_Iout_pdfs[vno]);  // alloced in GCABalloc
          gcabs->h_grad_pdfs[vno] = HISTOreadFrom(fp);
          gcabs->h_Iin_pdfs[vno] = HISTOreadFrom(fp);
          gcabs->h_Iout_pdfs[vno] = HISTOreadFrom(fp);
          for (i1 = 0; i1 < gcab->nintensity_bins; i1++)
            for (i2 = 0; i2 < gcab->nintensity_bins; i2++) gcabs->pdfs[vno][i1][i2] = freadFloat(fp);
        }
      }
    }
  }

  fclose(fp);
  return (gcab);
}
static int gcabUpdateDistributions(GCAB *gcab,
                                   MRI *mri_int,
                                   MRI *mri_seg,
                                   MRI *mri_dist,
                                   TRANSFORM *transform,
                                   int x,
                                   int y,
                                   int z,
                                   float xn,
                                   float yn,
                                   float zn,
                                   int target_label)
{
  float nx, ny, nz, xi, yi, zi, xo, yo, zo, dx, dy, dz, d, xi0, yi0, zi0, xo0, yo0, zo0, x0, y0, z0;
  // float xd, yd, zd;
  int vno, label;
  // int xb, yb, zb, xu, yu, zu;
  double grad, mag, val1, val2, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2;

  vno = gcabComputeBorderNormal(gcab, mri_seg, transform, x, y, z, &nx, &ny, &nz, target_label);
  if (vno < 0) return (NO_ERROR);  // no well-defined normal - skip it

  label = (int)MRIgetVoxVal(mri_seg, x, y, z, 0);
  if (label == target_label)  // current voxel is on inside border
  {
    xi = x;
    yi = y;
    zi = z;
    xo = x + nx;
    yo = y + ny;
    zo = z + nz;
  }
  else  // current voxel is on outside border
  {
    xo = x;
    yo = y;
    zo = z;
    xi = x - nx;
    yi = y - ny;
    zi = z - nz;
  }
  dx = xo - xi;
  dy = yo - yi;
  dz = zo - zi;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(d)) DiagBreak();
  dx /= d;
  dy /= d;
  dz /= d;

  // sample values and gradient across border
  x0 = (xi + xo) / 2;
  y0 = (yi + yo) / 2;
  z0 = (zi + zo) / 2;
  xi = x0 - SDIST * dx;
  yi = y0 - SDIST * dy;
  zi = z0 - SDIST * dz;
  xo = x0 + SDIST * dx;
  yo = y0 + SDIST * dy;
  zo = z0 + SDIST * dz;
  MRIsampleVolumeDerivativeScale(mri_int, x0, y0, z0, nx, ny, nz, &grad, GRAD_SIGMA);
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  if (grad < gcab->min_grad) grad = gcab->min_grad;

  MRIsampleVolume(mri_int, xi, yi, zi, &val1);
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  MRIsampleVolume(mri_int, xo, yo, zo, &val2);
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;

  gcabUpdateNode(gcab, nint(xn), nint(yn), nint(zn), val1, val2, grad, vno, 1);

  // build two tangent vectors for sampling exterior values
  e1x = -ny;
  e1y = -nz;
  e1z = nx;  // non-parallel vector
  CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z);
  CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z);
  mag = sqrt(e1x * e1x + e1y * e1y + e1z * e1z);
  if (FZERO(mag)) mag = 1e-5;
  e1x /= mag;
  e1y /= mag;
  e1z /= mag;
  mag = sqrt(e2x * e2x + e2y * e2y + e2z * e2z);
  if (FZERO(mag)) mag = 1e-5;
  e2x /= mag;
  e2y /= mag;
  e2z /= mag;

  xo0 = xo;
  yo0 = yo;
  zo0 = zo;
  for (d1 = -1; d1 <= 1.0; d1++)
    for (d2 = -1.0; d2 <= 1.0; d2++) {
      xo = xo0 + d1 * e1x + d2 * e2x;
      yo = yo0 + d1 * e1y + d2 * e2y;
      zo = zo0 + d1 * e1z + d2 * e2z;
      MRIsampleVolume(mri_int, xo, yo, zo, &val2);
      if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
      if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
      gcabUpdateNodeIout(gcab, nint(xn), nint(yn), nint(zn), val2, vno, 1);
    }

  // update all  nearest bins
  // xb = (int)xn;
  // yb = (int)yn;
  // zb = (int)zn;
  // xd = xn - xb;
  // yd = yn - yb;
  // zd = zn - zb;
  // xu = xb + 1;
  // yu = yb + 1;
  // zu = zb + 1;
  xi0 = xi;
  yi0 = yi;
  zi0 = zi;
  for (d = 0; d > -10; d -= 0.5)  // start 1/2mm inside border
  {
    xi = xi0 + d * dx;
    yi = yi0 + d * dy;
    zi = zi0 + d * dz;
    MRIsampleVolumeDerivativeScale(mri_dist, xi, yi, zi, dx, dy, dz, &mag, 0);
    if (mag < 0)  // closer to the opposite side of the surface
      break;
    MRIsampleVolume(mri_int, xi, yi, zi, &val1);
    if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
    if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
#if 0
    gcabUpdateNode(gcab, xb, yb, zb, val1, val2, grad,vno,(1-xd)*(1-yd)*(1-zd));
    gcabUpdateNode(gcab, xb, yb, zu, val1, val2, grad,vno,(1-xd)*(1-yd)*zd);
    gcabUpdateNode(gcab, xb, yu, zb, val1, val2, grad,vno,(1-xd)*yd    *(1-zd));
    gcabUpdateNode(gcab, xb, yu, zu, val1, val2, grad,vno,(1-xd)*yd    *zd);
    gcabUpdateNode(gcab, xu, yb, zb, val1, val2, grad,vno,xd    *(1-yd)*(1-zd));
    gcabUpdateNode(gcab, xu, yb, zu, val1, val2, grad,vno,xd    *(1-yd)*zd);
    gcabUpdateNode(gcab, xu, yu, zb, val1, val2, grad,vno,xd    *yd    *(1-zd));
    gcabUpdateNode(gcab, xu, yu, zu, val1, val2, grad,vno,xd    *yd    *zd);
#else
    gcabUpdateNodeIin(gcab, nint(xn), nint(yn), nint(zn), val1, vno, 1);
#endif
  }

  return (NO_ERROR);
}
static int gcabUpdateNodeIout(GCAB *gcab, int xn, int yn, int zn, double val2, int vno, float wt)
{
  float i1f, i2f, gbinf;
  // float i2d;
  // int i2b, i2u;
  GCABS *gcabs;

  if (xn < 0 || xn >= gcab->width || yn < 0 || yn >= gcab->height || zn < 0 || zn >= gcab->depth) return (NO_ERROR);

  gcabValsToBins(gcab, 0, val2, 0, &i1f, &i2f, &gbinf);

  // do linear interpolation
  // i2b = (int)i2f;
  // i2u = i2b + 1;
  // i2d = i2f - i2b;

  gcabs = &gcab->bs[xn][yn][zn];

  if (xn == Gx && yn == Gy && zn == Gz) {
    if (vno == 5) DiagBreak();
    DiagBreak();
  }
#if 0
  gcabs->pdfs[vno][i1b][i2b] += wt * (1.0-i1d) * (1.0-i2d);
  gcabs->pdfs[vno][i1b][i2u] += wt * (1.0-i1d) * i2d;
  gcabs->pdfs[vno][i1u][i2b] += wt * i1d       * (1.0-i2d);
  gcabs->pdfs[vno][i1u][i2u] += wt * i1d       * i2d;

  gcabs->h_grad_pdfs[vno]->counts[gbinb] += wt * (1-gbind) ;
  gcabs->h_grad_pdfs[vno]->counts[gbinu] += wt * (gbind) ;
#else
  gcabs->h_Iout_pdfs[vno]->counts[nint(i2f)]++;
#endif

  gcabs->ntraining[vno] += wt;
  return (NO_ERROR);
}

static int gcabUpdateNodeIin(GCAB *gcab, int xn, int yn, int zn, double val1, int vno, float wt)
{
  float i1f, i2f, gbinf;
  // float i1d;
  // int i1b, i1u;
  GCABS *gcabs;

  if (xn < 0 || xn >= gcab->width || yn < 0 || yn >= gcab->height || zn < 0 || zn >= gcab->depth) return (NO_ERROR);

  gcabValsToBins(gcab, val1, 0, 0, &i1f, &i2f, &gbinf);

  // do linear interpolation
  // i1b = (int)i2f;
  // i1u = i1b + 1;
  // i1d = i1f - i1b;

  gcabs = &gcab->bs[xn][yn][zn];

  if (xn == Gx && yn == Gy && zn == Gz) {
    if (vno == 5) DiagBreak();
    DiagBreak();
  }
#if 0
  gcabs->pdfs[vno][i1b][i2b] += wt * (1.0-i1d) * (1.0-i2d);
  gcabs->pdfs[vno][i1b][i2u] += wt * (1.0-i1d) * i2d;
  gcabs->pdfs[vno][i1u][i2b] += wt * i1d       * (1.0-i2d);
  gcabs->pdfs[vno][i1u][i2u] += wt * i1d       * i2d;

  gcabs->h_grad_pdfs[vno]->counts[gbinb] += wt * (1-gbind) ;
  gcabs->h_grad_pdfs[vno]->counts[gbinu] += wt * (gbind) ;
#else
  gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]++;
#endif

  gcabs->ntraining[vno] += wt;
  return (NO_ERROR);
}

static int gcabUpdateNode(GCAB *gcab, int xn, int yn, int zn, double val1, double val2, float grad, int vno, float wt)
{
  float i1f, i2f,  gbinf;
  // float i1d, i2d, gbind;
  // int i1b, i2b;
  // int i1u, i2u, gbinu, gbinb;
  GCABS *gcabs;

  if (xn < 0 || xn >= gcab->width || yn < 0 || yn >= gcab->height || zn < 0 || zn >= gcab->depth) return (NO_ERROR);

  gcabValsToBins(gcab, val1, val2, grad, &i1f, &i2f, &gbinf);

  // do linear interpolation
  // gbinb = (int)gbinf;
  // gbinu = gbinb + 1;
  // i1b = (int)i1f;
  // i1u = i1b + 1;
  // i2b = (int)i2f;
  // i2u = i2b + 1;
  // i1d = i1f - i1b;
  // i2d = i2f - i2b;
  // gbind = gbinf - gbinb;

  gcabs = &gcab->bs[xn][yn][zn];

  if (xn == Gx && yn == Gy && zn == Gz) {
    if (vno == 5) DiagBreak();
    DiagBreak();
  }
#if 0
  gcabs->pdfs[vno][i1b][i2b] += wt * (1.0-i1d) * (1.0-i2d);
  gcabs->pdfs[vno][i1b][i2u] += wt * (1.0-i1d) * i2d;
  gcabs->pdfs[vno][i1u][i2b] += wt * i1d       * (1.0-i2d);
  gcabs->pdfs[vno][i1u][i2u] += wt * i1d       * i2d;

  gcabs->h_grad_pdfs[vno]->counts[gbinb] += wt * (1-gbind) ;
  gcabs->h_grad_pdfs[vno]->counts[gbinu] += wt * (gbind) ;
#else
  if (vno == Gdiag_no && xn == Gx && yn == Gy && zn == Gz) {
    printf("grad (%d, %d, %d) = %2.2f\n", xn, yn, zn, grad);
    DiagBreak();
  }
  gcabs->h_grad_pdfs[vno]->counts[nint(gbinf)]++;
  gcabs->h_Iin_pdfs[vno]->counts[nint(i1f)]++;
  gcabs->h_Iout_pdfs[vno]->counts[nint(i2f)]++;
  gcabs->pdfs[vno][nint(i1f)][nint(i1f)]++;
#endif

  gcabs->ntraining[vno] += wt;
  return (NO_ERROR);
}

static int gcabComputeBorderNormal(GCAB *gcab,
                                   MRI *mri_seg,
                                   TRANSFORM *transform,
                                   int x0,
                                   int y0,
                                   int z0,
                                   float *pnx,
                                   float *pny,
                                   float *pnz,
                                   int target_label)
{
  int olabel, x1, y1, z1, xk, yk, zk, label, num, max_n;
  float nx, ny, nz, mag;

  olabel = (int)MRIgetVoxVal(mri_seg, x0, y0, z0, 0);

  for (nx = ny = nz = 0.0, num = 0, xk = -1; xk <= 1; xk++) {
    x1 = mri_seg->xi[x0 + xk];
    for (yk = -1; yk <= 1; yk++) {
      y1 = mri_seg->yi[y0 + yk];
      for (zk = -1; zk <= 1; zk++) {
        z1 = mri_seg->zi[z0 + zk];
        if (fabs(xk) + fabs(yk) + fabs(zk) > 1) continue;  // only 6-connected
        label = (int)MRIgetVoxVal(mri_seg, x1, y1, z1, 0);
        if ((label == target_label && olabel != target_label) || (label != target_label && olabel == target_label)) {
          nx += xk;
          ny += yk;
          nz += zk;
          num++;
        }
      }
    }
  }

  if (!FZERO(num)) {
    nx /= num;
    ny /= num;
    nz /= num;
  }
  else
    DiagBreak();

  if (olabel != target_label)  // make normal point outwards
  {
    nx *= -1;
    ny *= -1;
    nz *= -1;
  }

  mag = sqrt(nx * nx + ny * ny + nz * nz);
  if (!FZERO(mag)) {
    nx /= mag;
    ny /= mag;
    nz /= mag;
  }
  *pnx = nx;
  *pny = ny;
  *pnz = nz;  // return in image coords

  max_n = gcabFindClosestIcoVertex(gcab, transform, x0, y0, z0, nx, ny, nz);

  return (max_n);
}

static int gcabFindClosestIcoVertex(
    GCAB *gcab, TRANSFORM *transform, float x0, float y0, float z0, float nx, float ny, float nz)
{
  int n, max_n;
  float xa0, ya0, za0, xa1, ya1, za1, mag, dot, max_dot;

  // convert it to atlas coords
  TransformSampleReal(transform, x0, y0, z0, &xa0, &ya0, &za0);
  TransformSampleReal(transform, x0 + 2 * nx, y0 + 2 * ny, z0 + 2 * nz, &xa1, &ya1, &za1);
  nx = xa1 - xa0;
  ny = ya1 - ya0;
  nz = za1 - za0;
  mag = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= mag;
  ny /= mag;
  nz /= mag;

  max_dot = -1;
  max_n = -1;
  for (n = 0; n < gcab->nvertices; n++) {
    dot = nx * gcab->ico->vertices[n].x + ny * gcab->ico->vertices[n].y + nz * gcab->ico->vertices[n].z;
    if (dot > max_dot) {
      max_dot = dot;
      max_n = n;
    }
  }
  return (max_n);
}

static int gcabValsToBins(GCAB *gcab, double val1, double val2, float grad, float *pi1f, float *pi2f, float *pgf)
{
  if (val1 < gcab->min_intensity) val1 = gcab->min_intensity;
  if (val1 > gcab->max_intensity) val1 = gcab->max_intensity;
  if (val2 < gcab->min_intensity) val2 = gcab->min_intensity;
  if (val2 > gcab->max_intensity) val2 = gcab->max_intensity;
  if (grad < gcab->min_grad) grad = gcab->min_grad;
  if (grad > gcab->max_grad) grad = gcab->max_grad;
  *pi1f = gcab->iscale * (val1 - gcab->min_intensity);
  *pi2f = gcab->iscale * (val2 - gcab->min_intensity);
  *pgf = gcab->gscale * (grad - gcab->min_grad);
  return (NO_ERROR);
}

static MRI *gcabWritePDFToMRI(GCAB *gcab, MRI *mri, int x, int y, int z)
{
  int vno, i1, i2;
  GCABS *gcabs;

  if (mri == NULL) mri = MRIalloc(gcab->nvertices, gcab->nintensity_bins, gcab->nintensity_bins, MRI_FLOAT);

  gcabs = &gcab->bs[x][y][z];
  for (vno = 0; vno < gcab->nvertices; vno++)
    for (i1 = 0; i1 < gcab->nintensity_bins; i1++)
      for (i2 = 0; i2 < gcab->nintensity_bins; i2++) MRIsetVoxVal(mri, vno, i1, i2, 0, gcabs->pdfs[vno][i1][i2]);
  return (mri);
}
static MRI *gcabWriteMRIToPDF(GCAB *gcab, MRI *mri, int x, int y, int z)
{
  int vno, i1, i2;
  GCABS *gcabs;

  gcabs = &gcab->bs[x][y][z];
  for (vno = 0; vno < gcab->nvertices; vno++)
    for (i1 = 0; i1 < gcab->nintensity_bins; i1++)
      for (i2 = 0; i2 < gcab->nintensity_bins; i2++) gcabs->pdfs[vno][i1][i2] = MRIgetVoxVal(mri, vno, i1, i2, 0);
  return (NO_ERROR);
}
int GCABsmoothPDFs(GCAB *gcab, float sigma)
{
  int x, y, z;
  MRI *mri_kernel;

  mri_kernel = MRIgaussian1d(sigma, -1);

  for (x = 0; x < gcab->width; x++) {
    for (y = 0; y < gcab->height; y++) {
      for (z = 0; z < gcab->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        gcabSmoothPDF(gcab, x, y, z, mri_kernel);
      }
    }
  }
  MRIfree(&mri_kernel);

  return (NO_ERROR);
}
static int gcabSmoothPDF(GCAB *gcab, int x, int y, int z, MRI *mri_kernel)
{
  MRI *mri_pdf, *mri_tmp;
  // int width, height, depth
  int klen;
  float *kernel;

  kernel = &MRIFvox(mri_kernel, 0, 0, 0);
  klen = mri_kernel->width;

  mri_pdf = gcabWritePDFToMRI(gcab, NULL, x, y, z);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_pdf, "pdf_before.mgz");
  // width = mri_pdf->width;
  // height = mri_pdf->height;
  // depth = mri_pdf->depth;
  mri_tmp = MRIclone(mri_pdf, NULL);
  MRIconvolve1d(mri_pdf, mri_tmp, kernel, klen, MRI_HEIGHT, 0, 0);
  MRIconvolve1d(mri_tmp, mri_pdf, kernel, klen, MRI_DEPTH, 0, 0);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_pdf, "pdf_after.mgz");
  gcabWriteMRIToPDF(gcab, mri_pdf, x, y, z);

  MRIfree(&mri_pdf);
  return (NO_ERROR);
}

static double gcabSamplePDF(GCAB *gcab, GCABS *gcabs, int which, int vno, double val)
{
  double p, p1, p2;
  float index_f = 0, index_d, dummy;
  int index_b, index_u;
  HISTOGRAM *h;

  h = NULL;
  switch (which) {
    default:
    case GCAB_IIN:
      gcabValsToBins(gcab, val, 0, 0, &index_f, &dummy, &dummy);
      break;
    case GCAB_IOUT:
      gcabValsToBins(gcab, 0, val, 0, &dummy, &index_f, &dummy);
      break;
    case GCAB_GRAD:
      gcabValsToBins(gcab, 0, 0, val, &dummy, &dummy, &index_f);
      break;
  }
  index_b = (int)index_f;
  index_d = index_f - index_b;
  index_u = index_b + 1;
  switch (which) {
    case GCAB_IIN:
      h = gcabs->h_Iin_pdfs[vno];
      break;
    case GCAB_IOUT:
      h = gcabs->h_Iout_pdfs[vno];
      break;
    case GCAB_GRAD:
      h = gcabs->h_grad_pdfs[vno];
      break;
  }
  p1 = h->counts[index_b];
  p2 = h->counts[index_u];
  p = (1 - index_d) * p1 + index_d * p2;
  if (DZERO(p)) p = 1e-10;
  return (p);
}
