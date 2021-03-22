/**
 * @brief utilities for MRI data structure histograms
 *
 */
/*
 * Original Author: Bruce Fischl (1/8/97)
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

#include "romp_support.h"

#include "box.h"
#include "diag.h"
#include "error.h"
#include "filter.h"
#include "label.h"
#include "macros.h"
#include "minc.h"
#include "mri.h"
#include "mrinorm.h"
#include "proto.h"
#include "region.h"

/*-----------------------------------------------------
  MACROS AND CONSTANTS
  -------------------------------------------------------*/

#define DEBUG_POINT(x, y, z) (((x) == 15) && ((y) == 6) && ((z) == 15))

/*-----------------------------------------------------
  STATIC DATA
  -------------------------------------------------------*/

/*-----------------------------------------------------
  STATIC PROTOTYPES
  -------------------------------------------------------*/
static HISTOGRAM *mriHistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region);
static HISTOGRAM *mriHistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region);
static HISTOGRAM *mriHistogramLabel(MRI *mri, int nbins, HISTOGRAM *histo, LABEL *label);
/*-----------------------------------------------------
  GLOBAL FUNCTIONS
  -------------------------------------------------------*/
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIapplyHistogram(MRI *mri_src, MRI *mri_dst, HISTOGRAM *histo)
{
  int width, height, depth, x, y, z;
  BUFTYPE *psrc, *pdst, sval, dval;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        sval = *psrc++;
        if (sval >= histo->nbins) sval = histo->nbins - 1;
        dval = histo->counts[sval];
        *pdst++ = dval;
      }
    }
  }
  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIhistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region)
{
  int overlap;
  float fmin, fmax, bin_size;
  // BUFTYPE bmin, bmax;
  static MRI *mri_prev = NULL;
  static HISTOGRAM *h_prev;
  static MRI_REGION reg_prev;

#if 0
  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,"MRIhistogramRegion: must by type UCHAR"));
#endif

  MRIvalRangeRegion(mri, &fmin, &fmax, region);
  // bmin = (BUFTYPE)fmin;
  // bmax = (BUFTYPE)fmax;
  if (nbins <= 0 && FEQUAL(fmin, fmax))
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIhistogramRegion: input mri %s is blank", mri->fname));
  if (nbins <= 0) nbins = nint(fmax - fmin + 1.0);
  if (nbins <= 0) nbins = 255;

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    // histo->nbins = nbins;
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }
  HISTOclear(histo, histo);
  bin_size = (fmax - fmin + 1) / (float)nbins;
  histo->bin_size = bin_size;

  if (!mri_prev) /* first invocation, initialize state machine */
  {
    h_prev = HISTOcopy(histo, NULL);
    HISTOclear(h_prev, h_prev);
    REGIONclear(&reg_prev);
  }

  if (h_prev->nbins != histo->nbins) {
    HISTOfree(&h_prev);
    h_prev = HISTOcopy(histo, NULL);
  }
  /*
    note that the overlap only works with subsequent windows advancing only
    in the x direction.
  */
  /* check to see if regions overlap */
  overlap = ((mri == mri_prev) && ((region->x - region->dx) > reg_prev.x) && (region->y == reg_prev.y) &&
             (region->z == reg_prev.z));

  if (0 && overlap) {  // take advantage of overlapping regions
    MRI_REGION reg_left, reg_right;
    HISTOGRAM *histo_left, *histo_right;

    REGIONcopy(&reg_prev, &reg_left);
    reg_left.dx = region->x - reg_left.x;
    histo_left = mriHistogramRegion(mri, 0, NULL, &reg_left);
    REGIONcopy(&reg_prev, &reg_right);
    reg_right.x = reg_prev.x + reg_prev.dx;
    reg_right.dx = region->x + region->dx - reg_right.x;
    histo_right = mriHistogramRegion(mri, 0, NULL, &reg_right);

    HISTOsubtract(h_prev, histo_left, histo);
    HISTOadd(histo, histo_right, histo);
    HISTOfree(&histo_left);
    HISTOfree(&histo_right);
  } else {
    mriHistogramRegion(mri, nbins, histo, region);
  }

  mri_prev = mri;
  HISTOcopy(histo, h_prev);
  REGIONcopy(region, &reg_prev);
  return (histo);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIhistogramLabelStruct(MRI *mri, int nbins, HISTOGRAM *histo, LABEL *label)
{
  float fmin, fmax, bin_size;

  MRIvalRange(mri, &fmin, &fmax);
  if (!nbins) nbins = nint((fmax - fmin) + 1.0);

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    // histo->nbins = nbins ;
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }
  HISTOclear(histo, histo);
  if (nbins == 256) histo->bin_size = bin_size = 1;

  histo->bin_size = bin_size = (fmax - fmin + 1) / (float)nbins;

  mriHistogramLabel(mri, nbins, histo, label);

  return (histo);
}

static HISTOGRAM *mriHistogramLabel(MRI *mri, int nbins, HISTOGRAM *histo, LABEL *label)
{
  // int width, height, depth;
  int x, y, z, bin_no, i;
  float fmin, fmax, bin_size, val;
  double xv, yv, zv;

  if (mri->type == MRI_UCHAR) {
    fmin = 0;
    fmax = 255;
  }
  else
    MRIvalRange(mri, &fmin, &fmax);

  if (!nbins) nbins = nint(fmax - fmin + 1.0);

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    // histo->nbins = nbins ;
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }
  HISTOclear(histo, histo);

  if (nbins == 256) {
    fmax = 255;
    fmin = 0;
    bin_size = 1;
  }
  else
    bin_size = histo->bin_size = (fmax - fmin + 1) / (float)nbins;
  // width = mri->width;
  // height = mri->height;
  // depth = mri->depth;

  for (bin_no = 0; bin_no < nbins; bin_no++) histo->bins[bin_no] = (bin_no)*bin_size + fmin;

  switch (mri->type) {
    case MRI_UCHAR:
      for (i = 0; i < label->n_points; i++) {
        // MRIworldToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z,
        //                 &xv, &yv, &zv) ;
        MRIsurfaceRASToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z, &xv, &yv, &zv);
        x = nint(xv);
        y = nint(yv);
        z = nint(zv);
        val = (float)MRIvox(mri, x, y, z);
        bin_no = nint((float)(val - fmin) / (float)bin_size);
        histo->counts[bin_no]++;
      }
      break;
    case MRI_SHORT:
      for (i = 0; i < label->n_points; i++) {
        // MRIworldToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z,
        //                &xv, &yv, &zv) ;
        MRIsurfaceRASToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z, &xv, &yv, &zv);
        x = nint(xv);
        y = nint(yv);
        z = nint(zv);
        val = (float)MRISvox(mri, x, y, z);
        bin_no = nint((float)(val - fmin) / (float)bin_size);
        histo->counts[bin_no]++;
      }
      break;
    case MRI_FLOAT:
      for (i = 0; i < label->n_points; i++) {
        // MRIworldToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z,
        //                &xv, &yv, &zv) ;
        MRIsurfaceRASToVoxel(mri, label->lv[i].x, label->lv[i].y, label->lv[i].z, &xv, &yv, &zv);
        x = nint(xv);
        y = nint(yv);
        z = nint(zv);
        val = (float)MRIFvox(mri, x, y, z);
        bin_no = nint((float)(val - fmin) / (float)bin_size);
        histo->counts[bin_no]++;
      }
      break;
    default:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "mriHistogramLabel: unsupported mri type %d", mri->type));
      break;
  }

  return (histo);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static HISTOGRAM *mriHistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region)
{
  int width, height, depth, x, y, z, bin_no, x0, y0, z0;
  float fmin, fmax;
  BUFTYPE val, *psrc;
  // BUFTYPE bmin, bmax;
  short *spsrc;
  float *fpsrc;

#if 0
  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,"mriHistogramRegion: must by type UCHAR"));
#endif

  // WHY DIFFER BY MRI->TYPE??
  /*  if (mri->type == MRI_UCHAR)
      {
      fmin = 0 ; fmax = 255 ;
      }
      else
  */

  MRIvalRangeRegion(mri, &fmin, &fmax, region);

  // bmin = (BUFTYPE)fmin;
  // bmax = (BUFTYPE)fmax;

  if (!nbins) nbins = nint(fmax - fmin + 1.0);

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }
  HISTOclear(histo, histo);
  HISTOinit(histo, nbins, fmin, fmax);
  //  histo->bin_size = bin_size = (fmax - fmin + 1) / (float)nbins ;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  width = region->x + region->dx;
  if (width > mri->width) width = mri->width;
  height = region->y + region->dy;
  if (height > mri->height) height = mri->height;
  depth = region->z + region->dz;
  if (depth > mri->depth) depth = mri->depth;
  x0 = region->x;
  if (x0 < 0) x0 = 0;
  y0 = region->y;
  if (y0 < 0) y0 = 0;
  z0 = region->z;
  if (z0 < 0) z0 = 0;

  for (bin_no = 0; bin_no < nbins; bin_no++) histo->bins[bin_no] = (bin_no)*histo->bin_size + fmin;

  switch (mri->type) {
    case MRI_UCHAR:
      for (z = z0; z < depth; z++) {
        for (y = y0; y < height; y++) {
          psrc = &MRIvox(mri, x0, y, z);
          for (x = x0; x < width; x++) {
            val = *psrc++;
            bin_no = nint((float)(val - fmin) / (float)histo->bin_size);
            if (bin_no >= nbins) bin_no = nbins - 1;
            histo->counts[bin_no]++;
          }
        }
      }
      break;
    case MRI_SHORT:
      for (z = z0; z < depth; z++) {
        for (y = y0; y < height; y++) {
          spsrc = &MRISvox(mri, x0, y, z);
          for (x = x0; x < width; x++) {
            bin_no = nint((float)(*spsrc++ - fmin) / (float)histo->bin_size);
            if (bin_no < 0) bin_no = 0;
            if (bin_no >= histo->nbins) bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
          }
        }
      }
      break;
    case MRI_FLOAT:
      for (z = z0; z < depth; z++) {
        for (y = y0; y < height; y++) {
          fpsrc = &MRIFvox(mri, x0, y, z);
          for (x = x0; x < width; x++) {
            bin_no = nint((float)(*fpsrc++ - fmin) / (float)histo->bin_size);
            if (bin_no < 0) bin_no = 0;
            if (bin_no >= histo->nbins) bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
          }
        }
      }
      break;
    default:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "mriHistogramRegion: unsupported mri type %d", mri->type));
      break;
  }

  return (histo);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIhistoEqualizeRegion(MRI *mri_src, MRI *mri_dst, int low, MRI_REGION *region)
{
  HISTOGRAM *histo_eq;

  if (region == NULL) {
    MRI_REGION tmp;
    region = &tmp;
    region->x = region->y = region->z = 0;
    region->dx = mri_src->width;
    region->dy = mri_src->height;
    region->dz = mri_src->depth;
  }
  histo_eq = HISTOalloc(256);
  MRIgetEqualizeHistoRegion(mri_src, histo_eq, low, region, 1);
  mri_dst = MRIapplyHistogramToRegion(mri_src, mri_dst, histo_eq, region);
  HISTOfree(&histo_eq);
  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIapplyHistogramToRegion(MRI *mri_src, MRI *mri_dst, HISTOGRAM *histo, MRI_REGION *region)
{
  int width, height, depth, x, y, z, x0, y0, z0;
  BUFTYPE *psrc, *pdst, sval, dval;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = region->x + region->dx;
  if (width > mri_src->width) width = mri_src->width;
  height = region->y + region->dy;
  if (height > mri_src->height) height = mri_src->height;
  depth = region->z + region->dz;
  if (depth > mri_src->depth) depth = mri_src->depth;
  x0 = region->x;
  if (x0 < 0) x0 = 0;
  y0 = region->y;
  if (y0 < 0) y0 = 0;
  z0 = region->z;
  if (z0 < 0) z0 = 0;

  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      psrc = &MRIvox(mri_src, x0, y, z);
      pdst = &MRIvox(mri_dst, x0, y, z);
      for (x = x0; x < width; x++) {
        sval = *psrc++;
        if (sval >= histo->nbins) sval = histo->nbins - 1;
        dval = histo->counts[sval];
        *pdst++ = dval;
      }
    }
  }
  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIgetEqualizeHistoRegion(MRI *mri, HISTOGRAM *histo_eq, int low, MRI_REGION *region, int norm)
{
  int i, nbins;
  HISTOGRAM *histo;
  float *pc, *pdst, total, total_pix;

  if (region == NULL) {
    MRI_REGION tmp;
    region = &tmp;
    region->x = region->y = region->z = 0;
    region->dx = mri->width;
    region->dy = mri->height;
    region->dz = mri->depth;
  }
  if (mri->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIgetEqualizeHistoRegion: unsupported type %d", mri->type));
  histo = MRIhistogramRegion(mri, 0, NULL, region);
  nbins = histo->nbins;
  if (!histo_eq)
    histo_eq = HISTOalloc(nbins);
  else {
    // histo_eq->nbins = nbins ;
    if (histo_eq->nbins < nbins)
      HISTOrealloc(histo_eq, nbins);
    else
      histo_eq->nbins = nbins;
    HISTOclear(histo_eq, histo_eq);
  }

  for (pc = &histo->counts[0], total_pix = 0, i = low; i < nbins; i++) total_pix += *pc++;

  if (total_pix) {
    pc = &histo->counts[0];
    pdst = &histo_eq->counts[0];

    for (total = 0, i = low; i < nbins; i++) {
      total += *pc++;
      if (norm)
        *pdst++ = nint(255.0f * (float)total / (float)total_pix);
      else
        *pdst++ = total;
    }
  }

  HISTOfree(&histo);
  return (histo_eq);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIgetEqualizeHisto(MRI *mri, HISTOGRAM *histo_eq, int low, int high, int norm)
{
  int i, total, total_pix;
  HISTOGRAM *histo;

  histo = MRIhistogram(mri, 0);
  if (!histo_eq) histo_eq = HISTOalloc(histo->nbins);

  if (high > histo->max) high = histo->max;

  for (total_pix = 0, i = low; i <= high; i++) total_pix += histo->counts[i];

  for (i = 0; i < low; i++) {
    histo_eq->counts[i] = 0;
    histo_eq->bins[i] = i;
  }

  for (total = 0, i = low; i <= high; i++) {
    total += histo->counts[i];
    if (norm)
      histo_eq->counts[i] = nint(255.0f * (float)total / (float)total_pix);
    else
      histo_eq->counts[i] = total;
    histo->bins[i] = i;
  }
  for (i = high + 1; i < histo->nbins; i++) {
    histo_eq->counts[i] = histo_eq->counts[high] + i;
    histo_eq->bins[i] = i;
  }

  histo_eq->nbins = histo->nbins;
  HISTOfree(&histo);
  return (histo_eq);
}
MRI *MRIhistogramNormalize(MRI *mri_src, MRI *mri_template, MRI *mri_dst)
{
  float scale, offset;
  HISTOGRAM *h_src, *h_template;

  h_src = MRIhistogram(mri_src, 0);
  h_template = MRIhistogram(mri_template, 0);

  HISTOfindLinearFit(h_src, h_template, .5, 1, -20, 20, &scale, &offset);
  printf("scaling image intensities by %2.3f + %2.3f\n", scale, offset);
  mri_dst = MRIscaleIntensities(mri_src, mri_dst, scale, offset);
  HISTOfree(&h_src);
  HISTOfree(&h_template);
  return (mri_dst);
}
MRI *MRIscaleIntensities(MRI *mri_src, MRI *mri_dst, float scale, float offset)
{
  int x, y, z, f;
  float val;

  if (mri_dst == NULL) mri_dst = MRIclone(mri_src, NULL);

  for (f = 0; f < mri_src->nframes; f++)
    for (x = 0; x < mri_src->width; x++)
      for (y = 0; y < mri_src->height; y++)
        for (z = 0; z < mri_src->depth; z++) {
          val = MRIgetVoxVal(mri_src, x, y, z, f);
          val = val * scale + offset;
          if (mri_dst->type == MRI_UCHAR && val > 255) val = 255;
          MRIsetVoxVal(mri_dst, x, y, z, f, val);
        }
  return (mri_dst);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIhistoEqualize(MRI *mri_src, MRI *mri_template, MRI *mri_dst, int low, int high)
{
  HISTOGRAM *histo_src, *histo_template;
  int width, height, depth, x, y, z, index, sval, dval;
  float pct;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  histo_src = MRIgetEqualizeHisto(mri_src, NULL, low, high, 1);
  histo_template = MRIgetEqualizeHisto(mri_template, NULL, low, high, 1);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(histo_src, "cum_src.plt");
    HISTOplot(histo_template, "cum_template.plt");
  }
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        sval = MRIgetVoxVal(mri_src, x, y, z, 0);
        dval = sval;
        if (sval > 80 && sval < 95) DiagBreak();
        if (sval >= 83 && sval <= 85) DiagBreak();

        if ((sval >= low) && (sval <= high)) {
          pct = histo_src->counts[sval];
          for (index = low; index < histo_template->nbins; index++) {
            if (histo_template->counts[index] >= pct) {
              dval = index;
              break;
            }
          }
        }
        if (sval >= 83 && sval <= 85) DiagBreak();
        if (dval > 246) DiagBreak();
        MRIsetVoxVal(mri_dst, x, y, z, 0, dval);
      }
    }
  }
  HISTOfree(&histo_template);
  HISTOfree(&histo_src);
  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *MRIcrunch(MRI *mri_src, MRI *mri_dst)
{
  HISTOGRAM *histo;
  int b, deleted;

  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIcrunch: unsupported type %d", mri_src->type));
  histo = MRIhistogram(mri_src, 0);

  deleted = 0;
  for (b = 0; b < histo->nbins; b++) {
    if (!histo->counts[b]) deleted++;
    histo->counts[b] = b - deleted;
  }

  histo->nbins -= deleted;
  mri_dst = MRIapplyHistogram(mri_src, mri_dst, histo);
  HISTOfree(&histo);
  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:
     if nbins==0 then compute the max range
  Returns value:
     the computed histogram
  Description
     compute a histogram of all frames in the input image
  ------------------------------------------------------*/
HISTOGRAM *MRIhistogram(MRI *mri, int nbins)
{
  int width, height, depth, x, y, z, bin_no, frame;
  HISTOGRAM *histo;
  float fmin, fmax;
  double val;

  MRIvalRange(mri, &fmin, &fmax);  // fmin = is wrong!
  if (!nbins) nbins = nint(fmax - fmin + 1.0);
  if (fmin == fmax) fmax = 255;
  histo = HISTOalloc(nbins);
  HISTOinit(histo, nbins, fmin, fmax);
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  for (frame = 0; frame < mri->nframes; frame++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          val = MRIgetVoxVal(mri, x, y, z, frame);
          bin_no = nint((float)(val - fmin) / (float)histo->bin_size);
          histo->counts[bin_no]++;
        }
      }
    }
  }
  return (histo);
}
HISTOGRAM *MRIhistogramLabelRegion(MRI *mri, MRI *mri_labeled, MRI_REGION *region, int label, int nbins)
{
  int width, height, depth, x, y, z, bin_no, x0, x1, y0, y1, z0, z1;
  HISTOGRAM *histo;
  float fmin, fmax;
  BUFTYPE *psrc;
  int val, bmin;
  // int bmax;
  float fval;

  MRIvalRangeRegion(mri, &fmin, &fmax, region);
  bmin = (int)fmin;
  // bmax = (int)fmax;
  if (nbins <= 0) nbins = nint(fmax - fmin + 1.0);

  histo = HISTOalloc(nbins);
  HISTOinit(histo, nbins, fmin, fmax);

  if (FEQUAL(fmin, fmax)) {
    histo->bin_size = 1;
    ErrorReturn(histo, (ERROR_BADPARM, "MRIhistogramLabelRegion: constant image"));
  }
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  x0 = MAX(0, region->x);
  y0 = MAX(0, region->y);
  z0 = MAX(0, region->z);
  x1 = MIN(width, region->x + region->dx);
  y1 = MIN(height, region->y + region->dy);
  z1 = MIN(depth, region->z + region->dz);
  for (z = z0; z < z1; z++) {
    for (y = y0; y < y1; y++) {
      for (x = x0; x < x1; x++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        if (nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) != label) continue;
        switch (mri->type) {
          case MRI_UCHAR:
            /* 0 -> x */
            psrc = &MRIvox(mri, x, y, z);
            val = *psrc++;
            bin_no = nint((float)(val - bmin) / (float)histo->bin_size);
            histo->counts[bin_no]++;
            break;
          case MRI_SHORT:
            val = MRISvox(mri, x, y, z);
            bin_no = nint((float)(val - bmin) / (float)histo->bin_size);
            histo->counts[bin_no]++;
            break;
          case MRI_FLOAT:
            fval = MRIFvox(mri, x, y, z);
            bin_no = nint((fval - fmin) / (float)histo->bin_size);
            histo->counts[bin_no]++;
            break;

          default:
            ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIhistogramLabelRegion: unsupported type %d", mri->type));
            break;
        }
      }
    }
  }
  return (histo);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIhistogramLabel(MRI *mri, MRI *mri_labeled, int label, int nbins)
{
  int width, height, depth, x, y, z, bin_no;
  HISTOGRAM *histo;
  float fmin, fmax;
  int val, bmin;
  // int bmax;
  float fval;

  MRIlabelValRange(mri, mri_labeled, label, &fmin, &fmax);
  bmin = (int)fmin;
  // bmax = (int)fmax;
  if (nbins <= 0) nbins = nint(fmax - fmin + 1.0);

  histo = HISTOalloc(nbins);
  HISTOinit(histo, nbins, fmin, fmax);

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (MRIgetVoxVal(mri_labeled, x, y, z, 0) != label) continue;
        switch (mri->type) {
          case MRI_UCHAR:
            /* 0 -> x */
            val = MRIvox(mri, x, y, z);
            bin_no = nint((float)(val - bmin) / (float)histo->bin_size);
            if (bin_no == 0) DiagBreak();
            if (bin_no < 0)
              bin_no = 0;
            else if (bin_no >= histo->nbins)
              bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
            break;
          case MRI_SHORT:
            val = MRISvox(mri, x, y, z);
            bin_no = nint((float)(val - bmin) / (float)histo->bin_size);
            if (bin_no < 0)
              bin_no = 0;
            else if (bin_no >= histo->nbins)
              bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
            break;
          case MRI_FLOAT:
            fval = MRIFvox(mri, x, y, z);
            bin_no = nint((fval - fmin) / (float)histo->bin_size);
            if (bin_no < 0 || bin_no >= histo->nbins) DiagBreak();
            if (bin_no < 0)
              bin_no = 0;
            else if (bin_no >= histo->nbins)
              bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
            break;

          default:
            fval = MRIgetVoxVal(mri, x, y, z, 0);
            bin_no = nint((fval - fmin) / (float)histo->bin_size);
            if (bin_no < 0 || bin_no >= histo->nbins) DiagBreak();
            if (bin_no < 0)
              bin_no = 0;
            else if (bin_no >= histo->nbins)
              bin_no = histo->nbins - 1;
            histo->counts[bin_no]++;
            break;
        }
      }
    }
  }
  return (histo);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Remove small inconsistincies in the labelling of a volume.
  ------------------------------------------------------*/
#define MIN_VOXELS_PCT 0.15 /* peak must average 15% of max count */
#define PEAK_SEPARATION 11  /* min width of peak */
#define VALLEY_WIDTH 3
#define GRAY_LOW 65 /* lowest interesting gray matter */

#define SMOOTH_SIGMA 4.0f

#define X_DB2 98
#define Y_DB2 147
#define Z_DB2 93

#define X_DB1 106
#define Y_DB1 158
#define Z_DB1 137

#define X_DB 157
#define Y_DB 153
#define Z_DB 128

MRI *MRIhistoSegment(MRI *mri_src, MRI *mri_labeled, int wm_low, int wm_hi, int gray_hi, int wsize, float sigma)
{
  int width, height, depth, x, y, z, whalf, in_val, label, nvox, valley, wm_peak, gray_peak, nlabeled;
  // int thresh;
  BUFTYPE *pdst, *psrc;
  MRI_REGION region;
  HISTOGRAM *histo, *hsmooth;

  if (sigma < 0) sigma = SMOOTH_SIGMA;

  // thresh = (gray_hi + wm_low) / 2;
  histo = hsmooth = NULL;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_labeled) mri_labeled = MRIclone(mri_src, NULL);

  whalf = wsize / 2;
  region.dx = region.dy = region.dz = wsize;

#if 0
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_labeled, "label.mnc") ;
    MRIwrite(mri_src, "src.mnc") ;
  }
#endif
  for (nlabeled = nvox = z = 0; z < depth; z++) {
    DiagShowPctDone((float)z / (float)(depth - 1), 5);
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_labeled, 0, y, z);
      psrc = &MRIvox(mri_src, 0, y, z);
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();

        label = *pdst;
        in_val = *psrc++;
        if (x == X_DB && y == Y_DB && z == Z_DB) DiagBreak();
        if (label != MRI_AMBIGUOUS) {
          pdst++;
          continue;
        }
        nvox++;
        region.x = x - whalf;
        region.y = y - whalf;
        region.z = z - whalf;
        histo = mriHistogramRegion(mri_src, 0, histo, &region);

        /*        HISTOclearBins(histo, histo, in_val-1, in_val+1) ;*/
        hsmooth = HISTOsmooth(histo, hsmooth, sigma);
#if 0
        HISTOclearBins(histo, histo, 0, GRAY_LOW-5) ;
        HISTOclearBins(histo, histo, wm_hi+5, 255) ;
#endif
        wm_peak = HISTOfindLastPeakInRegion(
            hsmooth, PEAK_SEPARATION, MIN_VOXELS_PCT, wm_low /*GRAY_LOW*/, wm_hi - PEAK_SEPARATION / 2 - 2);
        if (wm_peak >= 0) wm_peak = hsmooth->bins[wm_peak];  // convert it to an image intensity
        gray_peak = HISTOfindLastPeakInRegion(hsmooth,
                                              PEAK_SEPARATION,
                                              MIN_VOXELS_PCT,
                                              GRAY_LOW + PEAK_SEPARATION / 2 + 2,
                                              wm_peak - PEAK_SEPARATION + 1);
        if (gray_peak >= 0) gray_peak = hsmooth->bins[gray_peak];  // convert it to an image intensity
        if (gray_peak >= 0 && wm_peak >= 0) {
          while (gray_peak > gray_hi) /* white matter is bimodal */
          {
            wm_peak = gray_peak;
            gray_peak = HISTOfindLastPeakInRegion(hsmooth,
                                                  PEAK_SEPARATION,
                                                  MIN_VOXELS_PCT,
                                                  GRAY_LOW + PEAK_SEPARATION / 2 + 2,
                                                  wm_peak - PEAK_SEPARATION + 1);
            if (gray_peak >= 0) gray_peak = hsmooth->bins[gray_peak];  // convert it to an image intensity
          }
        }

        if ((wm_peak < 0) || (gray_peak < 0)) /* unimodal - take best guess */
          valley = -1 /*thresh*/;
        else /* bimodal, find min between peaks and use it as descriminant */
          valley = HISTOfindValley(hsmooth, VALLEY_WIDTH, gray_peak + VALLEY_WIDTH - 1, wm_peak - VALLEY_WIDTH + 1);
        if (valley >= 0) valley = hsmooth->bins[valley];  // convert it to an image intensity
        if (valley > gray_hi)                             /* can't be proper discriminant */
          valley = -1;
#if 0
        if (x == X_DB && y == Y_DB && z == Z_DB)
        {
          FILE *fp = fopen("histo.dat", "w") ;
          fprintf(fp, "histogram at (%d, %d, %d) = %d\n", x, y, z, in_val) ;
          fprintf(fp, "wm peak = %d, gray peak = %d, valley = %d\n",
                  wm_peak, gray_peak, valley) ;
          HISTOdump(histo, fp) ;
          fprintf(fp, "smooth histo:\n") ;
          HISTOdump(hsmooth, fp) ;
          fclose(fp) ;
          HISTOplot(hsmooth, "hsmooth.plt") ;
          HISTOplot(histo, "histo.plt") ;
        }
        if (x == X_DB1 && y == Y_DB1 && z == Z_DB1)
        {
          FILE *fp = fopen("histo1.dat", "w") ;
          fprintf(fp, "histogram at (%d, %d, %d) = %d\n", x, y, z, in_val) ;
          fprintf(fp, "wm peak = %d, gray peak = %d, valley = %d\n",
                  wm_peak, gray_peak, valley) ;
          HISTOdump(histo, fp) ;
          fprintf(fp, "smooth histo:\n") ;
          HISTOdump(hsmooth, fp) ;
          fclose(fp) ;
          HISTOplot(hsmooth, "hsmooth1.plt") ;
          HISTOplot(histo, "histo1.plt") ;
        }
#endif
#if 0
        if (valley < 0)
          valley = (wm_peak + gray_peak) / 2 ;  /* assume equal sigmas */
#endif
        if (valley < 0 || abs(in_val - valley) <= 2)
          *pdst++ = MRI_AMBIGUOUS;
        else if (valley >= gray_hi)
          *pdst++ = MRI_AMBIGUOUS;
        else if (in_val >= valley) {
          nlabeled++;
          *pdst++ = MRI_WHITE;
        }
        else {
          nlabeled++;
          *pdst++ = MRI_NOT_WHITE;
        }
      }
    }
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr,
            "              %8d voxels processed (%2.2f%%)\n",
            nvox,
            100.0f * (float)nvox / (float)(width * height * depth));
    fprintf(stderr,
            "              %8d voxels labeled (%2.2f%%)\n",
            nlabeled,
            100.0f * (float)nlabeled / (float)(width * height * depth));
  }
  if (histo) HISTOfree(&histo);
  if (hsmooth) HISTOfree(&hsmooth);
  return (mri_labeled);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Remove small inconsistincies in the labelling of a volume.
  ------------------------------------------------------*/
MRI *MRIhistoSegmentVoxel(MRI *mri_src,
                          MRI *mri_labeled,
                          int wm_low,
                          int wm_hi,
                          int gray_hi,
                          int wsize,
                          int x,
                          int y,
                          int z,
                          HISTOGRAM *histo,
                          HISTOGRAM *hsmooth,
                          float sigma)
{
  int whalf, in_val, valley, wm_peak, gray_peak /*, thresh */;
  MRI_REGION region;
  float sig;

  // thresh = (gray_hi + wm_low) / 2;

  whalf = wsize / 2;
  region.dx = region.dy = region.dz = wsize;

  in_val = MRIvox(mri_src, x, y, z);
  if (x == X_DB && y == Y_DB && z == Z_DB) DiagBreak();

  region.x = x - whalf;
  region.y = y - whalf;
  region.z = z - whalf;
  histo = mriHistogramRegion(mri_src, 0, histo, &region);
  hsmooth = HISTOsmooth(histo, hsmooth, sigma);
/*        HISTOclearBins(histo, histo, in_val-1, in_val+1) ;*/

#if 1
  sig = sigma;
#else
  int npeaks, peaks[300];
  sig = 0.0f;
  do {
    hsmooth = HISTOsmooth(histo, hsmooth, sigma);
    npeaks = HISTOcountPeaksInRegion(
        hsmooth, PEAK_SEPARATION, MIN_VOXELS_PCT, peaks, 5, GRAY_LOW, wm_hi - PEAK_SEPARATION / 2 - 2);
    if (npeaks != 2) {
      sig += 0.5f;
      if (sig > 8.0f) {
        fprintf(stderr, "could not find two peaks\n");
        return (NULL);
      }
    }
  } while (npeaks != 2);
#endif

  HISTOclearBins(histo, histo, 0, GRAY_LOW - 5);
  HISTOclearBins(histo, histo, wm_hi + 5, 255);
  wm_peak =
      HISTOfindLastPeakInRegion(hsmooth, PEAK_SEPARATION, MIN_VOXELS_PCT, GRAY_LOW, wm_hi - PEAK_SEPARATION / 2 - 2);
  if (wm_peak >= 0) wm_peak = hsmooth->bins[wm_peak];  // convert it to an intensity
  gray_peak = HISTOfindLastPeakInRegion(
      hsmooth, PEAK_SEPARATION, MIN_VOXELS_PCT, GRAY_LOW + PEAK_SEPARATION / 2 + 2, wm_peak - PEAK_SEPARATION + 1);
  if (gray_peak >= 0) gray_peak = hsmooth->bins[gray_peak];  // convert it to an intensity
  while (gray_peak > gray_hi)                                /* white matter is bimodal */
  {
    wm_peak = gray_peak;
    gray_peak = HISTOfindLastPeakInRegion(
        hsmooth, PEAK_SEPARATION, MIN_VOXELS_PCT, GRAY_LOW + PEAK_SEPARATION / 2 + 2, wm_peak - PEAK_SEPARATION + 1);
    if (gray_peak >= 0) gray_peak = hsmooth->bins[gray_peak];  // convert it to an intensity
  }

  if ((wm_peak < 0) || (gray_peak < 0)) /* unimodal - take best guess */
    valley = -1;                        /* was thresh */
  else                                  /* bimodal, find min between peaks and use it as descriminant */
    valley = HISTOfindValley(hsmooth, VALLEY_WIDTH, gray_peak, wm_peak);
  if (valley >= 0) valley = hsmooth->bins[valley];  // convert it to an intensity

  {
    FILE *fp = fopen("histo.dat", "w");
    fprintf(fp, "histogram at (%d, %d, %d) = %d\n", x, y, z, in_val);
    fprintf(fp, "wm peak = %d, gray peak = %d, valley = %d, sigma = %2.2f\n", wm_peak, gray_peak, valley, sig);
    printf("histogram at (%d, %d, %d) = %d\n", x, y, z, in_val);
    printf("wm peak = %d, gray peak = %d, valley = %d\n", wm_peak, gray_peak, valley);
    HISTOdump(histo, fp);
    fprintf(fp, "smooth histo:\n");
    HISTOdump(hsmooth, fp);
    fclose(fp);
    HISTOplot(hsmooth, "hsmooth.plt");
    HISTOplot(histo, "histo.plt");
  }
  if (valley < 0) valley = -1; /*(wm_peak + gray_peak) / 2 ;*/ /* assume equal sigmas */
  if (valley < 0)
    MRIvox(mri_labeled, x, y, z) = MRI_AMBIGUOUS;
  else if (in_val >= valley)
    MRIvox(mri_labeled, x, y, z) = MRI_WHITE;
  else
    MRIvox(mri_labeled, x, y, z) = MRI_NOT_WHITE;
  return (mri_labeled);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
HISTOGRAM *MRIhistogramVoxel(
    MRI *mri, int nbins, HISTOGRAM *histo, int x0, int y0, int z0, int wsize, MRI *mri_thresh, float thresh)
{
  int whalf;
  float fmin, fmax;
  MRI_REGION region;

  whalf = (wsize - 1) / 2;
  region.x = x0 - whalf;
  region.dx = wsize;
  region.y = y0 - whalf;
  region.dy = wsize;
  region.z = z0 - whalf;
  region.dz = wsize;
  MRIvalRangeRegion(mri, &fmin, &fmax, &region);
  if (nbins <= 0) nbins = nint(fmax - fmin + 1.0);
  if (nbins <= 0) nbins = 255;

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }
  HISTOclear(histo, histo);
  HISTOinit(histo, nbins, fmin, fmax);

  if (mri_thresh == NULL)
    mriHistogramRegion(mri, nbins, histo, &region);
  else
    MRIhistogramRegionWithThreshold(mri, nbins, histo, &region, mri_thresh, thresh, 0);

  HISTOsoapBubbleZeros(histo, histo, 100);
  return (histo);
}
HISTOGRAM *MRIhistogramRegionWithThreshold(
    MRI *mri, int nbins, HISTOGRAM *histo, MRI_REGION *region, MRI *mri_thresh, float thresh, int frame)
{
  int width, height, depth, z, x0, y0, z0;
  float fmin, fmax;
#ifdef HAVE_OPENMP
  int tid;
  HISTOGRAM *histos[_MAX_FS_THREADS];
#else
  HISTOGRAM *histos[1];
#endif

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  width = region->x + region->dx;
  if (width > mri->width) width = mri->width;
  height = region->y + region->dy;
  if (height > mri->height) height = mri->height;
  depth = region->z + region->dz;
  if (depth > mri->depth) depth = mri->depth;
  x0 = region->x;
  if (x0 < 0) x0 = 0;
  y0 = region->y;
  if (y0 < 0) y0 = 0;
  z0 = region->z;
  if (z0 < 0) z0 = 0;

  fmin = 1e10;
  fmax = -fmin;
  for (z = z0; z < depth; z++) {
    int y, x;
    float val;
    for (y = y0; y < height; y++) {
      for (x = x0; x < width; x++) {
        val = MRIgetVoxVal(mri_thresh, x, y, z, frame);
        if (val < thresh) continue;
        val = MRIgetVoxVal(mri, x, y, z, frame);
        if (val < fmin) fmin = val;
        if (val > fmax) fmax = val;
      }
    }
  }

  if (!nbins) nbins = nint(fmax - fmin + 1.0);

  if (!histo)
    histo = HISTOalloc(nbins);
  else {
    if (histo->nbins < nbins)
      HISTOrealloc(histo, nbins);
    else
      histo->nbins = nbins;
  }

#ifdef HAVE_OPENMP
  for (tid = 0; tid < _MAX_FS_THREADS; tid++) {
    histos[tid] = HISTOalloc(nbins);
    HISTOclear(histos[tid], histos[tid]);
    HISTOinit(histos[tid], nbins, fmin, fmax);
  }
#else
  histos[0] = histo;
#endif
  HISTOclear(histo, histo);
  HISTOinit(histo, nbins, fmin, fmax);

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) shared(histos, width, height, depth, fmin, fmax, frame, x0, y0, z0, histo)
#endif
  for (z = z0; z < depth; z++) {
    ROMP_PFLB_begin
    
    int y, x, tid;
    float val;
    for (y = y0; y < height; y++) {
      for (x = x0; x < width; x++) {
        val = MRIgetVoxVal(mri_thresh, x, y, z, frame);
        if (val < thresh) continue;
        val = MRIgetVoxVal(mri, x, y, z, frame);
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        HISTOaddSample(histos[tid], val, fmin, fmax);
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

#ifdef HAVE_OPENMP
  for (tid = 0; tid < _MAX_FS_THREADS; tid++) {
    HISTOadd(histos[tid], histo, histo);
    HISTOfree(&histos[tid]);
  }
#endif

  return (histo);
}
#define NBINS 1000
double MRIfindPercentile(MRI *mri, double percentile, int frame)
{
  int x, y, z, val, bin;
  HISTO *histo, *hcdf;
  float min_val, max_val;

  MRIvalRange(mri, &min_val, &max_val);
  histo = HISTOinit(NULL, NBINS, min_val, max_val);

  for (x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        val = MRIgetVoxVal(mri, x, y, z, frame);
        if (FZERO(val)) continue;
        HISTOaddSample(histo, val, 0, 0);
      }

  hcdf = HISTOmakeCDF(histo, NULL);
  bin = HISTOfindBinWithCount(hcdf, percentile);
  val = hcdf->bins[bin];
  HISTOfree(&hcdf);
  HISTOfree(&histo);

  return (val);
}
