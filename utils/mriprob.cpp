/*
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

#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "proto.h"

#define DEBUG_POINT(x, y, z) (((x == 140) && (y == 74)) && ((z) == 174))

MRI *MRIapplyBayesLaw(MRI *mri_priors, MRI *mri_p1, MRI *mri_p2, MRI *mri_dst)
{
  int x, y, z, width, height, depth;
  BUFTYPE *ppriors, *pdst;
  float p, p1, p2, prior, *pp1, *pp2;

  if (!mri_dst) mri_dst = MRIclone(mri_priors, NULL);

  width = mri_dst->width;
  height = mri_dst->height;
  depth = mri_dst->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      ppriors = &MRIvox(mri_priors, 0, y, z);
      pp1 = &MRIFvox(mri_p1, 0, y, z);
      pp2 = &MRIFvox(mri_p2, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        p1 = (float)*pp1++;
        p2 = (float)*pp2++;
        prior = (float)*ppriors++;
        p1 /= 100.0f;
        p2 /= 100.0f;
        prior /= 100.0f;
        p1 *= prior;
        p2 *= (1.0f - prior);
        p = p1 / (p1 + p2);
        *pdst++ = (BUFTYPE)nint(p * 100.0f);
      }
    }
  }
  return (mri_dst);
}
MRI *MRIcomputeConditionalProbabilities(MRI *mri_T1, MRI *mri_mean, MRI *mri_std, MRI *mri_dst)
{
  int x, y, z, width, height, depth;
  BUFTYPE *pT1, *pmean, *pstd;
  float p, mean, std, val, n, *pdst;

  width = mri_T1->width;
  height = mri_T1->height;
  depth = mri_T1->depth;

  if (!mri_dst) mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pT1 = &MRIvox(mri_T1, 0, y, z);
      pmean = &MRIvox(mri_mean, 0, y, z);
      pstd = &MRIvox(mri_std, 0, y, z);
      pdst = &MRIFvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        val = (float)*pT1++;
        mean = (float)*pmean++;
        std = (float)*pstd++;
        if (FZERO(std)) std = 1.0;
#if 0
        if (std < 10.0)  /* hack!!!!! - not enough observations */
          std = 10.0 ;
#endif
        n = 1 / (std * sqrt(2.0 * M_PI));
        p = n * exp(-SQR(val - mean) / (2.0f * SQR(std)));
        *pdst++ = p * 100.0f;
      }
    }
  }
  return (mri_dst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
MRI *MRIprobabilityThresholdNeighborhoodOff(MRI *mri_src, MRI *mri_prob, MRI *mri_dst, float threshold, int nsize)
{
  BUFTYPE *pprob, out_val;
  int width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi, xmin, xmax, ymin, ymax, zmin, zmax, nchanged;
  float nvox;

  if (mri_prob->type != MRI_UCHAR) ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRI3Dthreshold: prob must be MRI_UCHAR"));

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  /* now apply the inverse morph to build an average wm representation
     of the input volume
     */

  xmin = width;
  ymin = height;
  zmin = depth;
  xmax = ymax = zmax = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pprob = &MRIvox(mri_prob, 0, y, z);
      for (x = 0; x < width; x++) {
        if (*pprob++ > threshold) {
          if (x < xmin) xmin = x;
          if (x > xmax) xmax = x;
          if (y < ymin) ymin = y;
          if (y > ymax) ymax = y;
          if (z < zmin) zmin = z;
          if (z > zmax) zmax = z;
        }
      }
    }
  }
  xmin = MAX(xmin - nsize, 0);
  ymin = MAX(ymin - nsize, 0);
  zmin = MAX(zmin - nsize, 0);
  xmax = MIN(xmax + nsize, width - 1);
  ymax = MIN(ymax + nsize, height - 1);
  zmax = MIN(zmax + nsize, depth - 1);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "bounding box = (%d:%d, %d:%d, %d:%d).\n", xmin, xmax, ymin, ymax, zmin, zmax);

  /* remove stuff outside bounding box */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if ((x < xmin || x > xmax) || (y < ymin || y > ymax) || (z < zmin || z > zmax)) MRIvox(mri_dst, x, y, z) = 0;
      }
    }
  }

  nchanged = 0;
  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        out_val = 0;
        for (z1 = -nsize; z1 <= nsize; z1++) {
          zi = mri_src->zi[z + z1];
          for (y1 = -nsize; y1 <= nsize; y1++) {
            yi = mri_src->yi[y + y1];
            for (x1 = -nsize; x1 <= nsize; x1++) {
              xi = mri_src->xi[x + x1];
              if (MRIvox(mri_prob, xi, yi, zi) > threshold) {
                out_val = MRIvox(mri_src, x, y, z);
                break;
              }
            }
            if (out_val > 0) break;
          }
          if (out_val > 0) break;
        }
        if (out_val == 0 && MRIvox(mri_src, x, y, z) > 0) nchanged++;
        MRIvox(mri_dst, x, y, z) = out_val;
      }
    }
  }

  nvox = (float)(width * height * depth);
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n", nchanged, (int)nvox, 100.0f * (float)nchanged / nvox);

  return (mri_dst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
MRI *MRIprobabilityThresholdNeighborhoodOn(
    MRI *mri_src, MRI *mri_prob, MRI *mri_dst, float threshold, int nsize, int out_label)
{
  BUFTYPE *pprob, out_val;
  int width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi, xmin, xmax, ymin, ymax, zmin, zmax, nchanged;
  float nvox;

  if (mri_prob->type != MRI_UCHAR) ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRI3Dthreshold: prob must be MRI_UCHAR"));

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  /* now apply the inverse morph to build an average wm representation
     of the input volume
     */

  xmin = width;
  ymin = height;
  zmin = depth;
  xmax = ymax = zmax = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pprob = &MRIvox(mri_prob, 0, y, z);
      for (x = 0; x < width; x++) {
        if (*pprob++ > threshold) {
          if (x < xmin) xmin = x;
          if (x > xmax) xmax = x;
          if (y < ymin) ymin = y;
          if (y > ymax) ymax = y;
          if (z < zmin) zmin = z;
          if (z > zmax) zmax = z;
        }
      }
    }
  }
  xmin = MAX(xmin - nsize, 0);
  ymin = MAX(ymin - nsize, 0);
  zmin = MAX(zmin - nsize, 0);
  xmax = MIN(xmax + nsize, width - 1);
  ymax = MIN(ymax + nsize, height - 1);
  zmax = MIN(zmax + nsize, depth - 1);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "bounding box = (%d:%d, %d:%d, %d:%d).\n", xmin, xmax, ymin, ymax, zmin, zmax);

#if 0
  xmin = ymin = zmin = 0 ;
  xmax = width-1 ;
  ymax = height-1 ;
  zmax = depth-1 ;
#endif

  nchanged = 0;
  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        out_val = out_label;
        for (z1 = -nsize; z1 <= nsize; z1++) {
          zi = mri_src->zi[z + z1];
          for (y1 = -nsize; y1 <= nsize; y1++) {
            yi = mri_src->yi[y + y1];
            for (x1 = -nsize; x1 <= nsize; x1++) {
              xi = mri_src->xi[x + x1];
              if (MRIvox(mri_prob, xi, yi, zi) < threshold) {
                out_val = MRIvox(mri_src, x, y, z);
                break;
              }
            }
          }
        }
        if (out_val == out_label && MRIvox(mri_src, x, y, z) != out_label) nchanged++;
        MRIvox(mri_dst, x, y, z) = out_val;
      }
    }
  }

  nvox = (float)(width * height * depth);
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n", nchanged, (int)nvox, 100.0f * (float)nchanged / nvox);

  return (mri_dst);
}

MRI *MRIprobabilityThreshold(MRI *mri_src, MRI *mri_prob, MRI *mri_dst, float threshold, int out_label)
{
  BUFTYPE *pprob, *pdst, *psrc, out_val, prob, in_val;
  int width, height, depth, x, y, z, nchanged, noff, non;
  float nvox;

  if (mri_prob->type != MRI_UCHAR) ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRI3Dthreshold: prob must be MRI_UCHAR"));

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  /* now apply the inverse morph to build an average wm representation
     of the input volume
     */

  noff = non = nchanged = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pprob = &MRIvox(mri_prob, 0, y, z);
      psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        out_val = 0;
        prob = *pprob++; /* value from inverse morphed volume */
        in_val = *psrc++;
        if (in_val < WM_MIN_VAL && prob >= threshold) /* probably on */
          out_val = out_label;
        else /* not sure, use original val */
          out_val = in_val;

        if (out_val != in_val) {
          if (out_val)
            non++;
          else
            noff++;
          nchanged++;
        }
        *pdst++ = out_val;
      }
    }
  }

  nvox = (float)(width * height * depth);
  fprintf(stderr, "%8d of %8d voxels changed - %2.1f%%.\n", nchanged, (int)nvox, 100.0f * (float)nchanged / nvox);
  fprintf(stderr, "%8d of %8d voxels off     - %2.1f%%.\n", noff, (int)nvox, 100.0f * (float)noff / nvox);
  fprintf(stderr, "%8d of %8d voxels on      - %2.1f%%.\n", nchanged, (int)nvox, 100.0f * (float)non / nvox);
  return (mri_dst);
}
