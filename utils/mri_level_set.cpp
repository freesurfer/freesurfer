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

/*----------------------------------------------------------------------
           File Name:

           Author:

           Description:

           Conventions:

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "diag.h"
#include "macros.h"
#include "mri.h"
#include "proto.h"

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

/*----------------------------------------------------------------------
                              FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#define DT 0.5
#define E_DIFFUSION 5e-8 /* increasing this gives a smoother solution */
#define E_CURV 0.4       /* increasing this gives less curvature to front */
#define MAG_SCALE 1.0    /* increasing this increases the effect of grad */
#define MAG_POW 4.0
#define WRITE_ITER 1000
#define MIN_DIST 2.0

/* solve ut + f |grad(u)| = e uxx */

MRI *MRIregionGrow(MRI *mri_src, MRI *mri_distance, float x0, float y0, float z0, int niter)
{
  int width, height, depth, t, x, y, z, xp1, xm1, yp1, ym1, zp1, zm1, write_iter;
  double dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz, dist_xp1, dist_xm1, dist_yp1, dist_ym1, dist_zp1, dist_zm1, dist, km,
      denom, mag, dist_xp1yp1, dist_xp1zp1, dist_yp1zp1, grad, F, laplacian, dx_f, dx_b, dy_f, dy_b, dz_f, dz_b, sdx,
      sdy, sdz, orig_mag;
  double mag_pow, mag_scale, dt, e_diffusion, e_curv;
  float f;
  BUFTYPE src, src_xp1, src_xm1, src_yp1, src_ym1, src_zp1, src_zm1, min_orig_src, max_orig_src;
  char *cp;

  cp = getenv("WRITE_ITER");
  write_iter = WRITE_ITER;
  if (cp) sscanf(cp, "%d", &write_iter);
  cp = getenv("MAG_POW");
  mag_pow = MAG_POW;
  if (cp) sscanf(cp, "%lf", &mag_pow);
  cp = getenv("MAG_SCALE");
  mag_scale = MAG_SCALE;
  if (cp) sscanf(cp, "%lf", &mag_scale);
  cp = getenv("E_DIFFUSION");
  e_diffusion = E_DIFFUSION;
  if (cp) sscanf(cp, "%lf", &e_diffusion);
  cp = getenv("E_CURV");
  e_curv = E_CURV;
  if (cp) sscanf(cp, "%lf", &e_curv);
  cp = getenv("DT");
  dt = DT;
  if (cp) sscanf(cp, "%lf", &dt);

  fprintf(stderr,
          "mag_scale = %2.1f^%2.1f, e_curv = %2.4f, e_diffusion = "
          "%2.2e, dt=%2.2f\n",
          mag_scale,
          mag_pow,
          e_curv,
          e_diffusion,
          dt);
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (x0 < 0.0f) x0 = (float)(width - 1) / 2.0f;
  if (y0 < 0.0f) y0 = (float)(height - 1) / 2.0f;
  if (z0 < 0.0f) z0 = (float)(depth - 1) / 2.0f;

  if (!mri_distance) mri_distance = MRIbuildDistanceMap(mri_src, NULL, x0, y0, z0, 2.0f);

  min_orig_src = 255;
  max_orig_src = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        src = MRIvox(mri_src, x, y, z);
        dist = MRIFvox(mri_distance, x, y, z);
        if (dist <= 0.0f) {
          if (src > max_orig_src) max_orig_src = src;
          if (src < min_orig_src) min_orig_src = src;
        }
      }
    }
  }

  fprintf(stderr, "min_src = %d, max_src = %d\n", min_orig_src, max_orig_src);
  for (t = 0; t < niter; t++) {
    if (!(t % 10)) fprintf(stderr, "\r%4.4d of %d     ", t, niter);
    if (!(t % write_iter) && (Gdiag & DIAG_WRITE)) {
      char fname[100];
      MRI *mri_interior;

      if (!(t % (write_iter * 10))) {
        sprintf(fname, "dist%d.mnc", t / write_iter);
        MRIwrite(mri_distance, fname);
      }
      sprintf(fname, "front%d.mnc", t / write_iter);
      mri_interior = MRIextractInterior(mri_src, mri_distance, NULL);
      MRIwrite(mri_interior, fname);
      MRIfree(&mri_interior);
    }
    for (z = 0; z < depth; z++) {
      zp1 = mri_src->zi[z + 1];
      zm1 = mri_src->zi[z - 1];
      for (y = 0; y < height; y++) {
        yp1 = mri_src->yi[y + 1];
        ym1 = mri_src->yi[y - 1];
        for (x = 0; x < width; x++) {
          f = MRIFvox(mri_distance, x, y, z);
          dist = (double)f;
          if (dist > MIN_DIST) continue;

          if (x == 23 && y == 9 && z == 24) DiagBreak();
          xp1 = mri_src->xi[x + 1];
          xm1 = mri_src->xi[x - 1];
          dist_xp1 = MRIFvox(mri_distance, xp1, y, z);
          dist_xm1 = MRIFvox(mri_distance, xm1, y, z);
          dist_yp1 = MRIFvox(mri_distance, x, yp1, z);
          dist_ym1 = MRIFvox(mri_distance, x, ym1, z);
          dist_zp1 = MRIFvox(mri_distance, x, y, zp1);
          dist_zm1 = MRIFvox(mri_distance, x, y, zm1);

          dist_xp1yp1 = MRIFvox(mri_distance, xp1, yp1, z);
          dist_xp1zp1 = MRIFvox(mri_distance, xp1, y, zp1);
          dist_yp1zp1 = MRIFvox(mri_distance, x, yp1, zp1);

          /* forward difference */
          dx_f = (dist_xp1 - dist);
          dy_f = (dist_yp1 - dist);
          dz_f = (dist_zp1 - dist);

          /* backward difference */
          dx_b = (dist - dist_xm1);
          dy_b = (dist - dist_ym1);
          dz_b = (dist - dist_zm1);

          /*
             use 'upwind' approximation of the derivatives. That is, calculate
             the derivatives using information from the direction of the front,
             not away from it.
          */
          dx = MAX(dx_b, 0) + MIN(dx_f, 0);
          dy = MAX(dy_b, 0) + MIN(dy_f, 0);
          dz = MAX(dz_b, 0) + MIN(dz_f, 0);

          /* mixed partials */
          dxy = (dist_xp1yp1 + dist - dist_xp1 - dist_yp1) / 2.0f;
          dxz = (dist_xp1zp1 + dist - dist_xp1 - dist_zp1) / 2.0f;
          dyz = (dist_yp1zp1 + dist - dist_yp1 - dist_zp1) / 2.0f;

          /* second order derivatives */
          dxx = dist_xp1 + dist_xm1 - 2.0f * dist;
          dyy = dist_yp1 + dist_ym1 - 2.0f * dist;
          dzz = dist_zp1 + dist_zm1 - 2.0f * dist;
          laplacian = dxx + dyy + dzz;

          denom = dx * dx + dy * dy + dz * dz;
          grad = sqrt(denom);
          denom = pow(denom, 3.0 / 2.0);
          if (DZERO(denom))
            km = 0.0;
          else
            km = ((dyy + dzz) * SQR(dx) + (dxx + dzz) * SQR(dy) + (dxx + dyy) * SQR(dz) -
                  2.0 * (dx * dy * dxy + dx * dz * dxz + dy * dz * dyz)) /
                 denom;
          if (!std::isfinite(km)) km = 0.0;

          /* speed function F based on normal and curvature */
          F = 1.0 - km * e_curv;

          /* modify it by gradient magnitude in source image */
          src = MRIvox(mri_src, x, y, z);
          src_xp1 = MRIvox(mri_src, xp1, y, z);
          src_xm1 = MRIvox(mri_src, xm1, y, z);
          src_yp1 = MRIvox(mri_src, x, yp1, z);
          src_ym1 = MRIvox(mri_src, x, ym1, z);
          src_zp1 = MRIvox(mri_src, x, y, zp1);
          src_zm1 = MRIvox(mri_src, x, y, zm1);

          /* forward difference */
          dx_f = (double)(src_xp1 - src);
          dy_f = (double)(src_yp1 - src);
          dz_f = (double)(src_zp1 - src);

          /* backward difference */
          dx_b = (double)(src - src_xm1);
          dy_b = (double)(src - src_ym1);
          dz_b = (double)(src - src_zm1);

          /*
             use 'upwind' approximation of the derivatives. That is, calculate
             the derivatives using information from the direction of the front,
             not away from it.
          */
          if (dx > 0.0)
            sdx = dx_b;
          else
            sdx = dx_f;
          if (dy > 0.0)
            sdy = dy_b;
          else
            sdy = dy_f;
          if (dz > 0.0) /* information coming from behind us */
            sdz = dz_b;
          else /* information coming from in front of us */
            sdz = dz_f;
          mag = sqrt(sdx * sdx + sdy * sdy + sdz * sdz);
          if (x == 23 && y == 9 && z == 24) {
            static double min_mag = 100000.0;
            if (mag < min_mag) {
              min_mag = mag;
              fprintf(stderr, "min mag set to %2.3f\n", min_mag);
              DiagBreak();
            }
          }
          mag *= mag_scale;
          mag = pow(mag, mag_pow);
          F *= 1.0 / (1.0 + mag);

          /* modify by intensity distance from original seed */
          if (src < min_orig_src)
            orig_mag = 0.0f /* (float)(min_orig_src - src)/5.0f */;
          else if (src > max_orig_src)
            orig_mag = (float)(src - max_orig_src);
          else
            orig_mag = 0.0f;
          orig_mag = pow(orig_mag, 2.0); /* was 8.0 */
          F *= 1.0 / (1.0 + orig_mag);

          dist += dt * (e_diffusion * laplacian - F * grad);
          MRIFvox(mri_distance, x, y, z) = (float)dist;
        }
      }
    }
    MRIupdateDistanceMap(mri_distance);
  }

  if (Gdiag & DIAG_WRITE) {
    char fname[100];
    MRI *mri_interior;

    sprintf(fname, "dist%d.mnc", t / write_iter);
    MRIwrite(mri_distance, fname);
    sprintf(fname, "front%d.mnc", t / write_iter);
    mri_interior = MRIextractInterior(mri_src, mri_distance, NULL);
    MRIwrite(mri_interior, fname);
    MRIfree(&mri_interior);
  }
  fprintf(stderr, "\n");
  return (mri_distance);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
MRI *MRIbuildDistanceMap(MRI *mri_src, MRI *mri_distance, float x0, float y0, float z0, float r)
{
  int width, height, depth, x, y, z;
  float dist, xdist, ydist, zdist, norm;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_distance) {
    mri_distance = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_distance);
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        xdist = (float)x - x0;
        ydist = (float)y - y0;
        zdist = (float)z - z0;
        norm = sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
        if (FZERO(norm))
          dist = -r;
        else
          dist = 1.0f - r / norm;
        MRIFvox(mri_distance, x, y, z) = dist;
      }
    }
  }

  return (mri_distance);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
MRI *MRIextractInterior(MRI *mri_src, MRI *mri_distance, MRI *mri_dst)
{
  int width, height, depth, x, y, z;
  float *pdist, dist;
  BUFTYPE *pdst, *psrc;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdist = &MRIFvox(mri_distance, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      psrc = &MRIvox(mri_src, 0, y, z);
      for (x = 0; x < width; x++, psrc++) {
        dist = *pdist++;
        if (dist <= 0.0f)
          *pdst++ = 140;
        else {
          *pdst++ = MIN(*psrc, 128);
        }
      }
    }
  }

  return (mri_dst);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
MRI *MRIupdateDistanceMap(MRI *mri_distance)
{
  int width, height, depth, x, y, z, xk, yk, zk, x0, y0, z0, xmin, ymin, zmin;
  float dist, min_dist, xdist, ydist, zdist;

  width = mri_distance->width;
  height = mri_distance->height;
  depth = mri_distance->depth;

  xmin = ymin = zmin = 1000; /* to supress compiler warning */
  for (z0 = 0; z0 < depth; z0++) {
    for (y0 = 0; y0 < height; y0++) {
      for (x0 = 0; x0 < width; x0++) {
        dist = MRIFvox(mri_distance, x0, y0, z0);
        if (dist >= MIN_DIST) {
          min_dist = 1000.0f;
          for (zk = -1; zk <= 1; zk++) {
            z = mri_distance->zi[z0 + zk];
            for (yk = -1; yk <= 1; yk++) {
              y = mri_distance->yi[y0 + yk];
              for (xk = -1; xk <= 1; xk++) {
                x = mri_distance->xi[x0 + xk];
                dist = MRIFvox(mri_distance, x, y, z);
                if (dist < min_dist) {
                  xmin = x0 + xk;
                  ymin = y0 + yk;
                  zmin = z0 + zk;
                  min_dist = dist;
                }
              }
            }
          }
          if (min_dist <= MIN_DIST) {
            if (!x0 && !y0 && !z0) DiagBreak();
            xdist = x0 - xmin;
            ydist = y0 - ymin;
            zdist = z0 - zmin;
            dist = min_dist + sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
            MRIFvox(mri_distance, x0, y0, z0) = dist;
          }
        }
      }
    }
  }

  return (mri_distance);
}
