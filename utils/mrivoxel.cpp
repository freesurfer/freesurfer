/**
 * @brief voxel utilities for MRI data structure
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
#include "macros.h"
#include "minc.h"
#include "mri.h"
#include "proto.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

/*!
  \fn float MRIvoxelDx(MRI *mri, int x, int y, int z)
  \brief This computes the intensity gradient in the column direction
  by computing the difference between the face and edge neighbors both
  before and after the target voxel within the slice-plane. The result
  of this funciton will be dependent on the volume orientation and
  volume resolution.
*/
float MRIvoxelDx(MRI *mri, int x, int y, int z)
{
  float left, right, dx;
  int xm1, xp1, ym1, yp1;

  x = mri->xi[x];
  xm1 = mri->xi[x - 1];
  xp1 = mri->xi[x + 1];
  y = mri->yi[y];
  ym1 = mri->yi[y - 1];
  yp1 = mri->yi[y + 1];
  z = mri->zi[z];

  // Sum the values within the slice-plane of the voxels
  // to the "left" of the target voxel, giving the
  // face neighbor twice the weight as the edge neighbors
  left  = 2*MRIgetVoxVal(mri, xm1, y,   z, 0);
  left +=   MRIgetVoxVal(mri, xm1, ym1, z, 0);
  left +=   MRIgetVoxVal(mri, xm1, yp1, z, 0);

  // Same for the voxels to the "right"
  right  = 2*MRIgetVoxVal(mri, xp1, y, z, 0);
  right +=   MRIgetVoxVal(mri, xp1, ym1, z, 0);
  right +=   MRIgetVoxVal(mri, xp1, yp1, z, 0);

  // Take the difference, divide by 8 because of 6 voxels, with two
  // weighted by x2
  dx = ((float)right - (float)left) / 8.0f;
  return (dx);
}
/*!
  \fn float MRIvoxelDy(MRI *mri, int x, int y, int z)
  \brief This computes the intensity gradient in the row direction
  by computing the difference between the face and edge neighbors
  both before and after the target voxel within the slice-plane.
*/
float MRIvoxelDy(MRI *mri, int x, int y, int z)
{
  float top, bottom, dy;
  int xm1, xp1, ym1, yp1;

  x = mri->xi[x];
  xm1 = mri->xi[x - 1];
  xp1 = mri->xi[x + 1];
  y = mri->yi[y];
  ym1 = mri->yi[y - 1];
  yp1 = mri->yi[y + 1];
  z = mri->zi[z];

  top = 2 * MRIgetVoxVal(mri, xm1, ym1, z, 0);
  top += MRIgetVoxVal(mri, x, ym1, z, 0);
  top += MRIgetVoxVal(mri, xp1, ym1, z, 0);

  bottom = 2 * MRIgetVoxVal(mri, xm1, yp1, z, 0);
  bottom += MRIgetVoxVal(mri, x, yp1, z, 0);
  bottom += MRIgetVoxVal(mri, xp1, yp1, z, 0);

  dy = ((float)bottom - (float)top) / 8.0f;
  return (dy);
}
/*!
  \fn float MRIvoxelDz(MRI *mri, int x, int y, int z)
  \brief This computes the intensity gradient in the slice direction
  by computing the difference between the face and edge neighbors
  both before and after the target voxel within the row-plane.
*/
float MRIvoxelDz(MRI *mri, int x, int y, int z)
{
  float dz, top, bottom;
  // int width, depth;
  int xm1, xp1, zm1, zp1;

  // width = mri->width;
  // depth = mri->depth;

  x = mri->xi[x];
  xm1 = mri->xi[x - 1];
  xp1 = mri->xi[x + 1];
  z = mri->zi[z];
  zm1 = mri->zi[z - 1];
  zp1 = mri->zi[z + 1];
  y = mri->yi[y];

  top = 2 * MRIgetVoxVal(mri, xm1, y, zm1, 0);
  top += MRIgetVoxVal(mri, x, y, zm1, 0);
  top += MRIgetVoxVal(mri, xp1, y, zm1, 0);

  bottom = 2 * MRIgetVoxVal(mri, xm1, y, zp1, 0);
  bottom += MRIgetVoxVal(mri, x, y, zp1, 0);
  bottom += MRIgetVoxVal(mri, xp1, y, zp1, 0);

  dz = ((float)bottom - (float)top) / 8.0f;
  return (dz);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIvoxelGradient(MRI *mri, int x, int y, int z, float *pdx, float *pdy, float *pdz)
{
  float mag, dx, dy, dz;

  *pdx = dx = MRIvoxelDx(mri, x, y, z);
  *pdy = dy = MRIvoxelDy(mri, x, y, z);
  *pdz = dz = MRIvoxelDz(mri, x, y, z);
  mag = sqrt(dx * dx + dy * dy + dz * dz);

  return (mag);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIvoxelMean(const MRI *mri, int x0, int y0, int z0, int wsize, int frame)
{
  float mean, total;
  int whalf, width, height, depth, x, y, z, npix, xmin, xmax, ymin, ymax, zmin, zmax;

  whalf = wsize / 2;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  total = 0;
  zmin = MAX(0, z0 - whalf);
  zmax = MIN(depth - 1, z0 + whalf);
  ymin = MAX(0, y0 - whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmin = MAX(0, x0 - whalf);
  xmax = MIN(width - 1, x0 + whalf);
  npix = (zmax - zmin + 1) * (ymax - ymin + 1) * (xmax - xmin + 1);

  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) total += MRIgetVoxVal(mri, x, y, z, frame);
    }
  }
  if (npix > 0)
    mean = (float)total / (float)npix;
  else
    mean = 0.0f;
  return (mean);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIvoxelMin(MRI *mri, int x0, int y0, int z0, int wsize)
{
  float min_val, val;
  int whalf, width, height, depth, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;

  whalf = wsize / 2;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  zmin = MAX(0, z0 - whalf);
  zmax = MIN(depth - 1, z0 + whalf);
  ymin = MAX(0, y0 - whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmin = MAX(0, x0 - whalf);
  xmax = MIN(width - 1, x0 + whalf);

  min_val = MRIgetVoxVal(mri, x0, y0, z0, 0);
  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (val < min_val) min_val = val;
      }
    }
  }
  return (min_val);
}
float MRIvoxelMax(MRI *mri, int x0, int y0, int z0, int wsize)
{
  float max_val, val;
  int whalf, width, height, depth, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;

  whalf = wsize / 2;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  zmin = MAX(0, z0 - whalf);
  zmax = MIN(depth - 1, z0 + whalf);
  ymin = MAX(0, y0 - whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmin = MAX(0, x0 - whalf);
  xmax = MIN(width - 1, x0 + whalf);

  max_val = MRIgetVoxVal(mri, x0, y0, z0, 0);
  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (val > max_val) max_val = val;
      }
    }
  }
  return (max_val);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIvoxelStd(MRI *mri, int x0, int y0, int z0, float mean, int wsize)
{
  float std, total, var;
  int whalf, width, height, depth, x, y, z, npix, xmin, xmax, ymin, ymax, zmin, zmax;
  BUFTYPE *psrc;

  whalf = wsize / 2;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  zmin = MAX(0, z0 - whalf);
  zmax = MIN(depth - 1, z0 + whalf);
  ymin = MAX(0, y0 - whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmin = MAX(0, x0 - whalf);
  xmax = MIN(width - 1, x0 + whalf);
  npix = (zmax - zmin + 1) * (ymax - ymin + 1) * (xmax - xmin + 1);

  total = 0.0;
  for (z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      psrc = &MRIvox(mri, xmin, y, z);
      for (x = xmin; x <= xmax; x++) {
        var = mean - (float)(*psrc++);
        total += var * var;
      }
    }
  }
  if (npix > 0)
    std = sqrt(total) / (float)npix;
  else
    std = 0.0f;

  return (std);
}
/*!
  \fn float MRIvoxelDirection(MRI *mri, int x0, int y0, int z0, int wsize)
  \param wsize - window size, typically 3 (ie, nearest face, edge, and corner neighbors)
*/
#define MAX_WINDOW 7
#define MAX_LEN (MAX_WINDOW * MAX_WINDOW * MAX_WINDOW)
float MRIvoxelDirection(MRI *mri, int x0, int y0, int z0, int wsize)
{
  int whalf, width, height, depth, x, y, z, xmax, ymax, zmax;
  float dx_win[MAX_LEN], dy_win[MAX_LEN], dz_win[MAX_LEN], dx, dy, dz, odx, ody, odz, *pdx, *pdy, *pdz, dir;

  if(wsize > MAX_WINDOW) 
    ErrorReturn(0.0f, (ERROR_BADPARM, "MRIvoxelDirection: window size %d too big", wsize));

  whalf = wsize / 2;
  width  = mri->width;
  height = mri->height;
  depth  = mri->depth;

  zmax = MIN(depth  - 1, z0 + whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmax = MIN(width  - 1, x0 + whalf);

  /* should do something smarter than this about border conditions */
  #if 0
  memset(dx_win, 0, MAX_LEN*sizeof(float)) ;
  memset(dy_win, 0, MAX_LEN*sizeof(float)) ;
  memset(dz_win, 0, MAX_LEN*sizeof(float)) ;
  #endif

  // For each voxel in the window, compute the intensity gradient in each direction
  pdx = dx_win;
  pdy = dy_win;
  pdz = dz_win;
  for (z = z0 - whalf; z <= zmax; z++) {
    for (y = y0 - whalf; y <= ymax; y++) {
      for (x = x0 - whalf; x <= xmax; x++) {
        *pdx++ = MRIvoxelDx(mri, x, y, z);
        *pdy++ = MRIvoxelDy(mri, x, y, z);
        *pdz++ = MRIvoxelDz(mri, x, y, z);
      }
    }
  }

  // Compute the intensity gradient for the target voxel
  odx = MRIvoxelDx(mri, x0, y0, z0);
  ody = MRIvoxelDy(mri, x0, y0, z0);
  odz = MRIvoxelDz(mri, x0, y0, z0);

  // In this loop, xyz are actually delta-xyz whereas they
  // are actually xyz in the loop above.
  pdx = dx_win;
  pdy = dy_win;
  pdz = dz_win;
  dir = 0.0f;
  for (z = -whalf; z <= whalf; z++) {
    for (y = -whalf; y <= whalf; y++) {
      for (x = -whalf; x <= whalf; x++) {
        dx = *pdx++;
        dy = *pdy++;
        dz = *pdz++;
	/* x*odx is dx*di0/dx, where dx=x is the signed distance from
  	   the center voxel to the current voxel and di0/dx=odx is the
  	   itensity gradient WRT x at the center voxel */
	/* dx*odx is (di/dx)*(di0/dx), the product of the intensity
	   gradients WRT x at the current voxel and the center voxel.  */
	/* Not sure what this is actually doing. Something like checking whether
	   the gradients of neibhboring voxels point in the same direction.*/
        dir += (x * odx + y * ody + z * odz) * (dx * odx + dy * ody + dz * odz);
      }
    }
  }
#define FSCALE 1000.0f
  return (dir * (1.0f / (FSCALE * 27.0f)));
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIvoxelGradientDir2ndDerivative(MRI *mri, int x0, int y0, int z0, int wsize)
{
  int whalf;
  // int width, height, depth;
  float odx, ody, odz, len, d2I_dg2, xf, yf, zf, d;
  double val;

  whalf = wsize / 2;
  // width = mri->width;
  // height = mri->height;
  // depth = mri->depth;

  odx = MRIvoxelDx(mri, x0, y0, z0);
  ody = MRIvoxelDy(mri, x0, y0, z0);
  odz = MRIvoxelDz(mri, x0, y0, z0);
  len = sqrt(odx * odx + ody * ody + odz * odz);
  if (FZERO(len)) return (0.0f);
  odx /= len;
  ody /= len;
  odz /= len;

  /* second derivative of intensity in gradient direction */
  MRIsampleVolume(mri, (double)x0, (double)y0, (double)z0, &val);
  d2I_dg2 = (wsize - 1) * val;
  for (d = 1; d <= whalf; d++) {
    xf = x0 + d * odx;
    yf = y0 + d * ody;
    zf = z0 + d * odz;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    d2I_dg2 -= val;

    xf = x0 - d * odx;
    yf = y0 - d * ody;
    zf = z0 - d * odz;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    d2I_dg2 -= val;
  }

  d2I_dg2 /= (float)(wsize);
  return (d2I_dg2);
}

float MRIvoxelZscore(MRI *mri, int x, int y, int z, int wsize)
{
  float mean, std, zscore, src;

  src = (float)MRIvox(mri, x, y, z);
  mean = MRIvoxelMean(mri, x, y, z, wsize, 0);
  std = MRIvoxelStd(mri, x, y, z, mean, wsize);
  if (!FZERO(std))
    zscore = (src - mean) / std;
  else
    zscore = 0.0f;

  return (zscore);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int compare_sort_array(const void *pc1, const void *pc2);
static int compare_sort_array(const void *pc1, const void *pc2)
{
  float c1, c2;

  c1 = *(float *)pc1;
  c2 = *(float *)pc2;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return (1);
  else if (c1 < c2)
    return (-1);

  return (0);
}
float MRIvoxelMedian(MRI *mri, int x0, int y0, int z0, int wsize)
{
  float median, val;
  int whalf, width, height, depth, x, y, z, xmin, xmax, median_index, ymin, ymax, zmin, zmax, wcubed;
  // int npix;
  static float *sort_array = NULL;
  float *sptr;

  whalf = wsize / 2;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  zmin = MAX(0, z0 - whalf);
  zmax = MIN(depth - 1, z0 + whalf);
  ymin = MAX(0, y0 - whalf);
  ymax = MIN(height - 1, y0 + whalf);
  xmin = MAX(0, x0 - whalf);
  xmax = MIN(width - 1, x0 + whalf);
  // npix = (zmax - zmin + 1) * (ymax - ymin + 1) * (xmax - xmin + 1);

  for (sptr = sort_array, z = zmin; z <= zmax; z++) {
    for (y = ymin; y <= ymax; y++) {
      for (x = xmin; x <= xmax; x++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        *sptr++ = val;
      }
    }
  }
  wcubed = (wsize * wsize * wsize);
  median_index = wcubed / 2;
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
#endif
  qsort(sort_array, wcubed, sizeof(float), compare_sort_array);
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
  median = sort_array[median_index];
  return (median);
}
