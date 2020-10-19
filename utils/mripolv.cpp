/**
 * @brief utils for calcing/filtring MRI data based on planes of least variance
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Inquiries: freesurfer@nmr.mgh.harvard.edu
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

#include "box.h"
#include "diag.h"
#include "error.h"
#include "filter.h"
#include "macros.h"
#include "minc.h"
#include "mri.h"
#include "mrimorph.h"
#include "mrisegment.h"
#include "proto.h"
#include "region.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define MAXLEN 256

/*-----------------------------------------------------
                    GLOBAL DATA
-------------------------------------------------------*/

/* vertices of an icosohedron (sp?) */
float ic_x_vertices[NVERTICES] = {0,       0.8944, 0.2764,  -0.7236, -0.7236, 0.2764, -0.4253, -0.8507,
                                  -0.4253, 0.1625, -0.2629, 0.5257,  0.6882,  0.1625, 0.6882,  -0.2629,
                                  -0.5878, 0,      -0.9511, -0.9511, -0.5878, 1.00};

float ic_y_vertices[NVERTICES] = {0,      0,       0.8507,  0.5257,  -0.5257, -0.8507, -0.3090, 0,
                                  0.3090, -0.5000, -0.8090, 0,       -0.5000, 0.5000,  0.5000,  0.8090,
                                  0.8090, 1.0000,  0.3090,  -0.3090, -0.8090, 0.00};

float ic_z_vertices[NVERTICES] = {1.0000, 0.4472, 0.4472, 0.4472, 0.4472, 0.4472, 0.8507, 0.5257,
                                  0.8507, 0.8507, 0.5257, 0.8507, 0.5257, 0.8507, 0.5257, 0.5257,
                                  0,      0,      0,      0,      0,      0.00};

static int vertices_initialized = 0;
static float e1_x_v[NVERTICES];
static float e1_y_v[NVERTICES];
static float e1_z_v[NVERTICES];
static float e2_x_v[NVERTICES];
static float e2_y_v[NVERTICES];
static float e2_z_v[NVERTICES];

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int compare_sort_array(const void *pc1, const void *pc2);
static int compare_sort_farray(const void *pc1, const void *pc2);
static void init_basis_vectors(void);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIpolvMean(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase, total;
  BUFTYPE *pdst, *pptr;

  init_basis_vectors();
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        vertex = *pptr++;
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the mean in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        /* now find the values in this plane */
        total = 0;
        for (yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            total += (float)MRIvox(mri_src, xi, yi, zi);
          }
        }
        *pdst++ = total / n;
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
MRI *MRIpolvNormalCurvature(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int width, height, depth, x, y, z, whalf, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi;
  float nx, ny, nz, mean, var, val, std, *pdst;
  BUFTYPE *pptr;

  init_basis_vectors();
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, MRI_FLOAT, mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst);
  }

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize - 1; /* excludes central point */
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = whalf; y < height; y++) {
      pdst = &MRIFvox(mri_dst, 0, y, z); /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        vertex = *pptr++;
        nx = ic_x_vertices[vertex]; /* normal vector */
        ny = ic_y_vertices[vertex];
        nz = ic_z_vertices[vertex];

        /*
           calculate the mean in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        /* now find the values in the normal direction */
        mean = var = 0.0f;
        for (yk = -whalf; yk <= whalf; yk++) {
          if (!yk) continue; /* skip central point */
          xi = pxi[nint((float)x + (float)yk * nx)];
          yi = pyi[nint((float)y + (float)yk * ny)];
          zi = pzi[nint((float)z + (float)yk * nz)];
          val = (float)MRIvox(mri_src, xi, yi, zi);
          var += val * val;
          mean += val;
        }
        mean /= (float)n;
        val = (float)MRIvox(mri_src, x, y, z);
        std = sqrt(var / (float)n - mean * mean);
        if (!FZERO(std))
          *pdst++ = (val - mean) / std;
        else
          *pdst++ = 0.0f;
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
MRI *MRIpolvZscore(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int width, height, depth, x, y, z, whalf, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi;
  float nx, ny, nz, mean, var, val, std, *pdst, xf, yf, zf;
  BUFTYPE *pptr;

  init_basis_vectors();
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, MRI_FLOAT, mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst);
  }

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = whalf; y < height; y++) {
      pdst = &MRIFvox(mri_dst, 0, y, z); /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        vertex = *pptr++;
        nx = ic_x_vertices[vertex]; /* normal vector */
        ny = ic_y_vertices[vertex];
        nz = ic_z_vertices[vertex];

        /*
           calculate the mean in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        /* now find the values in the normal direction */
        mean = var = 0.0f;
        xf = (float)x - whalf * nx;
        yf = (float)y - whalf * ny;
        zf = (float)z - whalf * nz;
        for (yk = -whalf; yk <= whalf; yk++) {
          xi = pxi[nint(xf)];
          yi = pyi[nint(yf)];
          zi = pzi[nint(zf)];
          val = (float)MRIvox(mri_src, xi, yi, zi);
          var += val * val;
          mean += val;
          xf += nx;
          yf += ny;
          zf += nz;
        }
        mean /= (float)n;
        val = (float)MRIvox(mri_src, x, y, z);
        std = sqrt(var / (float)n - mean * mean);
        if (!FZERO(std))
          *pdst++ = (val - mean) / std;
        else
          *pdst++ = 0.0f;
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
MRI *MRIpolvMedian(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  BUFTYPE *pdst, *pptr, plane_vals[MAXLEN], *pvals;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        vertex = *pptr++;
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the median in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        pvals = plane_vals;

        /* now find the values in this plane */
        for (yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            *pvals++ = (float)MRIvox(mri_src, xi, yi, zi);
          }
        }
        qsort(plane_vals, n, sizeof(BUFTYPE), compare_sort_array);
        *pdst++ = plane_vals[n / 2];
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
MRI *MRIpolvOrder(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, int thresh)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi, order;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  BUFTYPE *pdst, *pptr, plane_vals[MAXLEN], *pvals;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);  /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        vertex = *pptr++;
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the median in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        pvals = plane_vals;

        /* now find the values in this plane */
        for (yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            *pvals++ = (float)MRIvox(mri_src, xi, yi, zi);
          }
        }
        qsort(plane_vals, n, sizeof(BUFTYPE), compare_sort_array);

        /* find the 1st supra-threshold value in the array */
        pvals = plane_vals;
        for (order = 0; order < n; order++)
          if (*pvals++ > thresh) break;
        *pdst++ = (BUFTYPE)order;
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
MRI *MRIpolvMeanRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, MRI_REGION *region)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, x0, y0, z0, xi, yi, zi, *pxi, *pyi, *pzi;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase, total;
  BUFTYPE *pdst, *pptr;

  init_basis_vectors();
  if (mri_src->type != MRI_UCHAR) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIpolvMeanRegion: unsupported src type"));

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
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
  whalf = (wsize - 1) / 2;

  if (mri_dst && mri_dst->type != MRI_UCHAR)
    ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIpolvMeanRegion: unsupported dst type"));

  if (!mri_dst) {
    int w, h, d;

    w = width - region->x;
    h = height - region->y;
    d = depth - region->z;
    mri_dst = MRIalloc(w, h, d, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize;
  }

  n = wsize * wsize;
  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y - y0, z - z0); /* ptr to destination */
      pptr = &MRIvox(mri_polv, x0, y, z);         /* ptr to normal vectors */
      for (x = x0; x < width; x++) {
        vertex = *pptr++;
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the mean in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        /* now find the values in this plane */
        total = 0;
        for (yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            total += (float)MRIvox(mri_src, xi, yi, zi);
          }
        }
        *pdst++ = (BUFTYPE)((float)total / (float)n);
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
MRI *MRIpolvMedianRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, MRI_REGION *region)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, x0, y0, z0, xi, yi, zi, *pxi, *pyi, *pzi, median_index;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  BUFTYPE *pdst, *pptr, plane_vals[MAXLEN], *pvals;

  init_basis_vectors();
  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIpolvMedianRegion: unsupported src type"));

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
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
  whalf = (wsize - 1) / 2;

  if (mri_dst && (mri_dst->type != MRI_UCHAR))
    ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIpolvMedianRegion: unsupported dst type"));
  if (!mri_dst) {
    int w, h, d;

    w = width - region->x;
    h = height - region->y;
    d = depth - region->z;
    mri_dst = MRIalloc(w, h, d, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize;
  }

  n = wsize * wsize;
  median_index = n / 2;
  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y - y0, z - z0); /* ptr to destination */
      pptr = &MRIvox(mri_polv, x0, y, z);         /* ptr to normal vectors */
      for (x = x0; x < width; x++) {
        vertex = *pptr++;
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the median in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */
        pvals = plane_vals;
        for (yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            *pvals++ = (int)MRIvox(mri_src, xi, yi, zi);
          }
        }
        qsort(plane_vals, n, sizeof(BUFTYPE), compare_sort_array);
        *pdst++ = plane_vals[median_index];
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
MRI *MRIextractCpolvCoords(MRI *mri_src, int *px, int *py, int *pz, MRI *mri_polv, int x, int y, int z, int wsize)
{
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int vertex, whalf, xk, yk, xi, yi, zi;

  init_basis_vectors();
  whalf = (wsize - 1) / 2;

  vertex = (int)MRIvox(mri_polv, x, y, z);
  e1_x = e1_x_v[vertex]; /* get basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk * e1_x)];
      yi = mri_src->yi[nint(ybase + xk * e1_y)];
      zi = mri_src->zi[nint(zbase + xk * e1_z)];
      *px++ = xi;
      *py++ = yi;
      *pz++ = zi;
    }
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIextractCpolv(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int x, int y, int z, int wsize)
{
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int vertex, whalf, xk, yk, xi, yi, zi;

  init_basis_vectors();
  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = x - whalf * mri_dst->xsize;
    mri_dst->ystart = y - whalf * mri_dst->ysize;
    mri_dst->zstart = z - whalf * mri_dst->zsize;
    mri_dst->xend = mri_dst->xstart + wsize * mri_dst->xsize;
    mri_dst->yend = mri_dst->ystart + wsize * mri_dst->ysize;
    mri_dst->zend = mri_dst->zstart + wsize * mri_dst->zsize;
    mri_dst->imnr0 = z + mri_src->imnr0;
    mri_dst->imnr1 = mri_dst->imnr0;
  }

  vertex = (int)MRIvox(mri_polv, x, y, z);
  e1_x = e1_x_v[vertex]; /* get basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk * e1_x)];
      yi = mri_src->yi[nint(ybase + xk * e1_y)];
      zi = mri_src->zi[nint(zbase + xk * e1_z)];
      MRIvox(mri_dst, xk + whalf, yk + whalf, 0) = MRIvox(mri_src, xi, yi, zi);
    }
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIextractPolvPlane(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int x, int y, int z, int wsize)
{
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int vertex, whalf, xk, yk, xi, yi, zi;

  init_basis_vectors();
  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = x - whalf * mri_dst->xsize;
    mri_dst->ystart = y - whalf * mri_dst->ysize;
    mri_dst->zstart = z - whalf * mri_dst->zsize;
    mri_dst->xend = mri_dst->xstart + wsize * mri_dst->xsize;
    mri_dst->yend = mri_dst->ystart + wsize * mri_dst->ysize;
    mri_dst->zend = mri_dst->zstart + wsize * mri_dst->zsize;
    mri_dst->imnr0 = z + mri_src->imnr0;
    mri_dst->imnr1 = mri_dst->imnr0;
  }

  vertex = (int)MRIvox(mri_polv, x, y, z);
  e1_x = e1_x_v[vertex]; /* get basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk * e1_x)];
      yi = mri_src->yi[nint(ybase + xk * e1_y)];
      zi = mri_src->zi[nint(zbase + xk * e1_z)];
      MRIvox(mri_dst, xk + whalf, yk + whalf, 0) = MRIvox(mri_src, xi, yi, zi);
    }
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIcountPlanarAboveThreshold(MRI *mri_src, int vertex, int x, int y, int z, int wsize, int lo_lim, int hi_lim)
{
  int whalf, xk, yk, n, xi, yi, zi, *pxi, *pyi, *pzi, count;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  BUFTYPE val;

  whalf = (wsize - 1) / 2;

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  e1_x = e1_x_v[vertex]; /* basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  /*
     calculate the median in the plane orthogonal to (a,b,c),
     through the current origin (x,y,z).
     */

  /* now find the values in this plane */
  for (count = 0, yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = nint(xbase + xk * e1_x);
      xi = pxi[xi];
      yi = nint(ybase + xk * e1_y);
      yi = pyi[yi];
      zi = nint(zbase + xk * e1_z);
      zi = pzi[zi];
      val = MRIvox(mri_src, xi, yi, zi);
      if (val >= lo_lim && val <= hi_lim) count++;
    }
  }

  return (count);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

MRI *MRIplaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, vertex, xk, yk, zk, pno, mini, maxi;
  float min_var, max_var, a, b, c, total[MAXLEN], total_sq[MAXLEN], nv[MAXLEN], varv[MAXLEN], avgv[MAXLEN], val,
      total_var;
  BUFTYPE *psrc, *pdst;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  depth -= whalf; /* don't do outer ring of pixels, so we don't have */
  width -= whalf; /* to deal with boundary conditions */
  height -= whalf;
  for (z = whalf; z < depth; z++) {
    DiagHeartbeat((float)(z - whalf) / (float)(depth - whalf - 1));
    for (y = whalf; y < height; y++) {
      pdst = &MRIvox(mri_dst, whalf, y, z);
      for (x = whalf; x < width; x++) {
        /*
          for this point (x,y,z), go through a set of directions on the unit
          sphere. For each direction, find all the planes orthogonal to that
          dir. within our window, and pick the direction in which the variance
          of all the planes is smallest. This will hopefully be the normal to
          the cortical surface.
          */
        maxi = mini = -1;
        min_var = 100000.0f;  /* minimum variance of set of planes */
        max_var = -100000.0f; /* maximum variance of set of planes */
        if (MRIvox(mri_src, x, y, z) < 50) continue;

        for (vertex = 0; vertex < NVERTICES; vertex++) {
          a = ic_x_vertices[vertex]; /* vector on unit sphere */
          b = ic_y_vertices[vertex];
          c = ic_z_vertices[vertex];

          memset(total, 0, wsize * sizeof(float));
          memset(total_sq, 0, wsize * sizeof(float));
          memset(nv, 0, wsize * sizeof(float));

          /* calculate the variance of all planes orthogonal to (a,b,c) */
          for (zk = -whalf; zk <= whalf; zk++) {
            for (yk = -whalf; yk <= whalf; yk++) {
              psrc = &MRIvox(mri_src, x - whalf, y + yk, z + zk);
              for (xk = -whalf; xk <= whalf; xk++) {
                /* pno is the index (#) of the plane we are in */

                /* find the index of this plane (projection on normal) */
                pno = whalf + nint(zk * c + yk * b + xk * a);
                pno = (pno < 0) ? 0 : (pno > wsize - 1) ? wsize - 1 : pno;
                val = (float)*psrc++;
                total[pno] += val; /* sum of all values in this plane */
                total_sq[pno] += val * val;
                nv[pno]++; /* # of values in this plane */
              }
            }
          }
          total_var = 0.0f;
          for (pno = 0; pno < wsize; pno++) {
            if (!nv[pno])
              varv[pno] = avgv[pno] = 0.0f;
            else {
              avgv[pno] = total[pno] / nv[pno]; /* mean of this plane */
              varv[pno] = total_sq[pno] / nv[pno] - avgv[pno] * avgv[pno];
            }
            total_var += varv[pno]; /* total variance of set of planes */
          }
          total_var /= (float)wsize;
          if (total_var > max_var) {
            max_var = total_var;
            maxi = vertex;
          }
          if (total_var < min_var) {
            min_var = total_var;
            mini = vertex;
          }
        }
        /* done - put vector components into output */
        *pdst++ = (BUFTYPE)mini;
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
MRI *MRIcentralPlaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, vertex, xk, yk, mini, maxi, xi, yi, zi, *pxi, *pyi, *pzi, x1, y1, z1, x0,
      y0, z0;
  float min_mean, min_var, max_var, total, total_sq, nv, varv, avgv, val, background_val, fmax;
  BUFTYPE *pdst, max_val;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  if (getenv("USE_CACHED_CPOLV")) {
    char fname[100];
    MRI *mri_tmp;

    /* try and read previously computed CPOLV file from disk */
    int req = snprintf(fname, 100, "%s/cpolv.mnc", mri_src->fname);
    if( req >= 100 ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    mri_tmp = MRIread(fname);
    if (mri_tmp) {
      if (Gdiag & DIAG_SHOW) fprintf(stderr, "reading previously calculated cpolv %s\n", fname);
      if (mri_tmp->width == mri_src->width && mri_tmp->height == mri_src->height && mri_tmp->depth == mri_src->depth) {
        mri_dst = MRIcopy(mri_tmp, NULL);
        MRIfree(&mri_tmp);
        return (mri_dst);
      }
      MRIfree(&mri_tmp);
    }
  }

  init_basis_vectors();

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  MRIvalRange(mri_src, &background_val, &fmax);
  background_val *= 0.2f; /* anything smaller than 20% of peak is bg */

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

#if 0
  depth -= whalf ;   /* don't do outer ring of pixels, so we don't have */
  width -= whalf ;   /* to deal with boundary conditions */
  height -= whalf ;
#endif

  if (mri_src->roi.dx > 0) {
    x0 = MAX(0, mri_src->roi.x);
    y0 = MAX(0, mri_src->roi.y);
    z0 = MAX(0, mri_src->roi.z);
    x1 = MIN(x0 + mri_src->roi.dx - 1, width - 1);
    y1 = MIN(y0 + mri_src->roi.dy - 1, height - 1);
    z1 = MIN(z0 + mri_src->roi.dz - 1, depth - 1);
  }
  else {
    x0 = y0 = z0 = 0;
    x1 = width - 1;
    y1 = height - 1;
    z1 = depth - 1;
  }

  for (z = z0; z <= z1; z++) {
    DiagHeartbeat((float)(z - z0) / (float)z1);
    for (y = y0; y <= y1; y++) {
      pdst = &MRIvox(mri_dst, x0, y, z);
      for (x = x0; x <= x1; x++) {
        if (MRIvox(mri_src, x, y, z) < 20) /* speed things up */
        {
          *pdst++ = 0;
          continue;
        }
        /*
          for this point (x,y,z), go through a set of directions on the unit
          sphere. For each direction, find all the planes orthogonal to that
          dir. within our window, and pick the direction in which the variance
          of all the planes is smallest. This will hopefully be the normal to
          the cortical surface.
          */
        maxi = mini = -1;
        min_mean = 1000.0f;   /* mean of minimum variance plane */
        min_var = 100000.0f;  /* minimum variance of central planes */
        max_var = -100000.0f; /* maximum variance of central planes */
        pe1_x = e1_x_v;
        pe1_y = e1_y_v;
        pe1_z = e1_z_v;
        pe2_x = e2_x_v;
        pe2_y = e2_y_v;
        pe2_z = e2_z_v;
        for (vertex = 0; vertex < NVERTICES; vertex++) {
          e1_x = *pe1_x++; /* first in-plane basis vector */
          e1_y = *pe1_y++;
          e1_z = *pe1_z++;
          e2_x = *pe2_x++; /* second in-plane basis vector */
          e2_y = *pe2_y++;
          e2_z = *pe2_z++;

          total = total_sq = nv = 0.0f;
          max_val = 0;
          /* now find the values in this plane */
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              val = (float)MRIvox(mri_src, xi, yi, zi);
              total += val;          /* sum of all values in this plane */
              total_sq += val * val; /* sum of squared values in this plane */
              nv++;                  /* # of values in this plane */
            }
          }

          if (!nv)
            varv = avgv = 0.0f;
          else {
            avgv = total / nv; /* mean of this plane */
            varv = total_sq / nv - avgv * avgv;
          }

          if (varv > max_var) {
            max_var = varv;
            maxi = vertex;
          }
          if (varv < min_var) {
            min_var = varv;
            mini = vertex;
            min_mean = avgv;
          }
          if (FZERO(varv)) /* zero variance - won't find anything less */
            break;
        }
#if 0
        /* done - put vector components into output */
        MRIseq_vox(mri_dst, x, y, z, 1) = (BUFTYPE)nint(min_mean) ;
#endif
        *pdst++ = (BUFTYPE)mini;
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
MRI *MRIcentralPlaneOfLeastVarianceNormalMarked(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, vertex, xk, yk, mini, maxi, xi, yi, zi, *pxi, *pyi, *pzi, x1, y1, z1, x0,
      y0, z0;
  float min_mean, min_var, max_var, total, total_sq, nv, varv, avgv, val;

  BUFTYPE *pdst, max_val, *pmask;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  if (mri_src->roi.dx > 0) {
    x0 = MAX(0, mri_src->roi.x);
    y0 = MAX(0, mri_src->roi.y);
    z0 = MAX(0, mri_src->roi.z);
    x1 = MIN(x0 + mri_src->roi.dx - 1, width - 1);
    y1 = MIN(y0 + mri_src->roi.dy - 1, height - 1);
    z1 = MIN(z0 + mri_src->roi.dz - 1, depth - 1);
  }
  else {
    x0 = y0 = z0 = 0;
    x1 = width - 1;
    y1 = height - 1;
    z1 = depth - 1;
  }

  for (z = z0; z <= z1; z++) {
    for (y = y0; y <= y1; y++) {
      pdst = &MRIvox(mri_dst, x0, y, z);
      pmask = &MRIvox(mri_mask, x0, y, z);
      for (x = x0; x <= x1; x++) {
        if (*pmask++ == 0) /* don't do this one */
        {
          *pdst++ = 0;
          continue;
        }
        /*
          for this point (x,y,z), go through a set of directions on the unit
          sphere. For each direction, find all the planes orthogonal to that
          dir. within our window, and pick the direction in which the variance
          of all the planes is smallest. This will hopefully be the normal to
          the cortical surface.
          */
        maxi = mini = -1;
        min_mean = 1000.0f;   /* mean of minimum variance plane */
        min_var = 100000.0f;  /* minimum variance of central planes */
        max_var = -100000.0f; /* maximum variance of central planes */
        pe1_x = e1_x_v;
        pe1_y = e1_y_v;
        pe1_z = e1_z_v;
        pe2_x = e2_x_v;
        pe2_y = e2_y_v;
        pe2_z = e2_z_v;
        for (vertex = 0; vertex < NVERTICES; vertex++) {
          e1_x = *pe1_x++; /* first in-plane basis vector */
          e1_y = *pe1_y++;
          e1_z = *pe1_z++;
          e2_x = *pe2_x++; /* second in-plane basis vector */
          e2_y = *pe2_y++;
          e2_z = *pe2_z++;

          total = total_sq = nv = 0.0f;
          max_val = 0;
          /* now find the values in this plane */
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              val = (float)MRIvox(mri_src, xi, yi, zi);
              total += val;          /* sum of all values in this plane */
              total_sq += val * val; /* sum of squared values in this plane */
              nv++;                  /* # of values in this plane */
            }
          }

          if (!nv)
            varv = avgv = 0.0f;
          else {
            avgv = total / nv; /* mean of this plane */
            varv = total_sq / nv - avgv * avgv;
          }

          if (varv > max_var) {
            max_var = varv;
            maxi = vertex;
          }
          if (varv < min_var) {
            min_var = varv;
            mini = vertex;
            min_mean = avgv;
          }
          if (FZERO(varv)) /* zero variance - won't find anything less */
            break;
        }
        *pdst++ = (BUFTYPE)mini;
      }
    }
  }
  return (mri_dst);
}
static int compare_sort_farray(const void *pc1, const void *pc2)
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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int compare_sort_array(const void *pc1, const void *pc2)
{
  BUFTYPE c1, c2;

  c1 = *(BUFTYPE *)pc1;
  c2 = *(BUFTYPE *)pc2;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return (1);
  else if (c1 < c2)
    return (-1);

  return (0);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
             1st take the cross-product of ez with a random vector
             to obtain one vector in the plane, then take the cross
             product of that vector with the normal (ez) to obtain
             the 2nd in-plane basis vector.
------------------------------------------------------*/
static void init_basis_vectors(void)
{
  float vx, vy, vz, *px, *py, *pz, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e3_x, e3_y, e3_z, e1_x, e1_y, e1_z,
      e2_x, e2_y, e2_z, len;
  int vertex;

  if (vertices_initialized) return;

  px = ic_x_vertices;
  py = ic_y_vertices;
  pz = ic_z_vertices;
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    e3_x = *px; /* vector on unit sphere */
    e3_y = *py;
    e3_z = *pz;
    /*
       now we must scale the length of the vector so that it reaches the
       border of the next voxel. Thus, 'diagonal' vectors must be extended
       by sqrt(2) relative to those which lie along the cardinal axes.
    */
    vx = fabs(e3_x);
    vy = fabs(e3_y);
    vz = fabs(e3_z);            /* use symmetry */
    if ((vx > vy) && (vx > vz)) /* scale using x component */
      len = 1.0f / vx;
    else if (vy > vz) /* scale using y component */
      len = 1.0f / vy;
    else /* scale using z component */
      len = 1.0f / vz;
    *px++ = e3_x * len;
    *py++ = e3_y * len;
    *pz++ = e3_z * len;

    /* pick some other unit non-linearly dependent vector */
    vx = e3_y;
    vy = e3_z;
    vz = e3_x;

    /* don't care about sign (right-hand rule) */
    e1_x = vy * e3_z - vz * e3_y;
    e1_y = vz * e3_x - vx * e3_z;
    e1_z = vx * e3_y - vy * e3_x;

    vx = fabs(e1_x);
    vy = fabs(e1_y);
    vz = fabs(e1_z);            /* use symmetry */
    if ((vx > vy) && (vx > vz)) /* scale using x component */
      len = 1.0f / vx;
    else if (vy > vz) /* scale using y component */
      len = 1.0f / vy;
    else /* scale using z component */
      len = 1.0f / vz;

    e1_x = *pe1_x++ = e1_x * len;
    e1_y = *pe1_y++ = e1_y * len;
    e1_z = *pe1_z++ = e1_z * len;
    len = sqrt(e1_x * e1_x + e1_y * e1_y + e1_z * e1_z);

    e2_x = e1_y * e3_z - e1_z * e3_y;
    e2_y = e1_x * e3_z - e1_z * e3_x;
    e2_z = e1_y * e3_x - e1_x * e3_y;
    vx = fabs(e2_x);
    vy = fabs(e2_y);
    vz = fabs(e2_z);            /* use symmetry */
    if ((vx > vy) && (vx > vz)) /* scale using x component */
      len = 1.0f / vx;
    else if (vy > vz) /* scale using y component */
      len = 1.0f / vy;
    else /* scale using z component */
      len = 1.0f / vz;
    e2_x = *pe2_x++ = e2_x * len;
    e2_y = *pe2_y++ = e2_y * len;
    e2_z = *pe2_z++ = e2_z * len;
#if 0
    DiagFprintf(0L,
                "vertex %d: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f) "
                "x (%2.2f, %2.2f, %2.2f)\n",
                vertex, ic_x_vertices[vertex], ic_y_vertices[vertex],
                ic_z_vertices[vertex], e1_x_v[vertex], e1_y_v[vertex],
                e1_z_v[vertex],e2_x_v[vertex],e2_y_v[vertex],e2_z_v[vertex]) ;
    DiagFprintf(0L, "lengths: %2.3f, %2.3f, %2.3f\n",
                sqrt(e3_x*e3_x+e3_y*e3_y+e3_z*e3_z),
                sqrt(e1_x*e1_x+e1_y*e1_y+e1_z*e1_z),
                sqrt(e2_x*e2_x+e2_y*e2_y+e2_z*e2_z)) ;
#endif
  }
  vertices_initialized = 1;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIfindThinWMStrands(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi, thin, was_white, black_white,
      val;
  float nx, ny, nz, xf, yf, zf;
  BUFTYPE *pdst, *psrc;

  printf("MRIfindThinWMStrands() wsize=%d\n",wsize);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z); /* ptr to destination */
      pdst = &MRIvox(mri_dst, 0, y, z); /* ptr to destination */
      for (x = 0; x < width; x++) {
        thin = 0;
        if (*psrc++)
          for (vertex = 0; !thin && vertex < NVERTICES; vertex++) {
            was_white = -1;
            black_white = -wsize - 1;
            nx = ic_x_vertices[vertex]; /* normal vector */
            ny = ic_y_vertices[vertex];
            nz = ic_z_vertices[vertex];

            xf = (float)x - wsize * nx;
            yf = (float)y - wsize * ny;
            zf = (float)z - wsize * nz;
            for (yk = -wsize; yk <= wsize; yk++) {
              xi = pxi[nint(xf)];
              yi = pyi[nint(yf)];
              zi = pzi[nint(zf)];
              val = (float)MRIvox(mri_src, xi, yi, zi);
              if ((was_white > 0) && !val) /* white to black transition */
              {
                /* check to see if we previously had a b-to-w transition */
                if (black_white >= -wsize) {
                  thin = ((yk - black_white) <= wsize);
                  break;
                }
              }
              else if (!was_white && val > 0) /* black to white transition */
                black_white = yk;
              was_white = val > 0;

              xf += nx;
              yf += ny;
              zf += nz;
            }
          }

        *pdst++ = thin;
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
#define MAX_LABELS 10000
MRI *MRIresegmentThinWMStrands(MRI *mri_src, MRI *mri_dst, int thickness)
{
  int width, height, depth, x, y, z, vertex, thin, i;
  float nx, ny, nz, nd;
  BUFTYPE *pdst, *psrc;
  double val, xf, yf, zf, max_dist, up_dist, down_dist;
  MRI *mri_label;
  MRI_SEGMENTATION *mriseg;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  max_dist = thickness / 2 + 1; /* allow search to extend into non-white */

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  for (vertex = 0; vertex < NVERTICES; vertex++) {
    nx = ic_x_vertices[vertex]; /* normal vector */
    ny = ic_y_vertices[vertex];
    nz = ic_z_vertices[vertex];
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        psrc = &MRIvox(mri_src, 0, y, z); /* ptr to source */
        pdst = &MRIvox(mri_dst, 0, y, z); /* ptr to destination */
        for (x = 0; x < width; x++) {
          if (x == 153 && y == 151 && z == 149) DiagBreak();
          thin = 0;
          if (*psrc++) {
            /* first search 'downwards' to see if we can find non-white */
            down_dist = thickness + 1;
            for (nd = 0; nd <= thickness + 1; nd++) {
              xf = (double)x - nd * nx;
              yf = (double)y - nd * ny;
              zf = (double)z - nd * nz;
              MRIsampleVolume(mri_src, xf, yf, zf, &val);
              if (val < WM_MIN_VAL) {
                down_dist = nd - 1;
                break;
              }
            }

            /* now search 'upwards' to see if we can find non-white */
            up_dist = thickness + 1;
            for (nd = 0; nd <= thickness + 1; nd++) {
              xf = (double)x + nd * nx;
              yf = (double)y + nd * ny;
              zf = (double)z + nd * nz;
              MRIsampleVolume(mri_src, xf, yf, zf, &val);
              if (val < WM_MIN_VAL) {
                up_dist = nd - 1;
                break;
              }
            }
            thin = ((up_dist + down_dist) <= thickness);
          }

          if (thin)
            *pdst++ = thin;
          else
            pdst++;
        }
      }
    }

#if 0
    MRIopen(mri_dst, mri_dst) ;
    MRIclose(mri_dst, mri_dst) ;
    MRIwrite(mri_dst, "thin.mgh") ;
#endif
    fprintf(stderr, "segmenting thin strands...\n");
    mriseg = MRIsegment(mri_dst, 1, 255);

    fprintf(stderr, "%d segments found\n", mriseg->nsegments);
    {
      float max_area = 0.0;
      int max_i = -1;

      mri_label = NULL;
      for (i = 0; i < mriseg->max_segments; i++) {
        if (mriseg->segments[i].area > max_area) {
          max_area = mriseg->segments[i].area;
          max_i = i;
          if (mriseg->segments[i].nvoxels > 100) {
            fprintf(stderr,
                    "segment %3d, area = %2.2f, nvox=%d\n",
                    i,
                    mriseg->segments[i].area,
                    mriseg->segments[i].nvoxels);
            mri_label = MRIsegmentToImage(mri_src, mri_label, mriseg, i);
          }
        }
      }
      if (max_i >= 0) {
        if (!mri_label) mri_label = MRIsegmentToImage(mri_src, NULL, mriseg, i);
        MRIwrite(mri_label, "max_label.mgh");
        MRIfree(&mri_label);
      }
    }

    break;
  }

  return (mri_dst);
}
MRI *MRIfillPlanarHoles(MRI *mri_src, MRI *mri_segment, MRI *mri_dst, MRI_SEGMENT *mseg);
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define MAX_LABELS 10000
MRI *
MRIthickenThinWMStrands(MRI *mri_T1, MRI *mri_src,MRI *mri_dst,int thickness,
                        int nsegments, float wm_hi)
{
  int      width, height, depth, x, y, z, vertex, thin, i, total_filled,
  nfilled, nseg ;
  float    nx, ny, nz, nd ;
  double     val, xf, yf, zf, max_dist, up_dist, down_dist ;
  MRI_SEGMENTATION *mriseg ;
  MRI              *mri_thin ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  max_dist = thickness/2+1 ;  /* allow search to extend into non-white */

  init_basis_vectors() ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  mri_thin = MRIclone(mri_src, NULL) ;

  MRIcopy(mri_src, mri_dst) ;
  total_filled = 0 ;
  for (vertex = 0 /*17*/ ; vertex < NVERTICES ; vertex++)
  {
    MRIclear(mri_thin) ;
    nx = ic_x_vertices[vertex] ;  /* normal vector */
    ny = ic_y_vertices[vertex] ;
    nz = ic_z_vertices[vertex] ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 160 && y == 147 && z == 116)
            DiagBreak() ;
          thin = 0 ;
          if (*psrc++)
          {
            /* first search 'downwards' to see if we can find non-white */
            down_dist = thickness+1 ;
            for (nd = 0 ; nd <= thickness+1 ; nd++)
            {
              xf = (double)x - nd*nx ;
              yf = (double)y - nd*ny ;
              zf = (double)z - nd*nz ;
              MRIsampleVolume(mri_src, xf, yf, zf, &val) ;
              if (val < WM_MIN_VAL)
              {
                down_dist = nd-1 ;
                break ;
              }
            }

            /* now search 'upwards' to see if we can find non-white */
            up_dist = thickness+1 ;
            for (nd = 0 ; nd <= thickness+1 ; nd++)
            {
              xf = (double)x + nd*nx ;
              yf = (double)y + nd*ny ;
              zf = (double)z + nd*nz ;
              MRIsampleVolume(mri_src, xf, yf, zf, &val) ;
              if (val < WM_MIN_VAL)
              {
                up_dist = nd-1 ;
                break ;
              }
            }
            thin =  ((up_dist+down_dist) <= thickness) ;
          }

          if (thin)
            MRIsetVoxVal(mri_thin, x, y, z, 0,  thin) ;
        }
      }
    }

    mriseg = MRIsegment(mri_thin, 1, 255) ;
    nfilled = 0 ;
    for (nseg = 0 ; nseg < nsegments ; nseg++)
    {
      MRI_SEGMENT *mseg ;
      int         v, xd, yd, zd ;

      i = MRIsegmentMax(mriseg) ;
      if (i < 0)
        break ;
      mseg = &mriseg->segments[i] ;
      for (v = 0 ; v < mseg->nvoxels ; v++)
      {
        x = mseg->voxels[v].x ;
        y = mseg->voxels[v].y ;
        z = mseg->voxels[v].z;
        for (nd = -(thickness) ; nd <= thickness ; nd++)
        {
          xd = nint((double)x + nd*nx) ;
          yd = nint((double)y + nd*ny) ;
          zd = nint((double)z + nd*nz) ;
          if (xd == 110 && yd == 125 && zd == 172)  /* T1=148, wm=THICKEN */
            DiagBreak() ;
          if (xd == 126 && yd == 69 && zd == 127)
            DiagBreak() ;
          if ((MRIgetVoxVal(mri_src, xd, yd, zd, 0) <= wm_hi) &&
              (MRIgetVoxVal(mri_dst, xd, yd, zd, 0) == 0) &&
              MRIneighborsOn(mri_dst, xd, yd, zd, WM_MIN_VAL) >= 1)
          {
            nfilled++ ;
            MRIsetVoxVal(mri_dst, xd, yd, zd, 0, THICKEN_FILL) ;
          }
        }
      }
    }
    fprintf(stderr, "orientation %2d of %2d: %2d segments, %d filled\n",
            vertex, NVERTICES-1, nseg, nfilled) ;
    total_filled += nfilled ;
    MRIsegmentFree(&mriseg) ;
  }

  MRIfree(&mri_thin) ;
  fprintf(stderr, "%d voxels filled\n", total_filled) ;
  return(mri_dst) ;
}
#else
/*-----------------------------------------------------
Parameters:

Returns value:

Description
------------------------------------------------------*/
#define MAX_LABELS 10000
#define TOO_THIN 2
MRI *MRIthickenThinWMStrands(MRI *mri_T1, MRI *mri_src, MRI *mri_dst, int thickness, int nsegments, float wm_hi)
{
  int width, height, depth, x, y, z, thin, i, dont_fill, up_added, down_added, total_filled, nfilled, nseg, nx, ny, nz,
      xv, yv, zv, v;
  float nd;
  double val, xf, yf, zf, max_dist, up_dist, down_dist /*, xt, yt, zt*/;
  MRI_SEGMENTATION *mriseg;
  MRI *mri_thin, *mri_tmp;
  printf("MRIthickenThinWMStrands(): thickness=%d, nsegments=%d\n",thickness,nsegments);
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  max_dist = thickness / 2 + 1; /* allow search to extend into non-white */

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);
  mri_thin = MRIclone(mri_src, NULL);
  mri_tmp = MRIremoveIslands(mri_src, NULL, 3, 27-3);
  MRIclose(mri_tmp, mri_tmp);

  MRIcopy(mri_src, mri_dst);
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        thin = 0;
        if (MRIgetVoxVal(mri_src, x, y, z, 0) > 0) /* this is in closed volume to prevent fragmentation */
        {
          for (nz = 0; !thin && nz <= 1; nz++) {
            for (ny = 0; !thin && ny <= 1; ny++) {
              for (nx = 0; !thin && nx <= 1; nx++) {
                if (!nx && !ny && !nz) continue;
                if ((fabs(nx) + fabs(ny) + fabs(nz)) > 1) continue;
                /* hack - only allow superior-inferior thickening */
                if (!ny || nx || nz) continue;

                /* first search 'downwards' to see if we can find non-white */
                down_dist = thickness + 1;
                for (nd = 0; nd <= thickness + 1; nd++) {
                  xf = (double)x - nd * nx;
                  yf = (double)y - nd * ny;
                  zf = (double)z - nd * nz;
                  MRIsampleVolume(mri_tmp, xf, yf, zf, &val);
                  if (val < WM_MIN_VAL) {
                    down_dist = nd - 0.5;
                    break;
                  }
                }

                /* now search 'upwards' to see if we can find non-white */
                up_dist = thickness + 1;
                for (nd = 0; nd <= thickness + 1; nd++) {
                  xf = (double)x + nd * nx;
                  yf = (double)y + nd * ny;
                  zf = (double)z + nd * nz;
                  MRIsampleVolume(mri_tmp, xf, yf, zf, &val);
                  if (val < WM_MIN_VAL) {
                    up_dist = nd - 0.5;
                    break;
                  }
                }
                thin = ((up_dist + down_dist) <= thickness);
              }
            }
          }
        }
        if (thin) MRIsetVoxVal(mri_thin, x, y, z, 0, thin);
      }
    }
  }

  /* now thicken the strand, being careful not to connect with other  strands.  */
  total_filled = 0;
  mriseg = MRIsegment(mri_thin, 1, 255);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[200];
    MRI_SEGMENT *mseg;
    MRI *mri_tmp = NULL;

    for (nseg = 0; nseg < nsegments; nseg++) {
      i = MRIsegmentMax(mriseg); /* find largest remaining segment */
      if (i < 0) break;
      mseg = &mriseg->segments[i];
      mri_tmp = MRIsegmentToImage(mri_src, mri_tmp, mriseg, i);
      sprintf(fname, "seg%d.mgh", i);
      MRIwrite(mri_tmp, fname);
    }
    MRIfree(&mri_tmp);
  }

  int NDILATIONS=5;
  for (i = 0; i < NDILATIONS; i++) MRIsegmentDilate(mriseg, mri_src);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[200];
    MRI *mri_tmp = NULL;
    MRI_SEGMENT *mseg;

    for (nseg = 0; nseg < nsegments; nseg++) {
      i = MRIsegmentMax(mriseg); /* find largest remaining segment */
      if (i < 0) break;
      mseg = &mriseg->segments[i];
      mri_tmp = MRIsegmentToImage(mri_src, mri_tmp, mriseg, i);
      sprintf(fname, "dilated_seg%d.mgh", i);
      MRIwrite(mri_tmp, fname);
    }
    MRIfree(&mri_tmp);
  }

  for (nseg = 0; nseg < nsegments; nseg++) {
    MRI_SEGMENT *mseg;
    /*    int         v, xd, yd, zd ;*/

    i = MRIsegmentMax(mriseg); /* find largest remaining segment */
    if (i < 0) break;          /* no more segments to process */
    mseg = &mriseg->segments[i];
    MRIclear(mri_thin);                              /* mri_thin will now contain the segment */
    MRIsegmentToImage(mri_tmp, mri_thin, mriseg, i); /* build image of seg */

    nfilled = 0;

    for (v = 0; v < mseg->nvoxels; v++) /* for each voxel in segment */
    {
      x = mseg->voxels[v].x;
      y = mseg->voxels[v].y;
      z = mseg->voxels[v].z;
      if (x == 96 && y == 134 && z == 97) DiagBreak();
      if (MRIgetVoxVal(mri_src, x, y, z, 0) < WM_MIN_VAL) continue;

      /* find the direction in which this voxel is thin and thicken it */
      for (nz = 0; nz <= 1; nz++) {
        for (ny = 0; ny <= 1; ny++) {
          for (nx = 0; nx <= 1; nx++) {
            if (!nx && !ny && !nz) continue;
            if ((fabs(nx) + fabs(ny) + fabs(nz)) > 1) continue;
            if (!ny || nx || nz) continue; /* hack - only allow superior-inferior thickening */

            /* first search 'downwards' to see if we can find non-white */
            down_dist = TOO_THIN + 1;
            for (nd = 0; nd <= TOO_THIN + 1; nd++) {
              xf = (double)x - nd * nx;
              yf = (double)y - nd * ny;
              zf = (double)z - nd * nz;
              MRIsampleVolume(mri_dst, xf, yf, zf, &val);
              if (val < WM_MIN_VAL) {
                down_dist = nd - 0.5;
                break;
              }
            }

            /* now search 'upwards' to see if we can find non-white */
            up_dist = TOO_THIN + 1;
            for (nd = 0; nd <= TOO_THIN + 1; nd++) {
              xf = (double)x + nd * nx;
              yf = (double)y + nd * ny;
              zf = (double)z + nd * nz;
              MRIsampleVolume(mri_dst, xf, yf, zf, &val);
              if (val < WM_MIN_VAL) {
                up_dist = nd - 0.5;
                break;
              }
            }
            thin = ((up_dist + down_dist) <= TOO_THIN);
            while (thin > 0) {
              down_added = up_added = 0;
              if (up_dist <= down_dist) {
                xv = nint((double)x + (up_dist + 1.5) * nx);
                yv = nint((double)y + (up_dist + 1.5) * ny);
                zv = nint((double)z + (up_dist + 1.5) * nz);
                // bounds-checking
                if (xv >= mri_dst->width || xv < 0 || yv >= mri_dst->height || yv < 0 || zv >= mri_dst->depth ||
                    zv < 0) {
                  ErrorPrintf(ERROR_BADPARM,
                              "ERROR: MRIthickenThinWMStrands: "
                              "invalid x,y,z!\n");
                  break;
                }
                if ((MRIindexNotInVolume(mri_T1, xv, yv, zv) == 0) && !MRIgetVoxVal(mri_dst, xv, yv, zv, 0)) {
                  up_added = 1;
                  xv = nint((double)x + (up_dist + .5) * nx);
                  yv = nint((double)y + (up_dist + .5) * ny);
                  zv = nint((double)z + (up_dist + .5) * nz);

                  if (MRIgetVoxVal(mri_T1, xv, yv, zv, 0) < wm_hi) {
                    if (xv == 110 && yv == 125 && zv == 172) DiagBreak(); /* T1=148, wm=THICKEN */
                    MRIsetVoxVal(mri_dst, xv, yv, zv, 0, THICKEN_FILL);
                    MRIsetVoxVal(mri_thin, xv, yv, zv, 0, THICKEN_FILL);
                    nfilled++;
                  }
                }
              }
              if (up_dist >= down_dist) {
                xv = nint((double)x - (down_dist + 1.5) * nx);
                yv = nint((double)y - (down_dist + 1.5) * ny);
                zv = nint((double)z - (down_dist + 1.5) * nz);
                // bounds-checking
                if (xv >= mri_dst->width || xv < 0 || yv >= mri_dst->height || yv < 0 || zv >= mri_dst->depth ||
                    zv < 0) {
                  ErrorPrintf(ERROR_BADPARM,
                              "ERROR: MRIthickenThinWMStrands: "
                              "invalid x,y,z!\n");
                  break;
                }
                if (!MRIgetVoxVal(mri_dst, xv, yv, zv, 0)) {
                  down_added = 1;
                  xv = nint((double)x - (down_dist + .5) * nx);
                  yv = nint((double)y - (down_dist + .5) * ny);
                  zv = nint((double)z - (down_dist + .5) * nz);
                  if (MRIgetVoxVal(mri_T1, xv, yv, zv, 0) < wm_hi) {
                    if (xv == 110 && yv == 125 && zv == 172) DiagBreak(); /* T1=148, wm=THICKEN */
                    MRIsetVoxVal(mri_dst, xv, yv, zv, 0, THICKEN_FILL);
                    MRIsetVoxVal(mri_thin, xv, yv, zv, 0, THICKEN_FILL);
                    nfilled++;
                  }
                }
              }
              if (up_added) up_dist += 1.0f;
              if (down_added) down_dist += 1.0f;
              thin = ((up_dist + down_dist) <= TOO_THIN);
              if (!up_added && !down_added) break;
            }
          }
        }
      }
    }

    total_filled += nfilled;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr,
              "%d: segment %d, %d voxels, filled = %d, total = %d\n",
              nseg,
              i,
              mseg->nvoxels,
              nfilled,
              total_filled);

    mseg->ignore = 1;

    /* now that the strand has been thickened some apply a neighborhood
    filter to try to fill some holes.
    */
    do {
      nfilled = 0;
      for (z = 1; z < depth - 1; z++) {
        for (y = 1; y < height - 1; y++) {
          for (x = 1; x < width - 1; x++) {
            if (x == 160 && y == 140 && z == 121) DiagBreak();

            /* check if it should be filled */
            if (!MRIgetVoxVal(mri_dst, x, y, z, 0) && (MRIgetVoxVal(mri_T1, x, y, z, 0) < wm_hi)) {
              if ((MRIgetVoxVal(mri_thin, x + 1, y, z, 0) && MRIgetVoxVal(mri_thin, x - 1, y, z, 0)) ||
                  (MRIgetVoxVal(mri_thin, x, y + 1, z, 0) && MRIgetVoxVal(mri_thin, x, y - 1, z, 0)) ||
                  (MRIgetVoxVal(mri_thin, x, y, z + 1, 0) && MRIgetVoxVal(mri_thin, x, y, z - 1, 0))) {
                dont_fill = 0;
                for (nz = -1; !dont_fill && nz <= 1; nz++) {
                  for (ny = -1; !dont_fill && ny <= 1; ny++) {
                    for (nx = -1; !dont_fill && nx <= 1; nx++) {
                      /* if any neighboring voxels are on the src image
                      that are not part of this segment, don't fill it.
                      */
                      if (MRIgetVoxVal(mri_dst, x + nx, y + ny, z + nz, 0) &&
                          !MRIgetVoxVal(mri_thin, x + nx, y + ny, z + nz, 0))
                        dont_fill = 1;
                    }
                  }
                }

                if (!dont_fill) {
                  if (x == 160 && y == 140 && z == 122) DiagBreak();
                  MRIsetVoxVal(mri_dst, x, y, z, 0, NBHD_FILL);
                  MRIsetVoxVal(mri_thin, x, y, z, 0, NBHD_FILL);
                  nfilled++;
                }
              }
            }
          }
        }
      }
      total_filled += nfilled;
    } while (nfilled > 0);

    MRIfillPlanarHoles(mri_dst, mri_thin, mri_dst, mseg);
  }
  printf("  %2d segments, %d filled\n", nseg, total_filled);
  MRIsegmentClearIgnoreFlags(mriseg);

  MRIsegmentFree(&mriseg);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_thin, "thin.mgh");
  MRIfree(&mri_thin);
  MRIfree(&mri_tmp);
  return (mri_dst);
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIorderThreshold(MRI *mri_src, MRI *mri_dst, MRI *mri_order, int num)
{
  int width, height, depth, x, y, z, frame;
  BUFTYPE *psrc, *pdst, *porder, val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  if (mri_src->type != MRI_UCHAR || mri_dst->type != MRI_UCHAR || mri_order->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIorderThresh: all input MRs must be byte"));

  for (frame = 0; frame < mri_src->nframes; frame++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        psrc = &MRIseq_vox(mri_src, 0, y, z, frame);
        pdst = &MRIseq_vox(mri_dst, 0, y, z, frame);
        porder = &MRIseq_vox(mri_order, 0, y, z, frame);
        for (x = 0; x < width; x++, psrc++) {
          if (*porder++ >= num)
            val = *psrc;
          else
            val = 0;
          *pdst++ = val;
        }
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
MRI *MRIpolvCount(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize, int lo_lim, int hi_lim)
{
  int width, height, depth, x, y, z, whalf, xk, yk, n, vertex, xi, yi, zi, *pxi, *pyi, *pzi, order;
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  float val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;
  n = wsize * wsize;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        vertex = (int)MRIgetVoxVal(mri_polv, x, y, z, 0);
        e1_x = e1_x_v[vertex]; /* basis vectors for plane */
        e1_y = e1_y_v[vertex];
        e1_z = e1_z_v[vertex];
        e2_x = e2_x_v[vertex];
        e2_y = e2_y_v[vertex];
        e2_z = e2_z_v[vertex];

        /*
           calculate the median in the plane orthogonal to (a,b,c),
           through the current origin (x,y,z).
           */

        /* now find the values in this plane */
        for (order = 0, yk = -whalf; yk <= whalf; yk++) {
          xbase = (float)x + (float)yk * e2_x;
          ybase = (float)y + (float)yk * e2_y;
          zbase = (float)z + (float)yk * e2_z;
          for (xk = -whalf; xk <= whalf; xk++) {
            /* in-plane vect. is linear combination of scaled basis vects */
            xi = nint(xbase + xk * e1_x);
            xi = pxi[xi];
            yi = nint(ybase + xk * e1_y);
            yi = pyi[yi];
            zi = nint(zbase + xk * e1_z);
            zi = pzi[zi];
            val = MRIvox(mri_src, xi, yi, zi);
            if (val >= lo_lim && val <= hi_lim) order++;
          }
        }

        MRIsetVoxVal(mri_dst, x, y, z, 0, order);
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

#define WHITE_LOW 90
#define GRAY_HI 95
#define WHITE_HI 130
#define WSIZE 9

#define LABEL_NOTWHITE 0
#define LABEL_UNKNOWN 255
#define LABEL_AMBIGUOUS LABEL_UNKNOWN

#define PSLOPE 1.0
#define NSLOPE 1.0

#define AMBIGUOUS_PCT 0.8
#define N1_SIZE 27
#define AMBIGUOUS_THRESH (N1_SIZE * AMBIGUOUS_PCT)
#define N2_SIZE (5 * 5)
#define REVERSE_PCT 0.6
#define REVERSE_THRESH (N2_SIZE * REVERSE_PCT)

#define DEBUG_POINT(x, y, z) (((x) == 75) && ((y) == 96) && ((z) == 127))
MRI *MRIwmfilter(MRI *mri_src, MRI *mri_polv, MRI *mri_dst, float nslope, float pslope)
{
  int width, height, depth, x, y, z, whalf, vertex, xi, yi, zi, xo, yo, zo, *pxi, *pyi, *pzi, i, nwhite, nblack, count;
  float nx, ny, nz, dx, dy, dz, curv;
  BUFTYPE *pdst, *pptr, val0, /* *psrc, */ gray_hi, white_low /*,mean, *pmean*/, *plabel, *psrc, l;
  MRI *mri_curv, *mri_label /*, *mri_tmp*/;

  if (FZERO(pslope)) pslope = PSLOPE;
  if (FZERO(nslope)) nslope = NSLOPE;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (DIAG_VERBOSE_ON) mri_curv = MRIalloc(width, height, depth, MRI_FLOAT);
  mri_label = MRIalloc(width, height, depth, MRI_UCHAR);
  /*  mri_tmp = MRIalloc(width, height, depth, MRI_UCHAR) ;*/
  whalf = (WSIZE - 1) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  /*
     first do a preliminary intensity-based segmentation and put
     it into the mri_label volume.
     */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      plabel = &MRIvox(mri_label, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      pptr = &MRIvox(mri_polv, 0, y, z); /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();

        /* calculate through-polv curvature */
        vertex = *pptr++;
        nx = ic_x_vertices[vertex];
        ny = ic_y_vertices[vertex];
        nz = ic_z_vertices[vertex];
        val0 = *psrc++;

#if 1
        /* now compute the curvature in the normal direction */
        for (curv = 0.0f, i = -whalf; i <= whalf; i++) {
          double rval;
          if (!i) continue;
          dx = (float)i * nx;
          dy = (float)i * ny;
          dz = (float)i * nz;
          MRIsampleVolume(mri_src, x + dx, y + dy, z + dz, &rval);
          if (rval > val0)
            curv += 1.0f;
          else if (rval < val0)
            curv -= 1.0f;
        }
        curv /= (WSIZE - 1);
#else
        /* now compute the curvature in the normal direction */
        for (curv = 0.0f, i = -whalf; i <= whalf; i++) {
          if (!i) continue;
          dx = (float)i * nx;
          dy = (float)i * ny;
          dz = (float)i * nz;
#if 0
          xi = pxi[x+nint(dx)] ;
          yi = pyi[y+nint(dy)] ;
          zi = pzi[z+nint(dz)] ;
          val = (float)MRIvox(mri_src, xi, yi, zi) ;
#else
          {
            double rval;
            MRIsampleVolume(mri_src, x + dx, y + dy, z + dz, &rval);
            val = (float)rval;
            if (val < 70) val = 70;
          }
#endif
          dsq = sqrt(dx * dx + dy * dy + dz * dz);
          val -= (float)val0;
          curv += val / dsq;
        }
        curv /= (WSIZE - 1);
#endif
        if (DIAG_VERBOSE_ON) MRIFvox(mri_curv, x, y, z) = curv;

        white_low = WHITE_LOW;
        gray_hi = GRAY_HI;
        if (curv < 0.0f) /* gyrus */
        {
          white_low += nint(nslope * curv);
          gray_hi += nint(nslope * curv);
        }
        else if (curv > 0.0f) {
          gray_hi += nint(pslope * curv);
          white_low += nint(pslope * curv);
        }

        /* mark it as ambiguous or not in the mri_dst buffer */
        if (val0 >= white_low && val0 <= gray_hi)
          *pdst++ = LABEL_AMBIGUOUS;
        else if (val0 > WHITE_HI || val0 < white_low)
          *pdst++ = LABEL_NOTWHITE;
        else
          *pdst++ = val0;

        if (val0 > WHITE_HI) /* too high to be white matter */
          val0 = 0;
        else if (val0 < white_low) /* too low to be white matter */
          val0 = 0;
#if 0
        else if (val0 < gray_hi)    /* ambiguous */
        {
          if ((gray_hi - val0) > (val0 - white_low))
            val0 = 0 ;
        }
#endif

        *plabel++ = val0;
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_label, "label.mnc");

  /* now go through the tentatively labeled volume and definitively label
     any points which lie in an unambiguous neighborhood and put them into
     the mri_dst volume. Ambiguous voxels will be labeled as LABEL_AMBIGUOUS.
     */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z); /* ptr to destination */
      psrc = &MRIvox(mri_src, 0, y, z);
      plabel = &MRIvox(mri_label, 0, y, z);
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        l = *pdst;
        if (l != LABEL_AMBIGUOUS) {
          psrc++;
          pdst++;
          plabel++;
          continue;
        }
        nwhite = nblack = 0;
        for (zo = -1; zo <= 1; zo++) {
          zi = pzi[z + zo];
          for (yo = -1; yo <= 1; yo++) {
            yi = pyi[y + yo];
            for (xo = -1; xo <= 1; xo++) {
              xi = pxi[x + xo];
              l = MRIvox(mri_label, xi, yi, zi);
              if (l)
                nwhite++;
              else
                nblack++;
            }
          }
        }
        val0 = *psrc++;
        l = *plabel++;
        if (l && nwhite > AMBIGUOUS_THRESH)
          l = val0;
        else if (!l && nblack > AMBIGUOUS_THRESH)
          l = LABEL_NOTWHITE;
        else
          l = LABEL_AMBIGUOUS;
        *pdst++ = l;
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_dst, "dst1.mnc");
  /*
    now go through again, and for any voxel which was labeled ambiguous,
    look at the number of in-plane labels. If enough of them are different
    from the prior labeling, change the label.
    */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);     /* ptr to destination */
      plabel = &MRIvox(mri_label, 0, y, z); /* ptr to destination */
      pptr = &MRIvox(mri_polv, 0, y, z);    /* ptr to normal vectors */
      psrc = &MRIvox(mri_src, 0, y, z);     /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        l = *pdst;
        if (l != LABEL_AMBIGUOUS) {
          pdst++;
          pptr++;
          plabel++;
          psrc++;
          continue;
        }

        vertex = *pptr++;
        val0 = *plabel++;
        nx = ic_x_vertices[vertex];
        ny = ic_y_vertices[vertex];
        nz = ic_z_vertices[vertex];
        count = MRIcountPlanarAboveThreshold(mri_label, vertex, x, y, z, 5, 1, 255);
        if (val0 && count < N2_SIZE - REVERSE_THRESH)
          val0 = 0;
        else if (!val0 && count > REVERSE_THRESH)
          val0 = 0;
        *pdst++ = val0;
      }
    }
  }

#if 0
  /*
    now go through again, and for any voxel which is labeled differently
    than a majority of the other in-plane labels, use the curvature
    to reclassify it.
    */
  for (j = 0 ; j < 5 ; j++)
  {
    int examined, changed ;

    changed = examined = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pdst = &MRIvox(mri_dst, 0, y, z) ;  /* ptr to destination */
        pptr = &MRIvox(mri_polv, 0, y, z) ; /* ptr to normal vectors */
        psrc = &MRIvox(mri_src, 0, y, z) ; /* ptr to normal vectors */
        for (x = 0 ; x < width ; x++)
        {
          if (DEBUG_POINT(x,y,z))
            DiagBreak() ;
          vertex = *pptr++ ;
          if (*psrc++ < 30)
          {
            pdst++ ;
            continue ;
          }
          val0 = *pdst ;

          nx = ic_x_vertices[vertex] ;
          ny = ic_y_vertices[vertex] ;
          nz = ic_z_vertices[vertex] ;
          count = MRIcountPlanarAboveThreshold(mri_dst, vertex,x,y,z,5,
                                               3,255) ;
          if (!val0)
            count = 5*5 - count ;  /* number off */
          if (count < 5*5*.6)  /* a significant # have a differing opinion */
          {
            examined++ ;

            /* now compute the curvature in the normal direction */
            for (curv = 0.0f, i = -whalf ; i <= whalf ; i++)
            {
              if (!i)
                continue ;
              dx = (float)i * nx ;
              xi = pxi[x+nint(dx)] ;
              dy = (float)i * ny ;
              yi = pyi[y+nint(dy)] ;
              dz = (float)i * nz ;
              zi = pzi[z+nint(dz)] ;
              dsq = sqrt(dx*dx + dy*dy + dz*dz) ;
              val = (float)MRIvox(mri_src, xi, yi, zi) - (float)val0 ;
              curv += val / dsq ;
            }
            curv /= (WSIZE-1) ;
            if (!val0 && curv < 0)
              val0 = 110, changed++ ;
            else if (val0 && curv > 0)
              val0 = 0, changed++ ;

          }

          *pdst++ = val0 ;
        }
      }
    }
    if (Gdiag & DIAG_SHOW && examined > 0)
      fprintf(stderr, "%d examined (%2.1f%%), %d changed (%2.2f%%)\n",
              examined, 100.0f*(float)examined/(float)(width*height*depth),
              changed, 100.0f*(float)changed/(float)examined) ;
  }
#endif

  if (DIAG_VERBOSE_ON) {
    MRIwrite(mri_curv, "curv.mnc");
    MRIfree(&mri_curv);
  }
  MRIfree(&mri_label);
  return (mri_dst);
}
MRI *MRIwmfilterMarked(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int wsize, float pct, int onoff)
{
  int width, height, depth, x, y, z, whalf, vertex, num_on, num_off, num, changed, low_thresh;
  float thresh;
  BUFTYPE *pdst, *pmask, *psrc, l;

  thresh = wsize * wsize * pct;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;
  low_thresh = (wsize * wsize) / 2;

  init_basis_vectors();
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  /*
    now go through again, and for any voxel which was labeled ambiguous,
    look at the number of in-plane labels. If enough of them are different
    from the prior labeling, change the label.
    */
  for (num = changed = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);   /* ptr to destination */
      pmask = &MRIvox(mri_mask, 0, y, z); /* ptr to destination */
      psrc = &MRIvox(mri_src, 0, y, z);   /* ptr to normal vectors */
      for (x = 0; x < width; x++) {
        if (DEBUG_POINT(x, y, z)) DiagBreak();
        l = *pmask++;
        if (l == 0) /* don't check it */
        {
          *pdst++ = *psrc++;
          continue;
        }

        vertex = MRIcountCpolvAtVoxel(mri_src, x, y, z, wsize, &num_on, l);
        if (!vertex) DiagBreak();
        vertex = MRIcountCpolvAtVoxel(mri_src, x, y, z, wsize, &num_off, 0);
        if (!vertex) DiagBreak();

        /*
          we are only considering voxels that are off (which may represent
          white matter or non-white) which are in a neighborhood that is
          predominantly on. Only change it to be on if there are more on
          then off and the number of on exceeds a threshold.
          */
        num++;
        if ((num_on > thresh) && (num_off < low_thresh)) {
          changed++;
          *pdst++ = l;
          psrc++;
        }
        else
          *pdst++ = *psrc++;
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d holes examined, %d changed (%2.1f%%)\n", num, changed, 100 * (float)changed / (float)num);
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIcountCpolvAtVoxel(MRI *mri_src, int x, int y, int z, int wsize, int *pnum, int label_to_check)
{
  int whalf, vertex, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex, max_count, num;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  max_count = peak_vertex = 0;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    num = 0;
    e1_x = *pe1_x++; /* first in-plane basis vector */
    e1_y = *pe1_y++;
    e1_z = *pe1_z++;
    e2_x = *pe2_x++; /* second in-plane basis vector */
    e2_y = *pe2_y++;
    e2_z = *pe2_z++;

    /* now find the values in this plane */
    for (yk = -whalf; yk <= whalf; yk++) {
      xbase = (float)x + (float)yk * e2_x;
      ybase = (float)y + (float)yk * e2_y;
      zbase = (float)z + (float)yk * e2_z;
      for (xk = -whalf; xk <= whalf; xk++) {
        /* in-plane vect. is linear combination of scaled basis vects */
        xi = nint(xbase + xk * e1_x);
        xi = pxi[xi];
        yi = nint(ybase + xk * e1_y);
        yi = pyi[yi];
        zi = nint(zbase + xk * e1_z);
        zi = pzi[zi];
        label = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
        if (label == label_to_check) num++;
      }
    }
    if (num >= max_count) {
      peak_vertex = vertex;
      max_count = num;
    }
  }

  if (pnum) *pnum = max_count;

  return (peak_vertex);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIcountCpolvOnAtVoxel(MRI *mri_src, int x, int y, int z, int wsize, int *pnum)
{
  int whalf, vertex, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex, max_count, num;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  max_count = peak_vertex = 0;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    num = 0;
    e1_x = *pe1_x++; /* first in-plane basis vector */
    e1_y = *pe1_y++;
    e1_z = *pe1_z++;
    e2_x = *pe2_x++; /* second in-plane basis vector */
    e2_y = *pe2_y++;
    e2_z = *pe2_z++;

    /* now find the values in this plane */
    for (yk = -whalf; yk <= whalf; yk++) {
      xbase = (float)x + (float)yk * e2_x;
      ybase = (float)y + (float)yk * e2_y;
      zbase = (float)z + (float)yk * e2_z;
      for (xk = -whalf; xk <= whalf; xk++) {
        /* in-plane vect. is linear combination of scaled basis vects */
        xi = nint(xbase + xk * e1_x);
        xi = pxi[xi];
        yi = nint(ybase + xk * e1_y);
        yi = pyi[yi];
        zi = nint(zbase + xk * e1_z);
        zi = pzi[zi];
        label = MRIvox(mri_src, xi, yi, zi);
        if (label) num++;
      }
    }
    if (num >= max_count) {
      peak_vertex = vertex;
      max_count = num;
    }
  }

  if (pnum) *pnum = max_count;

  return (peak_vertex);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIcountCpolvOffAtVoxel(MRI *mri_src, int x, int y, int z, int wsize, int *pnum)
{
  int whalf, vertex, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex, max_count, num;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  peak_vertex = max_count = 0;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    num = 0;
    e1_x = *pe1_x++; /* first in-plane basis vector */
    e1_y = *pe1_y++;
    e1_z = *pe1_z++;
    e2_x = *pe2_x++; /* second in-plane basis vector */
    e2_y = *pe2_y++;
    e2_z = *pe2_z++;

    /* now find the values in this plane */
    for (yk = -whalf; yk <= whalf; yk++) {
      xbase = (float)x + (float)yk * e2_x;
      ybase = (float)y + (float)yk * e2_y;
      zbase = (float)z + (float)yk * e2_z;
      for (xk = -whalf; xk <= whalf; xk++) {
        /* in-plane vect. is linear combination of scaled basis vects */
        xi = nint(xbase + xk * e1_x);
        xi = pxi[xi];
        yi = nint(ybase + xk * e1_y);
        yi = pyi[yi];
        zi = nint(zbase + xk * e1_z);
        zi = pzi[zi];
        label = MRIvox(mri_src, xi, yi, zi);
        if (!label) num++;
      }
    }
    if (num >= max_count) {
      peak_vertex = vertex;
      max_count = num;
    }
  }

  if (pnum) *pnum = max_count;

  return (peak_vertex);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIcpolvThreshold(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int wm_low, int gray_hi, int wsize)
{
  int width, height, depth, x, y, z, whalf, vertex, xk, yk, xi, yi, zi, *pxi, *pyi, *pzi, label, nlabeled, num_white,
      num_ambiguous, num_non_white;
  BUFTYPE *pdst, src, thresh;
  float xbase, ybase, zbase, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  pxi = mri_labeled->xi;
  pyi = mri_labeled->yi;
  pzi = mri_labeled->zi;

  thresh = (wm_low + gray_hi) / 2;
  width = mri_labeled->width;
  height = mri_labeled->height;
  depth = mri_labeled->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_labeled, NULL);

  nlabeled = 0;
  for (z = 0; z < depth; z++) {
    DiagHeartbeat((float)(z) / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        label = MRIvox(mri_labeled, x, y, z);
        if (label == MRI_AMBIGUOUS) {
          nlabeled++;
          if (!FZERO(e2_z_v[0])) DiagBreak();
#if 0
          vertex = MRIneighborhoodPlanarDirection(mri_labeled,x,y,z,3,5);
#else
          vertex = MRIcentralPlaneOfLeastVarianceNormalVoxel(mri_src, wsize, x, y, z);
#endif
          e1_x = e1_x_v[vertex]; /* first in-plane basis vector */
          e1_y = e1_y_v[vertex];
          e1_z = e1_z_v[vertex];
          e2_x = e2_x_v[vertex]; /* second in-plane basis vector */
          e2_y = e2_y_v[vertex];
          e2_z = e2_z_v[vertex];

          /* now vote in this plane */
          num_white = num_ambiguous = num_non_white = 0;
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /*
                in-plane vect. is linear combination of scaled basis vects
                */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              src = MRIvox(mri_src, xi, yi, zi);
              if (src >= thresh)
                num_white++;
              else
                num_non_white++;
            }
          }

          if (num_white >= num_non_white)
            label = MRI_WHITE;
          else
            label = MRI_NOT_WHITE;
        }
        *pdst++ = label;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,
            "              %8d voxels processed (%2.2f%%)\n",
            nlabeled,
            100.0f * (float)nlabeled / (float)(width * height * depth));
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIcpolvVote(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int wsize, int niter, int use_all)
{
  int i, width, height, depth, x, y, z, whalf, vertex, xk, yk, xi, yi, zi, *pxi, *pyi, *pzi, label, nvox, nlabeled,
      num_white, num_ambiguous, num_non_white;
  BUFTYPE *pdst;
  float xbase, ybase, zbase, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;
  MRI *mri_tmp = NULL;

  init_basis_vectors();

  pxi = mri_labeled->xi;
  pyi = mri_labeled->yi;
  pzi = mri_labeled->zi;

  width = mri_labeled->width;
  height = mri_labeled->height;
  depth = mri_labeled->depth;
  whalf = (wsize - 1) / 2;

  if (!mri_dst) mri_dst = MRIclone(mri_labeled, NULL);
  if (niter > 1) {
    mri_tmp = MRIcopy(mri_labeled, NULL);
    mri_labeled = mri_tmp;
  }

  for (i = 0; i < niter; i++) {
    nlabeled = nvox = 0;
    for (z = 0; z < depth; z++) {
      DiagHeartbeat((float)(z) / (float)(depth - 1));
      for (y = 0; y < height; y++) {
        pdst = &MRIvox(mri_dst, 0, y, z);
        for (x = 0; x < width; x++) {
          label = MRIvox(mri_labeled, x, y, z);
          if (label == MRI_WHITE) DiagBreak();
          if (label == MRI_AMBIGUOUS) {
            nvox++;
            if (!FZERO(e2_z_v[0])) DiagBreak();
#if 0
            vertex = MRIneighborhoodPlanarDirection(mri_labeled,x,y,z,3,5);
#else
            vertex = MRIcentralPlaneOfLeastVarianceNormalVoxel(mri_src, wsize, x, y, z);
#endif
            e1_x = e1_x_v[vertex]; /* first in-plane basis vector */
            e1_y = e1_y_v[vertex];
            e1_z = e1_z_v[vertex];
            e2_x = e2_x_v[vertex]; /* second in-plane basis vector */
            e2_y = e2_y_v[vertex];
            e2_z = e2_z_v[vertex];

            /* now vote in this plane */
            num_white = num_ambiguous = num_non_white = 0;
            for (yk = -whalf; yk <= whalf; yk++) {
              xbase = (float)x + (float)yk * e2_x;
              ybase = (float)y + (float)yk * e2_y;
              zbase = (float)z + (float)yk * e2_z;
              for (xk = -whalf; xk <= whalf; xk++) {
                /*in-plane vect. is linear combination of scaled basis vects */
                xi = nint(xbase + xk * e1_x);
                xi = pxi[xi];
                yi = nint(ybase + xk * e1_y);
                yi = pyi[yi];
                zi = nint(zbase + xk * e1_z);
                zi = pzi[zi];
                label = MRIvox(mri_labeled, xi, yi, zi);
                switch (label) {
                  default:
                  case MRI_WHITE:
                    num_white++;
                    break;
                  case MRI_AMBIGUOUS:
                    num_ambiguous++;
                    break;
                  case MRI_NOT_WHITE:
                    num_non_white++;
                    break;
                }
              }
            }

            if (!use_all) {
              if (num_white >= num_non_white && num_white > num_ambiguous)
                label = MRI_WHITE;
              else if (num_non_white > num_ambiguous)
                label = MRI_NOT_WHITE;
              else
                label = MRI_AMBIGUOUS;
            }
            else {
              if (num_white >= num_non_white)
                label = MRI_WHITE;
              else
                label = MRI_NOT_WHITE;
            }
            if (label != MRI_AMBIGUOUS) nlabeled++;
          }
          *pdst++ = label;
        }
      }
    }
    if (i < niter - 1) MRIcopy(mri_dst, mri_labeled);
    if (Gdiag & DIAG_SHOW) {
      fprintf(stderr,
              "                %8d voxels processed (%2.2f%%)\n",
              nvox,
              100.0f * (float)nvox / (float)(width * height * depth));
      fprintf(stderr,
              "              %8d voxels labeled (%2.2f%%)\n",
              nlabeled,
              100.0f * (float)nlabeled / (float)(width * height * depth));
    }
  }
  if (mri_tmp) MRIfree(&mri_tmp);
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIcpolvSmooth(MRI *mri_orig, MRI *mri_src, MRI *mri_dst, int wsize, int low_val, int hi_val, int niter)
{
  int width, height, depth, x, y, z, label, num_white, total_vox, white_vertex, black_vertex, x0, y0, z0, xi, yi, zi,
      nskipped, dst_label, num_black, nwhite_to_black, nblack_to_white, i, skip;
  BUFTYPE *pdst, *psrc, *porig, orig;
  float thresh, hi_thresh, sthresh;
  MRI *mri_tmp;

  init_basis_vectors();

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  thresh = (float)(wsize * wsize) / 3.0f;
  hi_thresh = (float)(wsize * wsize) - thresh;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);
  if (niter > 1)
    mri_tmp = MRIcopy(mri_src, NULL);
  else
    mri_tmp = mri_src;
  mri_src = mri_tmp; /* don't overwrite input */

  sthresh = 3 * 3 * 3 * 0.6666f;
  for (i = 0; i < niter; i++) {
    nskipped = total_vox = nwhite_to_black = nblack_to_white = 0;
    for (z = 0; z < depth; z++) {
      DiagHeartbeat((float)(z) / (float)(depth - 1));
      for (y = 0; y < height; y++) {
        pdst = &MRIvox(mri_dst, 0, y, z);
        psrc = &MRIvox(mri_src, 0, y, z);
        porig = &MRIvox(mri_orig, 0, y, z);
        for (x = 0; x < width; x++) {
          if (x == 10 && y == 43 && z == 13) DiagBreak();
          if (x == 10 && y == 44 && z == 10) DiagBreak();
          if (x == 107 && y == 147 && z == 127) DiagBreak();
          dst_label = *psrc++;

          total_vox++;
          orig = *porig++;
          skip = (orig < low_val) || (orig > hi_val);
          num_black = num_white = 0;
          for (z0 = z - 1; !skip && z0 <= z + 1; z0++) {
            zi = mri_src->zi[z0];
            for (y0 = y - 1; !skip && y0 <= y + 1; y0++) {
              yi = mri_src->yi[y0];
              for (x0 = x - 1; !skip && x0 <= x + 1; x0++) {
                xi = mri_src->xi[x0];
                label = MRIvox(mri_src, xi, yi, zi);
                if (label == MRI_WHITE)
                  ++num_white;
                else
                  ++num_black;

                if (dst_label == MRI_WHITE) {
                  if ((float)num_white > sthresh) skip = 1;
                }
                else {
                  if ((float)num_black > sthresh) skip = 1;
                }
              }
            }
          }
          if (skip) /* not in an ambiguous region */
          {
            *pdst++ = dst_label;
            nskipped++;
            continue;
          }
          if (dst_label == MRI_WHITE) {
            white_vertex = MRIneighborhoodCpolv(mri_src, x, y, z, 3, wsize, NULL);
            num_white = MRIwhiteInPlane(mri_src, x, y, z, white_vertex, wsize);
            if (num_white < (wsize * wsize) / 2) {
              black_vertex = MRIneighborhoodBlackCpolv(mri_src, x, y, z, 3, wsize, NULL);
              num_black = wsize * wsize - MRIwhiteInPlane(mri_src, x, y, z, black_vertex, wsize);
              if ((num_white < num_black) && (((float)num_white < thresh) || (num_black >= hi_thresh))) {
                nwhite_to_black++;
                dst_label = MRI_NOT_WHITE;
              }
            }
          }
          else {
            black_vertex = MRIneighborhoodBlackCpolv(mri_src, x, y, z, 3, wsize, NULL);
            num_black = wsize * wsize - MRIwhiteInPlane(mri_src, x, y, z, black_vertex, wsize);
            if (num_black < wsize * wsize / 2) {
              white_vertex = MRIneighborhoodCpolv(mri_src, x, y, z, 3, wsize, NULL);
              num_white = MRIwhiteInPlane(mri_src, x, y, z, white_vertex, wsize);
              if ((num_white > num_black) && (((float)num_black < thresh) || (num_white >= hi_thresh))) {
                nblack_to_white++;
                dst_label = MRI_WHITE;
              }
            }
          }

          *pdst++ = dst_label;
        }
      }
    }

    if (Gdiag & DIAG_SHOW) {
      fprintf(stderr,
              "                %8d voxels deleted (%2.2f%%)\n",
              nwhite_to_black,
              100.0f * (float)nwhite_to_black / (float)total_vox);
      fprintf(stderr,
              "              %8d voxel added (%2.2f%%)\n",
              nblack_to_white,
              100.0f * (float)nblack_to_white / (float)total_vox);
      fprintf(stderr, "              %8d skipped (%%%2.2f)\n", nskipped, 100.0f * (float)nskipped / (float)total_vox);
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[100];
      sprintf(fname, "/tmp/smooth%d.mnc", i + 1);
      MRIwrite(mri_dst, fname);
    }
    if (i < niter - 1) MRIcopy(mri_dst, mri_src);
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIextractVertexCoords(MRI *mri_src, int *px, int *py, int *pz, int vertex, int x, int y, int z, int wsize)
{
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int whalf, xk, yk, xi, yi, zi;

  init_basis_vectors();
  whalf = (wsize - 1) / 2;

  e1_x = e1_x_v[vertex]; /* get basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk * e1_x)];
      yi = mri_src->yi[nint(ybase + xk * e1_y)];
      zi = mri_src->zi[nint(zbase + xk * e1_z)];
      *px++ = xi;
      *py++ = yi;
      *pz++ = zi;
    }
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIextractVertexPlane(MRI *mri_src, MRI *mri_dst, int vertex, int x0, int y0, int z0, int wsize)
{
  float e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase, x, y, z;
  int whalf, xk, yk;
  double val;

  init_basis_vectors();
  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = x0 - whalf * mri_dst->xsize;
    mri_dst->ystart = y0 - whalf * mri_dst->ysize;
    mri_dst->zstart = z0 - whalf * mri_dst->zsize;
    mri_dst->xend = mri_dst->xstart + wsize * mri_dst->xsize;
    mri_dst->yend = mri_dst->ystart + wsize * mri_dst->ysize;
    mri_dst->zend = mri_dst->zstart + wsize * mri_dst->zsize;
    mri_dst->imnr0 = z0 + mri_src->imnr0;
    mri_dst->imnr1 = mri_dst->imnr0;
  }

  e1_x = e1_x_v[vertex]; /* get basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x0 + (float)yk * e2_x;
    ybase = (float)y0 + (float)yk * e2_y;
    zbase = (float)z0 + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      x = xbase + xk * e1_x;
      y = ybase + xk * e1_y;
      z = zbase + xk * e1_z;
      MRIsampleVolume(mri_src, x, y, z, &val);
      MRIvox(mri_dst, xk + whalf, yk + whalf, 0) = (BUFTYPE)nint(val);
    }
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIneighborhoodBlackCpolv(MRI *mri_src, int xv, int yv, int zv, int nsize, int wsize, int *pnum_black)
{
  int best_plane[NVERTICES], whalf, nhalf, vertex, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex_index,
      max_count, x0, x1, y0, y1, z0, z1, x, y, z, num_black;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (xv == 166 && yv == 150 && zv == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;
  nhalf = (nsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  memset(best_plane, 0, sizeof(best_plane));
  x0 = MAX(0, xv - nhalf);
  x1 = MIN(mri_src->width - 1, xv + nhalf);
  y0 = MAX(0, yv - nhalf);
  y1 = MIN(mri_src->height - 1, yv + nhalf);
  z0 = MAX(0, zv - nhalf);
  z1 = MIN(mri_src->depth - 1, zv + nhalf);
  max_count = 0;
  for (z = z0; z <= z1; z++) {
    for (y = y0; y <= y1; y++) {
      for (x = x0; x <= x1; x++) {
        label = MRIvox(mri_src, x, y, z);
        if (label > MRI_NOT_WHITE) continue;
        pe1_x = e1_x_v;
        pe1_y = e1_y_v;
        pe1_z = e1_z_v;
        pe2_x = e2_x_v;
        pe2_y = e2_y_v;
        pe2_z = e2_z_v;
        for (vertex = 0; vertex < NVERTICES; vertex++) {
          num_black = 0;
          e1_x = *pe1_x++; /* first in-plane basis vector */
          e1_y = *pe1_y++;
          e1_z = *pe1_z++;
          e2_x = *pe2_x++; /* second in-plane basis vector */
          e2_y = *pe2_y++;
          e2_z = *pe2_z++;

          /* now find the values in this plane */
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              label = MRIvox(mri_src, xi, yi, zi);
              if (label < MRI_AMBIGUOUS) num_black++;
            }
          }
          if (num_black && (num_black >= max_count)) {
            if (num_black > max_count) /* new best-fit */
            {
              memset(best_plane, 0, sizeof(best_plane));
              max_count = num_black;
            }
            best_plane[vertex]++;
          }
        }
      }
    }
  }
  if (pnum_black) *pnum_black = max_count;
  peak_vertex_index = 0;

  /* now find the orientation with the most of the best-fitting planes */
  max_count = best_plane[0];
  for (vertex = 1; vertex < NVERTICES; vertex++) {
    if (best_plane[vertex] > max_count) {
      peak_vertex_index = vertex;
      max_count = best_plane[vertex];
    }
  }
  return (peak_vertex_index);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIneighborhoodCpolv(MRI *mri_src, int xv, int yv, int zv, int nsize, int wsize, int *pnum_white)
{
  int best_plane[NVERTICES], whalf, nhalf, vertex, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex_index,
      max_count, x0, x1, y0, y1, z0, z1, x, y, z, num_white;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (xv == 166 && yv == 150 && zv == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;
  nhalf = (nsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  memset(best_plane, 0, sizeof(best_plane));
  x0 = MAX(0, xv - nhalf);
  x1 = MIN(mri_src->width - 1, xv + nhalf);
  y0 = MAX(0, yv - nhalf);
  y1 = MIN(mri_src->height - 1, yv + nhalf);
  z0 = MAX(0, zv - nhalf);
  z1 = MIN(mri_src->depth - 1, zv + nhalf);
  max_count = 0;
  for (z = z0; z <= z1; z++) {
    for (y = y0; y <= y1; y++) {
      for (x = x0; x <= x1; x++) {
        label = MRIvox(mri_src, x, y, z);
        if (label != MRI_WHITE) continue;
        pe1_x = e1_x_v;
        pe1_y = e1_y_v;
        pe1_z = e1_z_v;
        pe2_x = e2_x_v;
        pe2_y = e2_y_v;
        pe2_z = e2_z_v;
        for (vertex = 0; vertex < NVERTICES; vertex++) {
          num_white = 0;
          e1_x = *pe1_x++; /* first in-plane basis vector */
          e1_y = *pe1_y++;
          e1_z = *pe1_z++;
          e2_x = *pe2_x++; /* second in-plane basis vector */
          e2_y = *pe2_y++;
          e2_z = *pe2_z++;

          /* now find the values in this plane */
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              label = MRIvox(mri_src, xi, yi, zi);
              if (label == MRI_WHITE) num_white++;
            }
          }
          if (num_white && (num_white >= max_count)) {
            if (num_white > max_count) /* new best-fit */
            {
              memset(best_plane, 0, sizeof(best_plane));
              max_count = num_white;
            }
            best_plane[vertex]++;
          }
        }
      }
    }
  }
  if (pnum_white) *pnum_white = max_count;
  peak_vertex_index = 0;

  /* now find the orientation with the most of the best-fitting planes */
  max_count = best_plane[0];
  for (vertex = 1; vertex < NVERTICES; vertex++) {
    if (best_plane[vertex] > max_count) {
      peak_vertex_index = vertex;
      max_count = best_plane[vertex];
    }
  }
  return (peak_vertex_index);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIneighborhoodPlanarDirection(MRI *mri_src, int xv, int yv, int zv, int nsize, int wsize)
{
  int best_white_plane[NVERTICES], best_black_plane[NVERTICES], whalf, nhalf, vertex, xk, yk, label, src_label,
      max_black_count, xi, yi, zi, *pxi, *pyi, *pzi, peak_vertex_index, max_white_count, x0, x1, y0, y1, z0, z1, x, y,
      z, num_white, num_black, max_count;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  if (xv == 166 && yv == 150 && zv == 127) DiagBreak(); /* 113 */

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;
  nhalf = (nsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  memset(best_white_plane, 0, sizeof(best_white_plane));
  x0 = MAX(0, xv - nhalf);
  x1 = MIN(mri_src->width - 1, xv + nhalf);
  y0 = MAX(0, yv - nhalf);
  y1 = MIN(mri_src->height - 1, yv + nhalf);
  z0 = MAX(0, zv - nhalf);
  z1 = MIN(mri_src->depth - 1, zv + nhalf);
  max_black_count = max_white_count = 0;
  for (z = z0; z <= z1; z++) {
    for (y = y0; y <= y1; y++) {
      for (x = x0; x <= x1; x++) {
        src_label = MRIvox(mri_src, x, y, z);
        if (src_label == MRI_AMBIGUOUS) continue;
        pe1_x = e1_x_v;
        pe1_y = e1_y_v;
        pe1_z = e1_z_v;
        pe2_x = e2_x_v;
        pe2_y = e2_y_v;
        pe2_z = e2_z_v;
        for (vertex = 0; vertex < NVERTICES; vertex++) {
          num_black = num_white = 0;
          e1_x = *pe1_x++; /* first in-plane basis vector */
          e1_y = *pe1_y++;
          e1_z = *pe1_z++;
          e2_x = *pe2_x++; /* second in-plane basis vector */
          e2_y = *pe2_y++;
          e2_z = *pe2_z++;

          /* now find the values in this plane */
          for (yk = -whalf; yk <= whalf; yk++) {
            xbase = (float)x + (float)yk * e2_x;
            ybase = (float)y + (float)yk * e2_y;
            zbase = (float)z + (float)yk * e2_z;
            for (xk = -whalf; xk <= whalf; xk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + xk * e1_x);
              xi = pxi[xi];
              yi = nint(ybase + xk * e1_y);
              yi = pyi[yi];
              zi = nint(zbase + xk * e1_z);
              zi = pzi[zi];
              label = MRIvox(mri_src, xi, yi, zi);
              if (label >= MRI_WHITE)
                num_white++;
              else if (label == MRI_NOT_WHITE)
                num_black++;
            }
          }
          if (src_label == MRI_NOT_WHITE && (num_black && (num_black >= max_black_count))) {
            if (num_black > max_black_count) /* new best-fit */
            {
              memset(best_black_plane, 0, sizeof(best_black_plane));
              max_black_count = num_black;
            }
            best_black_plane[vertex]++;
          }
          else if (src_label == MRI_WHITE && (num_white && (num_white >= max_white_count))) {
            if (num_white > max_white_count) /* new best-fit */
            {
              memset(best_white_plane, 0, sizeof(best_white_plane));
              max_white_count = num_white;
            }
            best_white_plane[vertex]++;
          }
        }
      }
    }
  }
  peak_vertex_index = 0;

  /* now find the orientation with the most of the best-fitting planes */
  max_count = best_white_plane[0] + best_black_plane[0];
  for (vertex = 1; vertex < NVERTICES; vertex++) {
    if (best_white_plane[vertex] + best_black_plane[vertex] > max_white_count) {
      peak_vertex_index = vertex;
      max_white_count = best_white_plane[vertex] + best_black_plane[vertex];
    }
  }
  return (peak_vertex_index);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIwhiteInPlane(MRI *mri_src, int x, int y, int z, int vertex, int wsize)
{
  int whalf, xk, yk, label, xi, yi, zi, *pxi, *pyi, *pzi, num_white;
  float xbase, ybase, zbase, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  e1_x = e1_x_v[vertex]; /* first in-plane basis vector */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex]; /* second in-plane basis vector */
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];
  num_white = 0;

  /* now find the values in this plane */
  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = nint(xbase + xk * e1_x);
      xi = pxi[xi];
      yi = nint(ybase + xk * e1_y);
      yi = pyi[yi];
      zi = nint(zbase + xk * e1_z);
      zi = pzi[zi];
      label = MRIvox(mri_src, xi, yi, zi);
      if (label >= MRI_WHITE) num_white++;
    }
  }
  return (num_white);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labeling of a volume.
------------------------------------------------------*/
MRI *MRIremoveHoles(MRI *mri_src, MRI *mri_dst, int wsize, float pct, int use_all)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, thresh, xi, yi, zi, num_on, num_off, in_val, nvox, nprocessed,
      nlabeled;
  BUFTYPE val, *pdst;
  float wcubed;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIcopy(mri_src, NULL);

  wcubed = (float)(wsize * wsize * wsize);
  thresh = nint((float)wcubed * pct);
  whalf = wsize / 2;

  nlabeled = nprocessed = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        if (x == 80 && y == 2 && z == 0) DiagBreak();
        in_val = MRIvox(mri_src, x, y, z);
        if (in_val != MRI_AMBIGUOUS) {
          *pdst++ = in_val;
          continue;
        }
        nprocessed++;
        for (nvox = num_on = num_off = 0, z0 = -whalf; z0 <= whalf; z0++) {
          zi = mri_src->zi[z + z0];
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = mri_src->xi[x + x0];
              val = MRIvox(mri_src, xi, yi, zi);
              if (val != MRI_AMBIGUOUS) {
                nvox++;
                if (val == MRI_NOT_WHITE)
                  num_off++;
                else if (val >= MRI_WHITE)
                  num_on++;
              }
            }
          }
        }
        if (!use_all) /* use pct of non-ambiguous voxels */
          thresh = nint((float)nvox * pct);
        if (num_off >= thresh)
          val = MRI_NOT_WHITE;
        else if (num_on >= thresh)
          val = MRI_WHITE;
        else
          val = MRI_AMBIGUOUS;
        if (val != MRI_AMBIGUOUS) nlabeled++;
        *pdst++ = val;
      }
    }
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr,
            "               %8d voxels processed (%2.2f%%)\n",
            nprocessed,
            100.0f * (float)nprocessed / (float)(width * height * depth));
    fprintf(stderr,
            "               %8d voxels labeled (%2.2f%%)\n",
            nlabeled,
            100.0f * (float)nlabeled / (float)(width * height * depth));
  }
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labeling of a volume.
------------------------------------------------------*/
MRI *MRImeanLabel(MRI *mri_src, MRI *mri_label, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi, in_val, label, total;
  BUFTYPE val, *pdst, *plabel, *psrc;
  float wcubed, mean;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      plabel = &MRIvox(mri_label, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      psrc = &MRIvox(mri_src, 0, y, z);
      for (x = 0; x < width; x++) {
        if (x == 80 && y == 2 && z == 0) DiagBreak();
        in_val = *psrc++;
        label = *plabel++;
        if (label != MRI_AMBIGUOUS) {
          *pdst++ = label;
          continue;
        }
        for (total = 0, z0 = -whalf; z0 <= whalf; z0++) {
          zi = mri_src->zi[z + z0];
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = mri_src->xi[x + x0];
              val = MRIvox(mri_src, xi, yi, zi);
              total += val;
            }
          }
        }
        mean = (float)total / wcubed;
        if (mean > in_val)
          val = MRI_WHITE;
        else
          val = MRI_NOT_WHITE;
        *pdst++ = val;
      }
    }
  }
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labeling of a volume.
------------------------------------------------------*/
MRI *MRIintensitySegmentation(MRI *mri_src, MRI *mri_labeled, float wm_low, float wm_hi, float gray_hi)
{
  int width, height, depth, x, y, z, nwhite, nblack, nambiguous;
  float val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_labeled) mri_labeled = MRIclone(mri_src, NULL);

  nwhite = nblack = nambiguous = 0;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == 105 && y == 75 && z == 127) DiagBreak();
        val = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (val < wm_low || val > wm_hi) {
          val = MRI_NOT_WHITE;
          nblack++;
        }
        else if (val <= gray_hi) {
          nambiguous++;
          val = MRI_AMBIGUOUS;
        }
        else {
          nwhite++;
          val = MRI_WHITE;
        }
        MRIsetVoxVal(mri_labeled, x, y, z, 0, val);
      }
    }
  }
  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr,
            "              %8d voxel white     (%.2f%%)\n",
            nwhite,
            100.0f * (float)nwhite / (float)(width * height * depth));
    fprintf(stderr,
            "              %8d voxel non white (%.2f%%)\n",
            nblack,
            100.0f * (float)nblack / (float)(width * height * depth));
    fprintf(stderr,
            "              %8d voxel ambiguous (%.2f%%)\n",
            nambiguous,
            100.0f * (float)nambiguous / (float)(width * height * depth));
  }
  return (mri_labeled);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labeling of a volume.
------------------------------------------------------*/
MRI *MRIthresholdLabel(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int wm_low)
{
  int width, height, depth, x, y, z;
  BUFTYPE val, *psrc, *plabel, *pdst, label;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      plabel = &MRIvox(mri_labeled, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      psrc = &MRIvox(mri_src, 0, y, z);
      for (x = 0; x < width; x++) {
        label = *plabel++;
        val = *psrc++;
        if (label == MRI_AMBIGUOUS) /* change label to white or non-white */
        {
          if (val < wm_low)
            label = MRI_NOT_WHITE;
          else
            label = MRI_WHITE;
        }

        *pdst++ = label;
      }
    }
  }
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Segment an MR image based on an order statistic.
------------------------------------------------------*/
MRI *MRIorderSegment(MRI *mri_src, MRI *mri_labeled, float thresh, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi, in_val, thresh_index, label, nvox;
  BUFTYPE val, *pdst, *sptr;
  float wcubed;
  static BUFTYPE *sort_array = NULL;
  static int sort_size = 0;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_labeled) mri_labeled = MRIclone(mri_src, NULL);

  wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;
  thresh_index = thresh * wcubed;

  if (sort_array && (wcubed != sort_size)) {
    free(sort_array);
    sort_array = NULL;
  }
  if (!sort_array) {
    sort_array = (BUFTYPE *)calloc(wcubed, sizeof(BUFTYPE));
    sort_size = wcubed;
  }
  for (nvox = z = 0; z < depth; z++) {
    DiagHeartbeat((float)z / (float)(depth - 1));
    for (y = 0; y < height; y++) {
      pdst = &MRIvox(mri_labeled, 0, y, z);
      for (x = 0; x < width; x++) {
        label = *pdst;
        in_val = MRIvox(mri_src, x, y, z);
        if (x == 28 && y == 10 && z == 63) DiagBreak();
        if (label != MRI_AMBIGUOUS) {
          pdst++;
          continue;
        }
        nvox++;
        for (sptr = sort_array, z0 = -whalf; z0 <= whalf; z0++) {
          zi = mri_src->zi[z + z0];
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = mri_src->xi[x + x0];
              *sptr++ = MRIvox(mri_src, xi, yi, zi);
            }
          }
        }

        qsort(sort_array, wcubed, sizeof(BUFTYPE), compare_sort_array);
        val = sort_array[thresh_index];
        if (val >= in_val)
          *pdst++ = MRI_NOT_WHITE;
        else
          *pdst++ = MRI_WHITE;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,
            "segmentation: %8d voxel processed (%2.2f%%)\n",
            nvox,
            100.0f * (float)nvox / (float)(width * height * depth));
  return (mri_labeled);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRImaskLabels(MRI *mri_src, MRI *mri_mask, MRI *mri_dst)
{
  int width, height, depth, x, y, z;
  float mask, src;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == 157 && y == 154 && z == 127) DiagBreak();
        mask = MRIgetVoxVal(mri_mask, x, y, z, 0);
        src = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (mask == MRI_AMBIGUOUS)
          MRIsetVoxVal(mri_dst, x, y, z, 0, MRI_AMBIGUOUS);
        else if (mask == MRI_WHITE)
          MRIsetVoxVal(mri_dst, x, y, z, 0, src);
        else
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);
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
int MRIcentralPlaneOfLeastVarianceNormalVoxel(MRI *mri_src, int wsize, int x, int y, int z)
{
  int width, height, depth, whalf, vertex, xk, yk, min_vertex, max_vertex, xi, yi, zi, *pxi, *pyi, *pzi;
  float min_mean, min_var, max_var, total, total_sq, nv, varv, avgv, val;
  BUFTYPE max_val;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;

  init_basis_vectors();

  pxi = mri_src->xi;
  pyi = mri_src->yi;
  pzi = mri_src->zi;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  max_vertex = min_vertex = -1;
  min_mean = 1000.0f;   /* mean of minimum variance plane */
  min_var = 100000.0f;  /* minimum variance of central planes */
  max_var = -100000.0f; /* maximum variance of central planes */
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    e1_x = *pe1_x++; /* first in-plane basis vector */
    e1_y = *pe1_y++;
    e1_z = *pe1_z++;
    e2_x = *pe2_x++; /* second in-plane basis vector */
    e2_y = *pe2_y++;
    e2_z = *pe2_z++;

    total = total_sq = nv = 0.0f;
    max_val = 0;
    /* now find the values in this plane */
    for (yk = -whalf; yk <= whalf; yk++) {
      xbase = (float)x + (float)yk * e2_x;
      ybase = (float)y + (float)yk * e2_y;
      zbase = (float)z + (float)yk * e2_z;
      for (xk = -whalf; xk <= whalf; xk++) {
        /* in-plane vect. is linear combination of scaled basis vects */
        xi = nint(xbase + xk * e1_x);
        xi = pxi[xi];
        yi = nint(ybase + xk * e1_y);
        yi = pyi[yi];
        zi = nint(zbase + xk * e1_z);
        zi = pzi[zi];
        val = (float)MRIvox(mri_src, xi, yi, zi);
        total += val;          /* sum of all values in this plane */
        total_sq += val * val; /* sum of squared values in this plane */
        nv++;                  /* # of values in this plane */
      }
    }

    if (!nv)
      varv = avgv = 0.0f;
    else {
      avgv = total / nv; /* mean of this plane */
      varv = total_sq / nv - avgv * avgv;
    }

    if (varv > max_var) {
      max_var = varv;
      max_vertex = vertex;
    }
    if (varv < min_var) {
      min_var = varv;
      min_vertex = vertex;
      min_mean = avgv;
    }
    if (FZERO(varv)) /* zero variance - won't find anything less */
      break;
  }
  return (min_vertex);
}
/*!
  \fn MRI *MRIcpolvMedianCurveSegment()
  \brief Reclassifies MRI_AMBIGUOUS voxels

*/
MRI *MRIcpolvMedianCurveSegment(
    MRI *mri, MRI *mri_labeled, MRI *mri_dst, int wsize, float len, float gray_hi, float wm_low)
{
  int x, y, z, width, height, depth, label, nlabeled, non, noff;

  printf("MRIcpolvMedianCurveSegment(): wsize=%d, len=%g, gmhi=%g, wmlow=%g\n",
	 wsize,len,gray_hi,wm_low);

  if (!mri_dst) mri_dst = MRIcopy(mri_labeled, NULL);

  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  for (non = noff = nlabeled = z = 0; z < depth; z++) {
    DiagShowPctDone((float)z / (float)(depth - 1), 5);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == 148 && y == 143 && z == 159) DiagBreak();
        label = MRIvox(mri_labeled, x, y, z);
        if (label != MRI_AMBIGUOUS) continue;
        nlabeled++;
        label = MRIcpolvMedianCurveVoxel(mri, mri_labeled, x, y, z, wsize, len, gray_hi, wm_low);
        if (label == MRI_WHITE)
          non++;
        else
          noff++;
        MRIvox(mri_dst, x, y, z) = label;
      }
    }
  }

  printf("  %8d voxels processed (%2.2f%%)\n", nlabeled,
            100.0f * (float)nlabeled / (float)(width * height * depth));
  printf("  %8d voxels white (%2.2f%%)\n", non,
            100.0f * (float)non / (float)(width * height * depth));
  printf("  %8d voxels non-white (%2.2f%%)\n",noff,
            100.0f * (float)noff / (float)(width * height * depth));

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIcpolvMedianCurveVoxel(
    MRI *mri, MRI *mri_labeled, int x0, int y0, int z0, int wsize, float len, float gray_hi, float wm_low)
{
  int vertex, what, i, maxi;
  FILE *fp = NULL;
  float dist, x, y, z, median, nx, ny, nz, white_dist, gray_dist, gray_val, white_val, val, white_grad, gray_grad,
      max_val, min_val, max_val_dist, min_val_dist;
  float medians[200], dists[200];

#if 0
  vertex = MRIneighborhoodCpolv(mri_labeled,x0, y0, z0, 3, wsize, NULL) ;
#else
  vertex = MRIcentralPlaneOfLeastVarianceNormalVoxel(mri, wsize, x0, y0, z0);
#endif
  nx = ic_x_vertices[vertex]; /* normal vector */
  ny = ic_y_vertices[vertex];
  nz = ic_z_vertices[vertex];
  dist = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= dist;
  ny /= dist;
  nz /= dist;

  white_dist = gray_dist = 1000.0f;
  what = MRI_AMBIGUOUS;
  gray_val = white_val = val = -1.0;
  min_val_dist = max_val_dist = 0.0;
  max_val = 0.0;
  min_val = 1000.0;
  for (i = 0, dist = -len; dist <= len; dist += 0.25f, i++) {
    x = x0 + nx * dist;
    y = y0 + ny * dist;
    z = z0 + nz * dist;
    median = MRIcpolvMedianAtVoxel(mri, vertex, x, y, z, wsize);
    medians[i] = median;
    dists[i] = dist;
    if ((median >= gray_hi) && (fabs(dist) < white_dist)) {
      white_val = median;
      white_dist = fabs(dist);
    }
    else if ((median <= wm_low) && (fabs(dist) < gray_dist)) {
      gray_val = median;
      gray_dist = fabs(dist);
    }
    if (median > max_val) {
      max_val_dist = fabs(dist);
      max_val = median;
    }
    if (median < min_val) {
      min_val_dist = fabs(dist);
      min_val = median;
    }

    if (FZERO(dist)) /* central vertex */
    {
      val = median;
      if (median >= gray_hi) {
        /*        fprintf(stderr, "median=%2.2f --> WHITE\n", median) ;*/
        what = MRI_WHITE;
      }
      else if (median <= wm_low) {
        what = MRI_NOT_WHITE;
        /*        fprintf(stderr, "median=%2.2f --> NOT WHITE\n", median) ;*/
      }
      if (!((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) && what != MRI_AMBIGUOUS) return (what);
    }
  }
#if 0
  if (white_val < 0)
  {
    white_dist  = max_val_dist ;
    white_val = max_val ;
  }
  if (gray_val < 0)
  {
    gray_dist = min_val_dist ;
    gray_val = min_val ;
  }
#endif

  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) {
    fp = fopen("med.plt", "w");
    maxi = i;
    for (i = 0; i < maxi; i++) fprintf(fp, "%2.2f  %2.3f\n", dists[i], medians[i]);
    fclose(fp);
  }

  if (white_dist > len) /* couldn't find white matter - look for peak */
  {
    maxi = i;
    for (i = 0; i < maxi; i++) {
      if ((medians[i] > white_val) && (medians[i] > wm_low)) {
        white_val = medians[i];
        white_dist = fabs(dists[i]);
      }
    }
    if (FZERO(white_dist)) what = MRI_WHITE; /* peak in this plane - must be white matter */
  }

  if (gray_dist > len) /* couldn't find gray matter - look for peak */
  {
    maxi = i;
    gray_val = 10000;
    for (i = 0; i < maxi; i++) {
      if ((medians[i] < gray_val) && (medians[i] < gray_hi)) {
        gray_val = medians[i];
        gray_dist = fabs(dists[i]);
      }
    }
    if (FZERO(gray_dist)) what = MRI_NOT_WHITE; /* peak in this plane - must be gray matter */
  }

  /* check to see if one distance is substantially bigger than the other */
  if (what == MRI_AMBIGUOUS) {
#define MAX_DIST_DIFF .75
    if (white_dist < (gray_dist - MAX_DIST_DIFF))
      what = MRI_WHITE;
    else if (gray_dist < (white_dist - MAX_DIST_DIFF))
      what = MRI_NOT_WHITE;
  }

  /* still don't know what it is - assign to region with smaller
     gradient magnitude.
  */
  if (what == MRI_AMBIGUOUS) {
    if (white_dist > 2.0) what = MRI_NOT_WHITE;
#if 1
    else {
#if 0
      white_grad = fabs(white_val-val) / white_dist ;
      gray_grad = fabs(gray_val-val) / gray_dist ;
#else
      white_grad = fabs(white_val - val);
      gray_grad = fabs(gray_val - val);
#endif
      if (white_grad > gray_grad)
        what = MRI_NOT_WHITE;
      else
        what = MRI_WHITE;
    }
#else
    else if (FEQUAL(gray_dist, white_dist)) {
      if (fabs(gray_val - val) < fabs(white_val - val))
        what = MRI_NOT_WHITE;
      else
        what = MRI_WHITE;
    }
    else if (gray_dist > white_dist)
      what = MRI_WHITE;
    else
      what = MRI_NOT_WHITE;
#endif
  }
  return (what);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float MRIcpolvMedianAtVoxel(MRI *mri_src, int vertex, float x, float y, float z, int wsize)
{
  int whalf, xk, yk;
  float xbase, ybase, zbase, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;
  float plane_vals[MAXLEN], *pvals;
  double rval, xr, yr, zr;

  init_basis_vectors();

  e1_x = e1_x_v[vertex]; /* basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  pvals = plane_vals;

  /* now find the values in this plane */
  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xr = xbase + xk * e1_x;
      yr = ybase + xk * e1_y;
      zr = zbase + xk * e1_z;
      MRIsampleVolume(mri_src, xr, yr, zr, &rval);
      *pvals++ = (float)rval;
    }
  }

  qsort(plane_vals, wsize * wsize, sizeof(float), compare_sort_farray);
  return (plane_vals[wsize * wsize / 2]);
}
int MRIvertexToVector(int vertex, float *pdx, float *pdy, float *pdz)
{
  if ((vertex < 0) || ((unsigned)vertex >= sizeof(ic_x_vertices) / sizeof(ic_x_vertices[0])))
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIvertexToVector(%d): index out of range", vertex));
  *pdx = ic_x_vertices[vertex];
  *pdy = ic_y_vertices[vertex];
  *pdz = ic_z_vertices[vertex];
  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labelling of a volume.
------------------------------------------------------*/
MRI *MRIremoveIslands(MRI *mri_src, MRI *mri_dst, int wsize, int thresh)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi, num_on, num_off;
  float val, out_val;
  float wcubed;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;

  for (z = whalf; z < depth; z++) {
    for (y = whalf; y < height; y++) {
      for (x = 0; x < width; x++) {
        for (out_val = num_on = num_off = 0, z0 = -whalf; z0 <= whalf; z0++) {
          zi = mri_src->zi[z + z0];
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = mri_src->xi[x + x0];
              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              if (!val)
                num_off++;
              else {
                if (!out_val) out_val = 0;
                num_on++;
              }
            }
          }
        }
        val = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (val && (num_off >= thresh))
          val = 0;
        else if (!val && (num_on >= thresh))
          val = THICKEN_FILL;
        MRIsetVoxVal(mri_dst, x, y, z, 0, val);
      }
    }
  }
  return (mri_dst);
}
MRI *MRIfillPlanarHoles(MRI *mri_src, MRI *mri_segment, MRI *mri_dst, MRI_SEGMENT *mseg)
{
  int width, height, depth, x, y, z, nfilled, total_filled, vertex;
  MRI *mri_binary_strand, *mri_strand_border;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  init_basis_vectors();
  mri_dst = MRIcopy(mri_src, mri_dst);

  /* make it 0-1 */
  mri_binary_strand = MRIbinarize(mri_segment, NULL, WM_MIN_VAL, 0, 1);

  /* mri_strand_border will be all voxels not in the segment, but within 2
     voxels of a segment voxel.
  */
  mri_strand_border = MRIdilate(mri_binary_strand, NULL);
  MRIdilate(mri_strand_border, mri_strand_border);
  MRIxor(mri_binary_strand, mri_strand_border, mri_strand_border, 1, 255);

  /*
    for each voxel that is off, but borders an 'on' voxel in the
    current segment, find max planar orientation at that point, and
    if a segment voxel is on in each of the 4 quadrants, and the current
    point doesn't neighbor a non-segment voxel, turn it on.
  */
  total_filled = 0;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    fprintf(stderr, "writing binary strand and border images...\n");
    MRIwrite(mri_binary_strand, "binary_strand.mgh");
    MRIwrite(mri_strand_border, "strand_border.mgh");
  }
  do {
    nfilled = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == 161 && y == 140 && z == 120) DiagBreak();

          /*
            find a voxel which is off but is close to the thin strand
            being considered.
          */
          if ((MRIvox(mri_dst, x, y, z) > WM_MIN_VAL) || (MRIvox(mri_strand_border, x, y, z) == 0)) continue;

          /* if this voxel neighbors an on voxel that is not in the
             strand then don't ever fill it.
          */
          /*
            check to see if it is in the plane of the strand,
            and if so fill it.
          */
          vertex = MRIcpolvMaxWhiteAtVoxel(mri_binary_strand, x, y, z, 5);
          if (MRIcpolvAllQuadrantsFilled(mri_binary_strand, x, y, z, vertex, 5)) {
            if (MRIneighborsOn(mri_binary_strand, x, y, z, 1) == MRIneighborsOn(mri_dst, x, y, z, WM_MIN_VAL)) {
              MRIsetVoxVal(mri_dst, x, y, z, 0, THICKEN_FILL);
              MRIsetVoxVal(mri_binary_strand, x, y, z, 0, 1);
              MRIsetVoxVal(mri_strand_border, x, y, z, 0, 0);
              nfilled++;
            }
          }
        }
      }
    }
    total_filled += nfilled;
  } while (nfilled > 0);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stderr, "%d planar holes filled\n", total_filled);
  MRIfree(&mri_binary_strand);
  MRIfree(&mri_strand_border);
  return (mri_dst);
}

int MRIcpolvMaxWhiteAtVoxel(MRI *mri, int x, int y, int z, int wsize)
{
  int whalf, vertex, xk, yk, peak_vertex, max_count, num, xi, yi, zi;
  float xbase, ybase, zbase, *pe1_x, *pe1_y, *pe1_z, xf, yf, zf, *pe2_x, *pe2_y, *pe2_z, e1_x, e1_y, e1_z, e2_x, e2_y,
      e2_z;
  double val;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  pe1_x = e1_x_v;
  pe1_y = e1_y_v;
  pe1_z = e1_z_v;
  pe2_x = e2_x_v;
  pe2_y = e2_y_v;
  pe2_z = e2_z_v;
  max_count = peak_vertex = 0;
  for (vertex = 0; vertex < NVERTICES; vertex++) {
    num = 0;
    e1_x = *pe1_x++; /* first in-plane basis vector */
    e1_y = *pe1_y++;
    e1_z = *pe1_z++;
    e2_x = *pe2_x++; /* second in-plane basis vector */
    e2_y = *pe2_y++;
    e2_z = *pe2_z++;

    /* now find the values in this plane */
    for (yk = -whalf; yk <= whalf; yk++) {
      xbase = (float)x + (float)yk * e2_x;
      ybase = (float)y + (float)yk * e2_y;
      zbase = (float)z + (float)yk * e2_z;
      for (xk = -whalf; xk <= whalf; xk++) {
        /* in-plane vect. is linear combination of scaled basis vects */
        xf = xbase + xk * e1_x;
        yf = ybase + xk * e1_y;
        zf = zbase + xk * e1_z;
        xi = mri->xi[nint(xf)];
        yi = mri->yi[nint(yf)];
        zi = mri->zi[nint(zf)];
#if 0
        MRIsampleVolume(mri, xi, yi, zi, &val) ;
#else
        val = MRIgetVoxVal(mri, xi, yi, zi, 0);
#endif
        if (val > 0.5) num++;
      }
    }
    if (num >= max_count) {
      peak_vertex = vertex;
      max_count = num;
    }
  }

  return (peak_vertex);
}

int MRIcpolvAllQuadrantsFilled(MRI *mri, int x, int y, int z, int vertex, int wsize)
{
#if 1
  int whalf;
  float xf, yf, zf, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, dist;
  double val;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  whalf = (wsize - 1) / 2;

  /*
    for this point (x,y,z), go through a set of directions on the unit
    sphere. For each direction, find all the planes orthogonal to that
    dir. within our window, and pick the direction in which the variance
    of all the planes is smallest. This will hopefully be the normal to
    the cortical surface.
    */
  e1_x = e1_x_v[vertex];
  e2_x = e2_x_v[vertex];
  e1_y = e1_y_v[vertex];
  e2_y = e2_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_z = e2_z_v[vertex];

#ifdef STEP_SIZE
#undef STEP_SIZE
#endif
#define STEP_SIZE 0.75

  /* positive 'x' direction */
  for (dist = STEP_SIZE; dist <= whalf; dist += STEP_SIZE) {
    xf = (float)x + dist * e1_x;
    yf = (float)y + dist * e1_y;
    zf = (float)z + dist * e1_z;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    if (val > 0.5) break;
  }
  if (val <= 0.5) return (0);

  /* negative 'y' direction */
  for (dist = STEP_SIZE; dist <= whalf; dist += STEP_SIZE) {
    xf = (float)x - dist * e1_x;
    yf = (float)y - dist * e1_y;
    zf = (float)z - dist * e1_z;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    if (val > 0.5) break;
  }
  if (val <= 0.5) return (0);

  /* positive 'y' direction */
  for (dist = STEP_SIZE; dist <= whalf; dist += STEP_SIZE) {
    xf = (float)x + dist * e2_x;
    yf = (float)y + dist * e2_y;
    zf = (float)z + dist * e2_z;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    if (val > 0.5) break;
  }
  if (val <= 0.5) return (0);

  /* negative 'y' direction */
  for (dist = STEP_SIZE; dist <= whalf; dist += STEP_SIZE) {
    xf = (float)x - dist * e2_x;
    yf = (float)y - dist * e2_y;
    zf = (float)z - dist * e2_z;
    MRIsampleVolume(mri, xf, yf, zf, &val);
    if (val > 0.5) break;
  }
  if (val <= 0.5) return (0);

  return (1);
#else
  int whalf, xk, yk, quads[2][2], xi, yi, i, j;
  float xbase, ybase, zbase, xf, yf, zf, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;
  double val;

  init_basis_vectors();

  if (x == 166 && y == 150 && z == 127) DiagBreak(); /* 113 */

  whalf = (wsize - 1) / 2;

  /*
  for this point (x,y,z), go through a set of directions on the unit
  sphere. For each direction, find all the planes orthogonal to that
  dir. within our window, and pick the direction in which the variance
  of all the planes is smallest. This will hopefully be the normal to
  the cortical surface.
  */
  e1_x = e1_x_v[vertex]; /* basis vectors for plane */
  e1_y = e1_y_v[vertex];
  e1_z = e1_z_v[vertex];
  e2_x = e2_x_v[vertex];
  e2_y = e2_y_v[vertex];
  e2_z = e2_z_v[vertex];
  memset(quads, 0, 4 * sizeof(int));
  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xf = xbase + xk * e1_x;
      yf = ybase + xk * e1_y;
      zf = zbase + xk * e1_z;
      MRIsampleVolume(mri, xf, yf, zf, &val);
      if (val > 0.5) {
        xi = xk < 0;
        yi = yk < 0;
        quads[xi][yi] = 1;
      }
    }
  }

  for (i = 0; i <= 1; i++)
    for (j = 0; j <= 1; j++)
      if (quads[i][j] == 0) return (0);

  return (1);
#endif
}
