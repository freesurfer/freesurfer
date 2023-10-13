/**
 * @brief filter routines on MRI data
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

#include "box.h"
#include "diag.h"
#include "error.h"
#include "fftutils.h"
#include "filter.h"
#include "macros.h"
#include "minc.h"
#include "mri.h"
#include "mrinorm.h"
#include "proto.h"
#include "region.h"
#include "talairachex.h"
#include "cma.h"

#include "romp_support.h"

#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*#define DEBUG_POINT(x,y,z)  (((x==8&&y==9) || (x==9&&y==8)) &&((z)==15))*/
/*#define DEBUG_POINT(x,y,z)  (((x==7&&y==9) || (x==9&&y==7)) &&((z)==15))*/
#define DEBUG_POINT(x, y, z) (((x == 21) && (y == 14)) && ((z) == 7))

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int compare_sort_array(const void *pc1, const void *pc2);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define NORM 1.0f

MRI *MRIdirectionMapUchar(MRI *mri_grad, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0;
  float dir, dot, *xpix, *ypix, *zpix, dx, dy, dz, ox, oy, oz, norm;
  BUFTYPE *dir_pix;

  whalf = wsize / 2;
  norm = (float)(wsize * wsize * wsize) * NORM;
  norm *= norm;
  width = mri_grad->width;
  height = mri_grad->height;
  depth = mri_grad->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    MRIcopyHeader(mri_grad, mri_dst);
  }

  for (z = whalf; z < depth - whalf; z++) {
    for (y = whalf; y < height - whalf; y++) {
      dir_pix = &MRIvox(mri_dst, whalf, y, z);
      for (x = whalf; x < width - whalf; x++) {
        ox = MRIFvox(mri_grad, x, y, z);
        oy = MRIFseq_vox(mri_grad, x, y, z, 1);
        oz = MRIFseq_vox(mri_grad, x, y, z, 2);
        for (dir = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          for (y0 = -whalf; y0 <= whalf; y0++) {
            xpix = &MRIFvox(mri_grad, x - whalf, y + y0, z + z0);
            ypix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 1);
            zpix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 2);
            for (x0 = -whalf; x0 <= whalf; x0++) {
              dx = *xpix++;
              dy = *ypix++;
              dz = *zpix++;
              dot = dx * ox + dy * oy + dz * oz;
              dir += (x0 * ox + y0 * oy + z0 * oz) * dot;
            }
          }
        }
        *dir_pix++ = (BUFTYPE)nint(127 * (1 + tanh(dir / norm)));
      }
    }
  }

  return (mri_dst);
}
MRI *MRIdirectionMap(MRI *mri_grad, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0;
  float dir, dot, *xpix, *ypix, *zpix, dx, dy, dz, ox, oy, oz, *dir_pix, norm;

  whalf = wsize / 2;
  norm = (float)(wsize * wsize * wsize) * NORM;
  norm *= norm;
  width = mri_grad->width;
  height = mri_grad->height;
  depth = mri_grad->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_grad, mri_dst);
  }

  for (z = whalf; z < depth - whalf; z++) {
    for (y = whalf; y < height - whalf; y++) {
      dir_pix = &MRIFvox(mri_dst, whalf, y, z);
      for (x = whalf; x < width - whalf; x++) {
        ox = MRIFvox(mri_grad, x, y, z);
        oy = MRIFseq_vox(mri_grad, x, y, z, 1);
        oz = MRIFseq_vox(mri_grad, x, y, z, 2);
        for (dir = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          for (y0 = -whalf; y0 <= whalf; y0++) {
            xpix = &MRIFvox(mri_grad, x - whalf, y + y0, z + z0);
            ypix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 1);
            zpix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 2);
            for (x0 = -whalf; x0 <= whalf; x0++) {
              dx = *xpix++;
              dy = *ypix++;
              dz = *zpix++;
              dot = dx * ox + dy * oy + dz * oz;
              dir += (x0 * ox + y0 * oy + z0 * oz) * dot;
            }
          }
        }
        *dir_pix++ = tanh(dir / norm);
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
#define SCALE_FACTOR 10.0f
MRI *MRIoffsetDirection(MRI *mri_grad, int wsize, MRI *mri_direction, MRI *mri_dir)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0;
  float dir, dot, *xpix, *ypix, *zpix, dx, dy, dz, ox, oy, oz, len, *dir_xpix, *dir_ypix, *dir_zpix;
  BUFTYPE *pdir;

  whalf = wsize / 2;
  width = mri_grad->width;
  height = mri_grad->height;
  depth = mri_grad->depth;

  if (!mri_direction) {
    mri_direction = MRIallocSequence(width, height, depth, MRI_FLOAT, 3);
    MRIcopyHeader(mri_grad, mri_direction);
  }

  for (z = whalf; z < depth - whalf; z++) {
    for (y = whalf; y < height - whalf; y++) {
      if (mri_dir) {
        pdir = &MRIvox(mri_dir, whalf, y, z);
      }
      else {
        pdir = NULL;
      }

      dir_xpix = &MRIFvox(mri_direction, whalf, y, z);
      dir_xpix = &MRIFvox(mri_direction, whalf, y, z);
      dir_ypix = &MRIFseq_vox(mri_direction, whalf, y, z, 1);
      dir_zpix = &MRIFseq_vox(mri_direction, whalf, y, z, 2);
      for (x = whalf; x < width - whalf; x++) {
#define DEBUG_OFFSET 0
#if DEBUG_OFFSET
        FILE *fp;
#endif
        ox = MRIFvox(mri_grad, x, y, z);
        oy = MRIFseq_vox(mri_grad, x, y, z, 1);
        oz = MRIFseq_vox(mri_grad, x, y, z, 2);
#if DEBUG_OFFSET
        if (DEBUG_POINT(x, y, z)) {
          fp = fopen("dir.log", "w");
          fprintf(fp, "direction (%d, %d, %d), O: (%2.1f, %2.1f, %2.1f)\n", x, y, z, ox, oy, oz);
          DiagBreak();
        }
#endif
        for (dir = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          for (y0 = -whalf; y0 <= whalf; y0++) {
            xpix = &MRIFvox(mri_grad, x - whalf, y + y0, z + z0);
            ypix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 1);
            zpix = &MRIFseq_vox(mri_grad, x - whalf, y + y0, z + z0, 2);
            for (x0 = -whalf; x0 <= whalf; x0++) {
              dx = *xpix++;
              dy = *ypix++;
              dz = *zpix++;
              dot = dx * ox + dy * oy + dz * oz;
#if 1
              if (dot < 0.0f) {
                dot = 0.0f;
              }
#endif
              dir += (x0 * ox + y0 * oy + z0 * oz) * dot;
#if DEBUG_OFFSET
              if (DEBUG_POINT(x, y, z))
                fprintf(fp,
                        "(%d, %d, %d): (%2.1f, %2.1f, %2.1f), "
                        "dot = %2.2f, prod = %2.1f (%2.1f)\n",
                        x + x0,
                        y + y0,
                        z + z0,
                        dx,
                        dy,
                        dz,
                        dot,
                        (x0 * ox + y0 * oy + z0 * oz),
                        dir);
#endif
            }
          }
        }
#if DEBUG_OFFSET
        if (DEBUG_POINT(x, y, z)) {
          fclose(fp);
        }
#endif
        if (ISSMALL(dir)) {
          ox = oy = oz = 0.0f;
        }
        else if (dir > 0.0f) {
          ox = -ox;
          oy = -oy;
          oz = -oz;
        }
        len = sqrt(ox * ox + oy * oy + oz * oz);
        if (FZERO(len)) {
          *dir_xpix++ = 0.0f;
          *dir_ypix++ = 0.0f;
          *dir_zpix++ = 0.0f;
        }
        else {
          *dir_xpix++ = ox;
          *dir_ypix++ = oy;
          *dir_zpix++ = oz;
        }

        if (pdir) {
          if (ISSMALL(dir)) {
            *pdir++ = OFFSET_ZERO;
          }
          else if (dir < 0.0f) {
            *pdir++ = OFFSET_GRADIENT_DIRECTION;
          }
          else {
            *pdir++ = OFFSET_NEGATIVE_GRADIENT_DIRECTION;
          }
        }
      }
    }
  }

  return (mri_direction);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIoffsetMagnitude(MRI *mri_src, MRI *mri_dst, int maxsteps)
{
  int width, height, depth, x, y, z, x1, y1, z1, steps;
  float odx, ody, odz, dx, dy, dz, *xpix, *ypix, *zpix, len, dot = 0.0f, adx, ady, adz;
  signed char *dst_xpix, *dst_ypix, *dst_zpix;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, MRI_UCHAR, 3);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      xpix = &MRIFvox(mri_src, 0, y, z);
      ypix = &MRIFseq_vox(mri_src, 0, y, z, 1);
      zpix = &MRIFseq_vox(mri_src, 0, y, z, 2);

      dst_xpix = &MRISCvox(mri_dst, 0, y, z);
      dst_ypix = &MRISCseq_vox(mri_dst, 0, y, z, 1);
      dst_zpix = &MRISCseq_vox(mri_dst, 0, y, z, 2);
      for (x = 0; x < width; x++) {
        odx = *xpix++;
        ody = *ypix++;
        odz = *zpix++;
        adx = fabs(odx);
        ady = fabs(ody);
        adz = fabs(odz);

        /* scale so that offset reaches next voxel at each step */
        if ((adx > ady) && (adx > adz)) /* scale using x component */
        {
          len = adx;
        }
        else if (ady > adz) /* scale using y component */
        {
          len = ady;
        }
        else /* scale using z component */
        {
          len = adz;
        }
        if (!FZERO(len)) /* a non-zero offset */
        {
          odx /= len;
          ody /= len;
          odz /= len;

          for (steps = 0; steps < maxsteps; steps++) {
            x1 = x + nint((float)steps * odx);
            y1 = y + nint((float)steps * ody);
            z1 = z + nint((float)steps * odz);
            if (x1 < 0 || x1 >= width || y1 < 0 || y1 >= height || z1 < 0 || z1 >= depth) {
              break;
            }
            dx = MRIFvox(mri_src, x1, y1, z1);
            dy = MRIFseq_vox(mri_src, x1, y1, z1, 1);
            dz = MRIFseq_vox(mri_src, x1, y1, z1, 2);
            dot = dx * odx + dy * ody + dz * odz;
            if (dot <= 0.0f) {
              break;
            }
          }
          if (!FZERO(dot)) {
            steps--;
          }
        }
        else {
          steps = 0;
        }
        x1 = x + nint((float)steps * odx);
        y1 = y + nint((float)steps * ody);
        z1 = z + nint((float)steps * odz);

        *dst_xpix++ = (x1 - x);
        *dst_ypix++ = (y1 - y);
        *dst_zpix++ = (z1 - z);
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
MRI *MRIapplyOffset(MRI *mri_src, MRI *mri_dst, MRI *mri_offset)
{
  int width, height, depth, x, y, z, dx, dy, dz;
  signed char *pdx, *pdy, *pdz;
  char *pdst;
  float *fdst;

  fdst = NULL;
  pdst = NULL; /* for compiler warning - I know it's silly */
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  if (mri_dst->type != mri_src->type)
    ErrorReturn(NULL, (ERROR_BADPARM, "MRIapplyOffset: source type %d != dst %d", mri_src->type, mri_dst->type));

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pdx = &MRISCvox(mri_offset, 0, y, z);
      pdy = &MRISCseq_vox(mri_offset, 0, y, z, 1);
      pdz = &MRISCseq_vox(mri_offset, 0, y, z, 2);
      switch (mri_dst->type) {
        case MRI_UCHAR:
          pdst = (char *)&MRIvox(mri_dst, 0, y, z);
          break;
        case MRI_FLOAT:
          fdst = &MRIFvox(mri_dst, 0, y, z);
          break;
        default:
          ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIapplyOffset: unsupported type %d", mri_dst->type));
          break;
      }
      for (x = 0; x < width; x++) {
        dx = *pdx++;
        dy = *pdy++;
        dz = *pdz++;

        switch (mri_src->type) {
          case MRI_UCHAR:
            *pdst++ = MRIvox(mri_src, x + dx, y + dy, z + dz);
            break;
          case MRI_FLOAT:
            *fdst++ = MRIFvox(mri_src, x + dx, y + dy, z + dz);
            break;
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
MRI *MRIdiffuse(MRI *mri_src, MRI *mri_dst, double k, int niter, int which, double slope)
{
  switch (which) {
    case FILTER_DIFFUSE_GRAD:
      mri_dst = MRIdiffusePerona(mri_src, mri_dst, k, niter, slope);
      break;
    case FILTER_DIFFUSE_CURV:
      mri_dst = MRIdiffuseCurvature(mri_src, mri_dst, k, niter, slope);
      break;
  }
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIdiffuseCurvature(MRI *mri_src, MRI *mri_dst, double A, int niter, double slope)
{
  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define C(grad, k) ((float)exp(-0.5f * SQR((float)fabs((grad)) / k)))

#define NUM_NEIGH 6
#define KERNEL_MUL (1.0f / ((float)NUM_NEIGH))

MRI *MRIdiffusePerona(MRI *mri_src, MRI *mri_dst, double k, int niter, double slope)
{
  int x, y, z, i, width, height, depth, *xE, *yN, *xW, *yS, *zD, *zU, ys, yn, xe, xw, zu, zd, ci;
  float c[NUM_NEIGH + 1], fvals[NUM_NEIGH + 1], dst_val;
  MRI *mri_tmp = NULL, *mri_mag = NULL, *mri_grad = NULL;
  ;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

#if 0
  if ((outImage->pixel_format != inImage->pixel_format) ||
      (outImage->pixel_format != PFFLOAT))
  {
    fprintf(stderr,
            "ImageDiffusePerona: input and output image format must both be "
            "in float format.\n");
    exit(2) ;
  }

  if ((outImage->rows != inImage->rows) ||
      (outImage->cols != inImage->cols))
  {
    fprintf(stderr,
            "ImageDiffusePerona: input and output image sizes must match.\n");
    exit(2) ;
  }

#endif

  if (!MRIcheckSize(mri_src, mri_mag, 0, 0, 0)) {
    if (mri_mag) {
      MRIfree(&mri_mag);
    }
    if (mri_tmp) {
      MRIfree(&mri_tmp);
    }
    mri_mag = MRIalloc(width, height, depth, MRI_FLOAT);
    mri_tmp = MRIalloc(width, height, depth, MRI_FLOAT);
  }

  /* build index tables for border */
  xE = (int *)calloc((unsigned int)width, sizeof(int));
  xW = (int *)calloc((unsigned int)width, sizeof(int));
  yN = (int *)calloc((unsigned int)height, sizeof(int));
  yS = (int *)calloc((unsigned int)height, sizeof(int));
  zU = (int *)calloc((unsigned int)depth, sizeof(int));
  zD = (int *)calloc((unsigned int)depth, sizeof(int));

  xW[0] = 0;
  for (x = 1; x < width; x++) {
    xW[x] = x - 1;
    xE[x - 1] = x;
  }
  xE[width - 1] = width - 1;
  yN[0] = 0;
  for (y = 1; y < height; y++) {
    yN[y] = y - 1;
    yS[y - 1] = y;
  }
  yS[height - 1] = height - 1;
  zD[0] = 0; /* up is positive z */
  for (z = 1; z < depth; z++) {
    zD[z] = z - 1;
    zU[z - 1] = z;
  }
  zU[depth - 1] = depth - 1;

  MRIcopy(mri_src, mri_tmp);

  memset(fvals, 0, NUM_NEIGH * sizeof(float));
  memset(c, 0, NUM_NEIGH * sizeof(float));
  for (i = 0; i < niter; i++) {
    mri_grad = MRIsobel(mri_tmp, mri_grad, mri_mag);

    for (x = 0; x < width; x++) {
      xe = xE[x];
      xw = xW[x];
      for (y = 0; y < height; y++) {
        yn = yN[y];
        ys = yS[y];

        for (z = 0; z < depth; z++) {
          zu = zU[z];
          zd = zD[z];
          /*
               |  C1  |
            ---------------
            C2 |  C0  | C3
            ---------------
               |  C4  |
          C6 (zU) is into the screen and c5 (zD) is out of the screen
          */

          fvals[0] = MRIFvox(mri_tmp, x, y, z);
          fvals[1] = MRIFvox(mri_tmp, x, yn, z);
          fvals[2] = MRIFvox(mri_tmp, xw, y, z);
          fvals[3] = MRIFvox(mri_tmp, xe, y, z);
          fvals[4] = MRIFvox(mri_tmp, x, ys, z);
          fvals[5] = MRIFvox(mri_tmp, x, y, zd);
          fvals[6] = MRIFvox(mri_tmp, x, y, zu);

          c[1] = KERNEL_MUL * C(MRIFvox(mri_mag, x, yn, z), k);
          c[2] = KERNEL_MUL * C(MRIFvox(mri_mag, xw, y, z), k);
          c[3] = KERNEL_MUL * C(MRIFvox(mri_mag, xe, y, z), k);
          c[4] = KERNEL_MUL * C(MRIFvox(mri_mag, x, ys, z), k);
          c[5] = KERNEL_MUL * C(MRIFvox(mri_mag, x, y, zd), k);
          c[6] = KERNEL_MUL * C(MRIFvox(mri_mag, x, y, zu), k);

          c[0] = 1.0f;
          for (ci = 1; ci <= NUM_NEIGH; ci++) {
            c[0] -= c[ci];
          }

          for (dst_val = 0.0f, ci = 0; ci <= NUM_NEIGH; ci++) {
            dst_val += fvals[ci] * c[ci];
          }
          MRIFvox(mri_dst, x, y, z) = dst_val;
        }
      }
    }

    MRIcopy(mri_dst, mri_tmp);
    k = k + k * slope;
  }

  free(xE);
  free(xW);
  free(yN);
  free(yS);
  free(zU);
  free(zD);

  return (mri_dst);
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Apply a sobel operator to a region of an image. The output images
          will be the size of the region, not the original image.
------------------------------------------------------*/
MRI *
MRIsobelRegion(MRI *mri_src, MRI *mri_grad, int domag, MRI_REGION *region)
{
  int      x, y, z, width, height, depth, w, h, d  ;
  int      x0, y0, z0, xi,  yi, zi ;
  float    *px, *py, *pz, xval,yval,zval, *pmag ;

  width = region->x + region->dx ;
  if (width > mri_src->width)
  {
    width = mri_src->width ;
  }
  height = region->y + region->dy ;
  if (height > mri_src->height)
  {
    height = mri_src->height ;
  }
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
  {
    depth = mri_src->depth ;
  }
  x0 = region->x ;
  if (x0 < 0)
  {
    x0 = 0 ;
  }
  y0 = region->y ;
  if (y0 < 0)
  {
    y0 = 0 ;
  }
  z0 = region->z ;
  if (z0 < 0)
  {
    z0 = 0 ;
  }

  w = width - region->x ;
  h = height - region->y ;
  d = depth - region->z ;

  if (!mri_grad)
  {
    mri_grad = MRIallocSequence(w, h, d, MRI_FLOAT, domag ? 4 : 3) ;
    MRIcopyHeader(mri_src, mri_grad) ;
    mri_grad->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_grad->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_grad->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_grad->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_grad->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_grad->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  MRIxSobelRegion(mri_src, mri_grad, 0, region) ;
  MRIySobelRegion(mri_src, mri_grad, 1, region) ;
  MRIzSobelRegion(mri_src, mri_grad, 2, region) ;

  if (domag)
  {
    for (z = 0 ; z < d ; z++)
    {
      for (y = 0 ; y < h ; y++)
      {
        px =   &MRIFvox(mri_grad, 0, y, z) ;
        py =   &MRIFseq_vox(mri_grad, 0, y, z, 1) ;
        pz =   &MRIFseq_vox(mri_grad, 0, y, z, 2) ;
        pmag = &MRIFseq_vox(mri_grad, 0, y, z, 3) ;

        for (x = 0 ; x < w ; x++)
        {
          xval = *px++ ;
          yval = *py++ ;
          zval = *pz++ ;
          *pmag++ = sqrt(xval*xval + yval*yval + zval*zval) ;
        }
      }
    }
  }

  return(mri_grad) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIxSobelRegion(MRI *mri_src, MRI *mri_x, int frame, MRI_REGION *region)
{
  BUFTYPE *tr_pix, *mr_pix, *br_pix ;
  int              x, y, z, width, height, depth, x0, y0, z0 ;
  float            *outPtr, left, middle, right ;

  width = region->x + region->dx ;
  if (width > mri_src->width)
  {
    width = mri_src->width ;
  }
  height = region->y + region->dy ;
  if (height > mri_src->height)
  {
    height = mri_src->height ;
  }
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
  {
    depth = mri_src->depth ;
  }
  x0 = region->x ;
  if (x0 < 0)
  {
    x0 = 0 ;
  }
  y0 = region->y ;
  if (y0 < 0)
  {
    y0 = 0 ;
  }
  z0 = region->z ;
  if (z0 < 0)
  {
    z0 = 0 ;
  }

  if (!mri_x)
  {
    int w, h, d ;
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_x = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_x) ;
    mri_x->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_x->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_x->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_x->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_x->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_x->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  /* applying sobel in x-y plane, don't worry about boundary conditions in z */
  /* don't apply sobel to outer ring to pixels to avoid border effects */
  width-- ;
  height-- ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 1 ; y < height ; y++)
    {
      tr_pix = &MRIvox(mri_src, 0, y-1, z) ;
      mr_pix = &MRIvox(mri_src, 0, y, z) ;
      br_pix = &MRIvox(mri_src, 0, y+1, z) ;
      outPtr = &MRIFseq_vox(mri_x, 1, y, z, frame) ;

      left =    (float)(*tr_pix++ + 2 * *mr_pix++ + *br_pix++) ;
      middle =  (float)(*tr_pix++ + 2 * *mr_pix++ + *br_pix++) ;


      for (x = 1 ; x < width ; x++)
      {
        right = (float)*tr_pix++ + 2.0f * (float)*mr_pix++ + (float)*br_pix++ ;
        *outPtr++ = (right - left) / 8.0f ;
        left = middle ;
        middle = right ;
      }
    }
  }

  return(mri_x) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIySobelRegion(MRI *mri_src, MRI *mri_y, int frame, MRI_REGION *region)
{
  BUFTYPE *tr_pix, *br_pix ;
  int              width, height, depth, x, y, z ;
  float            *outPtr, left, middle, right ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_y)
  {
    mri_y = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_y) ;
  }

  /* applying sobel in x-y plane, don't worry about boundary conditions in z */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ;
  height-- ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 1 ; y < height ; y++)
    {
      tr_pix = &MRIvox(mri_src, 0, y+1, z) ;
      br_pix = &MRIvox(mri_src, 0, y-1, z) ;
      outPtr = &MRIFseq_vox(mri_y, 1, y, z, frame) ;

      left = (float)*br_pix++ - (float)*tr_pix++ ;
      middle = (float)*br_pix++ - (float)*tr_pix++ ;

      for (x = 1 ; x < width ; x++)
      {
        right = (float)*br_pix++ - (float)*tr_pix++ ;
        *outPtr++ = (right + 2.0f * middle + left) / 8.0f ;
        left = middle ;
        middle = right ;
      }
    }
  }

  return(mri_y) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           calculate a z derivitive in the x-z plane
------------------------------------------------------*/
MRI *
MRIzSobelRegion(MRI *mri_src, MRI *mri_z, int frame, MRI_REGION *region)
{
  BUFTYPE *tr_pix, *br_pix ;
  int              width, height, depth, x, y, z ;
  float            *outPtr, left, middle, right ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_z)
  {
    mri_z = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_z) ;
  }

  /* applying sobel in x-z plane, don't worry about boundary conditions in y */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ;
  depth-- ;
  for (z = 1 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      tr_pix = &MRIvox(mri_src, 0, y, z-1) ;
      br_pix = &MRIvox(mri_src, 0, y, z+1) ;
      outPtr = &MRIFseq_vox(mri_z, 1, y, z, frame) ;

      left = (float)*br_pix++ - (float)*tr_pix++ ;
      middle = (float)*br_pix++ - (float)*tr_pix++ ;

      for (x = 1 ; x < width ; x++)
      {
        right = (float)*br_pix++ - (float)*tr_pix++ ;
        *outPtr++ = (right + 2.0f * middle + left) / 8.0f ;
        left = middle ;
        middle = right ;
      }
    }
  }

  return(mri_z) ;
}
#endif

MRI *MRIdivergence(MRI *mri_src, MRI *mri_divergence)
{
  MRI *mri_sobel;
  int x, y, z;
  double dU_dx, dU_dy, dU_dz;

  mri_divergence = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
  MRIcopyHeader(mri_src, mri_divergence);
  mri_sobel = MRIsobel(mri_src, NULL, NULL);
  for (x = 0; x < mri_src->width; x++)
    for (y = 0; y < mri_src->height; y++)
      for (z = 0; z < mri_src->depth; z++) {
        dU_dx = MRIgetVoxVal(mri_sobel, x, y, z, 0);
        dU_dy = MRIgetVoxVal(mri_sobel, x, y, z, 1);
        dU_dz = MRIgetVoxVal(mri_sobel, x, y, z, 2);
        MRIsetVoxVal(mri_divergence, x, y, z, 0, dU_dx + dU_dy + dU_dz);
      }

  MRIfree(&mri_sobel);

  return (mri_divergence);
}

// compute the laplacian of the input volume (6-connected)
MRI *MRIlaplacian(MRI *mri_src, MRI *mri_laplacian)
{
  int x, y, z, f;
  float lap;

  if (mri_laplacian == NULL) {
    mri_laplacian = MRIcloneDifferentType(mri_src, MRI_FLOAT);
  }

  for (f = 0; f < mri_src->nframes; f++)
    for (x = 0; x < mri_src->width; x++) {
      for (y = 0; y < mri_src->height; y++) {
        for (z = 0; z < mri_src->depth; z++) {
          lap = -6 * MRIgetVoxVal(mri_src, x, y, z, 0);
          lap += MRIgetVoxVal(mri_src, mri_src->xi[x - 1], y, z, f);
          lap += MRIgetVoxVal(mri_src, mri_src->xi[x + 1], y, z, f);
          lap += MRIgetVoxVal(mri_src, x, mri_src->yi[y - 1], z, f);
          lap += MRIgetVoxVal(mri_src, x, mri_src->yi[y + 1], z, f);
          lap += MRIgetVoxVal(mri_src, x, y, mri_src->zi[z - 1], f);
          lap += MRIgetVoxVal(mri_src, x, y, mri_src->zi[z + 1], f);
          MRIsetVoxVal(mri_laplacian, x, y, z, f, lap);
        }
      }
    }

  return (mri_laplacian);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:
          Two orthogonal basis vectors

        Description
          Compute two vectors orthogonal to each other and to a third
------------------------------------------------------*/
int compute_tangent_vectors(
    float nx, float ny, float nz, float *pe1x, float *pe1y, float *pe1z, float *pe2x, float *pe2y, float *pe2z)
{
  float vx, vy, vz, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, len;

  /* pick some other unit non-linearly dependent vector */
  if (!FEQUAL(nx, ny)) {
    vx = ny;
    vy = nx;
    vz = nz;
  }
  else {
    vx = ny;
    vy = nz;
    vz = nx;
  }

  /* don't care about sign (right-hand rule) */
  e1_x = vy * nz - vz * ny;
  e1_y = vz * nx - vx * nz;
  e1_z = vx * ny - vy * nx;

  len = sqrt(e1_x * e1_x + e1_y * e1_y + e1_z * e1_z);
  if (FZERO(len)) {
    fprintf(stderr, "zero tangent vector computed!!\n");
    DiagBreak();
    exit(0);
  }
  else {
    e1_x /= len;
    e1_y /= len;
    e1_z /= len;
  }

  e2_x = e1_y * nz - e1_z * ny;
  e2_y = e1_x * nz - e1_z * nx;
  e2_z = e1_y * nx - e1_x * ny;
  len = sqrt(e2_x * e2_x + e2_y * e2_y + e2_z * e2_z);
  e2_x /= len;
  e2_y /= len;
  e2_z /= len;

  *pe1x = e1_x;
  *pe1y = e1_y;
  *pe1z = e1_z;
  *pe2x = e2_x;
  *pe2y = e2_y;
  *pe2z = e2_z;

#if 0
  DiagFprintf(0L,
              "vertex %d: (%2.2f, %2.2f, %2.2f) --> (%2.2f, %2.2f, %2.2f) "
              "x (%2.2f, %2.2f, %2.2f)\n",
              vertex, ic_x_vertices[vertex], ic_y_vertices[vertex],
              ic_z_vertices[vertex], e1_x_v[vertex], e1_y_v[vertex],
              e1_z_v[vertex],e2_x_v[vertex], e2_y_v[vertex], e2_z_v[vertex]) ;
  DiagFprintf(0L, "lengths: %2.3f, %2.3f, %2.3f\n",
              sqrt(nx*nx+ny*ny+nz*nz),
              sqrt(e1_x*e1_x+e1_y*e1_y+e1_z*e1_z),
              sqrt(e2_x*e2_x+e2_y*e2_y+e2_z*e2_z)) ;
#endif
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:
          a volume with the 2nd derivative values

        Description
          compute the 2nd derivative in the specified direction
------------------------------------------------------*/
MRI *MRI2ndDirectionalDerivative(MRI *mri_src, MRI *mri_deriv, float nx, float ny, float nz)
{
  int x, y, z, dn, d1, d2;
  float deriv, xf, yf, zf, e1x, e1y, e1z, e2x, e2y, e2z, wt = 1;
  double val;

  if (mri_deriv == NULL) {
    mri_deriv = MRIcloneDifferentType(mri_src, MRI_FLOAT);
  }

  compute_tangent_vectors(nx, ny, nz, &e1x, &e1y, &e1z, &e2x, &e2y, &e2z);
  for (x = 0; x < mri_src->width; x++) {
    for (y = 0; y < mri_src->height; y++) {
      for (z = 0; z < mri_src->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        deriv = 0;
        for (dn = -1; dn <= 1; dn++)  // normal direction
        {
          switch (dn) {
            case 1:
            case -1:
              wt = 1;
              break;
            default:
            case 0:
              wt = -2;
              break;
          }
          for (d1 = -1; d1 <= 1; d1++)    // first tangent direction
            for (d2 = -1; d2 <= 1; d2++)  // secont tangent direction
            {
              xf = x + dn * nx + d1 * e1x + d2 * e2x;
              yf = y + dn * ny + d1 * e1y + d2 * e2y;
              zf = z + dn * nz + d1 * e1z + d2 * e2z;
              MRIsampleVolume(mri_src, xf, yf, zf, &val);
              deriv += (val * wt);
            }
        }
        MRIsetVoxVal(mri_deriv, x, y, z, 0, deriv);
      }
    }
  }

  return (mri_deriv);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIsobelFrame(MRI *mri_src, MRI *mri_grad, MRI *mri_mag, int frame)
{
  MRI *mri_frame;

  mri_frame = MRIcopyFrame(mri_src, NULL, frame, 0);
  mri_grad = MRIsobel(mri_frame, mri_grad, mri_mag);
  MRIfree(&mri_frame);
  return (mri_grad);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIsobel(MRI *mri_src, MRI *mri_grad, MRI *mri_mag)
{
  int x, y, z, width, height, depth;
  float *px, *py, *pz, xval, yval, zval, *pmag;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_grad) {
    mri_grad = MRIallocSequence(width, height, depth, MRI_FLOAT, 3);
    MRIcopyHeader(mri_src, mri_grad);
  }

  MRIxSobel(mri_src, mri_grad, 0);
  MRIySobel(mri_src, mri_grad, 1);
  MRIzSobel(mri_src, mri_grad, 2);

  if (mri_mag) {
    if (mri_mag->type == MRI_FLOAT && mri_grad->type == MRI_FLOAT) {
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          px = &MRIFvox(mri_grad, 0, y, z);
          py = &MRIFseq_vox(mri_grad, 0, y, z, 1);
          pz = &MRIFseq_vox(mri_grad, 0, y, z, 2);
          pmag = &MRIFvox(mri_mag, 0, y, z);

          for (x = 0; x < width; x++) {
            xval = *px++;
            yval = *py++;
            zval = *pz++;
            *pmag++ = sqrt(xval * xval + yval * yval + zval * zval);
          }
        }
      }
    }
    else {
      for (z = 0; z < depth; z++) {
        float mag;

        for (y = 0; y < height; y++) {
          px = &MRIFvox(mri_grad, 0, y, z);
          py = &MRIFseq_vox(mri_grad, 0, y, z, 1);
          pz = &MRIFseq_vox(mri_grad, 0, y, z, 2);
          pmag = &MRIFvox(mri_mag, 0, y, z);

          for (x = 0; x < width; x++) {
            xval = MRIgetVoxVal(mri_grad, x, y, z, 0);
            yval = MRIgetVoxVal(mri_grad, x, y, z, 1);
            zval = MRIgetVoxVal(mri_grad, x, y, z, 2);
            mag = sqrt(xval * xval + yval * yval + zval * zval);
            MRIsetVoxVal(mri_mag, x, y, z, 0, mag);
          }
        }
      }
    }
  }

  return (mri_grad);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIxSobel(MRI *mri_src, MRI *mri_x, int frame)
{
  BUFTYPE *tr_pix, *mr_pix, *br_pix;
  int x, y, z, width, height, depth;
  float *outPtr, left, middle, right;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (mri_src->type != MRI_UCHAR) {
    return MRIxSobelForAllTypes(mri_src, mri_x, frame);
  }

  if (!mri_x) {
    mri_x = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_x);
  }

  if (mri_src->type != MRI_UCHAR || mri_x->type != MRI_FLOAT)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIxSobel: unsupported src or dst type (%d or %d)", mri_src->type, mri_x->type));

  /* applying sobel in x-y plane, don't worry about boundary conditions in z */
  /* don't apply sobel to outer ring to pixels to avoid border effects */
  width--;
  height--;
  for (z = 0; z < depth; z++) {
    for (y = 1; y < height; y++) {
      tr_pix = &MRIvox(mri_src, 0, y - 1, z);
      mr_pix = &MRIvox(mri_src, 0, y, z);
      br_pix = &MRIvox(mri_src, 0, y + 1, z);
      outPtr = &MRIFseq_vox(mri_x, 1, y, z, frame);

      left = (float)(*tr_pix++ + 2 * *mr_pix++ + *br_pix++);
      middle = (float)(*tr_pix++ + 2 * *mr_pix++ + *br_pix++);

      for (x = 1; x < width; x++) {
        right = (float)*tr_pix++ + 2.0f * (float)*mr_pix++ + (float)*br_pix++;

        *outPtr++ = (right - left) / 8.0f;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_x);
}

MRI *MRIxSobelForAllTypes(MRI *mri_src, MRI *mri_x, int frame)
{
  int tr_pix[3], mr_pix[3], br_pix[3];
  float tr_pix_val, mr_pix_val, br_pix_val;
  int x, y, z, width, height, depth;
  float *outPtr;
  float left, middle, right, two, eight;

  two = 2.0;
  eight = 8.0;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_x) {
    mri_x = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_x);
  }

  /* applying sobel in x-y plane */
  width--;
  height--;
  for (z = 0; z < depth; z++) {
    for (y = 1; y < height; y++) {
      tr_pix[0] = 0;
      tr_pix[1] = y - 1;
      tr_pix[2] = z;
      mr_pix[0] = 0;
      mr_pix[1] = y;
      mr_pix[2] = z;
      br_pix[0] = 0;
      br_pix[1] = y + 1;
      br_pix[2] = z;
      outPtr = &MRIFseq_vox(mri_x, 1, y, z, frame);

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      mr_pix_val = MRIgetVoxVal(mri_src, mr_pix[0], mr_pix[1], mr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      left = (tr_pix_val + two * mr_pix_val + br_pix_val);
      tr_pix[0]++;
      mr_pix[0]++;
      br_pix[0]++;

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      mr_pix_val = MRIgetVoxVal(mri_src, mr_pix[0], mr_pix[1], mr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      middle = (tr_pix_val + two * mr_pix_val + br_pix_val);
      tr_pix[0]++;
      mr_pix[0]++;
      br_pix[0]++;

      for (x = 1; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
        mr_pix_val = MRIgetVoxVal(mri_src, mr_pix[0], mr_pix[1], mr_pix[2], 0);
        br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
        right = (tr_pix_val + two * mr_pix_val + br_pix_val);
        tr_pix[0]++;
        mr_pix[0]++;
        br_pix[0]++;

        *outPtr++ = (right - left) / eight;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_x);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIySobel(MRI *mri_src, MRI *mri_y, int frame)
{
  BUFTYPE *tr_pix, *br_pix;
  int width, height, depth, x, y, z;
  float *outPtr, left, middle, right;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (mri_src->type != MRI_UCHAR) {
    return MRIySobelForAllTypes(mri_src, mri_y, frame);
  }

  if (!mri_y) {
    mri_y = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_y);
  }

  if (mri_src->type != MRI_UCHAR || mri_y->type != MRI_FLOAT)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIySobel: unsupported src or dst type (%d or %d)", mri_src->type, mri_y->type));
  /* applying sobel in x-y plane, don't worry about boundary conditions in z */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width--;
  height--;
  for (z = 0; z < depth; z++) {
    for (y = 1; y < height; y++) {
      tr_pix = &MRIvox(mri_src, 0, y - 1, z);
      br_pix = &MRIvox(mri_src, 0, y + 1, z);
      outPtr = &MRIFseq_vox(mri_y, 1, y, z, frame);

      left = (float)*br_pix++ - (float)*tr_pix++;
      middle = (float)*br_pix++ - (float)*tr_pix++;

      for (x = 1; x < width; x++) {
        right = (float)*br_pix++ - (float)*tr_pix++;
        *outPtr++ = (right + 2.0f * middle + left) / 8.0f;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_y);
}

MRI *MRIySobelForAllTypes(MRI *mri_src, MRI *mri_y, int frame)
{
  int tr_pix[3], br_pix[3];
  float tr_pix_val, br_pix_val;
  int x, y, z, width, height, depth;
  float *outPtr;
  float left, middle, right, two, eight;

  two = 2.0;
  eight = 8.0;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_y) {
    mri_y = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_y);
  }

  /* applying sobel in x-y plane */
  width--;
  height--;
  for (z = 0; z < depth; z++) {
    for (y = 1; y < height; y++) {
      tr_pix[0] = 0;
      tr_pix[1] = y - 1;
      tr_pix[2] = z;
      br_pix[0] = 0;
      br_pix[1] = y + 1;
      br_pix[2] = z;
      outPtr = &MRIFseq_vox(mri_y, 1, y, z, frame);

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      left = (br_pix_val - tr_pix_val);
      tr_pix[0]++;
      br_pix[0]++;

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      middle = (br_pix_val - tr_pix_val);
      tr_pix[0]++;
      br_pix[0]++;

      for (x = 1; x < width; x++) {
        if (Gx == x && y == Gy && z == Gz) {
          DiagBreak();
        }
        tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
        br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
        right = (br_pix_val - tr_pix_val);
        tr_pix[0]++;
        br_pix[0]++;

        *outPtr++ = (right + two * middle + left) / eight;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_y);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           calculate a z derivitive in the x-z plane
------------------------------------------------------*/
MRI *MRIzSobel(MRI *mri_src, MRI *mri_z, int frame)
{
  BUFTYPE *tr_pix, *br_pix;
  int width, height, depth, x, y, z;
  float *outPtr, left, middle, right;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (mri_src->type != MRI_UCHAR) {
    return MRIzSobelForAllTypes(mri_src, mri_z, frame);
  }

  if (!mri_z) {
    mri_z = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_z);
  }

  if (mri_src->type != MRI_UCHAR || mri_z->type != MRI_FLOAT)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIzSobel: unsupported src or dst type (%d or %d)", mri_src->type, mri_z->type));

  /* applying sobel in x-z plane, don't worry about boundary conditions in y */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width--;
  depth--;
  for (z = 1; z < depth; z++) {
    for (y = 0; y < height; y++) {
      tr_pix = &MRIvox(mri_src, 0, y, z - 1);
      br_pix = &MRIvox(mri_src, 0, y, z + 1);
      outPtr = &MRIFseq_vox(mri_z, 1, y, z, frame);

      left = (float)*br_pix++ - (float)*tr_pix++;
      middle = (float)*br_pix++ - (float)*tr_pix++;

      for (x = 1; x < width; x++) {
        right = (float)*br_pix++ - (float)*tr_pix++;
        *outPtr++ = (right + 2.0f * middle + left) / 8.0f;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_z);
}

MRI *MRIzSobelForAllTypes(MRI *mri_src, MRI *mri_z, int frame)
{
  int tr_pix[3], br_pix[3];
  float tr_pix_val, br_pix_val;
  int x, y, z, width, height, depth;
  float *outPtr;
  float left, middle, right, two, eight;

  two = 2.0;
  eight = 8.0;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_z) {
    mri_z = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_z);
  }

  /* applying sobel in x-z plane */
  width--;
  depth--;
  for (z = 1; z < depth; z++) {
    for (y = 0; y < height; y++) {
      tr_pix[0] = 0;
      tr_pix[1] = y;
      tr_pix[2] = z - 1;
      br_pix[0] = 0;
      br_pix[1] = y;
      br_pix[2] = z + 1;
      outPtr = &MRIFseq_vox(mri_z, 1, y, z, frame);

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      left = (br_pix_val - tr_pix_val);
      tr_pix[0]++;
      br_pix[0]++;

      tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
      br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
      middle = (br_pix_val - tr_pix_val);
      tr_pix[0]++;
      br_pix[0]++;

      for (x = 1; x < width; x++) {
        tr_pix_val = MRIgetVoxVal(mri_src, tr_pix[0], tr_pix[1], tr_pix[2], 0);
        br_pix_val = MRIgetVoxVal(mri_src, br_pix[0], br_pix[1], br_pix[2], 0);
        right = (br_pix_val - tr_pix_val);
        tr_pix[0]++;
        br_pix[0]++;

        *outPtr++ = (right + two * middle + left) / eight;
        left = middle;
        middle = right;
      }
    }
  }

  return (mri_z);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          compute a zscore for an image
------------------------------------------------------*/
MRI *MRIzScore(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_std)
{
  MRI_REGION region;

  region.x = region.y = region.z = 0;
  region.dx = mri_src->width;
  region.dx = mri_src->width;
  region.dx = mri_src->width;
  return (MRIzScoreRegion(mri_src, mri_dst, mri_mean, mri_std, &region));
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          compute a zscore for an image in a region. Note that
          mri_std and mri_mean are presumed to have been computed
          with MRImeanRegion and MRIstdRegion and are therefore
          in the same coordinate system as the destination image,
          while mri_src is presumed to be the original (full size)
          image.
------------------------------------------------------*/
MRI *MRIzScoreRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_std, MRI_REGION *region)
{
  int width, height, depth, x, y, z, x0, y0, z0;
  BUFTYPE *psrc;
  float *pmean, *pstd, std, *pdst;

  width = region->x + region->dx;
  if (width > mri_src->width) {
    width = mri_src->width;
  }
  height = region->y + region->dy;
  if (height > mri_src->height) {
    height = mri_src->height;
  }
  depth = region->z + region->dz;
  if (depth > mri_src->depth) {
    depth = mri_src->depth;
  }
  x0 = region->x;
  if (x0 < 0) {
    x0 = 0;
  }
  y0 = region->y;
  if (y0 < 0) {
    y0 = 0;
  }
  z0 = region->z;
  if (z0 < 0) {
    z0 = 0;
  }

  if (!mri_dst) {
    int w, h, d;

    w = width - region->x;
    h = height - region->y;
    d = depth - region->z;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize;
  }

  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      pmean = &MRIFvox(mri_mean, 0, y - y0, z - z0);
      pstd = &MRIFvox(mri_std, 0, y - y0, z - z0);
      pdst = &MRIFvox(mri_dst, 0, y - y0, z - z0);
      psrc = &MRIvox(mri_src, x0, y, z);
      for (x = x0; x < width; x++) {
        std = *pstd++;
        if (!std) {
          *pdst++ = 0.0f;
          psrc++;
          pmean++;
        }
        else {
          *pdst++ = ((float)*psrc++ - *pmean++) / std;
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
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *MRImeanRegion(MRI *mri_src, MRI *mri_dst, int wsize, MRI_REGION *region)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, val, xk, yk, zk, xi, yi, zi;
  float *pdst, wcubed;

  wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;
  width = region->x + region->dx;
  if (width > mri_src->width) {
    width = mri_src->width;
  }
  height = region->y + region->dy;
  if (height > mri_src->height) {
    height = mri_src->height;
  }
  depth = region->z + region->dz;
  if (depth > mri_src->depth) {
    depth = mri_src->depth;
  }
  x0 = region->x;
  if (x0 < 0) {
    x0 = 0;
  }
  y0 = region->y;
  if (y0 < 0) {
    y0 = 0;
  }
  z0 = region->z;
  if (z0 < 0) {
    z0 = 0;
  }

  if (!mri_dst) {
    int w, h, d;

    w = width - region->x;
    h = height - region->y;
    d = depth - region->z;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize;
  }

  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      pdst = &MRIFvox(mri_dst, 0, y - y0, z - z0);
      for (x = x0; x < width; x++) {
        for (val = 0, zk = -whalf; zk <= whalf; zk++) {
          zi = mri_src->zi[z + zk];
          for (yk = -whalf; yk <= whalf; yk++) {
            yi = mri_src->yi[y + yk];
            for (xk = -whalf; xk <= whalf; xk++) {
              xi = mri_src->xi[x + xk];
              val += (int)MRIvox(mri_src, xi, yi, zi);
            }
          }
        }
        *pdst++ = (float)val / wcubed;
      }
    }
  }
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           calculate the standard deviation in a region. Note that
           mri_src is a full image, but mri_mean is presumed to have
           been computed with the same region, and is therefore in
           the same coordinate system as the destination image.
------------------------------------------------------*/
MRI *MRIstdRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, int wsize, MRI_REGION *region)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xk, yk, zk, xi, yi, zi;
  float *pdst, wcubed, mean, *pmean, variance, f;

  whalf = wsize / 2;
  width = region->x + region->dx;
  if (width > mri_src->width) {
    width = mri_src->width;
  }
  height = region->y + region->dy;
  if (height > mri_src->height) {
    height = mri_src->height;
  }
  depth = region->z + region->dz;
  if (depth > mri_src->depth) {
    depth = mri_src->depth;
  }
  x0 = region->x;
  if (x0 < 0) {
    x0 = 0;
  }
  y0 = region->y;
  if (y0 < 0) {
    y0 = 0;
  }
  z0 = region->z;
  if (z0 < 0) {
    z0 = 0;
  }

  if (!mri_dst) {
    int w, h, d;

    w = width - region->x;
    h = height - region->y;
    d = depth - region->z;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize;
  }

  for (z = z0; z < depth; z++) {
    for (y = y0; y < height; y++) {
      pdst = &MRIFvox(mri_dst, 0, y - y0, z - z0);
      pmean = &MRIFvox(mri_mean, 0, y - y0, z - z0);
      for (x = x0; x < width; x++) {
        mean = *pmean++;
        wcubed = 0;
        for (variance = 0.0f, zk = -whalf; zk <= whalf; zk++) {
          zi = z + zk;
          if (zi < 0 || zi >= depth) {
            continue;
          }
          for (yk = -whalf; yk <= whalf; yk++) {
            yi = y + yk;
            if (yi < 0 || yi >= height) {
              continue;
            }
            for (xk = -whalf; xk <= whalf; xk++) {
              xi = x + xk;
              if (xi < 0 || xi >= width) {
                continue;
              }
              wcubed++;
              f = (float)MRIvox(mri_src, xi, yi, zi) - mean;
              variance += (f * f);
            }
          }
        }
        if (FZERO(wcubed)) {
          *pdst = 0.0;
        }
        else {
          *pdst++ = sqrt(variance / wcubed);
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
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *MRImeanByte(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, val;
  BUFTYPE *psrc, *pdst;
  float wcubed;

  wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_dst->type != MRI_UCHAR) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRImeanByte: dst must be MRI_UCHAR"));

  for (z = whalf; z < depth - whalf; z++) {
    for (y = whalf; y < height - whalf; y++) {
      pdst = &MRIvox(mri_dst, whalf, y, z);
      for (x = whalf; x < width - whalf; x++) {
        for (val = 0, z0 = -whalf; z0 <= whalf; z0++) {
          for (y0 = -whalf; y0 <= whalf; y0++) {
            psrc = &MRIvox(mri_src, x - whalf, y + y0, z + z0);
            for (x0 = -whalf; x0 <= whalf; x0++) {
              val += (int)*psrc++;
            }
          }
        }
        *pdst++ = (BUFTYPE)nint((float)val / wcubed);
      }
    }
  }

  /* now copy information to borders from source image */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < whalf; z++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
      for (z = depth - whalf; z < depth; z++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
    }
  }
  for (x = 0; x < width; x++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < whalf; y++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
      for (y = height - whalf; y < height; y++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
    }
  }
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < whalf; x++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
      for (x = width - whalf; x < width; x++) {
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src, x, y, z);
      }
    }
  }
  return (mri_dst);
}

MRI *MRImeanInMask(MRI *mri_src, MRI *mri_dst, MRI *mri_mask, int wsize)
{
  int width, height, depth, whalf, num;
  // float wcubed;
  int x, y, z, x0, y0, z0;
  float val;

  // wcubed = (float)(wsize * wsize * wsize);
  whalf = wsize / 2;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0) {
          continue;
        }
        for (num = 0, val = 0.0, z0 = -whalf; z0 <= whalf; z0++) {
          if (z + z0 < 0 || z + z0 >= mri_src->depth) {
            continue;
          }
          for (y0 = -whalf; y0 <= whalf; y0++) {
            if (y + y0 < 0 || y + y0 >= mri_src->height) {
              continue;
            }
            for (x0 = -whalf; x0 <= whalf; x0++) {
              if (x + x0 < 0 || x + x0 >= mri_src->width) {
                continue;
              }
              if (MRIgetVoxVal(mri_mask, x + x0, y + y0, z + z0, 0) == 0) {
                continue;
              }
              val += MRIgetVoxVal(mri_src, x + x0, y + y0, z + z0, 0);
              num++;
            }
          }
        }
        if (FZERO(num) == 0) {
          val /= num;
          MRIsetVoxVal(mri_dst, x, y, z, 0, val);
        }
      }
    }
  }

  return (mri_dst);
}

MRI *MRIstdInMask(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_mask, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi;
  float wcubed, mean, variance, f;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_dst->type != MRI_FLOAT) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIstd: dst must be MRI_FLOAT"));

  whalf = wsize / 2;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0) {
          continue;
        }
        mean = MRIgetVoxVal(mri_mean, x, y, z, 0);
        wcubed = 0;
        for (variance = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          zi = z + z0;
          if (zi < 0 || zi >= depth) {
            continue;
          }
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = y + y0;
            if (yi < 0 || yi >= height) {
              continue;
            }
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = x + x0;
              if (xi < 0 || xi >= width) {
                continue;
              }
              if (MRIgetVoxVal(mri_mask, xi, yi, zi, 0) == 0) {
                continue;
              }
              wcubed++;
              f = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              f -= mean;
              variance += (f * f);
            }
          }
        }
        if (wcubed == 0) {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, sqrt(variance / wcubed));
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
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *MRImean(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int width, height, depth, whalf;
  // float wcubed;

  // wcubed = (float)(wsize * wsize * wsize);
  whalf = (wsize - 1) / 2;
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, MRI_FLOAT, mri_src->nframes);
    MRIcopyHeader(mri_src, mri_dst);
  }

  {
    int frame, z;

    for (frame = 0; frame < mri_src->nframes; frame++) {
      ROMP_PF_begin
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(experimental)
#endif
      for (z = 0; z < depth; z++) {
        ROMP_PFLB_begin
        int x, y, x0, y0, z0;
        float val, num;

        for (y = 0; y < height; y++) {
          for (x = 0; x < width; x++) {
            for (num = 0, val = 0.0, z0 = -whalf; z0 <= whalf; z0++) {
              if (z + z0 < 0 || z + z0 >= mri_src->depth) continue;

              for (y0 = -whalf; y0 <= whalf; y0++) {
                if (y + y0 < 0 || y + y0 >= mri_src->height) continue;

                for (x0 = -whalf; x0 <= whalf; x0++) {
                  if (x + x0 < 0 || x + x0 >= mri_src->width) continue;

                  val += MRIgetVoxVal(mri_src, x + x0, y + y0, z + z0, frame);
                  num++;
                }
              }
            }
            if (FZERO(num) == 0) {
              val /= num;
              MRIsetVoxVal(mri_dst, x, y, z, frame, val);
            }
          }
        }
        exec_progress_callback(frame * depth + z, mri_src->nframes * depth, 0, 1);
        ROMP_PFLB_end
      }
      ROMP_PF_end
    }
  }

  return (mri_dst);
}

MRI *MRIinvert(MRI *mri_src, MRI *mri_dst)
{
  int width, height, depth, x, y, z;
  float val, max_val;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  max_val = 0;

  // figure out what to do with zeros
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        val = MRIgetVoxVal(mri_src, x, y, z, 0);
        if (!DZERO(val)) {
          val = 1 / val;
          if (val > max_val) {
            max_val = val;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, val);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);  // set inverse of 0 to 0
        }
      }
    }
  }
#if 0
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z==Gz)
        {
          DiagBreak() ;
        }
        val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        if (DZERO(val))
        {
          val = 2*max_val ;
          MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
        }
      }
    }
  }
#endif
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a median filter on the input MRI
------------------------------------------------------*/
MRI *MRImedian(MRI *mri_src, MRI *mri_dst, int wsize, MRI_REGION *box)
{
  static float *sort_array = NULL;
  static int sort_size = 0;
  int width, height, depth, x, y, z, whalf, x0, y0, z0, frame, median_index, wcubed, yi, zi, xmin, xmax, ymin, ymax,
      zmin, zmax;
  float *sptr;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_src->depth == 1)  // do a 2D median instead of 3D
  {
    int wsquared;

    printf("performing 2D median filter...\n");

    wsquared = (wsize * wsize);
    median_index = wsquared / 2;
    whalf = wsize / 2;

    if (sort_array && (wsquared != sort_size)) {
      free(sort_array);
      sort_array = NULL;
    }
    if (!sort_array) {
      sort_array = (float *)calloc(wsquared, sizeof(float));
      sort_size = wsquared;
    }

    if (box) {
      xmin = box->x;
      ymin = box->y;
      zmin = box->z;
      xmax = box->x + box->dx - 1;
      ymax = box->y + box->dy - 1;
      zmax = box->z + box->dz - 1;
    }
    else {
      xmin = ymin = zmin = 0;
      xmax = width - 1;
      ymax = height - 1;
      zmax = depth - 1;
    }

    for (frame = 0; frame < mri_src->nframes; frame++) {
      for (y = ymin; y <= ymax; y++) {
        for (x = xmin; x <= xmax; x++) {
          for (sptr = sort_array, y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              *sptr++ = MRIgetVoxVal(mri_src, mri_src->xi[x + x0], yi, 0, frame);
            }
          }
          qsort(sort_array, sort_size, sizeof(float), compare_sort_array);
          MRIsetVoxVal(mri_dst, x, y, 0, frame, sort_array[median_index]);
        }
        exec_progress_callback(frame * depth + y, mri_src->nframes * height, 0, 1);
      }
    }
    return (mri_dst);
  }

  wcubed = (wsize * wsize * wsize);
  median_index = wcubed / 2;
  whalf = wsize / 2;
  printf("wsize %d wcubed %d  whalf %d\n",wsize,wcubed,whalf);

  if (sort_array && (wcubed != sort_size)) {
    free(sort_array);
    sort_array = NULL;
  }
  if (!sort_array) {
    sort_array = (float *)calloc(wcubed, sizeof(float));
    sort_size = wcubed;
  }
  if (box) {
    xmin = box->x;
    ymin = box->y;
    zmin = box->z;
    xmax = box->x + box->dx - 1;
    ymax = box->y + box->dy - 1;
    zmax = box->z + box->dz - 1;
  }
  else {
    xmin = ymin = zmin = 0;
    xmax = width - 1;
    ymax = height - 1;
    zmax = depth - 1;
  }

  for (frame = 0; frame < mri_src->nframes; frame++) {
    for (z = zmin; z <= zmax; z++) {
      for (y = ymin; y <= ymax; y++) {
        for (x = xmin; x <= xmax; x++) {
          for (sptr = sort_array, z0 = -whalf; z0 <= whalf; z0++) {
            zi = mri_src->zi[z + z0];
            for (y0 = -whalf; y0 <= whalf; y0++) {
              yi = mri_src->yi[y + y0];
              for (x0 = -whalf; x0 <= whalf; x0++) {
                *sptr++ = MRIgetVoxVal(mri_src, mri_src->xi[x + x0], yi, zi, frame);
              }
            }
          }

          qsort(sort_array, wcubed, sizeof(float), compare_sort_array);
          MRIsetVoxVal(mri_dst, x, y, z, frame, sort_array[median_index]);
        }
      }
      exec_progress_callback(frame * depth + z, mri_src->nframes * depth, 0, 1);
    }
  }
  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a generalized median filter on the input MRI
------------------------------------------------------*/
MRI *MRIorder(MRI *mri_src, MRI *mri_dst, int wsize, float pct)
{
  static float *sort_array = NULL;
  static int sort_size = 0;
  int width, height, depth, x, y, z, whalf, x0, y0, z0, order_index, wcubed, yi, zi;
  float *sptr;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_dst->type != MRI_UCHAR) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIorder: dst must be MRI_UCHAR"));

  wcubed = (wsize * wsize * wsize);
  order_index = pct * wcubed;
  whalf = wsize / 2;

  if (sort_array && (wcubed != sort_size)) {
    free(sort_array);
    sort_array = NULL;
  }
  if (!sort_array) {
    sort_array = (float *)calloc(wcubed, sizeof(float));
    sort_size = wcubed;
  }
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        for (sptr = sort_array, z0 = -whalf; z0 <= whalf; z0++) {
          zi = mri_src->zi[z + z0];
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = mri_src->yi[y + y0];
            for (x0 = -whalf; x0 <= whalf; x0++) {
              *sptr++ = MRIgetVoxVal(mri_src, mri_src->xi[x + x0], yi, zi, 0);
            }
          }
        }

        qsort(sort_array, wcubed, sizeof(float), compare_sort_array);
        MRIsetVoxVal(mri_dst, x, y, z, 0, sort_array[order_index]);
      }
    }
  }
  return (mri_dst);
}
static int compare_sort_array(const void *pc1, const void *pc2)
{
  float c1, c2;

  c1 = *(float *)pc1;
  c2 = *(float *)pc2;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2) {
    return (1);
  }
  else if (c1 < c2) {
    return (-1);
  }

  return (0);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a splined gaussian bump which tails to 0.
             Returns an image which is (8*sigma)+1
             (Nitzberg and Shiota, 1993)
----------------------------------------------------------------------*/
MRI *MRIgaussian1d(float sigma, int max_len)
{
  MRI *mri;
  float norm, two_sigma, fx, k;
  int x, half, len;

  // printf( "%s: sigma=%30.10f\n", __FUNCTION__, sigma );

  /* build the kernel in k */
  len = (int)nint(8.0f * sigma) + 1;
  if (ISEVEN(len)) /* ensure it's even */
  {
    len++;
  }
  if (max_len > 0 && (max_len < len)) {
    len = max_len;
  }
  half = len / 2;
  if (len < 1) {
    len = 1;
  }
  mri = MRIalloc(len, 1, 1, MRI_FLOAT);

  if (len <= 1) {
    MRIFvox(mri, 0, 0, 0) = 1.0;
    return (mri);
  }
  norm = 0.0f;
  two_sigma = 2.0f * sigma;

  for (norm = 0.0f, x = 0; x < len; x++) {
    fx = (float)(x - half);
    if (fabs(fx) <= two_sigma) {
      k = (float)exp((double)(-fx * fx / (two_sigma * sigma)));
    }
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f * sigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0f - fabs(fx) / (double)sigma, 4.0);
    else {
      k = 0;
    }

    MRIFvox(mri, x, 0, 0) = k;
    norm += k;
  }

  /* normalize kernel to sum to 1 */
  for (x = 0; x < len; x++) {
    MRIFvox(mri, x, 0, 0) /= norm;
  }

  return (mri);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           calculate the standard deviation of the MRI
           in a wsize window around each point
------------------------------------------------------*/
MRI *MRIstd(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi;
  float wcubed, mean, variance, f;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_dst->type != MRI_FLOAT) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIstd: dst must be MRI_FLOAT"));

  whalf = wsize / 2;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        mean = MRIgetVoxVal(mri_mean, x, y, z, 0);
        wcubed = 0;
        for (variance = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          zi = z + z0;
          if (zi < 0 || zi >= depth) {
            continue;
          }
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = y + y0;
            if (yi < 0 || yi >= height) {
              continue;
            }
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = x + x0;
              if (xi < 0 || xi >= width) {
                continue;
              }
              wcubed++;
              f = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              f -= mean;
              variance += (f * f);
            }
          }
        }
        if (wcubed == 0) {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, sqrt(variance / wcubed));
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
           calculate the standard deviation of the MRI
           in a wsize window around each point
------------------------------------------------------*/
MRI *MRIstdNonzero(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, int wsize)
{
  int width, height, depth, x, y, z, whalf, x0, y0, z0, xi, yi, zi;
  float wcubed, mean, variance, f;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  if (mri_dst->type != MRI_FLOAT) ErrorReturn(mri_dst, (ERROR_UNSUPPORTED, "MRIstd: dst must be MRI_FLOAT"));

  whalf = wsize / 2;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        mean = MRIgetVoxVal(mri_mean, x, y, z, 0);
        if (FZERO(mean)) {
          continue;
        }
        wcubed = 0;
        for (variance = 0.0f, z0 = -whalf; z0 <= whalf; z0++) {
          zi = z + z0;
          if (zi < 0 || zi >= depth) {
            continue;
          }
          for (y0 = -whalf; y0 <= whalf; y0++) {
            yi = y + y0;
            if (yi < 0 || yi >= height) {
              continue;
            }
            for (x0 = -whalf; x0 <= whalf; x0++) {
              xi = x + x0;
              if (xi < 0 || xi >= width) {
                continue;
              }
              f = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              if (FZERO(f)) {
                continue;
              }
              wcubed++;
              f -= mean;
              variance += (f * f);
            }
          }
        }
        if (wcubed == 0) {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, sqrt(variance / wcubed));
        }
      }
    }
  }
  return (mri_dst);
}

/*-----------------------------------------------------
MRIconvolveGaussian() - see also MRIgaussianSmooth();
------------------------------------------------------*/
MRI *MRIconvolveGaussian(MRI *mri_src, MRI *mri_dst, MRI *mri_gaussian)
{
  // int width, height, depth;
  int klen;
  int frame;
  MRI *mtmp1, *mri_tmp;
  float *kernel;

#if 0
  static unsigned int nCalls = 0;
  char fname[STRLEN];

  snprintf( fname, STRLEN-1, "mrigauss%05u.mgz", nCalls );
  fname[STRLEN-1] = '\0';
  MRIwrite( mri_src, fname );

  nCalls++;
#endif

  kernel = &MRIFvox(mri_gaussian, 0, 0, 0);
  klen = mri_gaussian->width;
  // width = mri_src->width;
  // height = mri_src->height;
  // depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  if (mri_src->type != mri_dst->type) {
    fprintf(stderr, "%s: Source and destination types differ!\n", __FUNCTION__);
    exit(EXIT_FAILURE);
  }

  if (mri_dst == mri_src) {
    mri_tmp = mri_dst = MRIclone(mri_src, NULL);
  }
  else {
    mri_tmp = NULL;
  }

  int nstart = global_progress_range[0];
  int nend = global_progress_range[1];
  int nstep = (nstart - nend) / mri_src->nframes;
  mtmp1 = NULL;
  for (frame = 0; frame < mri_src->nframes; frame++) {
    global_progress_range[1] = global_progress_range[0] + nstep / 3;
    mtmp1 = MRIcopyFrame(mri_src, mtmp1, frame, 0);
    MRIconvolve1d(mri_src, mtmp1, kernel, klen, MRI_WIDTH, frame, 0);
    global_progress_range[0] += nstep / 3;
    global_progress_range[1] += nstep / 3;
    MRIconvolve1d(mtmp1, mri_dst, kernel, klen, MRI_HEIGHT, 0, frame);
    global_progress_range[0] += nstep / 3;
    global_progress_range[1] += nstep / 3;
    MRIconvolve1d(mri_dst, mtmp1, kernel, klen, MRI_DEPTH, frame, 0);

    MRIcopyFrame(mtmp1, mri_dst, 0, frame); /* convert it back to UCHAR */
    global_progress_range[0] = global_progress_range[1];
  }

  MRIfree(&mtmp1);
  MRIcopyHeader(mri_src, mri_dst);

  if (mri_tmp)  // src and dst are the same
  {
    MRIcopy(mri_tmp, mri_src);
    mri_dst = mri_src;
    MRIfree(&mri_tmp);
  }

  return (mri_dst);
}

/*---------------------------------------------------------------------
  MRIgaussianSmooth() - performs isotropic gaussian spatial smoothing.
  The standard deviation of the gaussian is std. If norm must be set to 1.
  This function calls the non-isotropic gaussian smoother MRIgaussianSmoothNI().
  See also MRIconvolveGaussian() and MRImaskedGaussianSmooth().
  -------------------------------------------------------------------*/
MRI *MRIgaussianSmooth(MRI *src, double std, int norm, MRI *targ)
{
  if (!norm) {
    printf("ERROR: MRIgaussianSmooth: norm must be non-zero\n");
    return (NULL);
  }
  targ = MRIgaussianSmoothNI(src, std, std, std, targ);
  return (targ);
}

/*----------------------------------------------------------------------------
  MRImaskedGaussianSmooth() - smooths only within the binary mask.
  Source values outside the mask do not contribute to target values
  inside the mask. Target values outside the mask are set to 0. The
  target voxels are rescaled so that kernel at each voxel sums to 1.
  This only has an effect near the edge of the mask. In a noise data
  set, this will increase the variance near the edge, but the expected
  value will be correct. The mask *must* be binary (ie, 0 and 1). OK
  if the mask is NULL. As a test, you can input the mask as both
  the src and binmask. The targ should end up being the same as
  the mask. See also MRIgaussianSmooth().
  ---------------------------------------------------------------------------*/
MRI *MRImaskedGaussianSmooth(MRI *src, MRI *binmask, double std, MRI *targ)
{
  MRI *binmasksm, *srcmasked;
  int c, r, s, f;
  double m, v;

  // If the mask is null, just smooth unmasked source and return
  if (binmask == NULL) {
    targ = MRIgaussianSmooth(src, std, 1, targ);  // 1 means sum(g)=1
    return (targ);
  }

  // Mask the source so that values outside the mask are 0
  // so that they will not affect the values inside the mask
  // srcmasked = MRImask(src,binmask,NULL,0.0,0.0);
  // Can't use MRImask here because it might run
  // MRImaskDifferentGeometry() if the geometry is even slightly off
  // so just do it by hand here

  srcmasked = MRIcopy(src, NULL);
  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        m = MRIgetVoxVal(binmask, c, r, s, 0);
        if (m >= 0.5) continue;
        for (f = 0; f < src->nframes; f++) MRIsetVoxVal(srcmasked, c, r, s, f, 0.0);
      }
    }
  }

  // Smooth the masked source
  targ = MRIgaussianSmooth(srcmasked, std, 1, targ);  // 1 means sum(g)=1

  // At this point, the voxels inside the mask are smoothed, but the
  // voxels near the edge of the mask will be scaled incorrectly (ie,
  // the kernel weights will not sum to 1). Also, voxels outside the
  // mask will be non-zero. Both of these are fixed below.

  // Smooth the binary mask - this gives us the correction factor to
  // make the kernel sum to 1.
  binmasksm = MRIgaussianSmooth(binmask, std, 1, NULL);  // 1 makes sum(g)=1

  for (c = 0; c < src->width; c++) {
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        m = MRIgetVoxVal(binmask, c, r, s, 0);
        if (m < 0.5) {
          // If outside the mask, set target voxel to 0
          for (f = 0; f < src->nframes; f++) {
            MRIsetVoxVal(targ, c, r, s, f, 0.0);
          }
          continue;
        }
        else {
          // If inside the mask, divide target voxel by smoothed mask value
          m = MRIgetVoxVal(binmasksm, c, r, s, 0);
          for (f = 0; f < src->nframes; f++) {
            v = MRIgetVoxVal(targ, c, r, s, f);
            MRIsetVoxVal(targ, c, r, s, f, v / m);
          }
        }
      }
    }
  }
  MRIfree(&srcmasked);
  MRIfree(&binmasksm);
  return (targ);
}

/**
 * MRImodifySampledHeader
 *
 * dst volume is sampled by 2 and thus copying direction cosines
 * but cras and other info changed. modify cached transform also.
 *
 * @param src
 * @param dst
 */
static void MRImodifySampledHeader(MRI *mri_src, MRI *mri_dst)
{
  double c_r, c_a, c_s;

  if (!mri_dst) {
    ErrorExit(ERROR_BADPARM, "dst volume must be given");
  }

  MRIcopyHeader(mri_src, mri_dst);

  mri_dst->thick = mri_src->thick * 2.0f;
  mri_dst->ps = mri_src->ps * 2.0f;
  mri_dst->scale = mri_src->scale * 2;
  mri_dst->xsize = mri_src->xsize * 2.0f;
  mri_dst->ysize = mri_src->ysize * 2.0f;
  mri_dst->zsize = mri_src->zsize * 2.0f;

  // calculate c_ras value
  MRIcalcCRASforSampledVolume(mri_src, mri_dst, &c_r, &c_a, &c_s);
  mri_dst->c_a = c_a;
  mri_dst->c_s = c_s;
  mri_dst->c_r = c_r;

  // header modified update cached info
  MRIreInitCache(mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          reduce a special type of MR image which contains
          two frames. The first is a frame of means and
          can be reduced normally. The second is a frame
          of standard deviations, and must be turned into
          variances, reduced, then back to stds and appended
          to the mean image.
------------------------------------------------------*/
MRI *MRIreduceMeanAndStdByte(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_var, *mri_var_reduced, *mri_means_reduced, *mri_std_reduced;
  int nframes;

  if (mri_src->nframes < 2) /* just reduce it normally */
  {
    return (MRIreduceByte(mri_src, mri_dst));
  }

  /* this is a hack, but only want to reduce the first frame of mri_src,
     as the 2nd is stds which need to be convolved as variances. Not
     doing this would work but would take twice as much time and memory.
  */
  nframes = mri_src->nframes;
  mri_src->nframes = 1;
  mri_means_reduced = MRIreduceByte(mri_src, NULL);
  mri_src->nframes = nframes;
  mri_var = MRIstdsToVariances(mri_src, NULL, 1);
  mri_var_reduced = MRIreduceByte(mri_var, NULL);
  mri_std_reduced = MRIvariancesToStds(mri_var_reduced, NULL, 0);
  mri_dst = MRIconcatenateFrames(mri_means_reduced, mri_std_reduced, mri_dst);

  MRIfree(&mri_var);
  MRIfree(&mri_means_reduced);
  MRIfree(&mri_var_reduced);
  MRIfree(&mri_std_reduced);

  MRImodifySampledHeader(mri_src, mri_dst);

  mri_dst->mean = MRImeanFrame(mri_dst, 1);
#if 0
  {
    char fname[100] ;
    sprintf(fname, "means_and_stds%d.mnc", (int)mri_dst->thick) ;
    fprintf(stderr, "writing means and stds to %s\n", fname) ;
    MRIwrite(mri_dst, fname) ;
  }
#endif
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          reduce a special type of MR image which contains
          two frames. The first is a frame of means and
          can be reduced normally. The second is a frame
          of standard deviations, and must be turned into
          variances, reduced, then back to stds and appended
          to the mean image.
------------------------------------------------------*/
MRI *MRIreduceMeanAndStd(MRI *mri_src, MRI *mri_dst)
{
  MRI *mri_var, *mri_var_reduced, *mri_means_reduced, *mri_std_reduced;
  int nframes;

  if (mri_src->nframes < 2) /* just reduce it normally */
  {
    return (MRIreduce(mri_src, mri_dst));
  }

  /* this is a hack, but only want to reduce the first frame of mri_src,
     as the 2nd is stds which need to be convolved as variances. Not
     doing this would work but would take twice as much time and memory.
  */
  nframes = mri_src->nframes;
  mri_src->nframes = 1;
  mri_means_reduced = MRIreduce(mri_src, NULL);
  mri_src->nframes = nframes;
  mri_var = MRIstdsToVariances(mri_src, NULL, 1);
  mri_var_reduced = MRIreduce(mri_var, NULL);
  mri_std_reduced = MRIvariancesToStds(mri_var_reduced, NULL, 0);
  mri_dst = MRIconcatenateFrames(mri_means_reduced, mri_std_reduced, mri_dst);

  MRIfree(&mri_var);
  MRIfree(&mri_means_reduced);
  MRIfree(&mri_var_reduced);
  MRIfree(&mri_std_reduced);

  MRImodifySampledHeader(mri_src, mri_dst);

  mri_dst->mean = MRImeanFrame(mri_dst, 1);
#if 0
  {
    char fname[100] ;
    sprintf(fname, "means_and_stds%d.mnc", (int)mri_dst->thick) ;
    fprintf(stderr, "writing means and stds to %s\n", fname) ;
    MRIwrite(mri_dst, fname) ;
  }
#endif
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define KERNEL_SIZE 5
#define K_A 0.4f
static float kernel[KERNEL_SIZE] = {0.25f - K_A / 2.0f, .25f, K_A, 0.25f, 0.25f - K_A / 2.0f};
MRI *MRIreduce(MRI *mri_src, MRI *mri_dst)
{
  int width, height, depth;
  MRI *mtmp1, *mtmp2 = NULL;
  MATRIX *m_vox2ras, *m_scale, *m_tmp;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (width <= 1 || height <= 1 /* || depth <= 1*/)
    ErrorExit(ERROR_BADPARM, "MRIreduce: insufficient dimension (%d, %d, %d)", width, height, depth);

  mtmp1 = MRIconvolve1d(mri_src, NULL, kernel, KERNEL_SIZE, MRI_WIDTH, 0, 0);
  if (depth > 1) {
    // again MRI_FLOAT volume
    mtmp2 = MRIconvolve1d(mtmp1, NULL, kernel, KERNEL_SIZE, MRI_HEIGHT, 0, 0);
    // half dimensions of the original
    mri_dst = MRIreduce1d(mtmp2, mri_dst, kernel, KERNEL_SIZE, MRI_DEPTH);
    // null then create MRI_FLOAT volume
    MRIfree(&mtmp1);
  }
  else {
    mri_dst = MRIreduce1d(mtmp1, mri_dst, kernel, KERNEL_SIZE, MRI_HEIGHT);
  }

  MRIcopyHeader(mri_src, mri_dst);
  MRIsetResolution(mri_dst, mri_src->xsize * 2, mri_src->ysize * 2, depth == 1 ? mri_src->zsize : mri_src->zsize * 2);
  mri_dst->xstart = mri_src->xstart;
  mri_dst->ystart = mri_src->ystart;
  mri_dst->zstart = mri_src->zstart;
  mri_dst->xend = mri_src->xend;
  mri_dst->yend = mri_src->yend;
  mri_dst->zend = mri_src->zend;
  m_vox2ras = MRIgetVoxelToRasXform(mri_src);
  m_scale = MatrixIdentity(4, NULL);
  *MATRIX_RELT(m_scale, 1, 1) = 2.0;
  *MATRIX_RELT(m_scale, 2, 2) = 2.0;
  *MATRIX_RELT(m_scale, 3, 3) = depth == 1 ? 1.0 : 2.0;
  m_tmp = MatrixMultiply(m_vox2ras, m_scale, NULL);
  MatrixFree(&m_vox2ras);
  MatrixFree(&m_scale);
  m_vox2ras = m_tmp;
  MRIsetVoxelToRasXform(mri_dst, m_vox2ras);
  MatrixFree(&m_vox2ras);
  //  MRImodifySampledHeader(mri_src, mri_dst);

  if (mtmp2) MRIfree(&mtmp2);

  return (mri_dst);
}

MRI *MRIreduce2D(MRI *mri_src, MRI *mri_dst)
{
  int width, height;
  // int depth;
  MRI *mtmp1;

  width = mri_src->width;
  height = mri_src->height;
  // depth = mri_src->depth;

  if (width <= 1 || height <= 1)
    ErrorExit(ERROR_BADPARM, "MRIreduce2D: insufficient dimension (%d, %d)", width, height);

  mtmp1 = MRIconvolve1d(mri_src, NULL, kernel, KERNEL_SIZE, MRI_WIDTH, 0, 0);
  mri_dst = MRIreduceSlice(mtmp1, mri_dst, kernel, KERNEL_SIZE, MRI_HEIGHT);
  MRIfree(&mtmp1);

  MRImodifySampledHeader(mri_src, mri_dst);

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Sae as above, but produce a byte image as output.
------------------------------------------------------*/
MRI *MRIreduceByte(MRI *mri_src, MRI *mri_dst)
{
  int width, height, depth;
  MRI *mtmp1, *mtmp2;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (width <= 1 || height <= 1 || depth <= 1)
    ErrorExit(ERROR_BADPARM, "MRIreduce: insufficient dimension (%d, %d, %d)", width, height, depth);

  mtmp1 = MRIalloc(width, height, depth, mri_src->type);
  MRIconvolve1d(mri_src, mtmp1, kernel, KERNEL_SIZE, MRI_WIDTH, 0, 0);
  mtmp2 = MRIclone(mtmp1, NULL);
  MRIconvolve1d(mtmp1, mtmp2, kernel, KERNEL_SIZE, MRI_HEIGHT, 0, 0);
  MRIfree(&mtmp1);
  mri_dst = MRIreduce1dByte(mtmp2, mri_dst, kernel, KERNEL_SIZE, MRI_DEPTH);

  MRImodifySampledHeader(mri_src, mri_dst);

  MRIfree(&mtmp2);
  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis, int src_frame, int dst_frame)
{
  int width, height, depth;
  int x = 0, y = 0, z = 0, halflen, *xi, *yi, *zi;
  int i = 0;
  BUFTYPE *inBase = NULL;
  float *ki = NULL, total = 0, *inBase_f = NULL, *foutPix = NULL, val = 0;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  // if dimension in convolve direction is 1, skip convolving:
  if ((axis == MRI_WIDTH && width == 1) || (axis == MRI_HEIGHT && height == 1) || (axis == MRI_DEPTH && depth == 1)) {
    mri_dst = MRIcopy(mri_src, mri_dst);
    return mri_dst;
  }

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
  }

  if (mri_dst->type == MRI_UCHAR) {
    return (MRIconvolve1dByte(mri_src, mri_dst, k, len, axis, src_frame, dst_frame));
  }
  else if (mri_dst->type == MRI_SHORT) {
    return (MRIconvolve1dShort(mri_src, mri_dst, k, len, axis, src_frame, dst_frame));
  }
  else if (mri_dst->type == MRI_INT) {
    return (MRIconvolve1dInt(mri_src, mri_dst, k, len, axis, src_frame, dst_frame));
  }

  if (mri_dst->type != MRI_FLOAT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1d: unsupported dst pixel format %d", mri_dst->type));

  halflen = len / 2;

  xi = mri_src->xi;
  yi = mri_src->yi;
  zi = mri_src->zi;

  switch (mri_src->type) {
    case MRI_UCHAR:
      switch (axis) {
        case MRI_WIDTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, inBase, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_src, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              inBase = &MRIseq_vox(mri_src, 0, y, z, src_frame);
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  total += *ki++ * (float)(*(inBase + xi[x + i - halflen]));
                }

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
	  
          break;
        case MRI_HEIGHT:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                if (x == Gx && y == Gy && z == Gz) {
                  DiagBreak();
                }
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * (float)(MRIseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame));

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
        case MRI_DEPTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * (float)(MRIseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame));

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
      }
      break;
    case MRI_FLOAT:
      switch (axis) {
        case MRI_WIDTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(assume_reproducible) firstprivate(y, x, inBase_f, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              inBase_f = &MRIFseq_vox(mri_src, 0, y, z, src_frame);
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                if (Gx == x && Gy == y && Gz == z) {
                  DiagBreak();
                }
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  total += *ki++ * (*(inBase_f + xi[x + i - halflen]));
                }

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
        case MRI_HEIGHT:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(assume_reproducible) firstprivate(y, x, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                if (x == Gx && y == Gy && z == Gz) {
                  DiagBreak();
                }

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * MRIFseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame);

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
        case MRI_DEPTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(assume_reproducible) firstprivate(y, x, foutPix, ki, i, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * MRIFseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame);

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
      }
      break;
    default:
      switch (axis) {
        case MRI_WIDTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, foutPix, ki, i, val, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  val = MRIgetVoxVal(mri_src, x + i - halflen, y, z, src_frame);
                  total += *ki++ * val;
                }

                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
        case MRI_HEIGHT:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, foutPix, ki, i, val, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                if (x == Gx && y == Gy && z == Gz) {
                  DiagBreak();
                }
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  val = MRIgetVoxVal(mri_src, x, y + i - halflen, z, src_frame);
                  total += *ki++ * val;
                }
                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
        case MRI_DEPTH:
          ROMP_PF_begin
#ifdef HAVE_OPENMP
	  #pragma omp parallel for if_ROMP(experimental) firstprivate(y, x, foutPix, ki, i, val, total) \
    shared(depth, height, width, len, halflen, mri_dst, src_frame, dst_frame, xi, yi, zi) schedule(static, 1)
#endif
          for (z = 0; z < depth; z++) {
	    ROMP_PFLB_begin
	    
            for (y = 0; y < height; y++) {
              foutPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  val = MRIgetVoxVal(mri_src, x, y, z + i - halflen, src_frame);
                  total += *ki++ * val;
                }
                *foutPix++ = total;
              }
            }
            exec_progress_callback(z, depth, 0, 1);
	    ROMP_PFLB_end
          }
	  ROMP_PF_end
          break;
      }
      break;
  }

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIconvolve1dByte(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis, int src_frame, int dst_frame)
{
  int x, y, z, width, height, halflen, depth, *xi, *yi, *zi;
  int i;
  BUFTYPE *inBase;
  BUFTYPE *inBase_f, *outPix;
  float *ki, total;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
  }

  if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1dByte: unsupported dst pixel format %d", mri_dst->type));

  halflen = len / 2;

  xi = mri_src->xi;
  yi = mri_src->yi;
  zi = mri_src->zi;

  switch (mri_src->type) {
    case MRI_UCHAR:
      switch (axis) {
        case MRI_WIDTH:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              inBase = &MRIseq_vox(mri_src, 0, y, z, src_frame);
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  total += *ki++ * (float)(*(inBase + xi[x + i - halflen]));
                }

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
        case MRI_HEIGHT:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * (float)(MRIseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame));

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
        case MRI_DEPTH:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * (float)(MRIseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame));

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
      }
      break;
    case MRI_FLOAT:
      switch (axis) {
        case MRI_WIDTH:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              inBase_f = &MRIseq_vox(mri_src, 0, y, z, src_frame);
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++) {
                  total += *ki++ * (*(inBase_f + xi[x + i - halflen]));
                }

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
        case MRI_HEIGHT:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * MRIseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame);

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
        case MRI_DEPTH:
          for (z = 0; z < depth; z++) {
            for (y = 0; y < height; y++) {
              outPix = &MRIseq_vox(mri_dst, 0, y, z, dst_frame);
              for (x = 0; x < width; x++) {
                total = 0.0f;

                for (ki = k, i = 0; i < len; i++)
                  total += *ki++ * MRIseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame);

                *outPix++ = (BUFTYPE)nint(total);
              }
            }
          }
          break;
      }
      break;
    default:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1d: unsupported pixel format %d", mri_src->type));
  }

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIconvolve1dShort(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis, int src_frame, int dst_frame)
{
  int x, y, z, width, height, halflen, depth, *xi, *yi, *zi;
  int i;
  short *inBase;
  short *outPix;
  float *ki, total;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_SHORT);
  }

  if (mri_dst->type != MRI_SHORT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1dShort: unsupported dst pixel format %d", mri_dst->type));

  halflen = len / 2;

  xi = mri_src->xi;
  yi = mri_src->yi;
  zi = mri_src->zi;

  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          inBase = &MRISseq_vox(mri_src, 0, y, z, src_frame);
          outPix = &MRISseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++) {
              total += *ki++ * (float)(*(inBase + xi[x + i - halflen]));
            }

            *outPix++ = (short)nint(total);
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRISseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRISseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame));

            *outPix++ = (short)nint(total);
          }
        }
      }
      break;
    case MRI_DEPTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRISseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRISseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame));

            *outPix++ = (short)nint(total);
          }
        }
      }
      break;
  }

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIconvolve1dInt(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis, int src_frame, int dst_frame)
{
  int x, y, z, width, height, halflen, depth, *xi, *yi, *zi;
  int i;
  int *inBase;
  int *outPix;
  float *ki, total;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_INT);
  }

  if (mri_dst->type != MRI_INT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1dInt: unsupported dst pixel format %d", mri_dst->type));

  halflen = len / 2;

  xi = mri_src->xi;
  yi = mri_src->yi;
  zi = mri_src->zi;

  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          inBase = &MRIIseq_vox(mri_src, 0, y, z, src_frame);
          outPix = &MRIIseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++) {
              total += *ki++ * (float)(*(inBase + xi[x + i - halflen]));
            }

            *outPix++ = (int)nint(total);
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRIIseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRIIseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame));

            *outPix++ = (int)nint(total);
          }
        }
      }
      break;
    case MRI_DEPTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRIIseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRIIseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame));

            *outPix++ = (int)nint(total);
          }
        }
      }
      break;
  }

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIconvolve1dFloat(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis, int src_frame, int dst_frame)
{
  int x, y, z, width, height, halflen, depth, *xi, *yi, *zi;
  int i;
  float *inBase;
  float *outPix;
  float *ki, total;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT);
  }

  if (mri_dst->type != MRI_FLOAT)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIconvolve1dFloat: unsupported dst pixel format %d", mri_dst->type));

  halflen = len / 2;

  xi = mri_src->xi;
  yi = mri_src->yi;
  zi = mri_src->zi;

  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          inBase = &MRIFseq_vox(mri_src, 0, y, z, src_frame);
          outPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++) {
              total += *ki++ * (float)(*(inBase + xi[x + i - halflen]));
            }

            *outPix++ = total;
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRIFseq_vox(mri_src, x, yi[y + i - halflen], z, src_frame));

            *outPix++ = total;
          }
        }
      }
      break;
    case MRI_DEPTH:
      for (z = 0; z < depth; z++) {
        for (y = 0; y < height; y++) {
          outPix = &MRIFseq_vox(mri_dst, 0, y, z, dst_frame);
          for (x = 0; x < width; x++) {
            total = 0.0f;

            for (ki = k, i = 0; i < len; i++)
              total += *ki++ * (float)(MRIFseq_vox(mri_src, x, y, zi[z + i - halflen], src_frame));

            *outPix++ = total;
          }
        }
      }
      break;
  }

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIreduce1d(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis)
{
  int x, y, z, i, dheight, dwidth, ddepth, xi, yi, zi, halflen;
  int sheight, swidth, sdepth;
  float total, val;

  swidth = mri_src->width;
  sheight = mri_src->height;
  sdepth = mri_src->depth;

  if (!mri_dst) {
#if 0
    mri_dst = MRIalloc(swidth/2, sheight/2, sdepth/2, MRI_UCHAR) ;
#else
    mri_dst = MRIalloc(swidth / 2, sheight / 2, sdepth / 2, mri_src->type);
#endif
    mri_dst->imnr1 = mri_dst->imnr0 + mri_dst->depth - 1;
    MRIsetResolution(mri_dst, mri_src->xsize * 2, mri_src->ysize * 2, mri_src->zsize * 2);
  }

#if 0
  if (((mri_dst->type != MRI_UCHAR) &&
       (mri_dst->type != MRI_FLOAT)) ||
      (mri_src->type != MRI_FLOAT))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRIreduce1d: src %d or dst %d format unsupported",
                 mri_src->type, mri_dst->type)) ;
#endif

  dwidth = mri_dst->width;
  dheight = mri_dst->height;
  ddepth = mri_dst->depth;

  halflen = (len - 1) / 2;
  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          yi = 2 * y;
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              xi = 2 * x + i - halflen;
              if (xi < 0) {
                xi = 0;
              }
              else if (xi >= swidth) {
                xi = swidth - 1;
              }

              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              total = total + k[i] * val;
            }
            MRIsetVoxVal(mri_dst, x, y, z, 0, total);
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;
            xi = 2 * x;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              yi = 2 * y + i - halflen;
              if (yi < 0) {
                yi = 0;
              }
              else if (yi >= sheight) {
                yi = sheight - 1;
              }

              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              total = total + k[i] * val;
            }
            MRIsetVoxVal(mri_dst, x, y, z, 0, total);
          }
        }
      }
      break;
    case MRI_DEPTH:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          yi = 2 * y;
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;
            xi = 2 * x;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              zi = 2 * z + i - halflen;
              if (zi < 0) {
                zi = 0;
              }
              else if (zi >= sdepth) {
                zi = sdepth - 1;
              }

              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              total = total + k[i] * val;
            }
            MRIsetVoxVal(mri_dst, x, y, z, 0, total);
          }
        }
      }
      break;
  }

  /*  mri_dst = MRIcopy(mri_src, mri_dst) ;*/

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIreduceSlice(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis)
{
  int x, y, z, i, dheight, dwidth, ddepth, xi, yi, zi, halflen;
  int sheight, swidth, sdepth;
  float total, val;

  swidth = mri_src->width;
  sheight = mri_src->height;
  sdepth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(swidth / 2, sheight / 2, sdepth, MRI_UCHAR);
    mri_dst->imnr1 = mri_dst->imnr0 + mri_dst->depth - 1;
    MRIsetResolution(mri_dst, mri_src->xsize * 2, mri_src->ysize * 2, mri_src->zsize);
  }

#if 0
  if (((mri_dst->type != MRI_UCHAR) &&
       (mri_dst->type != MRI_FLOAT)) ||
      (mri_src->type != MRI_FLOAT))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRIreduce1d: src %d or dst %d format unsupported",
                 mri_src->type, mri_dst->type)) ;
#endif

  dwidth = mri_dst->width;
  dheight = mri_dst->height;
  ddepth = mri_dst->depth;

  halflen = (len - 1) / 2;
  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < ddepth; z++) {
        zi = z;
        for (y = 0; y < dheight; y++) {
          yi = 2 * y;
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              xi = 2 * x + i - halflen;
              if (xi < 0) {
                xi = 0;
              }
              else if (xi >= swidth) {
                xi = swidth - 1;
              }

              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              total = total + k[i] * val;
            }
            MRIsetVoxVal(mri_dst, x, y, z, 0, total);
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < ddepth; z++) {
        zi = z;
        for (y = 0; y < dheight; y++) {
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;
            xi = 2 * x;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              yi = 2 * y + i - halflen;
              if (yi < 0) {
                yi = 0;
              }
              else if (yi >= sheight) {
                yi = sheight - 1;
              }

              val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
              total = total + k[i] * val;
            }
            MRIsetVoxVal(mri_dst, x, y, z, 0, total);
          }
        }
      }
      break;
  }

  /*  mri_dst = MRIcopy(mri_src, mri_dst) ;*/

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *MRIreduce1dByte(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis)
{
  int x, y, z, i, dheight, dwidth, ddepth, xi, yi, zi, halflen;
  int sheight, swidth, sdepth;
  float total;

  swidth = mri_src->width;
  sheight = mri_src->height;
  sdepth = mri_src->depth;

  if (!mri_dst) {
    mri_dst = MRIalloc(swidth / 2, sheight / 2, sdepth / 2, MRI_UCHAR);
    mri_dst->imnr1 = mri_dst->imnr0 + mri_dst->depth - 1;
  }

  if ((mri_dst->type != MRI_UCHAR) || (mri_src->type != MRI_UCHAR))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, "MRIreduce1d: src %d or dst %d format unsupported", mri_src->type, mri_dst->type));

  dwidth = mri_dst->width;
  dheight = mri_dst->height;
  ddepth = mri_dst->depth;

  halflen = (len - 1) / 2;
  switch (axis) {
    case MRI_WIDTH:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          yi = 2 * y;
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              xi = 2 * x + i - halflen;
              if (xi < 0) {
                xi = 0;
              }
              else if (xi >= swidth) {
                xi = swidth - 1;
              }

              total = total + k[i] * (float)MRIvox(mri_src, xi, yi, zi);
            }
            MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total);
          }
        }
      }
      break;
    case MRI_HEIGHT:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;
            xi = 2 * x;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              yi = 2 * y + i - halflen;
              if (yi < 0) {
                yi = 0;
              }
              else if (yi >= sheight) {
                yi = sheight - 1;
              }

              total = total + k[i] * (float)MRIvox(mri_src, xi, yi, zi);
            }
            MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total);
          }
        }
      }
      break;
    case MRI_DEPTH:
      for (z = 0; z < ddepth; z++) {
        zi = 2 * z;
        for (y = 0; y < dheight; y++) {
          yi = 2 * y;
          for (x = 0; x < dwidth; x++) {
            total = 0.0f;
            xi = 2 * x;

            for (i = 0; i < len; i++) {
              /* Neumann boundary conditions */
              zi = 2 * z + i - halflen;
              if (zi < 0) {
                zi = 0;
              }
              else if (zi >= sdepth) {
                zi = sdepth - 1;
              }

              total = total + k[i] * (float)MRIvox(mri_src, xi, yi, zi);
            }
            MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total);
          }
        }
      }
      break;
  }

  /*  mri_dst = MRIcopy(mri_src, mri_dst) ;*/

  return (mri_dst);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate two MRIs, limitting the extent of the
          correlation to window_size
------------------------------------------------------*/
MRI *MRIxcorrWindow(MRI *mri_ref, MRI *mri_in, MRI *mri_dst, int window_size)
{
  int width, height, depth, xoff, yoff, zoff, x0, y0, z0, kx, ky, kz, whalf, x, y, z, dst_width, dst_height, dst_depth,
      kxend;
  BUFTYPE *pref, *pin, ref, in;
  float dst, *pdst;

  if (!mri_dst) {
    /* decide on approriate size for destination */
    dst_width = window_size /* MAX(mri_ref->width, window_size)*/;
    dst_height = window_size /*MAX(mri_ref->height, window_size)*/;
    dst_depth = window_size /*MAX(mri_ref->depth, window_size)*/;

    mri_dst = MRIalloc(dst_width, dst_height, dst_depth, MRI_FLOAT);
    MRIcopyHeader(mri_ref, mri_dst);
  }

  dst_width = mri_dst->width;
  dst_height = mri_dst->height;
  dst_depth = mri_dst->depth;

  height = mri_ref->height;
  width = mri_ref->width;
  depth = mri_ref->depth;

  if (window_size <= 0) {
    window_size = width;
  }
  whalf = (window_size - 1) / 2;

  x0 = (dst_width - 1) / 2;
  y0 = (dst_height - 1) / 2;
  z0 = (dst_depth - 1) / 2;

  /*
     for each element of the destination array, correlate the reference
     MRI at that offset with the (unshifted) input MRI
     */
  for (zoff = -whalf; zoff <= whalf; zoff++) {
    for (yoff = -whalf; yoff <= whalf; yoff++) {
      pdst = &MRIFvox(mri_dst, x0 - whalf, y0 + yoff, z0 + zoff);
      for (xoff = -whalf; xoff <= whalf; xoff++) {
        /* apply offset template (ref) to input */
        dst = 0.0f;

        /* now calculate value of cross-correlation at these offsets */

        /*
          x, y, and z are in 'input' coordinates, while kx, ky, kz are
          in 'reference' coordinates'.
        */
        for (z = 0, kz = zoff; kz < depth + zoff; kz++, z++) {
          if (kz < 0 || kz >= depth) {
            continue;
          }
          for (y = 0, ky = yoff; ky < height + yoff; ky++, y++) {
            if (ky < 0 || ky >= height) {
              continue;
            }
#if 1
            pref = &MRIvox(mri_ref, xoff, ky, kz);
            pin = &MRIvox(mri_in, 0, y, z);
#else
            pref = &mri_ref->slices[kz][ky][xoff];
            pin = mri_in->slices[z][y];
#endif
            kxend = MIN(width + xoff, width);
            for (x = 0, kx = xoff; kx < kxend; kx++, x++, pin++, pref++) {
              if (kx < 0 || kx >= width) {
                continue;
              }

              in = *pin;
              ref = *pref;
              dst += (float)in * (float)ref;
            }
          }
        }
        *pdst++ = dst;
      }
    }
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate two MRIs, limitting the extent of the
          correlation to window_size and normalizing by the energy
          in the window.
------------------------------------------------------*/
MRI *MRInxcorrWindow(MRI *mri_ref, MRI *mri_in, MRI *mri_dst, int window_size)
{
  int width, height, depth, xoff, yoff, zoff, x0, y0, z0, kx, ky, kz, kxend, whalf, x, y, z, dst_width, dst_height,
      dst_depth;
  BUFTYPE *pref, *pin, ref, in;
  float dst, norm;
  // float *pdst;

  if (!mri_dst) {
    /* decide on approriate size for destination */
    dst_width = window_size /* MAX(mri_ref->width, window_size)*/;
    dst_height = window_size /*MAX(mri_ref->height, window_size)*/;
    dst_depth = window_size /*MAX(mri_ref->depth, window_size)*/;

    mri_dst = MRIalloc(dst_width, dst_height, dst_depth, MRI_FLOAT);
    MRIcopyHeader(mri_ref, mri_dst);
  }

  dst_width = mri_dst->width;
  dst_height = mri_dst->height;
  dst_depth = mri_dst->depth;

  height = mri_ref->height;
  width = mri_ref->width;
  depth = mri_ref->depth;

  if (window_size <= 0) {
    window_size = width;
  }
  whalf = (window_size - 1) / 2;

  x0 = (dst_width - 1) / 2;
  y0 = (dst_height - 1) / 2;
  z0 = (dst_depth - 1) / 2;

  /*
     for each element of the destination array, correlate the reference
     MRI at that offset with the (unshifted) input MRI
     */
  for (zoff = -whalf; zoff <= whalf; zoff++) {
    for (yoff = -whalf; yoff <= whalf; yoff++) {
      // pdst = &MRIFvox(mri_dst, x0 - whalf, y0 + yoff, z0 + zoff);
      for (xoff = -whalf; xoff <= whalf; xoff++) {
        /* apply offset template (ref) to input */
        dst = 0.0f;

        /* now calculate value of cross-correlation at these offsets */

        /*
          x, y, and z are in 'input' coordinates, while kx, ky, kz are
          in 'reference' coordinates'.
        */

        norm = 0.0f;
        for (z = 0, kz = zoff; kz < depth + zoff; kz++, z++) {
          if (kz < 0 || kz >= depth) {
            continue;
          }
          for (y = 0, ky = yoff; ky < height + yoff; ky++, y++) {
            if (ky < 0 || ky >= height) {
              continue;
            }

            pref = &MRIvox(mri_ref, xoff, ky, kz);
            pin = &MRIvox(mri_in, 0, y, z);
            kxend = MIN(width + xoff, width);
            for (x = 0, kx = xoff; kx < kxend; kx++, x++, pin++, pref++) {
              if (kx < 0 || kx >= width) {
                continue;
              }

              in = *pin;
              ref = *pref;

              norm += fabs((float)in * (float)ref);
              dst += (float)in * (float)ref;
            }
          }
        }
        norm = 1.0f; /* disable normalization (for now) */
        if (!FZERO(norm)) {
          dst /= norm;
        }
        else {
          dst = 0.0f;
        }

        MRIFvox(mri_dst, x0 + xoff, y0 + yoff, z0 + zoff) = dst;
      }
    }
  }

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate to MRIs
------------------------------------------------------*/
long MRIcorrelate(MRI *mri_ref, MRI *mri_in, int xoff, int yoff, int zoff)
{
  int row, col, slice, width, height, depth;
  long total;
  BUFTYPE *pref, *pin, ref, in;

  depth = mri_ref->depth;
  height = mri_ref->height;
  width = mri_ref->width;

  total = 0L;
  for (slice = 0; slice < depth; slice++) {
    for (row = 0; row < height; row++) {
      pref = mri_ref->slices[slice][row];
      pin = mri_in->slices[slice][row];
      for (col = 0; col < width; col++) {
        in = *pin++;
        ref = *pref++;
        total += (long)in * (long)ref;
      }
    }
  }

  return (total);
}
#if 0
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Remove small inconsistincies in the labelling of a volume.
------------------------------------------------------*/
MRI *
MRIremoveHoles(MRI *mri_src, MRI*mri_dst, int wsize, float pct)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0, thresh,xi,yi,zi,
          num_on, num_off ;
  BUFTYPE val, *pdst, out_val ;
  float   wcubed ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

  wcubed = (float)(wsize*wsize*wsize) ;
  thresh = nint((float)wcubed*pct) ;
  whalf = wsize/2 ;

  for (z = whalf ; z < depth ; z++)
  {
    for (y = whalf ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        for (out_val=num_on = num_off = 0, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          zi = mri_src->zi[z+z0] ;
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            yi = mri_src->yi[y+y0] ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              xi = mri_src->xi[x+x0] ;
              val = MRIvox(mri_src, xi, yi, zi) ;
              if (!val)
              {
                num_off++ ;
              }
              else
              {
                if (!out_val)
                {
                  out_val = 0 ;
                }
                num_on++ ;
              }
            }
          }
        }
        val = MRIvox(mri_src, x, y, z) ;
        if (val && (num_off >= thresh))
        {
          val = 0 ;
        }
        else if (!val && (num_on >= thresh))
        {
          val = 255 ;
        }
        *pdst++ = val ;
      }
    }
  }
  return(mri_dst) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           turn a frame of standard deviations (byte) into
           a frame of variances (also byte, but scaled down
           by a factor of 20 to fit in a byte image (float
           images are just too big).
------------------------------------------------------*/
#define VARIANCE_SCALING 15
MRI *MRIstdsToVariances(MRI *mri_std, MRI *mri_var, int source_frame)
{
  int x, y, z, width, height, depth;
  float std, var;
  BUFTYPE *pstd, *pvar;

  if (mri_var == NULL) {
    mri_var = MRIclone(mri_std, NULL);
  }

  width = mri_var->width;
  height = mri_var->height;
  depth = mri_var->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pstd = &MRIseq_vox(mri_std, 0, y, z, source_frame);
      pvar = &MRIvox(mri_var, 0, y, z);
      for (x = 0; x < width; x++) {
        std = (float)*pstd++;
        var = (std * std) / VARIANCE_SCALING;
        if (var > 255) {
          var = 255; /* can't help it - some clipping will occur */
        }
        *pvar++ = (BUFTYPE)nint(var);
      }
    }
  }
  return (mri_var);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Inverse of above procedure.
------------------------------------------------------*/
MRI *MRIvariancesToStds(MRI *mri_var, MRI *mri_std, int dst_frame)
{
  int x, y, z, width, height, depth;
  float std, var;
  BUFTYPE *pstd, *pvar;

  if (mri_std == NULL) {
    mri_std = MRIclone(mri_var, NULL);
  }

  width = mri_var->width;
  height = mri_var->height;
  depth = mri_var->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pstd = &MRIseq_vox(mri_std, 0, y, z, dst_frame);
      pvar = &MRIvox(mri_var, 0, y, z);
      for (x = 0; x < width; x++) {
        var = (float)*pvar++;
        std = sqrt(var * VARIANCE_SCALING);
        if (std > 255) {
          std = 255; /* can't help it - some clipping will occur */
        }
        *pstd++ = (BUFTYPE)nint(std);
      }
    }
  }
  return (mri_std);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          convolve a special type of MR image which contains
          two frames. The first is a frame of means and
          can be convolved normally. The second is a frame
          of standard deviations, and must be turned into
          variances, convolved, then back to stds and appended
          to the mean image.
------------------------------------------------------*/
MRI *MRIconvolveGaussianMeanAndStdByte(MRI *mri_src, MRI *mri_dst, MRI *mri_gaussian)
{
  MRI *mri_var, *mri_var_blurred, *mri_means_blurred, *mri_std_blurred;
  int nframes;

  if (mri_src->nframes < 2) /* just convolve it normally */
  {
    return (MRIconvolveGaussian(mri_src, mri_dst, mri_gaussian));
  }

  /* this is a hack, but only want to convolve the first frame of mri_src,
     as the 2nd is stds which need to be convolved as variances. Not
     doing this would work but would take twice as much time and memory.
  */
  nframes = mri_src->nframes;
  mri_src->nframes = 1;
  mri_means_blurred = MRIconvolveGaussian(mri_src, NULL, mri_gaussian);
  mri_src->nframes = nframes;
  mri_var = MRIstdsToVariances(mri_src, NULL, 1);
  mri_var_blurred = MRIconvolveGaussian(mri_var, NULL, mri_gaussian);
  MRIfree(&mri_var);
  mri_std_blurred = MRIvariancesToStds(mri_var_blurred, NULL, 0);
  MRIfree(&mri_var_blurred);
  mri_dst = MRIconcatenateFrames(mri_means_blurred, mri_std_blurred, mri_dst);

  MRIfree(&mri_means_blurred);
  MRIfree(&mri_std_blurred);

  MRIcopyHeader(mri_src, mri_dst);

  mri_dst->mean = MRImeanFrame(mri_dst, 1);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON && 0) {
    char fname[100];
    sprintf(fname, "means_and_stds%d.mgh", (int)mri_dst->thick);
    fprintf(stderr, "writing means and stds to %s\n", fname);
    MRIwrite(mri_dst, fname);
  }

  return (mri_dst);
}
/*! -----------------------------------------------------
  \fn MRI *MRImarkBorderVoxels(MRI *mri_src, MRI *mri_dst)
  \param mri_src - binarized volume
  \parm mri_dst - trinerized output 
  \brief Trinerizes mri_src. If an input voxel is MRI_WHITE and one
  of its neighbors is MRI_NOT_WHITE, then the output set to MRI_WHITE.
  If an input voxel is MRI_NOT_WHITE and one of its neighbors is
  MRI_WHITE, then the output set to MRI_NOT_WHITE. All other voxels
  are set to MRI_AMBIGUOUS. Basically, it is marking the boundaries
  between white and not white. In practice (eg, mris_make_surfaces),
  mri_src is often wm.mgz, which could give inaccurate segmentation.
  MRI_AMBIGUOUS (128) - not a border voxel
  MRI_WHITE (255) - white source voxel with at least one non-white neighbor
  MRI_NOT_WHITE (1) - non-white source voxel with at least one white neighbor
------------------------------------------------------*/
MRI *MRImarkBorderVoxels(MRI *mri_src, MRI *mri_dst)
{
  int x, y, z, width, height, depth, xk, yk, zk, xi, yi, zi, slabel, dlabel, nw, ng;
  BUFTYPE *psrc, *pdst;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (mri_src->type == MRI_UCHAR) {
    for (nw = ng = z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        psrc = &MRIvox(mri_src, 0, y, z);
        pdst = &MRIvox(mri_dst, 0, y, z);
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          slabel = *psrc++;

          dlabel = MRI_AMBIGUOUS;
          if (slabel == MRI_WHITE) /* see if it is next to an MRI_NOT_WHITE */
          {
            for (zk = -1; zk <= 1; zk++) {
              zi = mri_src->zi[z + zk];
              for (yk = -1; yk <= 1; yk++) {
                yi = mri_src->yi[y + yk];
                for (xk = -1; xk <= 1; xk++) {
                  xi = mri_src->xi[x + xk];
                  if (MRIvox(mri_src, xi, yi, zi) == MRI_NOT_WHITE) {
                    dlabel = MRI_WHITE; /* mark it as a white border voxel */
                    break;
                  }
                }
                if (dlabel == MRI_WHITE) {
                  break;
                }
              }
              if (dlabel == MRI_WHITE) {
                break;
              }
            }
          }
          else if (slabel == MRI_NOT_WHITE) /*see if it is next
                                              to an MRI_WHITE */
          {
            for (zk = -1; zk <= 1; zk++) {
              zi = mri_src->zi[z + zk];
              for (yk = -1; yk <= 1; yk++) {
                yi = mri_src->yi[y + yk];
                for (xk = -1; xk <= 1; xk++) {
                  xi = mri_src->xi[x + xk];
                  if (MRIgetVoxVal(mri_src, xi, yi, zi, 0) == MRI_WHITE) {
                    dlabel = MRI_NOT_WHITE; /* mark it as a
                                               gray border voxel */
                    break;
                  }
                }
                if (dlabel == MRI_NOT_WHITE) {
                  break;
                }
              }
              if (dlabel == MRI_NOT_WHITE) {
                break;
              }
            }
          }

          if (dlabel == MRI_WHITE) {
            nw++;
          }
          else if (dlabel == MRI_NOT_WHITE) {
            ng++;
          }
          *pdst++ = dlabel; // assign value of voxel to dlabel
        }
      }
    }
  }
  else  // not UCHAR
  {
    for (nw = ng = z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          slabel = nint(MRIgetVoxVal(mri_src, x, y, z, 0));

          dlabel = MRI_AMBIGUOUS;
          if (slabel == MRI_WHITE) /* see if it is next to an MRI_NOT_WHITE */
          {
            for (zk = -1; zk <= 1; zk++) {
              zi = mri_src->zi[z + zk];
              for (yk = -1; yk <= 1; yk++) {
                yi = mri_src->yi[y + yk];
                for (xk = -1; xk <= 1; xk++) {
                  xi = mri_src->xi[x + xk];
                  if ((int)MRIgetVoxVal(mri_src, xi, yi, zi, 0) == MRI_NOT_WHITE) {
                    dlabel = MRI_WHITE; /* mark it as a
                                           white border voxel */
                    break;
                  }
                }
                if (dlabel == MRI_WHITE) {
                  break;
                }
              }
              if (dlabel == MRI_WHITE) {
                break;
              }
            }
          }
          else if (slabel == MRI_NOT_WHITE) /*see if it is next
                                              to an MRI_WHITE */
          {
            for (zk = -1; zk <= 1; zk++) {
              zi = mri_src->zi[z + zk];
              for (yk = -1; yk <= 1; yk++) {
                yi = mri_src->yi[y + yk];
                for (xk = -1; xk <= 1; xk++) {
                  xi = mri_src->xi[x + xk];
                  if ((int)MRIgetVoxVal(mri_src, xi, yi, zi, 0) == MRI_WHITE) {
                    dlabel = MRI_NOT_WHITE; /* mark it as a
                                               gray border voxel */
                    break;
                  }
                }
                if (dlabel == MRI_NOT_WHITE) {
                  break;
                }
              }
              if (dlabel == MRI_NOT_WHITE) {
                break;
              }
            }
          }

          if (dlabel == MRI_WHITE) {
            nw++;
          }
          else if (dlabel == MRI_NOT_WHITE) {
            ng++;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, dlabel);
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr, "border white:  %8d voxels (%2.2f%%)\n", nw, 100.0f * (float)nw / (float)(width * height * depth));
    fprintf(stderr, "border gray    %8d voxels (%2.2f%%)\n", ng, 100.0f * (float)ng / (float)(width * height * depth));
  }
  return (mri_dst);
}

/*!
  \fn int MRIborderClassifyVoxel(MRI *mri_src, MRI *mri_labeled, int wsize, int x, int y, int z, float *ppw, float *ppg)
  \brief Classifies voxel as either MRI_NOT_WHITE=1 or MRI_WHITE=255. Works
  by computing mean and stddev of each class (defined by mri_labeled)
  based on intensities in a wsize neighborhood around the given
  voxel. The probability of each class is computed from means and
  stddevs and the most likely is chosen. wsize is typically 11, meaning
  the neighborhood is +/-5 voxels which works out to 1331 voxels, but
  it might be much less than this if not all voxels in mri_labeled are
  actually labeled (eg, MRIreclassify() uses only border voxels)
 */
int MRIborderClassifyVoxel(MRI *mri_src, MRI *mri_labeled, int wsize, int x, int y, int z, float *ppw, float *ppg)
{
  int xk, yk, zk, xi, yi, zi, nw, ng, whalf, val;
  // int width, height, depth;
  float wvar, wmean, gvar, gmean, wsq, gsq, pg, pw, dist, ptotal;

  // width = mri_src->width;
  // height = mri_src->height;
  // depth = mri_src->depth;
  whalf = (wsize - 1) / 2;
  wvar = wmean = gvar = gmean = wsq = gsq = 0.0;
  nw = ng = 0;
  for (zk = -whalf; zk <= whalf; zk++) {
    zi = mri_src->zi[z + zk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri_src->yi[y + yk];
      for (xk = -whalf; xk <= whalf; xk++) {
        if (zk==0 && yk==0 && xk==0)
          continue; /* don't use central value as we are classifying it */
        xi = mri_src->xi[x + xk];
        val = MRIgetVoxVal(mri_src, xi, yi, zi, 0);
        if (MRIgetVoxVal(mri_labeled, xi, yi, zi, 0) == MRI_NOT_WHITE) {
          ng++;
          gsq += val * val;
          gmean += val;
        }
        else if (MRIgetVoxVal(mri_labeled, xi, yi, zi, 0) == MRI_WHITE) {
          nw++;
          wsq += val * val;
          wmean += val;
        }
      }
    }
  }
  if (nw) {
    wmean /= (float)nw;
    wvar = wsq / (float)nw - wmean * wmean;
  }
  else {
    wmean = wvar = 1.0;
  }

  if (ng) {
    gmean /= (float)ng;
    gvar = gsq / (float)ng - gmean * gmean;
  }
  else {
    gmean = gvar = 1.0;
  }

  val = MRIgetVoxVal(mri_src, x, y, z, 0);
  dist = (float)val - wmean;
  dist *= dist;
  pw = exp(-dist / (2 * wvar));

  val = MRIgetVoxVal(mri_src, x, y, z, 0);
  dist = (float)val - gmean;
  dist *= dist;
  pg = exp(-dist / (2 * wvar));

  ptotal = pg + pw;

  if (!FZERO(ptotal)){ 
    pg /= ptotal;
  }
  pw /= ptotal; // is this a bug?
  if(ppg) *ppg = pg;
  if(ppw) *ppw = pw;

  return (pg > pw ? MRI_NOT_WHITE : MRI_WHITE);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIreclassifyBorder(MRI *mri_src, MRI *mri_labeled, MRI *mri_border, MRI *mri_dst, int wsize)
{
  int x, y, z, width, height, depth, nw, ng, nchanged, ntested;
  float label, new_label, border;

  if (!mri_dst) {
    mri_dst = MRIcopy(mri_labeled, NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (ntested = nchanged = nw = ng = z = 0; z < depth; z++) {
    DiagShowPctDone((float)(z) / (float)(depth - 1), 5);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        label = MRIgetVoxVal(mri_labeled, x, y, z, 0);
        border = MRIgetVoxVal(mri_border, x, y, z, 0);
        if (border != MRI_AMBIGUOUS) /* a border voxel */
        {
          ntested++;
          new_label = MRIborderClassifyVoxel(mri_src, mri_border, wsize, x, y, z, NULL, NULL);
          if (new_label != label) {
            nchanged++;
          }
          label = new_label;
        }
        MRIsetVoxVal(mri_labeled, x, y, z, 0, label);
      }
    }
  }

  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr,
            "               %8d voxels tested (%2.2f%%)\n",
            ntested,
            100.0f * (float)ntested / (float)(width * height * depth));
    fprintf(stderr,
            "               %8d voxels changed (%2.2f%%)\n",
            nchanged,
            100.0f * (float)nchanged / (float)(width * height * depth));
  }

  return (NO_ERROR);
}
/*!
  \fn int MRIreclassify(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, float wm_low, float gray_hi, int wsize)

  \brief Reclassifies voxels where mri_src intensity is between wm_low
  and gray_hi using MRIborderClassifyVoxel(), which looks in a
  neighborhood around the voxel to compute means and stddevs of the
  two classes (WM, GM), then computes the probability that the given
  voxel came from each class, then classifies based on which is
  higher. If the prob is less than 0.7, then the size of the
  neighborhood is changed from 5 to 21 (step 4). The probs for WM and
  GM are computed at each scale, and the class with the highest sum of
  probs is chosen. Every voxel is considered.
*/
int MRIreclassify(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, float wm_low, float gray_hi, int wsize)
{
  int x, y, z, width, height, depth, nw, ng, nchanged, ntested, w, nmulti;
  int label, new_label, src;
  float pw, pg, gtotal, wtotal;
  MRI *mri_border;

  printf("MRIreclassify(): wm_low=%g, gray_hi=%g, wsize=%d\n",wm_low,gray_hi,wsize);
  fflush(stdout);

  // Get WM and non-WM labels along the border. mri_label
  // is not used again for reclassification, so it is ok 
  // that its contents are changed along the way.
  mri_border = MRImarkBorderVoxels(mri_labeled, NULL);

  if (!mri_dst)
    mri_dst = MRIcopy(mri_labeled, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  nmulti = ntested = nchanged = nw = ng = 0;
  for (z = 0; z < depth; z++) {
    DiagShowPctDone((float)(z) / (float)(depth - 1), 5);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        src = (int)MRIgetVoxVal(mri_src, x, y, z, 0);
        label = (int)MRIgetVoxVal(mri_labeled, x, y, z, 0);
	// Only reclassify voxels between wm_low and gray_hi
        if (src > wm_low && src < gray_hi) {
          ntested++;

	  // Get a new label based on stats of voxels near the border
          new_label = MRIborderClassifyVoxel(mri_src, mri_border, wsize, x, y, z, &pw, &pg);

	  // If the likelihood of the chosen label is sufficiently
	  // high, then just take it. If not, then do a multi-scale
	  // search. Where the neighborhood size changes 5:4:21. At
	  // each neighborhood, the prob of WM and GM are computed.
	  // The probabilities are summed over the sizes, and the
	  // class that has the highest sum is chosen.
          if (MAX(pw, pg) < .7) {
            nmulti++;
            gtotal = wtotal = 0.0;
            for (w = 5; w < 21; w += 4) {
              MRIborderClassifyVoxel(mri_src, mri_border, w, x, y, z, &pw, &pg);
              gtotal += pg;
              wtotal += pw;
            }
            if(gtotal > wtotal) 
              new_label = MRI_NOT_WHITE;
            else
              new_label = MRI_WHITE;
          }
          if(new_label != label)
            nchanged++;
          label = new_label;
        }
        MRIsetVoxVal(mri_dst, x, y, z, 0, label);
      }
    }
  }

  printf("  %8d voxels tested (%2.2f%%)\n",
            ntested, 100.0f * (float)ntested / (float)(width * height * depth));
  printf("  %8d voxels changed (%2.2f%%)\n",
            nchanged,100.0f * (float)nchanged / (float)(width * height * depth));
  printf("  %8d multi-scale searches  (%2.2f%%)\n",
            nmulti, 100.0f * (float)nmulti / (float)(width * height * depth));
  fflush(stdout);

  MRIfree(&mri_border);
  return (NO_ERROR);
}


/*!
  \fn MRI *MRIreclassifyWMCtxNonPar(MRI *norm, MRI *seg, int nitersmax, MRI *newseg)
  \brief Reclassifies WM and Cortex based on non-parametric
  probability intensity distributions.  The passed seg (eg,
  aseg.presurf.mgz) is used to determine the non-parametric PDFs for
  WM and cortex and an optimal threshold is chosen. Voxels that are on
  the border between WM and cortex are then examined and possibly
  recoded based on intensity. This process is iterated until either
  nitersmax is reached or no more voxels are recoded. This funciton
  itself could be iterated because the PDFs are computed from the
  segmentation; however, this iteration does not stablize (in one case
  it caused WM to take over all GM).
 */
MRI *MRIreclassifyWMCtxNonPar(MRI *norm, MRI *seg, int nitersmax, MRI *newseg)
{
  int wmsegs[] = {2,41,251,252,253,254,255},nwmsegs;
  int ctxsegs[] = {3,42},nctxsegs;
  int nwm,nctx;
  MRI *wmmask, *ctxmask;
  float wmmin,  wmmax,  wmrange,  wmmean,  wmstd;
  float ctxmin, ctxmax, ctxrange, ctxmean, ctxstd;
  int wmm, ctxm, nbrwmm, nbrctxm, isborder, c,r,s,dc, dr, ds,nrecode,nborder,nhits,segid;
  HISTOGRAM *hwm, *hctx;
  double val,min,thresh,hbinmin,hbinmax,d;
  int n,nmin=0,k;

  printf("MRIreclassifyWMCtxNonPar(): nitersmax = %d \n",nitersmax);fflush(stdout);
  if(seg==newseg){
    printf("ERROR: cannot be done in-place\n");
    return(NULL);
  }

  /* Create binary masks of WM and Cortex (note corpus callosum will not be
     in aseg.presurf.mgz, so 251-255 will not be relevant*/
  nwmsegs  = sizeof(wmsegs)/sizeof(int);
  wmmask = MRImatchSegs(seg,wmsegs, nwmsegs, 1, NULL);//1=output value
  nctxsegs = sizeof(ctxsegs)/sizeof(int);
  ctxmask = MRImatchSegs(seg,ctxsegs, nctxsegs, 1, NULL);//1=output value

  // Compute stats from WM
  nwm = MRIsegStats(wmmask, 1, norm, 0, &wmmin, &wmmax, &wmrange,&wmmean,&wmstd);
  if(nwm == 0){
    printf("ERROR: could not find any WM voxels\n");
    return(NULL);
  }
  printf("wmstats  n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nwm,wmmean,wmstd,wmmin,wmmax);

  // Compute stats from Cortex
  nctx = MRIsegStats(ctxmask, 1, norm, 0, &ctxmin, &ctxmax, &ctxrange,&ctxmean,&ctxstd);
  if(nctx == 0){
    printf("ERROR: could not find any CTX voxels\n");
    return(NULL);
  }
  printf("ctxstats n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nctx,ctxmean,ctxstd,ctxmin,ctxmax);
  fflush(stdout);

  // Get intensity limits for histogram
  hbinmin = MIN(wmmin,ctxmin);
  hbinmax = MIN(wmmax,ctxmax);
  printf("hbminmin = %g, hbinmax = %g\n",hbinmin,hbinmax);

  // Create histogram (PDF) of WM voxels, then CDF
  hwm = HISTOseg(wmmask, 1, norm, hbinmin,hbinmax,1);
  HISTOsumNorm(hwm);
  HISTOcumsum(hwm,hwm);

  // Create histogram (PDF) of Cortex voxels, then CDF
  hctx = HISTOseg(ctxmask, 1, norm, hbinmin,hbinmax,1);
  HISTOsumNorm(hctx);
  HISTOcumsum(hctx,hctx);

  //HISTOwriteTxt(hwm, "wmhist.dat");
  //HISTOwriteTxt(hctx, "ctxhist.dat");

  /* Find threshold where (1-CtxCDF) = WMCDF. This is the place where
     the probability that a voxel is WM equals that of Cortex.  This
     is a non-parametric way of computing it (ie, instead of computing
     means, stddevs, and z-scores). This assumes T1-weighting. */
  min = 10e10;
  for(n=0; n < hwm->nbins; n++){
    d = (1.0-hctx->counts[n]) - hwm->counts[n];
    if(fabs(d)<min){
      min = d;
      nmin = n;
    }
  }
  thresh = hwm->bins[nmin];
  printf("thresh = %g\n",thresh);
  HISTOfree(&hwm);
  HISTOfree(&hctx);

  /* Iteratively recode voxels based on the threshold. Only voxels
     that are face border voxels are eligible to be recoded in
     a given iteration. */
  for(k=0; k < nitersmax; k++){
    // Copy seg into newseg to keep changes separate. 
    // Note: seg = newseg from previous iteration
    newseg = MRIcopy(seg,newseg); 

    // Create new wm and cortex masks
    wmmask = MRImatchSegs(seg,wmsegs, nwmsegs, 1, wmmask);//1=output value
    ctxmask = MRImatchSegs(seg,ctxsegs, nctxsegs, 1, ctxmask);//1=output value
    
    nhits = 0;   // total number of wm or ctx voxels
    nborder = 0; // total number of face border voxels
    nrecode=0;   // number of voxels recoded on this iteration
    for(c=0; c < wmmask->width; c++){
      for(r=0; r < wmmask->height; r++){
	for(s=0; s < wmmask->depth; s++){
	  wmm  = MRIgetVoxVal(wmmask,c,r,s,0);
	  ctxm = MRIgetVoxVal(ctxmask,c,r,s,0);
	  // If this voxel is neither WM or Cortex, then skip
	  if(!wmm && !ctxm) continue;
	  nhits ++;
	  // Check whether it is a border voxel by determining
	  // whether any face neighbors have a different tissue class 
	  isborder = 0;
	  for(dc=-1; dc <= +1; dc++){
	    if(isborder) break;
	    for(dr=-1; dr <= +1; dr++){
	      if(isborder) break;
	      for(ds=-1; ds <= +1; ds++){
		if(isborder) break;
		if(abs(dc)+abs(dr)+abs(ds) != 1) continue; // face
		nbrwmm  = MRIgetVoxVal(wmmask,c+dc,r+dr,s+ds,0);
		nbrctxm = MRIgetVoxVal(ctxmask,c+dc,r+dr,s+ds,0);
		if((wmm && nbrctxm) || (ctxm && nbrwmm))
		  isborder = 1;
	      }
	    }
	  }
	  // If it is not a border voxel, skip it
	  if(!isborder) continue;
	  nborder++;

	  // Classify depending upon whether the value is above or
	  // below threshold. Note: corpus callosum will not be in
	  // aseg.presurf.mgz, so 251-255 will not be relevant
	  val = MRIgetVoxVal(norm,c,r,s,0);	
	  segid = MRIgetVoxVal(seg,c,r,s,0); // current seg
	  if(val > thresh){
	    // Recode cortex to WM
	    if(segid == Left_Cerebral_Cortex){
	      MRIsetVoxVal(newseg,c,r,s,0,Left_Cerebral_White_Matter);
	      nrecode++;
	    }
	    if(segid == Right_Cerebral_Cortex){
	      MRIsetVoxVal(newseg,c,r,s,0,Right_Cerebral_White_Matter);
	      nrecode++;
	    }
	  }
	  else {
	    // Recode WM to cortex
	    if(segid == Left_Cerebral_White_Matter){
	      MRIsetVoxVal(newseg,c,r,s,0,Left_Cerebral_Cortex);
	      nrecode++;
	    }
	    if(segid == Right_Cerebral_White_Matter){
	      MRIsetVoxVal(newseg,c,r,s,0,Right_Cerebral_Cortex);
	      nrecode++;
	    }
	  }
	}
      }
    }
    printf("iter %2d/%2d ",k,nitersmax);
    printf("nborder = %d, nrecode = %d\n",nborder,nrecode);
    fflush(stdout);
    if(nrecode==0) break;
    seg = newseg;
  }

  MRIfree(&wmmask);
  MRIfree(&ctxmask);

  return(newseg);
}

/*!
  \fn MRI *MRIclassifyAmbiguousNonPar(MRI *norm, MRI *seg, MRI *statseg, double Qwm, double Qctx, int NdilWM, int NdilCtx)
  \brief Relabels all voxels into MRI_AMBIGUOUS=128, MRI_WHITE=255, or
  MRI_NOT_WHITE=1. Candidate voxels must be within NdilWM voxels of WM
  and NdilCtx voxels of cortex as indicated by the seg
  volume. Ambiguous voxels are defined as being between low and hi
  thresholds. The low threshold is determined by the value at which
  Prob(Ctx) > Qctx*Prob(WM). The hi threshold is determined by the
  value at which Prob(WM) > Qwm*Prob(WM). Prob() is determined by
  computing the PDF of the histogram of WM and Cortex based on the
  statseg. This second seg is a bit of a hack because seg usually
  comes from MRIreclassifyWMCtxNonPar(), which has already used
  intensity to segment. Example parameters Qwm=Qctx=3.0, NdilWM=NdilCtx=2.
 */
MRI *MRIclassifyAmbiguousNonPar(MRI *norm, MRI *seg, MRI *statseg, double Qwm, double Qctx, int NdilWM, int NdilCtx)
{
  int wmsegs[] = {2,41,251,252,253,254,255},nwmsegs;
  int ctxsegs[] = {3,42},nctxsegs;
  int nwm,nctx;
  MRI *wmmask, *ctxmask, *wmmaskdil, *ctxmaskdil, *ambi;
  float wmmin,  wmmax,  wmrange,  wmmean,  wmstd;
  float ctxmin, ctxmax, ctxrange, ctxmean, ctxstd;
  int wmm, ctxm, c,r,s,nrecode,nhits;
  HISTOGRAM *hwm, *hctx, *wmcdf, *ctxcdf;
  double val,min,thresh,hbinmin,hbinmax,d,thresh_low,thresh_hi,q;
  int n,nmin=0;

  printf("MRIclassifyAmbiguousNonPar(): Qwm=%g, Qctx=%g, NdilWM=%d, NdilCtx=%d \n",Qwm,Qctx,NdilWM,NdilCtx);

  /* Create binary masks of WM and Cortex (note corpus callosum will not be
     in aseg.presurf.mgz, so 251-255 will not be relevant*/
  nwmsegs  = sizeof(wmsegs)/sizeof(int);
  nctxsegs = sizeof(ctxsegs)/sizeof(int);
  wmmask = MRImatchSegs(statseg,wmsegs, nwmsegs, 1, NULL);//1=output value
  ctxmask = MRImatchSegs(statseg,ctxsegs, nctxsegs, 1, NULL);//1=output value

  // Compute stats from WM
  nwm = MRIsegStats(wmmask, 1, norm, 0, &wmmin, &wmmax, &wmrange,&wmmean,&wmstd);
  if(nwm == 0){
    printf("ERROR: could not find any WM voxels\n");
    return(NULL);
  }
  printf("wmstats  n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nwm,wmmean,wmstd,wmmin,wmmax);

  // Compute stats from Cortex
  nctx = MRIsegStats(ctxmask, 1, norm, 0, &ctxmin, &ctxmax, &ctxrange,&ctxmean,&ctxstd);
  if(nctx == 0){
    printf("ERROR: could not find any CTX voxels\n");
    return(NULL);
  }
  printf("ctxstats n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nctx,ctxmean,ctxstd,ctxmin,ctxmax);
  fflush(stdout);

  // Get intensity limits for histogram
  hbinmin = MIN(wmmin,ctxmin);
  hbinmax = MIN(wmmax,ctxmax);
  printf("hbminmin = %g, hbinmax = %g\n",hbinmin,hbinmax);

  // Create histogram (PDF) of WM voxels, then CDF
  hwm = HISTOseg(wmmask, 1, norm, hbinmin,hbinmax,1);
  HISTOsumNorm(hwm);
  wmcdf = HISTOcumsum(hwm,NULL);

  // Create histogram (PDF) of Cortex voxels, then CDF
  hctx = HISTOseg(ctxmask, 1, norm, hbinmin,hbinmax,1);
  HISTOsumNorm(hctx);
  ctxcdf = HISTOcumsum(hctx,NULL);

  HISTOwriteTxt(hwm, "wm.pdf.dat");
  HISTOwriteTxt(hctx,"ctx.pdf.dat");
  //HISTOwriteTxt(wmcdf, "wm.cdf.dat");
  //HISTOwriteTxt(ctxcdf,"ctx.cdf.dat");

  /* Find threshold where (1-CtxCDF) = WMCDF. This is the place where
     the probability that a voxel is WM equals that of Cortex.  This
     is a non-parametric way of computing it (ie, instead of computing
     means, stddevs, and z-scores). This assumes T1-weighting. */
  min = 10e10;
  for(n=0; n < wmcdf->nbins; n++){
    d = (1.0-ctxcdf->counts[n]) - wmcdf->counts[n];
    if(fabs(d)<min){
      min = d;
      nmin = n;
    }
  }
  thresh = hwm->bins[nmin];
  // From this minimum, look back toward 0 to find where the prob of Ctx >= Qctx*prob of WM
  for(n=nmin; n>=0;  n--){
    q = hctx->counts[n]/(hwm->counts[n]+1e-6);
    if(q >= Qctx) break;
  }
  thresh_low = hctx->bins[n];
  // From this minimum, look ahead toward max to find where the prob of Cwm >= Qctx*prob of Ctx
  for(n=0; n < hwm->nbins; n++){
    q = hwm->counts[n]/(hctx->counts[n]+1e-6);
    if(q >= Qwm) break;
  }
  thresh_hi = hwm->bins[n];
  printf("thresh = (%g,%g,%g)\n",thresh,thresh_low,thresh_hi);

  HISTOfree(&hwm);
  HISTOfree(&wmcdf);
  HISTOfree(&hctx);
  HISTOfree(&ctxcdf);

  wmmask  = MRImatchSegs(seg,wmsegs,  nwmsegs,  1, wmmask);//1=output value
  ctxmask = MRImatchSegs(seg,ctxsegs, nctxsegs, 1, ctxmask);//1=output value
  wmmaskdil = NULL;
  for(n=0; n < NdilWM; n++)  wmmaskdil = MRIdilate(wmmask, wmmaskdil);
  ctxmaskdil = NULL;
  for(n=0; n < NdilCtx; n++)  ctxmaskdil = MRIdilate(ctxmask, ctxmaskdil);

  ambi = MRIalloc(seg->width,seg->height,seg->depth,MRI_UCHAR);
  MRIcopyHeader(seg,ambi);

  nhits = 0;   // total number of wm or ctx voxels
  nrecode=0;   // number of voxels recoded on this iteration
  for(c=0; c < wmmask->width; c++){
    for(r=0; r < wmmask->height; r++){
      for(s=0; s < wmmask->depth; s++){
	wmm  = MRIgetVoxVal(wmmask,c,r,s,0);
	if(wmm) MRIsetVoxVal(ambi,c,r,s,0,MRI_WHITE);
	else 	MRIsetVoxVal(ambi,c,r,s,0,MRI_NOT_WHITE);
	wmm  = MRIgetVoxVal(wmmaskdil,c,r,s,0);
	ctxm = MRIgetVoxVal(ctxmaskdil,c,r,s,0);
	// If this voxel is neither WM or Cortex in dilation, then skip
	if(!wmm && !ctxm) continue;
	nhits ++;
	val = MRIgetVoxVal(norm,c,r,s,0);	
	if(val > thresh_low && val < thresh_hi){
	  MRIsetVoxVal(ambi,c,r,s,0,MRI_AMBIGUOUS);
	  nrecode++;
	}
      }
    }
  }
  printf("nhits = %d, nrecode = %d\n",nhits,nrecode);
  fflush(stdout);

  MRIfree(&wmmask);
  MRIfree(&ctxmask);
  MRIfree(&wmmaskdil);
  MRIfree(&ctxmaskdil);

  return(ambi);
}



/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRIclassifyAmbiguous(MRI *mri_src, MRI *mri_labeled, MRI *mri_border, MRI *mri_dst, int wsize)
{
  int x, y, z, width, height, depth, nw, ng, nchanged, ntested;
  BUFTYPE *plabeled, *pdst, label, new_label;
  // BUFTYPE *psrc;

  if (!mri_dst) {
    mri_dst = MRIcopy(mri_labeled, NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (ntested = nchanged = nw = ng = z = 0; z < depth; z++) {
    DiagShowPctDone((float)(z) / (float)(depth - 1), 5);
    for (y = 0; y < height; y++) {
      // psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      plabeled = &MRIvox(mri_labeled, 0, y, z);
      for (x = 0; x < width; x++) {
        label = *plabeled++;
        if (label == MRI_AMBIGUOUS) /* classify it */
        {
          ntested++;
          new_label = MRIborderClassifyVoxel(mri_src, mri_border, wsize, x, y, z, NULL, NULL);
          if (new_label != label) {
            nchanged++;
          }
          label = new_label;
        }
        *pdst++ = label;
      }
    }
  }

  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr,
            "               %8d voxels tested (%2.2f%%)\n",
            ntested,
            100.0f * (float)ntested / (float)(width * height * depth));
    fprintf(stderr,
            "               %8d voxels changed (%2.2f%%)\n",
            nchanged,
            100.0f * (float)nchanged / (float)(width * height * depth));
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *MRIremoveBrightStuff(MRI *mri_src, MRI *mri_dst, int threshold)
{
  int x, y, z, width, height, depth;
  BUFTYPE *psrc, *pdst, val;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      pdst = &MRIvox(mri_dst, 0, y, z);
      for (x = 0; x < width; x++) {
        val = *psrc++;
        if (val > threshold) {
          val = 0;
        }
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
------------------------------------------------------*/
MRI *MRIwindow(MRI *mri_src, MRI *mri_dst, int which, float x0, float y0, float z0, float parm)
{
  int x, y, z, height, depth, width, val;
  double w, dist, dx, dy, dz;

  mri_dst = MRIcopy(mri_src, mri_dst);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        dx = (double)x - x0;
        dy = (double)y - y0;
        dz = (double)z - z0;
        dist = sqrt(dx * dx + dy * dy + dz * dz);
        switch (which) {
          default:
          case WINDOW_GAUSSIAN:
            ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRIwindow: unsupported window type %d", which));
            break;
          case WINDOW_HANNING:
#if 0
          w = .5*(1 - cos(2*M_PI*dist/parm)); /* symmetric */
#else
            w = 0.5 * cos(M_PI * dist / parm) + 0.5; /* one-sided */
#endif
            if (dist > parm) {
              w = 0.0;
            }
            break;
          case WINDOW_HAMMING:
#if 0
          w = .54 - .46*cos(2*M_PI*dist/parm);  /* symmetric */
#else
            w = .54 + .46 * cos(M_PI * dist / parm); /* one-sided */
#endif
            if (dist > parm) {
              w = 0.0;
            }
            break;
        }
        val = nint(w * MRIvox(mri_dst, x, y, z));
        MRIvox(mri_dst, x, y, z) = val;
      }
    }
  }
  return (mri_dst);
}
/*!
  \fn int MRIcomputeClassStatistics()
  \brief Computes the mean and stddev of WM and NOT_WM (ie, GM) classes
  \param mri_T1 - volume of T1 intensities in the range of 0-255
  \param mri_labeled - binary volume of WM voxels (wm.mgz)
  \param gray_low (eg,  30) - not necessarily for gray (see mris_make_surfaces)
  \param gray_hi  (eg, 110) - not necessarily for gray (see mris_make_surfaces)

  The mri_labeled is first reduced to boundary voxels that are either WM or NOT WM. 
  The segmentation is further refined with WM voxels needing to have an intensity
  greather than (gray_low+gray_hi)/2 (~70) and NOT WM voxels needing to be
  between gray_low (30) and gray_high (110). Seems like it would be better to
  do this with the aseg.

  Setting gray_low <0 disables the intensity range checking.  This is
  useful for images with contrast that is inverted with respect to
  normal T1-weighted in vivo scans.
  ------------------------------------------------------*/
int MRIcomputeClassStatistics(MRI *mri_T1,
                              MRI *mri_labeled,
                              float gray_low,
                              float gray_hi,
                              float *pmean_wm,
                              float *psigma_wm,
                              float *pmean_gm,
                              float *psigma_gm)     // BEVIN mris_make_surfaces
{
  MRI *mri_border;
  int x, y, z, width, height, depth, label, border_label, ngray, nwhite, nbins, bin, peak;
  double white_min, white_max, white_mean, white_std, gray_min, gray_max, gray_mean, gray_std, val;
  HISTOGRAM *h_white, *h_gray;
  float min_val, max_val, white_mode, gray_mode;

  MRIvalRange(mri_T1, &min_val, &max_val);
  min_val = 0;
  nbins = ceil(max_val - min_val) + 1;
  h_white = HISTOalloc(nbins);
  h_gray = HISTOalloc(nbins);
  for (bin = 0; bin < nbins; bin++) {
    h_white->bins[bin] = min_val + bin;
    h_gray->bins[bin] = min_val + bin;
  }
  HISTOinit(h_gray, h_gray->nbins, 0, 0);
  HISTOinit(h_white, h_white->nbins, 0, 0);

  // This labels two classes (WM and NOT_WM) along the border of WM.
  // All other voxels are labeled as ambiguous.
  mri_border = MRImarkBorderVoxels(mri_labeled, NULL);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    MRIwrite(mri_border, "border.mgz");
  }

  width = mri_T1->width;
  height = mri_T1->height;
  depth = mri_T1->depth;

  white_mean = gray_mean = white_std = gray_std = 0.0;
  white_min = gray_min = 100000;
  white_max = gray_max = -white_min;
  nwhite = ngray = 0;
  
  // Measurement showed this not worth improving
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(serial)
#endif
  for (z = 0; z < depth; z++) {
    ROMP_PFLB_begin
    
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak();

	// Ignore any voxels in the input marked as ambiguous
        label = MRIgetVoxVal(mri_labeled, x, y, z, 0);
        if (label == MRI_AMBIGUOUS)
          continue;

        val = MRIgetVoxVal(mri_T1, x, y, z, 0);
        border_label = MRIgetVoxVal(mri_border, x, y, z, 0);

        if (border_label == MRI_WHITE) {
          // This is a voxel labeled as WM in the input (and border)
          if ((val >= (gray_low + gray_hi) / 2) || (gray_low < 0)) {
            // The intensity is > the passed mid gray intensity
            // Record stats
            nwhite++;
            white_mean += (double)val;
            white_std += (double)(val * val);
            if (val < white_min) 
              white_min = (double)val;
            if (val > white_max) 
              white_max = val;
            bin = nint(val - min_val);
            if (bin < 0 || bin >= h_white->nbins)
              DiagBreak();
            h_white->counts[bin]++;
          }
        }
        else if (border_label == MRI_NOT_WHITE) 
        {
          /* Labeled as non-white (gray) in input (and bordering) */
          if ((val >= gray_low && val <= gray_hi) || gray_low < 0) {
            // value is between gray low and hi
            ngray++;
            gray_mean += (double)val;
            gray_std += (double)(val * val);
            if (val < gray_min) {
              gray_min = (double)val;
            }
            if (val > gray_max) {
              gray_max = val;
            }
            bin = nint(val - min_val);
            if (bin < 0 || bin >= h_gray->nbins) {
              DiagBreak();
            }
            h_gray->counts[bin]++;
          }
        }
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)  // debugging
  {
    double mn, st;
    HISTOrobustGaussianFit(h_gray, .1, &mn, &st);
    DiagBreak();
  }
  peak = HISTOfindHighestPeakInRegion(h_white, 0, h_white->nbins);
  white_mode = h_white->bins[peak];
  peak = HISTOfindHighestPeakInRegion(h_gray, 0, h_gray->nbins);
  gray_mode = h_gray->bins[peak];

  // use median instead of mode. More robust, but ignored below
  white_mode = HISTOfindMedian(h_white);
  gray_mode = HISTOfindMedian(h_gray);
  gray_mean /= (double)ngray;
  white_mean /= (double)nwhite;
  white_std = sqrt(white_std / (double)nwhite - white_mean * white_mean);
  gray_std = sqrt(gray_std / (double)ngray - gray_mean * gray_mean);
  printf("CCS WM (%2.1f): %2.1f +- %2.1f [%2.1f --> %2.1f]\n",
          white_mode, white_mean, white_std, white_min, white_max);
  printf("CCS GM (%2.1f) : %2.1f +- %2.1f [%2.1f --> %2.1f]\n", 
    gray_mode, gray_mean, gray_std, gray_min, gray_max);
  MRIfree(&mri_border);
  *pmean_wm = white_mode;
  *pmean_gm = gray_mode;

  *pmean_wm = white_mean;  // use means for now - works better on tutorial data
  *pmean_gm = gray_mean;
  *psigma_wm = white_std;
  *psigma_gm = gray_std;
  return (NO_ERROR);
}

/*!
  \fn int MRIcomputeClassStatisticsSeg(MRI *norm, MRI *seg, float *wmmean, float *wmstd, float *ctxmean, float *ctxstd)
  \brief Computes WM and cortex mean and stddev based on the segementation.
 */
int MRIcomputeClassStatisticsSeg(MRI *norm, MRI *seg, float *wmmean, float *wmstd, float *ctxmean, float *ctxstd)
{
  int wmsegs[] = {2,41,251,252,253,254,255},nwmsegs;
  int ctxsegs[] = {3,42},nctxsegs;
  int nwm,nctx;
  float wmmin,  wmmax,  wmrange;
  float ctxmin, ctxmax, ctxrange;
  MRI *wmmask, *ctxmask;

  /* Create binary masks of WM and Cortex */
  nwmsegs  = sizeof(wmsegs)/sizeof(int);
  nctxsegs = sizeof(ctxsegs)/sizeof(int);
  wmmask = MRImatchSegs(seg,wmsegs, nwmsegs, 1, NULL);//1=output value
  ctxmask = MRImatchSegs(seg,ctxsegs, nctxsegs, 1, NULL);//1=output value

  // Compute stats from WM
  nwm = MRIsegStats(wmmask, 1, norm, 0, &wmmin, &wmmax, &wmrange,wmmean,wmstd);
  if(nwm == 0){
    printf("ERROR: could not find any WM voxels\n");
    return(1);
  }
  printf("wmstats  n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nwm,*wmmean,*wmstd,wmmin,wmmax);

  // Compute stats from Cortex
  nctx = MRIsegStats(ctxmask, 1, norm, 0, &ctxmin, &ctxmax, &ctxrange,ctxmean,ctxstd);
  if(nctx == 0){
    printf("ERROR: could not find any CTX voxels\n");
    return(1);
  }
  printf("ctxstats n=%6d, mean=%g, std=%g, min=%g, max=%g\n",nctx,*ctxmean,*ctxstd,ctxmin,ctxmax);
  MRIfree(&wmmask);
  MRIfree(&ctxmask);
  fflush(stdout);
  return(0);
}

MRI *MRIthreshModeFilter(MRI *mri_src, MRI *mri_dst, int niter, float thresh)
{
  int x, y, z, n, width, height, depth, histo[256], xk, yk, zk, xi, yi, zi, val, i, max_histo, max_i, changed,
      orig_label;

  mri_src = MRIcopy(mri_src, NULL); /* make a copy we can modify */

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }
  MRIcopy(mri_src, mri_dst);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (n = 0; n < niter; n++) {
    changed = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          memset(histo, 0, sizeof(histo));
          orig_label = MRIvox(mri_src, x, y, z);
          for (zk = -1; zk <= 1; zk++) {
            zi = mri_src->zi[z + zk];
            for (yk = -1; yk <= 1; yk++) {
              yi = mri_src->yi[y + yk];
              for (xk = -1; xk <= 1; xk++) {
                xi = mri_src->xi[x + xk];
                val = MRIvox(mri_src, xi, yi, zi);
                histo[val]++;
              }
            }
          }
          for (max_histo = max_i = i = 0; i < 256; i++) {
            if (histo[i] > max_histo) {
              max_histo = histo[i];
              max_i = i;
            }
          }
          if (max_i == orig_label) {
            continue; /* current label is already mode of histogram */
          }

          if ((((float)max_histo / (3 * 3 * 3)) >= thresh) || ((float)histo[orig_label] / (3 * 3 * 3) < (1 - thresh))) {
            changed++;
            MRIvox(mri_dst, x, y, z) = max_i;
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr,
              "mode filter pass %d:  %8d voxels changed (%2.2f%%)\n",
              n,
              changed,
              100.0f * (float)changed / (float)(width * height * depth));
    if (n < niter - 1) {
      MRIcopy(mri_dst, mri_src);
    }
  }

  MRIfree(&mri_src);

  return (mri_dst);
}

MRI *MRImodeFilterWithControlPoints(MRI *mri_src, MRI *mri_ctrl, MRI *mri_dst, int niter)
{
  int x, y, z, n, width, height, depth, histo[256], xk, yk, zk, xi, yi, zi, val, i, max_histo, max_i, npts;
  MRI *mri_tmp;

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  mri_tmp = MRIcopy(mri_src, NULL);
  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (n = 0; n < niter; n++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          if (MRIvox(mri_ctrl, x, y, z) == 1) /* an original
                                                 control point
                                                 - don't change */
          {
            MRIvox(mri_dst, x, y, z) = MRIvox(mri_tmp, x, y, z);
            continue;
          }

          memset(histo, 0, sizeof(histo));
          for (npts = 0, zk = -1; zk <= 1; zk++) {
            zi = mri_src->zi[z + zk];
            for (yk = -1; yk <= 1; yk++) {
              yi = mri_src->yi[y + yk];
              for (xk = -1; xk <= 1; xk++) {
                xi = mri_src->xi[x + xk];
                if (MRIvox(mri_ctrl, xi, yi, zi) == 0) {
                  continue;
                }
                npts++;
                val = MRIvox(mri_tmp, xi, yi, zi);
                histo[val]++;
              }
            }
          }
          if (npts == 0) {
            MRIvox(mri_dst, x, y, z) = MRIvox(mri_tmp, x, y, z);
            continue;
          }
          for (max_histo = max_i = i = 0; i < 256; i++) {
            if (histo[i] > max_histo) {
              max_histo = histo[i];
              max_i = i;
            }
          }
          MRIvox(mri_dst, x, y, z) = max_i;
          MRIvox(mri_ctrl, x, y, z) = 2; /* comes from a control point
                                            - use it to spread values */
        }
      }
    }
    MRIcopy(mri_dst, mri_tmp);
  }
  MRIfree(&mri_tmp);
  return (mri_dst);
}

MRI *MRImodeFilter(MRI *mri_src, MRI *mri_dst, int niter)
{
  int x, y, z, n, width, height, depth, *histo, xk, yk, zk, xi, yi, zi, val, i, max_histo, max_i, max_val, start_val;
  MRI *mri_tmp;
  float fmin, fmax;

  mri_tmp = MRIcopy(mri_src, NULL);
  MRIvalRange(mri_src, &fmin, &fmax);
  if (DZERO(fmax)) {
    return (mri_tmp);
  }
  max_val = (int)ceil(fmax);
  histo = (int *)calloc(max_val + 1, sizeof(int));

  if (!mri_dst) {
    mri_dst = MRIclone(mri_src, NULL);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (n = 0; n < niter; n++) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          memset(histo, 0, sizeof(int) * (max_val + 1));
          start_val = MRIgetVoxVal(mri_tmp, x, y, z, 0);
          for (zk = -1; zk <= 1; zk++) {
            zi = mri_src->zi[z + zk];
            for (yk = -1; yk <= 1; yk++) {
              yi = mri_src->yi[y + yk];
              for (xk = -1; xk <= 1; xk++) {
                xi = mri_src->xi[x + xk];
                val = nint(MRIgetVoxVal(mri_tmp, xi, yi, zi, 0));
                histo[val]++;
              }
            }
          }
          max_histo = histo[start_val];
          max_i = start_val;
          for (i = 0; i <= max_val; i++) {
            if (histo[i] > max_histo) {
              max_histo = histo[i];
              max_i = i;
            }
          }
          if (max_i > 0) {
            DiagBreak();
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, max_i);
        }
      }
    }
    MRIcopy(mri_dst, mri_tmp);
  }
  free(histo);
  MRIfree(&mri_tmp);
  return (mri_dst);
}

MRI *MRIgradientDir2ndDerivative(MRI *mri_src, MRI *mri_dst, int wsize)
{
  float dI_d2;
  int x, y, z, width, height, depth;

  mri_dst = MRIclone(mri_src, mri_dst);

  width = mri_dst->width;
  height = mri_dst->height;
  depth = mri_dst->depth;

  if (mri_dst->type == MRI_UCHAR) {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          dI_d2 = MRIvoxelGradientDir2ndDerivative(mri_src, x, y, z, wsize);
          MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(127 * (1 + tanh(dI_d2 / 30)));
        }
      }
    }
  }
  else {
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          dI_d2 = MRIvoxelGradientDir2ndDerivative(mri_src, x, y, z, wsize);
          MRIsetVoxVal(mri_dst, x, y, z, 0, 127 * (1 + tanh(dI_d2 / 30)));
        }
      }
    }
  }

  return (mri_dst);
}

MRI *MRIsmoothLabel(MRI *mri_intensity, MRI *mri_label, MRI *mri_smooth, int niter, int label, float min_change)
{
  int x, y, z, n, xi, yi, zi, xk, yk, zk, i, l;
  float val, val_mean, change, max_change;
  MRI *mri_tmp;

  mri_tmp = MRIcopy(mri_intensity, NULL);
  mri_smooth = MRIcopy(mri_intensity, mri_smooth);

  for (i = 0; i < niter; i++) {
    max_change = 0.0;
    for (x = 0; x < mri_tmp->width; x++) {
      for (y = 0; y < mri_tmp->height; y++) {
        for (z = 0; z < mri_tmp->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();

          l = nint(MRIgetVoxVal(mri_label, x, y, z, 0));
          if (l != label) continue;

          val_mean = 0;
          n = 0;
          for (xk = -1; xk <= 1; xk++) {
            xi = mri_tmp->xi[xk + x];
            for (yk = -1; yk <= 1; yk++) {
              yi = mri_tmp->yi[yk + y];
              for (zk = -1; zk <= 1; zk++) {
                zi = mri_tmp->zi[zk + z];
                l = nint(MRIgetVoxVal(mri_label, xi, yi, zi, 0));
                if (l != label) continue;

                val = MRIgetVoxVal(mri_smooth, xi, yi, zi, 0);
                val_mean += val;
                n++;
              }
            }
          }
          if (n > 0) {
            val = MRIgetVoxVal(mri_smooth, x, y, z, 0);
            val_mean /= (float)n;
            change = fabs(val - val_mean);
            if (change > max_change) max_change = change;
            MRIsetVoxVal(mri_tmp, x, y, z, 0, val_mean);
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_smooth);
    if (max_change < min_change) break;
  }

  MRIfree(&mri_tmp);
  return (mri_smooth);
}

MRI *MRIsmoothLabel6Connected(
    MRI *mri_intensity, MRI *mri_label, MRI *mri_smooth, int niter, int label, int fixed_label, float min_change)
{
  int x, n, i, nthreads;
  float max_change[_MAX_FS_THREADS], max_change_overall;
  MRI *mri_tmp;

  mri_tmp = MRIcopy(mri_intensity, NULL);
  mri_smooth = MRIcopy(mri_intensity, mri_smooth);
#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif

  for (i = 0; i < niter; i++) {
    memset(max_change, 0, sizeof(max_change));
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate(max_change) shared(mri_tmp, mri_label, Gx, Gy, Gz) schedule(static, 1)
#endif
    for (x = 0; x < mri_tmp->width; x++) {
      ROMP_PFLB_begin
      
      int tid, y, z, xi, yi, zi, xk, yk, zk, l, n;
      float change, val, val_mean;
#ifdef HAVE_OPENMP
      tid = omp_get_thread_num();
#else
      tid = 0;
#endif
      for (y = 0; y < mri_tmp->height; y++) {
        for (z = 0; z < mri_tmp->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();

          l = nint(MRIgetVoxVal(mri_label, x, y, z, 0));
          if (l != label) continue;

          val_mean = 0;
          n = 0;
          for (xk = -1; xk <= 1; xk++) {
            xi = mri_tmp->xi[xk + x];
            for (yk = -1; yk <= 1; yk++) {
              yi = mri_tmp->yi[yk + y];
              for (zk = -1; zk <= 1; zk++) {
                if (abs(xk) + abs(yk) + abs(zk) > 1) continue;  // enforce 6-connectivity

                zi = mri_tmp->zi[zk + z];
                l = nint(MRIgetVoxVal(mri_label, xi, yi, zi, 0));
                if (l != label && l != fixed_label) continue;

                val = MRIgetVoxVal(mri_smooth, xi, yi, zi, 0);
                val_mean += val;
                n++;
              }
            }
          }
          if (n > 0) {
            val = MRIgetVoxVal(mri_smooth, x, y, z, 0);
            val_mean /= (float)n;
            change = fabs(val - val_mean);
            if (change > max_change[tid]) max_change[tid] = change;
            MRIsetVoxVal(mri_tmp, x, y, z, 0, val_mean);
          }
        }
      }
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    MRIcopy(mri_tmp, mri_smooth);
    for (n = 0, max_change_overall = 0.0; n < nthreads; n++)
      if (max_change[n] > max_change_overall) max_change_overall = max_change[n];

    if (max_change_overall < min_change) break;
  }

  MRIfree(&mri_tmp);
  return (mri_smooth);
}

MRI *MRIxDerivative(MRI *mri_src, MRI *mri_dx)
{
  int x, y, z;
  double dx;

  if (mri_dx == NULL) {
    mri_dx = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dx);
  }

  for (x = 0; x < mri_src->width; x++)
    for (y = 0; y < mri_src->height; y++)
      for (z = 0; z < mri_src->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
#if 0
        dx = MRIsampleXDerivative(mri_src, x, y, z,1) ;
#else
        dx = MRIvoxelDx(mri_src, x, y, z);
#endif
        MRIsetVoxVal(mri_dx, x, y, z, 0, dx);
      }
  return (mri_dx);
}

MRI *MRIyDerivative(MRI *mri_src, MRI *mri_dy)
{
  int x, y, z;
  double dy;

  if (mri_dy == NULL) {
    mri_dy = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dy);
  }

  for (x = 0; x < mri_src->width; x++)
    for (y = 0; y < mri_src->height; y++)
      for (z = 0; z < mri_src->depth; z++) {
#if 0
        dy = MRIsampleYDerivative(mri_src, x, y, z,1) ;
#else
        dy = MRIvoxelDy(mri_src, x, y, z);
#endif
        MRIsetVoxVal(mri_dy, x, y, z, 0, dy);
      }
  return (mri_dy);
}

MRI *MRIzDerivative(MRI *mri_src, MRI *mri_dz)
{
  int x, y, z;
  double dz;

  if (mri_dz == NULL) {
    mri_dz = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dz);
  }

  for (x = 0; x < mri_src->width; x++)
    for (y = 0; y < mri_src->height; y++)
      for (z = 0; z < mri_src->depth; z++) {
#if 0
        dz = MRIsampleZDerivative(mri_src, x, y, z,1) ;
#else
        dz = MRIvoxelDz(mri_src, x, y, z);
#endif
        MRIsetVoxVal(mri_dz, x, y, z, 0, dz);
      }
  return (mri_dz);
}

/*!
 \fn MRI *MRI_Gaussian(int len, float std, int norm, int xsize, int ysize, int zsize)
 \param len - dimension of the MRI image
 \param std - standard deviation of the gaussian filter
 \param norm - 1 for normalization of the gaussian field

  \brief creates an MRI image with a 3D gaussian distribution
  */
MRI *MRI_Gaussian(int len, float std, int norm, int xsize, int ysize, int zsize)
{
  int x, y, z;
  float val;
  float var = std * std;
  float f = sqrt(2 * M_PI) * std;

  MRI *g;
  g = MRIallocSequence(len, len, len, MRI_FLOAT, 1);
  for (y = 0; y < len; y++)
    for (x = 0; x < len; x++)
      for (z = 0; z < len; z++) {
        val = exp((2 - FFTdist(x, y, z, len)) / (2 * var)) / f;
        MRIsetVoxVal(g, x, y, z, 0, val);
      }
  return (g);
}
/*!
 \fn MRI* MRI_fft(MRI *mri_src, MRI* dst)
 \param mri_src - the MRI source (spatial)
 \param dst - the MRI destination (freq - 2 frames)

  \brief MRI_fft applies a FFT transform to
  the mri_src. The result is in 2 frames : the
  first is the modulus (whose quarters has been
  shifted), the second is the argument
*/
MRI *MRI_fft(MRI *mri_src, MRI *dst)
{
  int dimension, w, h, d;
  int x, y, z, k, j, f;
  float ***src, ***dst_r, ***dst_i;
  MRI *new_mri;

  w = mri_src->width;
  h = mri_src->height;
  d = mri_src->depth;

  // be sure it is a power of 2 size image
  if (!FFTisPowerOf2(mri_src->depth) || !FFTisPowerOf2(mri_src->height) || !FFTisPowerOf2(mri_src->width)) {
    if (d >= w) {
      dimension = (d >= h ? d : h);
    }
    else {
      dimension = (h >= w ? h : w);
    }
    dimension = FFTpow2(FFTlog2(dimension));
    printf(
        "The dimension is not a power of 2. "
        "Build a new image (size : %i)\n",
        dimension);

    new_mri = MRIallocSequence(dimension, dimension, dimension, MRI_FLOAT, mri_src->nframes);
    MRIcopyHeader(mri_src, new_mri);

    for (f = 0; f < mri_src->nframes; f++)
      for (z = 0; z < dimension; z++)
        for (y = 0; y < dimension; y++)
          for (x = 0; x < dimension; x++) {
            if (z >= mri_src->depth || y >= mri_src->height || x >= mri_src->width) {
              MRIsetVoxVal(new_mri, x, y, z, f, 0);
            }
            else
              MRIsetVoxVal(new_mri, x, y, z, f, MRIgetVoxVal(mri_src, x, y, z, f));
          }
    mri_src = new_mri;
  }
  else {
    dimension = mri_src->depth;
  }

  if (dst == NULL) {
    dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT, 2 * mri_src->nframes);
    if (dst == NULL) {
      printf("ERROR: MRI_fft: could not alloc\n");
      return (NULL);
    }
    MRIcopy(mri_src, dst);
  }
  else {
    if (mri_src->width != dst->width) {
      printf("ERROR: MRI_fft: width dimension mismatch\n");
      return (NULL);
    }
    if (mri_src->height != dst->height) {
      printf("ERROR: MRI_fft: height dimension mismatch\n");
      return (NULL);
    }
    if (mri_src->depth != dst->depth) {
      printf("ERROR: MRI_fft: depth dimension mismatch\n");
      return (NULL);
    }
    if (2 * mri_src->nframes != dst->nframes) {
      printf("ERROR: MRI_fft: frames dimension mismatch\n");
      return (NULL);
    }
    MRIcopyHeader(mri_src, dst);
  }

  src = (float ***)malloc(dimension * sizeof(float **));
  dst_r = (float ***)malloc(dimension * sizeof(float **));
  dst_i = (float ***)malloc(dimension * sizeof(float **));
  for (k = 0; k < dimension; k++) {
    src[k] = (float **)malloc(dimension * sizeof(float *));
    dst_r[k] = (float **)malloc(dimension * sizeof(float *));
    dst_i[k] = (float **)malloc(dimension * sizeof(float *));
    for (j = 0; j < dimension; j++) {
      src[k][j] = (float *)malloc(dimension * sizeof(float));
      dst_r[k][j] = (float *)malloc(dimension * sizeof(float));
      dst_i[k][j] = (float *)malloc(dimension * sizeof(float));
    }
  }

  for (f = 0; f < mri_src->nframes; f++) {
    // Initialize the image and intermediate MRI
    for (z = 0; z < dimension; z++)
      for (y = 0; y < dimension; y++)
        for (x = 0; x < dimension; x++) {
          src[x][y][z] = MRIgetVoxVal(mri_src, x, y, z, f);
        }

    // FFT
    if (f == 0) {
      printf("FFT : 1D -");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        RFFTforward(src[x][y], dimension, dst_r[x][y], dst_i[x][y]);
      }
    FFTswitch_with_z(dst_r, dimension, 0);
    FFTswitch_with_z(dst_i, dimension, 0);
    if (f == 0) {
      printf(" 2D -");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        CFFTforward(dst_r[x][y], dst_i[x][y], dimension);
      }
    FFTswitch_with_z(dst_r, dimension, 0);
    FFTswitch_with_z(dst_i, dimension, 0);
    FFTswitch_with_z(dst_r, dimension, 1);
    FFTswitch_with_z(dst_i, dimension, 1);
    if (f == 0) {
      printf(" 3D\n");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        CFFTforward(dst_r[x][y], dst_i[x][y], dimension);
      }
    FFTswitch_with_z(dst_r, dimension, 1);
    FFTswitch_with_z(dst_i, dimension, 1);

    // from (real, imaginary) to (mudulus, argument)
    FFTreim_to_modarg(dst_r, dst_i, dimension);
    dst_r = FFTinv_quarter(dst_r, dimension);
    // Write the result
    for (z = 0; z < dimension; z++)
      for (y = 0; y < dimension; y++)
        for (x = 0; x < dimension; x++) {
          MRIsetVoxVal(dst, x, y, z, 2 * f, dst_r[x][y][z]);
          MRIsetVoxVal(dst, x, y, z, 2 * f + 1, dst_i[x][y][z]);
        }
  }
  free(src);
  free(dst_r);
  free(dst_i);
  return (dst);
}
/*!
 \fn MRI *MRI_ifft(MRI *src, MRI *dst, int w, int h, int d)
 \param src - the MRI source (freq domain - 2 frames)
 \param dst - the MRI destination (back into spatial domain)
 \param w, h, d - dimension of the original image in the
 spatial domain

  \brief MRI_ifft applies an inverse FFT to the mri_src :
  the modulus on the first frame and the argument in the
  second. The output is a single frame MRI with the real
   inverse.w h d, specifies the output size, if different
   from the size of the image in the frequency domain.
*/

MRI *MRI_ifft(MRI *src, MRI *dst, int w, int h, int d)
{
  int dimension;
  int x, y, z, k, j, f;
  float ***dst_r, ***dst_i;

  dimension = src->width;

  if (src->nframes % 2 != 0) {
    printf("ERROR: MRI_ifft: frames of input\n");
  }

  if (dst == NULL) {
    dst = MRIallocSequence(w, h, d, MRI_FLOAT, src->nframes / 2);
    if (dst == NULL) {
      printf("ERROR: MRI_ifft: could not alloc\n");
      return (NULL);
    }
    MRIcopyHeader(src, dst);
  }
  else {
    if (w != dst->width) {
      printf("ERROR: MRI_ifft: width dimension mismatch\n");
      return (NULL);
    }
    if (h != dst->height) {
      printf("ERROR: MRI_ifft: height dimension mismatch\n");
      return (NULL);
    }
    if (d != dst->depth) {
      printf("ERROR: MRI_ifft: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != 2 * dst->nframes) {
      printf("ERROR: MRI_ifft: frames dimension mismatch\n");
      return (NULL);
    }
    MRIcopyHeader(src, dst);
  }

  dst_r = (float ***)malloc(dimension * sizeof(float **));
  dst_i = (float ***)malloc(dimension * sizeof(float **));
  for (k = 0; k < dimension; k++) {
    dst_r[k] = (float **)malloc(dimension * sizeof(float *));
    dst_i[k] = (float **)malloc(dimension * sizeof(float *));
    for (j = 0; j < dimension; j++) {
      dst_r[k][j] = (float *)malloc(dimension * sizeof(float));
      dst_i[k][j] = (float *)malloc(dimension * sizeof(float));
    }
  }

  for (f = 0; f < src->nframes / 2; f++) {
    // Initialize the image and intermediate MRI
    for (z = 0; z < dimension; z++)
      for (y = 0; y < dimension; y++)
        for (x = 0; x < dimension; x++) {
          dst_r[x][y][z] = MRIgetVoxVal(src, x, y, z, 2 * f);
          dst_i[x][y][z] = MRIgetVoxVal(src, x, y, z, 2 * f + 1);
        }

    // back to (real, imaginary)
    dst_r = FFTinv_quarter(dst_r, dimension);
    FFTmodarg_to_reim(dst_r, dst_i, dimension);

    // FFT inverse
    if (f == 0) {
      printf("inverse FFT : 1D -");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        CFFTbackward(dst_r[x][y], dst_i[x][y], dimension);
      }
    FFTswitch_with_z(dst_r, dimension, 0);
    FFTswitch_with_z(dst_i, dimension, 0);
    if (f == 0) {
      printf(" 2D -");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        CFFTbackward(dst_r[x][y], dst_i[x][y], dimension);
      }
    FFTswitch_with_z(dst_r, dimension, 0);
    FFTswitch_with_z(dst_i, dimension, 0);
    FFTswitch_with_z(dst_r, dimension, 1);
    FFTswitch_with_z(dst_i, dimension, 1);
    if (f == 0) {
      printf(" 3D\n");
    }
    for (x = 0; x < dimension; x++)
      for (y = 0; y < dimension; y++) {
        CFFTbackward(dst_r[x][y], dst_i[x][y], dimension);
      }
    FFTswitch_with_z(dst_r, dimension, 1);
    FFTswitch_with_z(dst_i, dimension, 1);

    // Write the result
    for (z = 0; z < d; z++)
      for (y = 0; y < h; y++)
        for (x = 0; x < w; x++) {
          MRIsetVoxVal(dst, x, y, z, f, dst_r[x][y][z]);
        }
  }

  free(dst_r);
  free(dst_i);
  return (dst);
}
/*!
 \fn MRI *MRI_fft_gaussian(MRI *src, MRI *dst, float std, int norm)
 \param src - the MRI source (spatial domain)
 \param dst - the MRI destination (spatial domain)
 \param std - standard deviation of the gaussian filter
 \param norm - normalizes the gaussian field

  \brief MRI_fft_smooth computes the fft transform of src, then
  apply a gaussian smoothing filter in the freq domain: ie multiply
  the shifted modulus by a gaussian MRI, and then computes the
  inverse FFT. The result is in dst.
*/
MRI *MRI_fft_gaussian(MRI *src, MRI *dst, float std, int norm)
{
  MRI *src_fft = NULL, *fft_sm = NULL, *g;
  MRI *g_fft, *temp;
  int dimension, f;

  if (dst == NULL) {
    dst = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (dst == NULL) {
      printf("ERROR: MRI_fft_gaussian: could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, dst);
  }
  else {
    if (src->width != dst->width) {
      printf("ERROR: MRI_fft_gaussian: width dimension mismatch\n");
      return (NULL);
    }
    if (src->height != dst->height) {
      printf("ERROR: MRI_fft_gaussian: height dimension mismatch\n");
      return (NULL);
    }
    if (src->depth != dst->depth) {
      printf("ERROR: MRI_fft_gaussian: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != dst->nframes) {
      printf("ERROR: MRI_fft_gaussian: frames dimension mismatch\n");
      return (NULL);
    }
    if (src != dst) {
      MRIcopy(src, dst);
    }
  }

  src_fft = MRI_fft(src, src_fft);
  dimension = src_fft->width;

  printf("Gaussian ");
  g = MRI_Gaussian(dimension, std, norm, src_fft->xsize, src_fft->ysize, src_fft->zsize);
  g_fft = MRI_fft(g, NULL);
  MRIfree(&g);

  temp = MRIallocSequence(dimension, dimension, dimension, MRI_FLOAT, 1);
  for (f = 0; f < src->nframes; f++) {
    temp = MRIcopyFrame(src_fft, temp, 2 * f, 0);
    fft_sm = MRImultiply(temp, g_fft, NULL);  // only frame 0
    src_fft = MRIcopyFrame(fft_sm, src_fft, 0, 2 * f);
  }
  MRI_ifft(src_fft, dst, src->width, src->height, src->depth);
  MRIfree(&src_fft);

  return (dst);
}

/*!
 \fn MRI *MRI_fft_lowpass(MRI *src, MRI *dst, int percent)
 \param src - the MRI source (spatial domain)
 \param dst - the MRI destination (spatial domain)
 \param percent - parameter of the lowpass filter

  \brief MRI_fft_smooth computes the fft transform of src, then
  apply a lowpass filter in the freq domain: ie keep only
  "percent" percent of the lowest freqs, and then computes the
  inverse FFT. The result is in dst.
*/

MRI *MRI_fft_lowpass(MRI *src, MRI *dst, int percent)
{
  MRI *src_fft;
  int x, y, z, f;

  if (dst == NULL) {
    dst = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (dst == NULL) {
      printf("ERROR: MRI_fft_lowpass: could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, dst);
  }
  else {
    if (src->width != dst->width) {
      printf("ERROR: MRI_fft_lowpass: width dimension mismatch\n");
      return (NULL);
    }
    if (src->height != dst->height) {
      printf("ERROR: MRI_fft_lowpass: height dimension mismatch\n");
      return (NULL);
    }
    if (src->depth != dst->depth) {
      printf("ERROR: MRI_fft_lowpass: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != dst->nframes) {
      printf("ERROR: MRI_fft_lowpass: frames dimension mismatch\n");
      return (NULL);
    }
    if (src != dst) {
      MRIcopy(src, dst);
    }
  }

  src_fft = MRI_fft(src, NULL);
  float threshold = (float)percent * src_fft->depth / 100;

  for (f = 0; f < src->nframes; f++)
    for (z = 0; z < src_fft->depth; z++)
      for (y = 0; y < src_fft->height; y++)
        for (x = 0; x < src_fft->width; x++) {
          if (FFTdist(x, y, z, src_fft->depth) > threshold * threshold) {
            MRIsetVoxVal(src_fft, x, y, z, 2 * f, 0);
          }
        }

  MRI_ifft(src_fft, dst, src->width, src->height, src->depth);
  MRIfree(&src_fft);
  return (dst);
}

/*!
 \fn MRI *MRI_fft_highpass(MRI *src, MRI *dst, int percent)
 \param src - the MRI source (spatial domain)
 \param dst - the MRI destination (spatial domain)
 \param percent - parameter of the highpass filter

  \brief MRI_fft_smooth computes the fft transform of src, then
  apply a lowpass filter in the freq domain: ie remove the
  "percent" percent of the lowest freqs, and then computes the
  inverse FFT. The result is in dst.
*/
MRI *MRI_fft_highpass(MRI *src, MRI *dst, int percent)
{
  MRI *src_fft;
  int x, y, z, f;

  if (dst == NULL) {
    dst = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (dst == NULL) {
      printf("ERROR: MRI_fft_highpass: could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, dst);
  }
  else {
    if (src->width != dst->width) {
      printf("ERROR: MRI_fft_highpass: width dimension mismatch\n");
      return (NULL);
    }
    if (src->height != dst->height) {
      printf("ERROR: MRI_fft_highpass: height dimension mismatch\n");
      return (NULL);
    }
    if (src->depth != dst->depth) {
      printf("ERROR: MRI_fft_highpass: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != dst->nframes) {
      printf("ERROR: MRI_fft_highpass: frames dimension mismatch\n");
      return (NULL);
    }
    if (src != dst) {
      MRIcopy(src, dst);
    }
  }
  src_fft = MRI_fft(src, NULL);
  float threshold = (float)percent * src_fft->depth / 100;

  for (f = 0; f < src->nframes; f++)
    for (z = 0; z < src_fft->depth; z++)
      for (y = 0; y < src_fft->height; y++)
        for (x = 0; x < src_fft->width; x++) {
          if (FFTdist(x, y, z, src_fft->depth) < threshold * threshold) {
            MRIsetVoxVal(src_fft, x, y, z, 2 * f, 0);
          }
        }

  MRI_ifft(src_fft, dst, src->width, src->height, src->depth);
  MRIfree(&src_fft);
  return (dst);
}
/*---------------------------------------------------------------------
  MRIgaussianSmoothNI() - performs non-isotropic gaussian spatial
  smoothing.  The standard deviation of the gaussian is std.  The mean
  is preserved (ie, sets the kernel integral to 1).  Can be done
  in-place. Handles multiple frames. See also MRIconvolveGaussian()
  and MRImaskedGaussianSmooth(). Has the capacity to do a 2-Gaussian
  mixture model using external variables. 
  -------------------------------------------------------------------*/
MRI *MRIgaussianSmoothNI(MRI *src, double cstd, double rstd, double sstd, MRI *targ)
{
  int c, r, s, cstop, rstop, sstop;
  MATRIX *G;
  MATRIX *vr = NULL, *vc = NULL, *vs = NULL;
  long double scale, vmf;
  // mixture model 
  extern float smni_cw1, smni_cstd2, smni_rw1, smni_rstd2, smni_sw1, smni_sstd2;

  if (targ == NULL) {
    targ = MRIallocSequence(src->width, src->height, src->depth, MRI_FLOAT, src->nframes);
    if (targ == NULL) {
      printf("ERROR: MRIgaussianSmoothNI: could not alloc\n");
      return (NULL);
    }
    MRIcopy(src, targ);
  }
  else {
    if (src->width != targ->width) {
      printf("ERROR: MRIgaussianSmoothNI: width dimension mismatch\n");
      return (NULL);
    }
    if (src->height != targ->height) {
      printf("ERROR: MRIgaussianSmoothNI: height dimension mismatch\n");
      return (NULL);
    }
    if (src->depth != targ->depth) {
      printf("ERROR: MRIgaussianSmoothNI: depth dimension mismatch\n");
      return (NULL);
    }
    if (src->nframes != targ->nframes) {
      printf("ERROR: MRIgaussianSmoothNI: frames dimension mismatch\n");
      return (NULL);
    }
    if (src != targ) MRIcopy(src, targ);
  }

#ifdef HAVE_OPENMP
  if (Gdiag_no > 0)
    printf("MRIgaussianSmoothNI(): %d avail. processors, running in %d threads\n", omp_get_num_procs(), omp_get_max_threads());
#endif

  if(Gdiag_no > 0 && (smni_cw1 != 1 || smni_rw1 != 1 || smni_sw1 != 1)){
    printf("MRIgaussianSmoothNI(): G2: %g %g %g  %g %g %g\n",
	   smni_cw1, smni_cstd2, smni_rw1, smni_rstd2, smni_sw1, smni_sstd2);
  }
  fflush(stdout);

  /* -----------------Smooth the columns -----------------------------*/
  if (cstd > 0) {
    if(smni_cw1 == 1) G = GaussianMatrix(src->width, cstd / src->xsize, 1, NULL);
    else  G = GaussianMatrix2(src->width, cstd/src->xsize, smni_cstd2/src->xsize, smni_cw1, 1, NULL);
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for (r = 0; r < src->height; r++) {
      ROMP_PFLB_begin
      int s, f, c;
      MATRIX *v = MatrixAlloc(src->width, 1, MATRIX_REAL);
      MATRIX *vg = MatrixAlloc(src->width, 1, MATRIX_REAL);
      if (Gdiag_no > 0 && DIAG_VERBOSE_ON) {
        printf("%d ", r);
        if (r % 10 == 9) printf("\n");
      }
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) {
          for (c = 0; c < src->width; c++) v->rptr[c + 1][1] = MRIgetVoxVal(targ, c, r, s, f);
          MatrixMultiply(G, v, vg);
          for (c = 0; c < src->width; c++) MRIsetVoxVal(targ, c, r, s, f, vg->rptr[c + 1][1]);
        }
      }
      MatrixFree(&v);
      MatrixFree(&vg);
      ROMP_PFLB_end
    }
    ROMP_PF_end
    if (Gdiag_no > 0 && DIAG_VERBOSE_ON) printf("\n");
    // This is for scaling
    vc = MatrixAlloc(src->width, 1, MATRIX_REAL);
    if (src->width > 1)
      for (c = 0; c < src->width; c++) vc->rptr[c + 1][1] = G->rptr[src->width / 2][c + 1];
    else
      vc->rptr[1][1] = 1;
    MatrixFree(&G);
  }

  /* -----------------Smooth the rows -----------------------------*/
  if (rstd > 0) {
    if (Gdiag_no > 0 && DIAG_VERBOSE_ON) printf("Smoothing rows\n");
    if(smni_rw1 == 1) G = GaussianMatrix(src->height, (double)rstd / src->ysize, 1, NULL);
    else  G = GaussianMatrix2(src->height, rstd/src->ysize, smni_rstd2/src->ysize, smni_rw1, 1, NULL);

    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for (c = 0; c < src->width; c++) {
      ROMP_PFLB_begin
      
      int s, f, r;
      MATRIX *v = MatrixAlloc(src->height, 1, MATRIX_REAL);
      MATRIX *vg = MatrixAlloc(src->height, 1, MATRIX_REAL);
      if (Gdiag_no > 0) {
        printf("%d ", c);
        if (c % 10 == 9) printf("\n");
      }
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) {
          for (r = 0; r < src->height; r++) v->rptr[r + 1][1] = MRIgetVoxVal(targ, c, r, s, f);
          MatrixMultiply(G, v, vg);
          for (r = 0; r < src->height; r++) MRIsetVoxVal(targ, c, r, s, f, vg->rptr[r + 1][1]);
        }
      }
      MatrixFree(&v);
      MatrixFree(&vg);
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    if (Gdiag_no > 0) printf("\n");

    // This is for scaling
    vr = MatrixAlloc(src->height, 1, MATRIX_REAL);
    if (src->height > 1)
      for (r = 0; r < src->height; r++) vr->rptr[r + 1][1] = G->rptr[src->height / 2][r + 1];
    else
      vr->rptr[1][1] = 1;
    MatrixFree(&G);
  }

  /* Smooth the slices */
  if (sstd > 0) {
    // printf("Smoothing slices by std=%g\n",sstd);
    
    if(smni_sw1 == 1) G = GaussianMatrix(src->depth, sstd / src->zsize, 1, NULL);
    else  G = GaussianMatrix2(src->depth, sstd/src->zsize, smni_sstd2/src->zsize, smni_sw1, 1, NULL);

    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for (c = 0; c < src->width; c++) {
      ROMP_PFLB_begin
      
      int r, f, s;
      MATRIX *v = MatrixAlloc(src->depth, 1, MATRIX_REAL);
      MATRIX *vg = MatrixAlloc(src->depth, 1, MATRIX_REAL);
      if (Gdiag_no > 0) {
        printf("%d ", c);
        if (c % 10 == 9) printf("\n");
      }
      for (r = 0; r < src->height; r++) {
        for (f = 0; f < src->nframes; f++) {
          for (s = 0; s < src->depth; s++) v->rptr[s + 1][1] = MRIgetVoxVal(targ, c, r, s, f);
          MatrixMultiply(G, v, vg);
          for (s = 0; s < src->depth; s++) MRIsetVoxVal(targ, c, r, s, f, vg->rptr[s + 1][1]);
        }
      }
      MatrixFree(&v);
      MatrixFree(&vg);
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    if (Gdiag_no > 0) printf("\n");
    // This is for scaling
    vs = MatrixAlloc(src->depth, 1, MATRIX_REAL);
    if (src->depth > 1)
      for (s = 0; s < src->depth; s++) vs->rptr[s + 1][1] = G->rptr[src->depth / 2][s + 1];
    else
      vs->rptr[1][1] = 1;
    MatrixFree(&G);
  }

  // Compute the sum of the kernel. Note: the expected variance
  // will be sum(k^2)
  scale = 0;
  vmf = 0;
  if (cstd > 0.0)
    cstop = src->width;
  else
    cstop = 1;
  if (rstd > 0.0)
    rstop = src->height;
  else
    rstop = 1;
  if (sstd > 0.0)
    sstop = src->depth;
  else
    sstop = 1;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : scale, vmf)
#endif
  for (c = 0; c < cstop; c++) {
    ROMP_PFLB_begin
    int r, s;
    double aa = 1, bb = 1, cc = 1, val;
    for (r = 0; r < rstop; r++) {
      for (s = 0; s < sstop; s++) {
        if (cstd > 0.0) aa = vc->rptr[c + 1][1];
        if (rstd > 0.0) bb = vr->rptr[r + 1][1];
        if (sstd > 0.0) cc = vs->rptr[s + 1][1];
        val = (aa * bb * cc);
        scale += (val);
        vmf += (val * val);
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  if (Gdiag_no > 0) {
    printf("MRIguassianSmoothNI(): scale = %Lg\n", scale);
    printf("MRIguassianSmoothNI(): VMF = %Lg, VRF = %Lg\n", vmf, 1.0 / vmf);
  }

/* With 'double scale, vmf', the following section of multi-threaded codes will cause small differences in voxel values 
 * when smoothing large images with fwhm >= 6 across different number of threads.
 * Changing 'double scale, vmf' to 'long double' resolves the differences.
 */
// Divide by the sum of the kernel so that a smoothed delta function
// will sum to one and so that a constant input yields const output.
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (c = 0; c < src->width; c++) {
    ROMP_PFLB_begin
    int r, s, f;
    double val;
    for (r = 0; r < src->height; r++) {
      for (s = 0; s < src->depth; s++) {
        for (f = 0; f < src->nframes; f++) {
          val = MRIgetVoxVal(targ, c, r, s, f);
          MRIsetVoxVal(targ, c, r, s, f, val / scale);
        }
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  if (vc) MatrixFree(&vc);
  if (vr) MatrixFree(&vr);
  if (vs) MatrixFree(&vs);

  return (targ);
}

MRI *MRIsegmentationSurfaceNormals(MRI *mri_seg, MRI *mri_normals, int label, MRI **pmri_ctrl)
{
  int x, y, z, border;
  float nx, ny, nz;
  MRI *mri_border;

  mri_border = MRImarkLabelBorderVoxels(mri_seg, NULL, label, 1, 0);
  mri_normals = MRIallocSequence(mri_seg->width, mri_seg->height, mri_seg->depth, MRI_FLOAT, 3);
  MRIcopyHeader(mri_seg, mri_normals);

  for (x = 0; x < mri_seg->width; x++) {
    for (y = 0; y < mri_seg->height; y++) {
      for (z = 0; z < mri_seg->height; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        border = nint(MRIgetVoxVal(mri_border, x, y, z, 0));
        if (border == 0) {
          continue;
        }
        if (MRIcomputeBorderNormalAtVoxel(mri_seg, x, y, z, &nx, &ny, &nz, label) > 0) {
          MRIsetVoxVal(mri_normals, x, y, z, 0, nx);
          MRIsetVoxVal(mri_normals, x, y, z, 1, ny);
          MRIsetVoxVal(mri_normals, x, y, z, 2, nz);
          if (pmri_ctrl && (!FZERO(nx) || !FZERO(ny) || !(FZERO(nz)))) {
            MRIsetVoxVal(*pmri_ctrl, x, y, z, 0, CONTROL_MARKED);
          }
        }
      }
    }
  }

  MRIfree(&mri_border);
  return (mri_normals);
}

int MRIcomputeBorderNormalAtVoxel(
    MRI *mri_seg, int x0, int y0, int z0, float *pnx, float *pny, float *pnz, int target_label)
{
  int olabel, x1, y1, z1, xk, yk, zk, label, num;
  float nx, ny, nz, mag;

  olabel = (int)MRIgetVoxVal(mri_seg, x0, y0, z0, 0);

  for (nx = ny = nz = 0.0, num = 0, xk = -1; xk <= 1; xk++) {
    x1 = mri_seg->xi[x0 + xk];
    for (yk = -1; yk <= 1; yk++) {
      y1 = mri_seg->yi[y0 + yk];
      for (zk = -1; zk <= 1; zk++) {
        z1 = mri_seg->zi[z0 + zk];
        if (fabs(xk) + fabs(yk) + fabs(zk) > 1) {
          continue;  // only 6-connected
        }
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
  else {
    DiagBreak();
  }

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

  return (num);
}

/*
  normalize the length of the sequence at each voxel to be unity.
*/
MRI *MRInormalizeFrameVectorLength(MRI *mri_src, MRI *mri_dst)
{
  int x, y, z, f;
  float val, mag;

  if (mri_src != mri_dst) {
    MRIcopy(mri_src, mri_dst);
  }

  for (x = 0; x < mri_dst->width; x++)
    for (y = 0; y < mri_dst->height; y++)
      for (z = 0; z < mri_dst->depth; z++) {
        for (mag = 0.0, f = 0; f < mri_dst->nframes; f++) {
          val = MRIgetVoxVal(mri_src, x, y, z, f);
          mag += (val * val);
        }
        mag = sqrt(mag / mri_src->nframes);
        if (!FZERO(mag)) {
          for (f = 0; f < mri_dst->nframes; f++) {
            val = MRIgetVoxVal(mri_src, x, y, z, f);
            val /= mag;
            MRIsetVoxVal(mri_dst, x, y, z, f, val);
          }
        }
      }

  return (mri_dst);
}
/*
  normalize the length of the sequence at each voxel to be unity.
*/
MRI *MRIcomputeFrameVectorLength(MRI *mri_src, MRI *mri_dst)
{
  int x, y, z, f;
  float val, mag;

  if (mri_dst == NULL) {
    mri_dst = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (x = 0; x < mri_dst->width; x++)
    for (y = 0; y < mri_dst->height; y++)
      for (z = 0; z < mri_dst->depth; z++) {
        for (mag = 0.0, f = 0; f < mri_src->nframes; f++) {
          val = MRIgetVoxVal(mri_src, x, y, z, f);
          mag += (val * val);
        }
        mag = sqrt(mag / mri_src->nframes);
        MRIsetVoxVal(mri_dst, x, y, z, 0, mag);
      }

  return (mri_dst);
}

MRI *MRIcomputeFrameVectorL1Length(MRI *mri_src, MRI *mri_dst)
{
  int x, y, z, f;
  float val, mag;

  if (mri_dst == NULL) {
    mri_dst = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT);
    MRIcopyHeader(mri_src, mri_dst);
  }

  for (x = 0; x < mri_dst->width; x++)
    for (y = 0; y < mri_dst->height; y++)
      for (z = 0; z < mri_dst->depth; z++) {
        for (mag = 0.0, f = 0; f < mri_src->nframes; f++) {
          val = MRIgetVoxVal(mri_src, x, y, z, f);
          mag += fabs(val);
        }
        mag /= mri_src->nframes;
        MRIsetVoxVal(mri_dst, x, y, z, 0, mag);
      }

  return (mri_dst);
}

MRI *MRInbrThresholdLabel(MRI *mri_src, MRI *mri_dst, int label, int out_label, int whalf, float thresh)
{
  int x, y, z, count, total;

  if (mri_dst != mri_src) mri_dst = MRIcopy(mri_src, mri_dst);

  for (total = x = 0; x < mri_dst->width; x++)
    for (y = 0; y < mri_dst->height; y++)
      for (z = 0; z < mri_dst->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        if (nint(MRIgetVoxVal(mri_dst, x, y, z, 0)) != label) continue;  // only filter this label
        count = MRIlabelsInNbhd(mri_dst, x, y, z, whalf, label);
        if (count < thresh) {
          total++;
          MRIsetVoxVal(mri_dst, x, y, z, 0, out_label);
        }
      }

  //  printf("%d voxels thresholded due to insufficient neighbors\n", total) ;
  return (mri_dst);
}

/*
  \fn MRI *MRImotionBlur2D(MRI *src, MB2D *mb, MRI *out)
  \brief Apply 2D motion (ie, radial) blurring to a volume.  Motion
  blurring is smoothing in 1D along a radius with the amount of
  smoothing changing with distance from the center. Each slice is
  blurred in the exact same way. The center does not need to be inside
  the volume. The parameters are controlled through the MB2D
  structure.  The basic implementation is to create a line passing
  through a given voxel and the center. A segment of this line is
  selected based on the FHWM at that point. The input image is sampled
  every DeltaD along that line and convolved with the Gaussian
  kernel. Cannot be done in-place. Offset is in mm. Slope is per
  centimeter (not mm).
 */
MRI *MRImotionBlur2D(MRI *src, MB2D *mb, MRI *out)
{
  int c;

  if (src == out) {
    printf("ERROR: MRImotionBlur2D(): cannot be done in-place\n");
    return (NULL);
  }

  if (mb->Interp != SAMPLE_NEAREST && mb->Interp != SAMPLE_TRILINEAR) {
    printf("ERROR: MRImotionBlur2D(): Interp = %d, must be %d or %d\n", mb->Interp, SAMPLE_NEAREST, SAMPLE_TRILINEAR);
    return (NULL);
  }

  if (out == NULL) {
    out = MRIcloneBySpace(src, MRI_FLOAT, -1);
    if (out == NULL) return (NULL);
  }

  if (Gdiag_no > 0) {
    printf("MB: type=%d cutoff=%lf, DeltaD=%lf, offset=%lf, slope=%lf %d c0=%d, r0=%d, cR=%d, rR=%d\n",
           mb->type,
           mb->cutoff,
           mb->DeltaD,
           mb->offset,
           mb->slope,
           mb->Interp,
           mb->c0,
           mb->r0,
           mb->cR,
           mb->rR);
    fflush(stdout);
  }

  // These are two structures to save slice-based parameters
  if (mb->d0) MRIfree(&mb->d0);
  if (mb->theta) MRIfree(&mb->theta);
  if (mb->fwhm) MRIfree(&mb->fwhm);
  if (mb->dmin) MRIfree(&mb->dmin);
  if (mb->nd) MRIfree(&mb->nd);
  mb->theta = MRIalloc(src->width, src->height, 1, MRI_FLOAT);
  mb->dmin = MRIalloc(src->width, src->height, 1, MRI_FLOAT);
  mb->fwhm = MRIalloc(src->width, src->height, 1, MRI_FLOAT);
  mb->nd = MRIalloc(src->width, src->height, 1, MRI_FLOAT);
  mb->d0 = MRIalloc(src->width, src->height, 1, MRI_FLOAT);
  MRIcopyHeader(src, mb->theta);
  MRIcopyHeader(src, mb->dmin);
  MRIcopyHeader(src, mb->fwhm);
  MRIcopyHeader(src, mb->nd);
  MRIcopyHeader(src, mb->d0);

// Fill the slice-based parameters
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < out->width; c++) {
    ROMP_PFLB_begin
    
    int r;
    for (r = 0; r < out->height; r++) {
      double theta, dmin, fwhm, dx, dy, d0, stddev, dlim;
      int dc, dr, nd, ndlim;
      dc = c + mb->cR - mb->c0;            // #cols from center
      dr = r + mb->rR - mb->r0;            // #rows from center
      theta = atan2(dr, dc);               // Angle relative to center
      dx = dc * src->xsize;                // mm distance in x from center
      dy = dr * src->ysize;                // mm distance in y from center
      d0 = sqrt(dx * dx + dy * dy);        // mm dist to cur vox
      // divide slope by 100 to make it more in the range of FWHM
      fwhm = mb->offset + d0*mb->slope/100; // compute FWHM based on distance
      stddev = fwhm / sqrt(log(256.0));
      ndlim = ceil(mb->cutoff * stddev / mb->DeltaD);
      nd = 2 * ndlim + 1;         // number of samples to integrate over
      dlim = ndlim * mb->DeltaD;  // sets the min and (implicitly) the max
      dmin = d0 - dlim;           // distance at which to start the integration
      MRIsetVoxVal(mb->theta, c, r, 0, 0, theta);
      MRIsetVoxVal(mb->fwhm, c, r, 0, 0, fwhm);
      MRIsetVoxVal(mb->dmin, c, r, 0, 0, dmin);
      MRIsetVoxVal(mb->nd, c, r, 0, 0, nd);
      MRIsetVoxVal(mb->d0, c, r, 0, 0, d0);
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

// Apply the smoothing
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (c = 0; c < out->width; c++) {
    ROMP_PFLB_begin
    
    int r;
    for (r = 0; r < out->height; r++) {
      double cd, rd, theta, dmin, fwhm, *kernel, ksum, d0, stddev, d = 0;
      double vsrc = 0, vdst, x, y, cosTheta = 0, sinTheta = 0, dd;
      double thetaarc = 0, thetadelta = 0, thetastart = 0;
      double coef11 = 0, coef12 = 0, coef21 = 0, coef22 = 0, cT, rT;
      int s, f, nd, nthd, cdi = 0, rdi = 0, c1 = 0, c2 = 0, r1 = 0, r2 = 0;
      theta = MRIgetVoxVal(mb->theta, c, r, 0, 0);
      fwhm = MRIgetVoxVal(mb->fwhm, c, r, 0, 0);
      dmin = MRIgetVoxVal(mb->dmin, c, r, 0, 0);
      nd = MRIgetVoxVal(mb->nd, c, r, 0, 0);
      d0 = MRIgetVoxVal(mb->d0, c, r, 0, 0);
      stddev = fwhm / sqrt(log(256.0));

      if (stddev == 0) {
        for (s = 0; s < src->depth; s++) {
          for (f = 0; f < src->nframes; f++) {
            vsrc = MRIgetVoxVal(src, c, r, s, f);
            MRIsetVoxVal(out, c, r, s, f, vsrc);
          }  // slice
        }    // frame
        continue;
      }

      // Compute the smoothing kernel
      kernel = (double *)calloc(nd, sizeof(double));
      ksum = 0;
      for (nthd = 0; nthd < nd; nthd++) {
        dd = dmin + nthd * mb->DeltaD - d0;  // mm dist to cur voxel
        kernel[nthd] = exp(-dd * dd / (2 * stddev * stddev));
        ksum += kernel[nthd];
      }
      for (nthd = 0; nthd < nd; nthd++) kernel[nthd] /= ksum;

      if (mb->type == MB_TANGENTIAL) {
        thetaarc = mb->DeltaD * nd / d0;      // tot angle of integration arc
        thetadelta = mb->DeltaD / d0;         // rad angle of each step
        thetastart = theta - thetaarc / 2.0;  // abs angle of integ start point
        d = d0;
      }
      if (mb->type == MB_RADIAL) {
        cosTheta = cos(theta);
        sinTheta = sin(theta);
      }

      for (nthd = 0; nthd < nd; nthd++) {
        if (mb->type == MB_RADIAL)
          d = dmin + nthd * mb->DeltaD;            // mm distance from center to this point
        else {                                     // Tangential smoothing
          theta = thetastart + nthd * thetadelta;  // distance along circle in rads
          cosTheta = cos(theta);
          sinTheta = sin(theta);
        }
        x = d * cosTheta;                       // mm x distance from center to this point
        y = d * sinTheta;                       // mm y distance from center to this point
        cd = mb->c0 + x / src->xsize - mb->cR;  // float col in REGION of this point
        rd = mb->r0 + y / src->ysize - mb->rR;  // float row in REGION of this point
        if (mb->Interp == SAMPLE_NEAREST) {
          // Compute nearest here so don't have to in slice/frame loop
          cdi = nint(cd);
          rdi = nint(rd);
          if (cdi < 0) continue;
          if (cdi >= src->width) continue;
          if (rdi < 0) continue;
          if (rdi >= src->height) continue;
        }
        if (mb->Interp == SAMPLE_TRILINEAR) {
          // Compute coefficients here so don't have to in slice/frame loop
          c1 = floor(cd);
          c2 = ceil(cd);
          r1 = floor(rd);
          r2 = ceil(rd);
          if (c1 < 0) continue;
          if (c1 >= src->width) continue;
          if (c2 < 0) continue;
          if (c2 >= src->width) continue;
          if (r1 < 0) continue;
          if (r1 >= src->height) continue;
          if (r2 < 0) continue;
          if (r2 >= src->height) continue;
          cT = cd - c1;
          rT = rd - r1;
          coef11 = (1 - cT) * (1 - rT);
          coef12 = (1 - cT) * rT;
          coef21 = cT * (1 - rT);
          coef22 = cT * rT;
        }
        // Loop through all the slices and frames
        for (s = 0; s < src->depth; s++) {
          for (f = 0; f < src->nframes; f++) {
            if (mb->Interp == SAMPLE_NEAREST) vsrc = MRIgetVoxVal(src, cdi, rdi, s, f);
            if (mb->Interp == SAMPLE_TRILINEAR)
              vsrc = coef11 * MRIgetVoxVal(src, c1, r1, s, f) + coef12 * MRIgetVoxVal(src, c1, r2, s, f) +
                     coef21 * MRIgetVoxVal(src, c2, r1, s, f) + coef22 * MRIgetVoxVal(src, c2, r2, s, f);
            vdst = MRIgetVoxVal(out, c, r, s, f);
            vdst += vsrc * kernel[nthd];
            MRIsetVoxVal(out, c, r, s, f, vdst);

          }  // slice
        }    // frame
      }      // nthd dist
      free(kernel);
    }  // row
    
    ROMP_PFLB_end
  }    // col
  ROMP_PF_end

  return (out);
}

int MB2Dfree(MB2D **pmb)
{
  MB2D *mb = *pmb;

  if (mb->d0) MRIfree(&mb->d0);
  if (mb->theta) MRIfree(&mb->theta);
  if (mb->fwhm) MRIfree(&mb->fwhm);
  if (mb->dmin) MRIfree(&mb->dmin);
  if (mb->nd) MRIfree(&mb->nd);
  free(mb);
  pmb = NULL;
  return (0);
}

MB2D *MB2Dcopy(MB2D *src, int CopyMRI, MB2D *copy)
{
  if (copy == NULL) copy = (MB2D *)calloc(sizeof(MB2D), 1);
  copy->type = src->type;
  copy->offset = src->offset;
  copy->slope = src->slope;
  copy->DeltaD = src->DeltaD;
  copy->c0 = src->c0;
  copy->r0 = src->r0;
  copy->cutoff = src->cutoff;
  copy->Interp = src->Interp;
  copy->cR = src->cR;
  copy->rR = src->rR;

  if (CopyMRI) {
    if (copy->d0) MRIfree(&copy->d0);
    if (copy->theta) MRIfree(&copy->theta);
    if (copy->fwhm) MRIfree(&copy->fwhm);
    if (copy->dmin) MRIfree(&copy->dmin);
    if (copy->nd) MRIfree(&copy->nd);
    copy->d0 = MRIcopy(src->d0, NULL);
    copy->theta = MRIcopy(src->theta, NULL);
    copy->fwhm = MRIcopy(src->fwhm, NULL);
    copy->dmin = MRIcopy(src->dmin, NULL);
    copy->nd = MRIcopy(src->nd, NULL);
  }
  return (copy);
}

/*
  \fn MRI *MB2Dgrid(MRI *mri_template, int skip, MRI *outvol)
  \brief Creates a grid of ones separated by skip. All the
  slices are the same. Good for testing MRImotionBlur2D().
*/
MRI *MB2Dgrid(MRI *mri_template, int skip, MRI *outvol)
{
  int c, r, s, c0, r0, dc, dr, f;
  double radmax;

  if (outvol == NULL) {
    outvol = MRIcloneBySpace(mri_template, MRI_FLOAT, 1);
    if (outvol == NULL) return (NULL);
  }
  // Make sure it is 0
  MRIconst(outvol->width, outvol->height, outvol->depth, 1, 0, outvol);
  c0 = nint(mri_template->width / 2.0);
  r0 = nint(mri_template->height / 2.0);
  radmax = sqrt(2) * mri_template->width / 2;
  for (dc = -radmax; dc < radmax; dc += skip) {
    c = nint(c0 + dc);
    if (c < 0) continue;
    if (c >= mri_template->width) continue;
    for (dr = -radmax; dr < radmax; dr += skip) {
      r = r0 + dr;
      if (r < 0) continue;
      if (r >= mri_template->width) continue;
      for (s = 0; s < outvol->depth; s++)
        for (f = 0; f < outvol->nframes; f++) MRIsetVoxVal(outvol, c, r, s, f, 1);
    }  // r
  }    // c
  return (outvol);
}


