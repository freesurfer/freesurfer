/*
 *       FILE NAME:   mri.c
 *
 *       DESCRIPTION: utilities for MRI  data structure
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        1/8/97
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>

#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "filter.h"
#include "box.h"
#include "region.h"
#include "nr.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*#define DEBUG_POINT(x,y,z)  (((x==8&&y==9) || (x==9&&y==8)) &&((z)==15))*/
/*#define DEBUG_POINT(x,y,z)  (((x==7&&y==9) || (x==9&&y==7)) &&((z)==15))*/
#define DEBUG_POINT(x,y,z)  (((x==21) && (y==14)) &&((z)==7))

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int compare_sort_array(const void *pc1, const void *pc2) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define NORM  1.0f

MRI *
MRIdirectionMap(MRI *mri_grad, MRI *mri_dst, int wsize)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0 ;
  float   dir, dot, *xpix, *ypix, *zpix, dx, dy, dz, ox, oy, oz, *dir_pix, 
          norm ;

  whalf = wsize/2 ;
  norm = (float)(wsize*wsize*wsize)*NORM ;
  norm *= norm ;
  width = mri_grad->width ;
  height = mri_grad->height ;
  depth = mri_grad->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_grad, mri_dst) ;
  }

  for (z = whalf ; z < depth-whalf ; z++)
  {
    for (y = whalf ; y < height-whalf ; y++)
    {
      dir_pix = &MRIFvox(mri_dst, whalf, y, z) ;
      for (x = whalf ; x < width-whalf ; x++)
      {
        ox = MRIFvox(mri_grad, x, y, z) ;
        oy = MRIFseq_vox(mri_grad, x, y, z, 1) ;
        oz = MRIFseq_vox(mri_grad, x, y, z, 2) ;
        for (dir = 0.0f, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            xpix = &MRIFvox(mri_grad, x-whalf, y+y0, z+z0) ;
            ypix = &MRIFseq_vox(mri_grad, x-whalf, y+y0, z+z0, 1) ;
            zpix = &MRIFseq_vox(mri_grad, x-whalf, y+y0, z+z0, 2) ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              dx = *xpix++ ;
              dy = *ypix++ ;
              dz = *zpix++ ;
              dot = dx*ox + dy*oy + dz*oz ;
              dir += (x0*ox + y0*oy + z0*oz) * dot ;
            }
          }
        }
        *dir_pix++ = tanh(dir/norm) ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define SCALE_FACTOR 10.0f
MRI *
MRIoffsetDirection(MRI *mri_grad, int wsize, MRI *mri_direction, MRI *mri_dir)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0 ;
  float   dir, dot, *xpix, *ypix, *zpix, dx, dy, dz, ox, oy, oz, len,
          *dir_xpix, *dir_ypix, *dir_zpix ;
  BUFTYPE *pdir ;

  whalf = wsize/2 ;
  width = mri_grad->width ;
  height = mri_grad->height ;
  depth = mri_grad->depth ;

  if (!mri_direction)
  {
    mri_direction = MRIallocSequence(width, height, depth, MRI_FLOAT, 3) ;
    MRIcopyHeader(mri_grad, mri_direction) ;
  }

  for (z = whalf ; z < depth-whalf ; z++)
  {
    for (y = whalf ; y < height-whalf ; y++)
    {
      if (mri_dir)
        pdir = &MRIvox(mri_dir, whalf, y, z) ;
      else
        pdir = NULL ;

      dir_xpix = &MRIFvox(mri_direction, whalf, y, z) ;
      dir_xpix = &MRIFvox(mri_direction, whalf, y, z) ;
      dir_ypix = &MRIFseq_vox(mri_direction, whalf, y, z, 1) ;
      dir_zpix = &MRIFseq_vox(mri_direction, whalf, y, z, 2) ;
      for (x = whalf ; x < width-whalf ; x++)
      {
#define DEBUG_OFFSET 0
#if DEBUG_OFFSET
        FILE *fp ;
#endif
        ox = MRIFvox(mri_grad, x, y, z) ;
        oy = MRIFseq_vox(mri_grad, x, y, z, 1) ;
        oz = MRIFseq_vox(mri_grad, x, y, z, 2) ;
#if DEBUG_OFFSET
        if (DEBUG_POINT(x,y,z))
        {
          fp = fopen("dir.log", "w") ;
          fprintf(fp, "direction (%d, %d, %d), O: (%2.1f, %2.1f, %2.1f)\n",
                  x, y, z, ox, oy, oz) ;
          DiagBreak() ;
        }
#endif
        for (dir = 0.0f, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            xpix = &MRIFvox(mri_grad, x-whalf, y+y0, z+z0) ;
            ypix = &MRIFseq_vox(mri_grad, x-whalf, y+y0, z+z0, 1) ;
            zpix = &MRIFseq_vox(mri_grad, x-whalf, y+y0, z+z0, 2) ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              dx = *xpix++ ;
              dy = *ypix++ ;
              dz = *zpix++ ;
              dot = dx*ox + dy*oy + dz*oz ;
#if 1
              if (dot < 0.0f)
                dot = 0.0f ;
#endif
              dir += (x0*ox + y0*oy + z0*oz) * dot ;
#if DEBUG_OFFSET
if (DEBUG_POINT(x,y,z))
  fprintf(fp, 
    "(%d, %d, %d): (%2.1f, %2.1f, %2.1f), dot = %2.2f, prod = %2.1f (%2.1f)\n",
    x+x0, y+y0, z+z0, dx, dy, dz, dot, (x0*ox + y0*oy + z0*oz), dir) ;
#endif
            }
          }
        }
#if DEBUG_OFFSET
if (DEBUG_POINT(x,y,z))
    fclose(fp) ;
#endif
        if (ISSMALL(dir))
          ox = oy = oz = 0.0f ;
        else if (dir > 0.0f)
        {
          ox = -ox ;
          oy = -oy ;
          oz = -oz ;
        }
        len = sqrt(ox*ox + oy*oy + oz*oz) ;
        if (FZERO(len))
        {
          *dir_xpix++ = 0.0f ;
          *dir_ypix++ = 0.0f ;
          *dir_zpix++ = 0.0f ;
        }
        else
        {
          *dir_xpix++ = ox ;
          *dir_ypix++ = oy ;
          *dir_zpix++ = oz ;
        }

        if (pdir)
        {
          if (ISSMALL(dir))
            *pdir++ = OFFSET_ZERO ;
          else if (dir < 0.0f)
            *pdir++ = OFFSET_GRADIENT_DIRECTION ;
          else
            *pdir++ = OFFSET_NEGATIVE_GRADIENT_DIRECTION ;
        }
      }
    }
  }

  return(mri_direction) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIoffsetMagnitude(MRI *mri_src, MRI *mri_dst, int maxsteps)
{
  int          width, height, depth, x, y, z, x1, y1, z1, steps ;
  float        odx, ody, odz, dx, dy, dz, *xpix, *ypix, *zpix, len, dot = 0.0f,
               adx, ady, adz ;
  signed char  *dst_xpix, *dst_ypix, *dst_zpix ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width, height, depth, MRI_UCHAR, 3) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      xpix = &MRIFvox(mri_src, 0, y, z) ;
      ypix = &MRIFseq_vox(mri_src, 0, y, z, 1) ;
      zpix = &MRIFseq_vox(mri_src, 0, y, z, 2) ;

      dst_xpix = &MRISCvox(mri_dst, 0, y, z) ;
      dst_ypix = &MRISCseq_vox(mri_dst, 0, y, z, 1) ;
      dst_zpix = &MRISCseq_vox(mri_dst, 0, y, z, 2) ;
      for (x = 0 ; x < width ; x++)
      {
        odx = *xpix++ ; ody = *ypix++ ; odz = *zpix++ ;
        adx = fabs(odx) ; ady = fabs(ody) ; adz = fabs(odz) ;

        /* scale so that offset reaches next voxel at each step */
        if ((adx > ady) && (adx > adz))  /* scale using x component */
          len = adx ;
        else if (ady > adz)              /* scale using y component */
          len = ady ;
        else                             /* scale using z component */
          len = adz ;
        if (!FZERO(len))   /* a non-zero offset */
        {
          odx /= len ; ody /= len ; odz /= len ;
          
          for (steps = 0 ; steps < maxsteps ; steps++)
          {
            x1 = x + nint((float)steps * odx) ;
            y1 = y + nint((float)steps * ody) ;
            z1 = z + nint((float)steps * odz) ;
            if (x1 < 0 || x1 >= width || y1 < 0 || y1 >= height ||
                z1 < 0 || z1 >= depth)
              break ;
            dx = MRIFvox(mri_src, x1, y1, z1) ;
            dy = MRIFseq_vox(mri_src, x1, y1, z1, 1) ;
            dz = MRIFseq_vox(mri_src, x1, y1, z1, 2) ;
            dot = dx * odx + dy * ody + dz * odz ;
            if (dot <= 0.0f)
              break ;
          }
          if (!FZERO(dot))
            steps-- ;
        }
        else
          steps = 0 ;
        x1 = x + nint((float)steps * odx) ;
        y1 = y + nint((float)steps * ody) ;
        z1 = z + nint((float)steps * odz) ;

        *dst_xpix++ = (x1 - x) ;
        *dst_ypix++ = (y1 - y) ;
        *dst_zpix++ = (z1 - z) ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIapplyOffset(MRI *mri_src, MRI *mri_dst, MRI *mri_offset)
{
  int     width, height, depth, x, y, z, dx, dy, dz ;
  signed char    *pdx, *pdy, *pdz ;
  char    *pdst ;
  float   *fdst ;

  fdst = NULL ;
  pdst = NULL ;   /* for compiler warning - I know it's silly */
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (mri_dst->type != mri_src->type)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIapplyOffset: source type %d != dst %d",
                 mri_src->type, mri_dst->type)) ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdx = &MRISCvox(mri_offset, 0, y, z) ;
      pdy = &MRISCseq_vox(mri_offset, 0, y, z, 1) ;
      pdz = &MRISCseq_vox(mri_offset, 0, y, z, 2) ;
      switch (mri_dst->type)
      {
      case MRI_UCHAR:
        pdst = &MRIvox(mri_dst, 0, y, z) ;
        break ;
      case MRI_FLOAT:
        fdst = &MRIFvox(mri_dst, 0, y, z) ;
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_UNSUPPORTED, "MRIapplyOffset: unsupported type %d",
                     mri_dst->type)) ;
        break ;
      }
      for (x = 0 ; x < width ; x++)
      {
        dx = *pdx++ ;
        dy = *pdy++ ;
        dz = *pdz++ ;

        switch (mri_src->type)
        {
        case MRI_UCHAR:
          *pdst++ = MRIvox(mri_src,x+dx,y+dy,z+dz)  ;
          break ;
        case MRI_FLOAT:
          *fdst++ = MRIFvox(mri_src,x+dx,y+dy,z+dz)  ;
          break ;
        }
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIdiffuse(MRI *mri_src, MRI *mri_dst, double k, 
                    int niter, int which, double slope)
{
  switch (which)
  {
  case FILTER_DIFFUSE_GRAD:
    mri_dst = MRIdiffusePerona(mri_src, mri_dst, k, niter, slope) ;
    break ;
  case FILTER_DIFFUSE_CURV:
    mri_dst = MRIdiffuseCurvature(mri_src, mri_dst, k, niter, slope) ;
    break ;
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIdiffuseCurvature(MRI *mri_src, MRI *mri_dst, 
                            double A,int niter, double slope)
{
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define C(grad,k)    ((float)exp(-0.5f * SQR((float)fabs((grad))/k)))

#define NUM_NEIGH    6
#define KERNEL_MUL   (1.0f / ((float)NUM_NEIGH))

MRI *
MRIdiffusePerona(MRI *mri_src, MRI *mri_dst, 
                             double k, int niter,double slope)
{
  int     x, y, z, i, width, height, depth, *xE, *yN, *xW, *yS, *zD, *zU,
          ys, yn, xe, xw, zu, zd, ci ;
  float   c[NUM_NEIGH+1], fvals[NUM_NEIGH+1], dst_val ;
  MRI     *mri_tmp = NULL, *mri_mag = NULL, *mri_grad = NULL ; ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
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

  if (!MRIcheckSize(mri_src, mri_mag, 0, 0, 0))
  {
    if (mri_mag)
      MRIfree(&mri_mag) ;
    if (mri_tmp)
      MRIfree(&mri_tmp) ;
    mri_mag = MRIalloc(width, height, depth, MRI_FLOAT) ;
    mri_tmp = MRIalloc(width, height, depth, MRI_FLOAT) ;
  }

  /* build index tables for border */
  xE = (int *)calloc((unsigned int)width, sizeof(int)) ;
  xW = (int *)calloc((unsigned int)width, sizeof(int)) ;
  yN = (int *)calloc((unsigned int)height, sizeof(int)) ;
  yS = (int *)calloc((unsigned int)height, sizeof(int)) ;
  zU = (int *)calloc((unsigned int)depth, sizeof(int)) ;
  zD = (int *)calloc((unsigned int)depth, sizeof(int)) ;

  xW[0] = 0 ;
  for (x = 1 ; x < width ; x++)
  {
    xW[x] = x - 1 ;
    xE[x-1] = x ;
  }
  xE[width-1] = width-1 ;
  yN[0] = 0 ;
  for (y = 1 ; y < height ; y++)
  {
    yN[y] = y - 1 ;
    yS[y-1] = y ;
  }
  yS[height-1] = height-1 ;
  zD[0] = 0 ;   /* up is positive z */
  for (z = 1 ; z < depth ; z++)
  {
    zD[z] = z - 1 ;
    zU[z-1] = z ;
  }
  zU[depth-1] = depth-1 ;

  MRIcopy(mri_src, mri_tmp) ;

  memset(fvals, 0, NUM_NEIGH*sizeof(float)) ;
  memset(c, 0, NUM_NEIGH*sizeof(float)) ;
  for (i = 0 ; i < niter ; i++)
  {
    mri_grad = MRIsobel(mri_tmp, mri_grad, mri_mag) ;

    for (x = 0 ; x < width ; x++)
    {
      xe = xE[x] ;
      xw = xW[x] ;
      for (y = 0 ; y < height ; y++)
      {
        yn = yN[y] ;
        ys = yS[y] ;
          
        for (z = 0 ; z < depth ; z++)
        {
          zu = zU[z] ;
          zd = zD[z] ;
/*
     |  C1  | 
  ---------------
  C2 |  C0  | C3
  ---------------
     |  C4  | 
C6 (zU) is into the screen and c5 (zD) is out of the screen
*/

          fvals[0] = MRIFvox(mri_tmp, x, y, z) ;
          fvals[1] = MRIFvox(mri_tmp, x, yn, z) ;
          fvals[2] = MRIFvox(mri_tmp, xw, y, z) ;
          fvals[3] = MRIFvox(mri_tmp, xe, y, z) ;
          fvals[4] = MRIFvox(mri_tmp, x, ys, z) ;
          fvals[5] = MRIFvox(mri_tmp, x, y, zd) ;
          fvals[6] = MRIFvox(mri_tmp, x, y, zu) ;
          
          c[1] = KERNEL_MUL * C(MRIFvox(mri_mag, x, yn,z), k) ;
          c[2] = KERNEL_MUL * C(MRIFvox(mri_mag, xw, y,z), k) ;
          c[3] = KERNEL_MUL * C(MRIFvox(mri_mag, xe, y,z), k) ;
          c[4] = KERNEL_MUL * C(MRIFvox(mri_mag, x, ys,z), k) ;
          c[5] = KERNEL_MUL * C(MRIFvox(mri_mag, x, y, zd), k) ;
          c[6] = KERNEL_MUL * C(MRIFvox(mri_mag, x, y, zu), k) ;
          
          c[0] = 1.0f ;
          for (ci = 1 ; ci <= NUM_NEIGH ; ci++)
            c[0] -= c[ci] ;
          
          for (dst_val = 0.0f, ci = 0 ; ci <= NUM_NEIGH ; ci++)
            dst_val += fvals[ci] * c[ci] ;
          MRIFvox(mri_dst, x, y, z) = dst_val ;
        }
      }
    }
    
    MRIcopy(mri_dst, mri_tmp) ;
    k = k + k * slope ;
  }

  free(xE) ;
  free(xW) ;
  free(yN) ;
  free(yS) ;
  free(zU) ;
  free(zD) ;


  return(mri_dst) ;
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
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

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
  register BUFTYPE *tr_pix, *mr_pix, *br_pix ;
  int              x, y, z, width, height, depth, x0, y0, z0 ;
  float            *outPtr, left, middle, right ;

  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

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

/* applying sobel in x-y plane, don't worry about boundary conditions in z  */
  /* don't apply sobel to outer ring to pixels to avoid border effects */
  width-- ; height-- ; 
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
  register BUFTYPE *tr_pix, *br_pix ;
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

/* applying sobel in x-y plane, don't worry about boundary conditions in z  */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ; height-- ; 
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
  register BUFTYPE *tr_pix, *br_pix ;
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

/* applying sobel in x-z plane, don't worry about boundary conditions in y  */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ; depth-- ;
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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI *
MRIsobel(MRI *mri_src, MRI *mri_grad, MRI *mri_mag)
{
  int      x, y, z, width, height, depth  ;
  float    *px, *py, *pz, xval,yval,zval, *pmag ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_grad)
  {
    mri_grad = MRIallocSequence(width, height, depth, MRI_FLOAT, 3) ;
    MRIcopyHeader(mri_src, mri_grad) ;
  }

  MRIxSobel(mri_src, mri_grad, 0) ;
  MRIySobel(mri_src, mri_grad, 1) ;
  MRIzSobel(mri_src, mri_grad, 2) ;

  if (mri_mag)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        px = &MRIFvox(mri_grad, 0, y, z) ;
        py = &MRIFseq_vox(mri_grad, 0, y, z, 1) ;
        pz = &MRIFseq_vox(mri_grad, 0, y, z, 2) ;
        pmag = &MRIFvox(mri_mag, 0, y, z) ;
        
        for (x = 0 ; x < width ; x++)
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
MRIxSobel(MRI *mri_src, MRI *mri_x, int frame)
{
  register BUFTYPE *tr_pix, *mr_pix, *br_pix ;
  int              x, y, z, width, height, depth ;
  float            *outPtr, left, middle, right ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_x)
  {
    mri_x = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_x) ;
  }

/* applying sobel in x-y plane, don't worry about boundary conditions in z  */
  /* don't apply sobel to outer ring to pixels to avoid border effects */
  width-- ; height-- ; 
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
MRIySobel(MRI *mri_src, MRI *mri_y, int frame)
{
  register BUFTYPE *tr_pix, *br_pix ;
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

/* applying sobel in x-y plane, don't worry about boundary conditions in z  */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ; height-- ; 
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 1 ; y < height ; y++)
    {
      tr_pix = &MRIvox(mri_src, 0, y-1, z) ;
      br_pix = &MRIvox(mri_src, 0, y+1, z) ;
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
MRIzSobel(MRI *mri_src, MRI *mri_z, int frame)
{
  register BUFTYPE *tr_pix, *br_pix ;
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

/* applying sobel in x-z plane, don't worry about boundary conditions in y  */
  /* don't apply sobel to outer ring of pixels to avoid border effects */
  width-- ; depth-- ;
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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          compute a zscore for an image
------------------------------------------------------*/
MRI *
MRIzScore(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_std)
{
  MRI_REGION region ;

  region.x = region.y = region.z = 0 ;
  region.dx = mri_src->width ;
  region.dx = mri_src->width ;
  region.dx = mri_src->width ;
  return(MRIzScoreRegion(mri_src, mri_dst, mri_mean, mri_std, &region)) ;
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
MRI *
MRIzScoreRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_std,
                MRI_REGION *region)
{
  int     width, height, depth, x, y, z, x0, y0, z0 ;
  BUFTYPE *psrc ;
  float   *pmean, *pstd, std, *pdst ;

  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  if (!mri_dst)
  {
    int  w, h, d ;
    
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      pmean = &MRIFvox(mri_mean, 0, y-y0, z-z0) ;
      pstd = &MRIFvox(mri_std,   0, y-y0, z-z0) ;
      pdst = &MRIFvox(mri_dst,   0, y-y0, z-z0) ;
      psrc = &MRIvox(mri_src, x0, y, z) ;
      for (x = x0 ; x < width ; x++)
      {
        std = *pstd++ ;
        if (!std)
        {
          *pdst++ = 0.0f ;
          psrc++ ;
          pmean++ ;
        }
        else
          *pdst++ = ((float)*psrc++ - *pmean++) / std ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *
MRImeanRegion(MRI *mri_src, MRI *mri_dst, int wsize,MRI_REGION *region)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0, val,
          xk, yk, zk, xi, yi, zi ;
  float   *pdst, wcubed ;

  wcubed = (float)(wsize*wsize*wsize) ;
  whalf = wsize/2 ;
  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  if (!mri_dst)
  {
    int  w, h, d ;
    
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      pdst = &MRIFvox(mri_dst, 0, y-y0, z-z0) ;
      for (x = x0 ; x < width ; x++)
      {
        for (val = 0, zk = -whalf ; zk <= whalf ; zk++)
        {
          zi = mri_src->zi[z+zk] ;
          for (yk = -whalf ; yk <= whalf ; yk++)
          {
            yi = mri_src->yi[y+yk] ;
            for (xk = -whalf ; xk <= whalf ; xk++)
            {
              xi = mri_src->xi[x+xk] ;
              val += (int)MRIvox(mri_src, xi, yi, zi) ;
            }
          }
        }
        *pdst++ = (float)val / wcubed ;
      }
    }
  }
  return(mri_dst) ;
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
MRI *
MRIstdRegion(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize,
             MRI_REGION *region)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0,
          xk, yk, zk, xi, yi, zi ;
  float   *pdst, wcubed, mean, *pmean, variance, f ;

  wcubed = (float)(wsize*wsize*wsize) ;
  whalf = wsize/2 ;
  width = region->x + region->dx ;
  if (width > mri_src->width)
    width = mri_src->width ;
  height = region->y + region->dy ;
  if (height > mri_src->height)
    height = mri_src->height ;
  depth = region->z + region->dz ;
  if (depth > mri_src->depth)
    depth = mri_src->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  if (!mri_dst)
  {
    int  w, h, d ;
    
    w = width - region->x ;
    h = height - region->y ;
    d = depth - region->z ;
    mri_dst = MRIalloc(w, h, d, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = mri_src->xstart + region->x * mri_src->xsize ;
    mri_dst->ystart = mri_src->ystart + region->y * mri_src->ysize ;
    mri_dst->zstart = mri_src->zstart + region->z * mri_src->zsize ;
    mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
    mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
    mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  }

  for (z = z0 ; z < depth ; z++)
  {
    for (y = y0 ; y < height ; y++)
    {
      pdst = &MRIFvox(mri_dst, 0, y-y0, z-z0) ;
      pmean = &MRIFvox(mri_mean, 0, y-y0, z-z0) ;
      for (x = x0 ; x < width ; x++)
      {
        mean = *pmean++ ;
        for (variance = 0.0f, zk = -whalf ; zk <= whalf ; zk++)
        {
          zi = mri_src->zi[z+zk] ;
          for (yk = -whalf ; yk <= whalf ; yk++)
          {
            yi = mri_src->yi[y+yk] ;
            for (xk = -whalf ; xk <= whalf ; xk++)
            {
              xi = mri_src->xi[x+xk] ;
              f = (float)MRIvox(mri_src, xi, yi, zi) - mean ;
              variance += (f * f) ;
            }
          }
        }
        *pdst++ = sqrt(variance / wcubed) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *
MRImeanByte(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0, val ;
  BUFTYPE *psrc, *pdst ;
  float   wcubed ;
  

  wcubed = (float)(wsize*wsize*wsize) ;
  whalf = wsize/2 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRImeanByte: dst must be MRI_UCHAR")) ;

  for (z = whalf ; z < depth-whalf ; z++)
  {
    for (y = whalf ; y < height-whalf ; y++)
    {
      pdst = &MRIvox(mri_dst, whalf, y, z) ;
      for (x = whalf ; x < width-whalf ; x++)
      {
        for (val = 0, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            psrc = &MRIvox(mri_src, x-whalf, y+y0, z+z0) ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              val += (int)*psrc++ ;
            }
          }
        }
        *pdst++ = (BUFTYPE)nint((float)val / wcubed) ;
      }
    }
  }

  /* now copy information to borders from source image */
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < whalf ; z++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
      for (z = depth-whalf ; z < depth ; z++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
    }
  }
  for (x = 0 ; x < width ; x++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < whalf ; y++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
      for (y = height-whalf ; y < height ; y++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
    }
  }
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < whalf ; x++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
      for (x = width-whalf ; x < width ; x++)
        MRIvox(mri_dst, x, y, z) = MRIvox(mri_src,x,y,z) ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a mean filter on the input MRI
------------------------------------------------------*/
MRI *
MRImean(MRI *mri_src, MRI *mri_dst, int wsize)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0, val ;
  BUFTYPE *psrc ;
  float   *pdst, wcubed ;

  wcubed = (float)(wsize*wsize*wsize) ;
  whalf = wsize/2 ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (mri_dst->type != MRI_FLOAT)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRImean: dst must be MRI_FLOAT")) ;

  for (z = whalf ; z < depth-whalf ; z++)
  {
    for (y = whalf ; y < height-whalf ; y++)
    {
      pdst = &MRIFvox(mri_dst, whalf, y, z) ;
      for (x = whalf ; x < width-whalf ; x++)
      {
        for (val = 0, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            psrc = &MRIvox(mri_src, x-whalf, y+y0, z+z0) ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              val += (int)*psrc++ ;
            }
          }
        }
        *pdst++ = (float)val / wcubed ;
      }
    }
  }

  /* now copy information to borders from source image */
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < whalf ; z++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
      for (z = depth-whalf ; z < depth ; z++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
    }
  }
  for (x = 0 ; x < width ; x++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < whalf ; y++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
      for (y = height-whalf ; y < height ; y++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
    }
  }
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < whalf ; x++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
      for (x = width-whalf ; x < width ; x++)
        MRIFvox(mri_dst, x, y, z) = (float)MRIvox(mri_src,x,y,z) ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a median filter on the input MRI
------------------------------------------------------*/
MRI *
MRImedian(MRI *mri_src, MRI *mri_dst, int wsize)
{
  static BUFTYPE *sort_array = NULL ;
  static int   sort_size = 0 ;
  int     width, height, depth, x, y, z, whalf, x0, y0, z0,
          median_index, wcubed, yi, zi ;
  BUFTYPE *sptr, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRImedian: dst must be MRI_UCHAR")) ;

  wcubed = (wsize*wsize*wsize) ;
  median_index = wcubed / 2 ;
  whalf = wsize/2 ;

  if (sort_array && (wcubed != sort_size))
  {
    free(sort_array) ;
    sort_array = NULL ;
  }
  if (!sort_array)
  {
    sort_array = (BUFTYPE *)calloc(wcubed, sizeof(BUFTYPE)) ;
    sort_size = wcubed ;
  }
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        for (sptr = sort_array, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          zi = mri_src->zi[z+z0] ;
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            yi = mri_src->yi[y+y0] ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
              *sptr++ = MRIvox(mri_src, mri_src->xi[x+x0], yi, zi) ;
          }
        }

        qsort(sort_array, wcubed, sizeof(BUFTYPE), compare_sort_array) ;
        *pdst++ = sort_array[median_index] ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform a generalized median filter on the input MRI
------------------------------------------------------*/
MRI *
MRIorder(MRI *mri_src, MRI *mri_dst, int wsize, float pct)
{
  static BUFTYPE *sort_array = NULL ;
  static int   sort_size = 0 ;
  int     width, height, depth, x, y, z, whalf, x0, y0, z0,
          order_index, wcubed, yi, zi ;
  BUFTYPE *sptr, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (mri_dst->type != MRI_UCHAR)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRIorder: dst must be MRI_UCHAR")) ;

  wcubed = (wsize*wsize*wsize) ;
  order_index = pct*wcubed ;
  whalf = wsize/2 ;

  if (sort_array && (wcubed != sort_size))
  {
    free(sort_array) ;
    sort_array = NULL ;
  }
  if (!sort_array)
  {
    sort_array = (BUFTYPE *)calloc(wcubed, sizeof(BUFTYPE)) ;
    sort_size = wcubed ;
  }
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        for (sptr = sort_array, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          zi = mri_src->zi[z+z0] ;
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            yi = mri_src->yi[y+y0] ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
              *sptr++ = MRIvox(mri_src, mri_src->xi[x+x0], yi, zi) ;
          }
        }

        qsort(sort_array, wcubed, sizeof(BUFTYPE), compare_sort_array) ;
        *pdst++ = sort_array[order_index] ;
      }
    }
  }
  return(mri_dst) ;
}
static int
compare_sort_array(const void *pc1, const void *pc2)
{
  register BUFTYPE c1, c2 ;

  c1 = *(BUFTYPE *)pc1 ;
  c2 = *(BUFTYPE *)pc2 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2)
    return(1) ;
  else if (c1 < c2)
    return(-1) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a splined gaussian bump which tails to 0. 
             Returns an image which is (8*sigma)+1 
             (Nitzberg and Shiota, 1993)
----------------------------------------------------------------------*/
MRI *
MRIgaussian1d(float sigma, int max_len)
{
  MRI       *mri ;
  float     norm, two_sigma, fx, k ;
  int       x, half, len ;

  /* build the kernel in k */
  len = (int)nint(8.0f * sigma)+1 ;
  if (ISEVEN(len))   /* ensure it's even */
    len++ ;
  if (max_len && (max_len < len))
    len = max_len ;
  half = len/2 ;
  mri = MRIalloc(len, 1, 1, MRI_FLOAT) ;

  norm = 0.0f ;
  two_sigma = 2.0f * sigma ;

  for (norm = 0.0f, x = 0 ; x < len ; x++)
  {
    fx = (float)(x-half) ;
    if (fabs(fx) <= two_sigma)
      k = (float)exp((double)(-fx*fx/(two_sigma*sigma))) ;
    else if (two_sigma < fabs(fx) && fabs(fx) <= 4.0f*sigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * 
        (float)pow(4.0f - fabs(fx)/(double)sigma, 4.0) ;
    else
      k = 0 ;

    MRIFvox(mri, x, 0, 0) = k ;
    norm += k ;
  }

  /* normalize kernel to sum to 1 */
  for (x = 0 ; x < len ; x++)
    MRIFvox(mri, x, 0, 0) /= norm ;

  return(mri) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           calculate the standard deviation of the MRI
           in a wsize window around each point
------------------------------------------------------*/
MRI *
MRIstd(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize)
{
  int     width, height, depth, x, y, z, whalf, x0, y0, z0 ;
  BUFTYPE *psrc ;
  float   *pdst, wcubed, mean, *pmean, variance, f ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (mri_dst->type != MRI_FLOAT)
    ErrorReturn(mri_dst, 
                (ERROR_UNSUPPORTED, "MRIstd: dst must be MRI_FLOAT")) ;

  wcubed = (float)(wsize*wsize*wsize) ;
  whalf = wsize/2 ;

  for (z = whalf ; z < depth-whalf ; z++)
  {
    for (y = whalf ; y < height-whalf ; y++)
    {
      pdst = &MRIFvox(mri_dst, whalf, y, z) ;
      pmean = &MRIFvox(mri_mean, whalf, y, z) ;
      for (x = whalf ; x < width-whalf ; x++)
      {
        mean = *pmean++ ;
        for (variance = 0.0f, z0 = -whalf ; z0 <= whalf ; z0++)
        {
          for (y0 = -whalf ; y0 <= whalf ; y0++)
          {
            psrc = &MRIvox(mri_src, x-whalf, y+y0, z+z0) ;
            for (x0 = -whalf ; x0 <= whalf ; x0++)
            {
              f = (float)*psrc++ - mean ;
              variance += (f * f) ;
            }
          }
        }
        *pdst++ = variance / wcubed ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIconvolveGaussian(MRI *mri_src, MRI *mri_dst, MRI *mri_gaussian)
{
  int  width, height, depth, klen ;
  MRI  *mtmp1, *mtmp2 ;
  float *kernel ;

  kernel = &MRIFvox(mri_gaussian, 0, 0, 0) ;
  klen = mri_gaussian->width ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (width <= 1 || height <= 1 || depth <= 1)
    ErrorExit(ERROR_BADPARM, 
              "MRIconvolveGaussian: insufficient dimension (%d, %d, %d)",
              width, height, depth) ;

  mtmp1 = MRIconvolve1d(mri_src, NULL, kernel, klen, MRI_WIDTH) ;
  mtmp2 = MRIconvolve1d(mtmp1, NULL, kernel, klen, MRI_HEIGHT) ;
  MRIconvolve1d(mtmp2, mtmp1, kernel, klen, MRI_DEPTH) ;
  MRIfree(&mtmp2) ;
  
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  MRIcopy(mtmp1, mri_dst) ;    /* convert it back to UCHAR */

  MRIcopyHeader(mri_src, mri_dst) ;
  
  MRIfree(&mtmp1) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define KERNEL_SIZE 5
#define K_A         0.4f
static float kernel[KERNEL_SIZE]  = 
{
  0.25f - K_A/2.0f, .25f, K_A, 0.25f, 0.25f-K_A/2.0f
} ;
MRI *
MRIreduce(MRI *mri_src, MRI *mri_dst)
{
  int  width, height, depth ;
  MRI  *mtmp1, *mtmp2 ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (width <= 1 || height <= 1 || depth <= 1)
    ErrorExit(ERROR_BADPARM, "MRIreduce: insufficient dimension (%d, %d, %d)",
              width, height, depth) ;

  mtmp1 = MRIconvolve1d(mri_src, NULL, kernel, KERNEL_SIZE, MRI_WIDTH) ;
  mtmp2 = MRIconvolve1d(mtmp1, NULL, kernel, KERNEL_SIZE, MRI_HEIGHT) ;
  mri_dst = MRIreduce1d(mtmp2, mri_dst, kernel, KERNEL_SIZE, MRI_DEPTH) ;

  MRIcopyHeader(mri_src, mri_dst) ;
  
  mri_dst->thick = mri_src->thick * 2.0f ;
  mri_dst->scale = mri_src->scale * 2 ;
  mri_dst->xsize = mri_src->xsize * 2.0f ;
  mri_dst->ysize = mri_src->ysize * 2.0f ;
  mri_dst->zsize = mri_src->zsize * 2.0f ;

  MRIfree(&mtmp1) ;
  MRIfree(&mtmp2) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis)
{
  int           x, y, z, width, height, halflen, depth, *xi, *yi, *zi ;
  register int  i ;
  BUFTYPE       *inBase ;
  float         *ki, total, *inBase_f, *outPix ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst = MRIalloc(width, height, depth, MRI_FLOAT) ;

  if (mri_dst->type != MRI_FLOAT)
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRIconvolve1d: unsupported dst pixel format %d",
                 mri_dst->type)) ;

  halflen = len/2 ;

  xi = mri_src->xi ; yi = mri_src->yi ; zi = mri_src->zi ;

  switch (mri_src->type)
  {
  case MRI_UCHAR:
    switch (axis)
    {
    case MRI_WIDTH:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          inBase = &MRIvox(mri_src, 0, y, z) ;
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;
            
            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * (float)(*(inBase + xi[x+i-halflen])) ;
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    case MRI_HEIGHT:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;

            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * (float)(MRIvox(mri_src, x,yi[y+i-halflen],z));
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    case MRI_DEPTH:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;
            
            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * (float)(MRIvox(mri_src, x,y,zi[z+i-halflen]));
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    }
    break ;
  case MRI_FLOAT:
    switch (axis)
    {
    case MRI_WIDTH:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          inBase_f = &MRIFvox(mri_src, 0, y, z) ;
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;
            
            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * (*(inBase_f + xi[x+i-halflen])) ;
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    case MRI_HEIGHT:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;

            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * MRIFvox(mri_src, x,yi[y+i-halflen],z);
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    case MRI_DEPTH:
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          outPix = &MRIFvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            total = 0.0f ;
            
            for (ki = k, i = 0 ; i < len ; i++)
              total += *ki++ * MRIFvox(mri_src, x,y,zi[z+i-halflen]);
            
            *outPix++ = total ;
          }
        }
      }
      break ;
    }
    break ;
  default:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, 
                 "MRIconvolve1d: unsupported pixel format %d",mri_src->type));
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
MRI *
MRIreduce1d(MRI *mri_src, MRI *mri_dst, float *k, int len, int axis)
{
  int    x, y, z, i, dheight, dwidth, ddepth, xi,yi, zi, halflen  ;
  int    sheight, swidth, sdepth ;
  float  total ;

  swidth = mri_src->width ;
  sheight = mri_src->height ;
  sdepth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(swidth/2, sheight/2, sdepth/2, MRI_UCHAR) ;
    mri_dst->imnr1 = mri_dst->depth-1 ;
  }

  if ((mri_dst->type != MRI_UCHAR) || (mri_src->type != MRI_FLOAT))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED, 
                 "MRIreduce1d: src %d or dst %d format unsupported",
                 mri_src->type, mri_dst->type)) ;

  dwidth = mri_dst->width ;
  dheight = mri_dst->height ;
  ddepth = mri_dst->depth ;


  halflen = (len-1)/2 ;
  switch (axis)
  {
  case MRI_WIDTH:
    for (z = 0 ; z < ddepth ; z++)
    {
      zi = 2*z ;
      for (y = 0 ; y < dheight ; y++)
      {
        yi = 2*y ;
        for (x = 0 ; x < dwidth ; x++)
        {
          total = 0.0f ;
          
          for (i = 0 ; i < len ; i++)
          {
            /* Neumann boundary conditions */
            xi = 2*x + i - halflen ;
            if (xi < 0)
              xi = 0 ;
            else if (xi >= swidth)
              xi = swidth - 1 ;
            
            total = total + k[i] * MRIFvox(mri_src, xi, yi, zi) ;
          }
          MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total) ;
        }
      }
    }
    break ;
  case MRI_HEIGHT:
    for (z = 0 ; z < ddepth ; z++)
    {
      zi = 2*z ;
      for (y = 0 ; y < dheight ; y++)
      {
        for (x = 0 ; x < dwidth ; x++)
        {
          total = 0.0f ;
          xi = 2*x ;
          
          for (i = 0 ; i < len ; i++)
          {
            /* Neumann boundary conditions */
            yi = 2*y + i - halflen ;
            if (yi < 0)
              yi = 0 ;
            else if (yi >= sheight)
              yi = sheight - 1 ;
            
            total = total + k[i] * MRIFvox(mri_src, xi, yi, zi) ;
          }
          MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total) ;
        }
      }
    }
    break ;
  case MRI_DEPTH:
    for (z = 0 ; z < ddepth ; z++)
    {
      zi = 2*z ;
      for (y = 0 ; y < dheight ; y++)
      {
        yi = 2*y ;
        for (x = 0 ; x < dwidth ; x++)
        {
          total = 0.0f ;
          xi = 2*x ;
          
          for (i = 0 ; i < len ; i++)
          {
            /* Neumann boundary conditions */
            zi = 2*z + i - halflen ;
            if (zi < 0)
              zi = 0 ;
            else if (zi >= sdepth)
              zi = sdepth - 1 ;
            
            total = total + k[i] * MRIFvox(mri_src, xi, yi, zi) ;
          }
          MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(total) ;
        }
      }
    }
    break ;
  }

/*  mri_dst = MRIcopy(mri_src, mri_dst) ;*/

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate two MRIs, limitting the extent of the 
          correlation to window_size
------------------------------------------------------*/
MRI *
MRIxcorrWindow(MRI *mri_ref, MRI *mri_in, MRI *mri_dst, int window_size)
{
  int     width, height, depth, xoff, yoff, zoff, x0, y0, z0, kx, ky, kz, 
          whalf, x, y, z, dst_width, dst_height, dst_depth, kxend ;
  BUFTYPE *pref, *pin, ref, in ;
  float   dst, *pdst ;

  if (!mri_dst)
  {
    /* decide on approriate size for destination */
    dst_width = window_size /* MAX(mri_ref->width, window_size)*/ ;
    dst_height = window_size /*MAX(mri_ref->height, window_size)*/ ;
    dst_depth = window_size /*MAX(mri_ref->depth, window_size)*/ ;

    mri_dst = MRIalloc(dst_width, dst_height, dst_depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_ref, mri_dst) ;
  }

  dst_width = mri_dst->width ;
  dst_height = mri_dst->height ;
  dst_depth = mri_dst->depth ;

  height = mri_ref->height ;
  width = mri_ref->width ;
  depth = mri_ref->depth ;

  if (window_size <= 0)
    window_size = width ;
  whalf = (window_size-1)/2 ;

  x0 = (dst_width-1) / 2 ;
  y0 = (dst_height-1) / 2 ;
  z0 = (dst_depth-1) / 2 ;

  /* 
     for each element of the destination array, correlate the reference
     MRI at that offset with the (unshifted) input MRI
     */
  for (zoff = -whalf ; zoff <= whalf ; zoff++)
  {
    for (yoff = -whalf ; yoff <= whalf ; yoff++)
    {
      pdst = &MRIFvox(mri_dst, x0-whalf, y0+yoff, z0+zoff) ;
      for (xoff = -whalf ; xoff <= whalf ; xoff++)
      {
        /* apply offset template (ref) to input */
        dst = 0.0f ;

        /* now calculate value of cross-correlation at these offsets */

        /* 
          x, y, and z are in 'input' coordinates, while kx, ky, kz are
          in 'reference' coordinates'.
        */
        for (z = 0, kz = zoff ; kz < depth+zoff ; kz++, z++)
        {
          if (kz < 0 || kz >= depth)
            continue ;
          for (y = 0, ky = yoff ; ky < height+yoff ; ky++, y++)
          {
            if (ky < 0 || ky >= height)
              continue ;

#if 1
            pref = &MRIvox(mri_ref,xoff,ky,kz) ;
            pin = &MRIvox(mri_in, 0, y, z) ;
#else
            pref = &mri_ref->slices[kz][ky][xoff] ;
            pin = mri_in->slices[z][y] ;
#endif
            kxend = MIN(width+xoff, width) ;
            for (x = 0, kx = xoff ; kx < kxend ; kx++, x++, pin++, pref++)
            {
              if (kx < 0 || kx >= width)
                continue ;

              in = *pin ;
              ref = *pref ;
              dst += (float)in * (float)ref ;
            }
          }
        }
        *pdst++ = dst ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate two MRIs, limitting the extent of the 
          correlation to window_size and normalizing by the energy
          in the window.
------------------------------------------------------*/
MRI *
MRInxcorrWindow(MRI *mri_ref, MRI *mri_in, MRI *mri_dst, int window_size)
{
  int     width, height, depth, xoff, yoff, zoff, x0, y0, z0, kx, ky, kz, kxend,
          whalf, x, y, z, dst_width, dst_height, dst_depth;
  BUFTYPE *pref, *pin, ref, in ;
  float   dst, norm, *pdst ;

  if (!mri_dst)
  {
    /* decide on approriate size for destination */
    dst_width = window_size /* MAX(mri_ref->width, window_size)*/ ;
    dst_height = window_size /*MAX(mri_ref->height, window_size)*/ ;
    dst_depth = window_size /*MAX(mri_ref->depth, window_size)*/ ;

    mri_dst = MRIalloc(dst_width, dst_height, dst_depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_ref, mri_dst) ;
  }

  dst_width = mri_dst->width ;
  dst_height = mri_dst->height ;
  dst_depth = mri_dst->depth ;

  height = mri_ref->height ;
  width = mri_ref->width ;
  depth = mri_ref->depth ;

  if (window_size <= 0)
    window_size = width ;
  whalf = (window_size-1)/2 ;

  x0 = (dst_width-1) / 2 ;
  y0 = (dst_height-1) / 2 ;
  z0 = (dst_depth-1) / 2 ;

  /* 
     for each element of the destination array, correlate the reference
     MRI at that offset with the (unshifted) input MRI
     */
  for (zoff = -whalf ; zoff <= whalf ; zoff++)
  {
    for (yoff = -whalf ; yoff <= whalf ; yoff++)
    {
      pdst = &MRIFvox(mri_dst, x0-whalf, y0+yoff, z0+zoff) ;
      for (xoff = -whalf ; xoff <= whalf ; xoff++)
      {
        /* apply offset template (ref) to input */
        dst = 0.0f ;

        /* now calculate value of cross-correlation at these offsets */

        /* 
          x, y, and z are in 'input' coordinates, while kx, ky, kz are
          in 'reference' coordinates'.
        */

        norm = 0.0f ;
        for (z = 0, kz = zoff ; kz < depth+zoff ; kz++, z++)
        {
          if (kz < 0 || kz >= depth)
            continue ;
          for (y = 0, ky = yoff ; ky < height+yoff ; ky++, y++)
          {
            if (ky < 0 || ky >= height)
              continue ;

            pref = &MRIvox(mri_ref,xoff,ky,kz) ;
            pin = &MRIvox(mri_in, 0, y, z) ;
            kxend = MIN(width+xoff, width) ;
            for (x = 0, kx = xoff ; kx < kxend ; kx++, x++, pin++, pref++)
            {
              if (kx < 0 || kx >= width)
                continue ;

              in = *pin ;
              ref = *pref ;

              norm += (float)(abs(in)*abs(ref)) ;
              dst += (float)in * (float)ref ;
            }
          }
        }
        norm = 1.0f ;    /* disable normalization (for now) */
        if (!FZERO(norm))
          dst /= norm ;
        else
          dst = 0.0f ;

        MRIFvox(mri_dst,x0+xoff,y0+yoff,z0+zoff) = dst ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Cross-correlate to MRIs
------------------------------------------------------*/
long
MRIcorrelate(MRI *mri_ref, MRI *mri_in, int xoff, int yoff, int zoff)
{
  int     row, col, slice, width, height, depth ;
  long    total ;
  BUFTYPE *pref, *pin, ref, in ;

  depth = mri_ref->depth ;
  height = mri_ref->height ;
  width = mri_ref->width ;

  total = 0L ;
  for (slice = 0 ; slice < depth ; slice++)
  {
    for (row = 0 ; row < height ; row++)
    {
      pref = mri_ref->slices[slice][row] ;
      pin = mri_in->slices[slice][row] ;
      for (col = 0 ; col < width ; col++)
      {
        in = *pin++ ;
        ref = *pref++ ;
        total += (long)in * (long)ref ;
      }
    }
  }

  return(total) ;
}
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
    mri_dst = MRIclone(mri_src, NULL) ;

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
                num_off++ ;
              else
              {
                if (!out_val)
                  out_val = 0 ;
                num_on++ ;
              }
            }
          }
        }
        val = MRIvox(mri_src, x, y, z) ;
        if (val && (num_off >= thresh))
          val = 0 ;
        else if (!val && (num_on >= thresh))
          val = 255 ;
        *pdst++ = val ;
      }
    }
  }
  return(mri_dst) ;
}
