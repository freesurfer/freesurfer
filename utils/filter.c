/*
 *       FILE NAME:   filter.c
 *
 *       DESCRIPTION: image processing filters
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        6/18/96
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>

#include <hipl_format.h>

#include "hips.h"

#include "image.h"
#include "error.h"
#include "utils.h"
#include "macros.h"
#include "machine.h"
#include "proto.h"
#include "diag.h"
#include "filter.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define NITSHI_TAU  2.0f

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

#define DEBUG_FILTER 0
IMAGE *Ifilter = NULL ;
#if DEBUG_FILTER
extern int Gx, Gy ;
#endif

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageNitShiFilter(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize, double sigma,
                  IMAGE *Idst)
{
  static IMAGE  *IE = NULL, *IF = NULL, *IG = NULL, *Iblur = NULL ;
  IMAGE  *Iin, *Iout ;
  int    rows, cols, x, y, whalf, xk, yk, xs, ys ;
  float  norm, total, *dpix, *spix, fmin, fmax,fval, Eval,Fval,Gval, sigma_sq;

#if DEBUG_FILTER
  if (Ifilter)
    ImageFree(&Ifilter) ;
  Ifilter = ImageAlloc(wsize, wsize, PFFLOAT, 1) ;
#endif

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (Isrc->pixel_format != PFFLOAT)
  {
    Iin = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    ImageCopy(Isrc, Iin) ;
    ImageValRange(Isrc, &fmin, &fmax) ;
  }
  else
    Iin = Isrc ;

  if (!Idst)
    Idst = ImageAlloc(rows, cols, PFFLOAT, 1) ;

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  else
    Iout = Idst ;

  if (!ImageCheckSize(Iin, IE, 0, 0, 0))
  {
    if (IE)
    {
      ImageFree(&IE) ;
      ImageFree(&IF) ;
      ImageFree(&IG) ;
    }
    IE = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    IF = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    IG = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  }

  if (!Iblur)
    Iblur = ImageGaussian1d(NITSHI_TAU, 0) ;
  whalf = (wsize-1)/2 ;

  /* calculate IE, IF, and IG */
  ImageMul(Ix, Ix, IE) ;
  ImageMul(Ix, Iy, IF) ;
  ImageMul(Iy, Iy, IG) ;
  ImageConvolveGaussian(IE, Iblur, IE, 0) ;
  ImageConvolveGaussian(IF, Iblur, IF, 0) ;
  ImageConvolveGaussian(IG, Iblur, IG, 0) ;
  
ImageWrite(Ix, "Ix.hipl") ;
ImageWrite(Iy, "Iy.hipl") ;
ImageWrite(IE, "IE.hipl");
ImageWrite(IF, "IF.hipl");
ImageWrite(IG, "IG.hipl");

  /* now apply actual filter */
  sigma_sq = 2.0f * (float)(sigma * sigma) ;
  dpix = IMAGEFpix(Iout, 0, 0) ;
  for (y = 0 ; y < rows ; y++)
  {
    for (x = 0 ; x < cols ; x++)
    {
      norm = total = 0.0f ;

      /* apply kernel */
      for (yk = -whalf ; yk <= whalf ; yk++)
      {
        ys = y + yk ;
        if (ys < 0)
          ys = 0 ;
        else if (ys >= rows)
          ys = rows - 1 ;

        for (xk = -whalf ; xk <= whalf ; xk++)
        {
          xs = x + xk ;
          if (xs < 0)
            xs = 0 ;
          else if (xs >= cols)
            xs = cols - 1 ;

          spix = IMAGEFpix(Iin, xs, ys) ;
          Eval = *IMAGEFpix(IE, xs, ys) ;
          Fval = *IMAGEFpix(IF, xs, ys) ;
          Gval = *IMAGEFpix(IG, xs, ys) ;

          fval = (float)exp(-(Eval*xk*xk + 2*Fval*xk*yk + Gval*yk*yk)/sigma_sq);
#if DEBUG_FILTER
          if (x == Gx && y == Gy)
            *IMAGEFpix(Ifilter, xk+whalf, yk+whalf) = fval ;
#endif          
          total += fval * *spix ;
          norm += fval ;
        }
      }

#if DEBUG_FILTER
      if (x == Gx && y == Gy)
        for (yk = -whalf ; yk <= whalf ; yk++)
        {
          for (xk = -whalf ; xk <= whalf ; xk++)
          {
            *IMAGEFpix(Ifilter, xk+whalf, yk+whalf) /= norm ;
          }
        }
#endif          
      *dpix++ = total / norm ;
    }
  }
  
  if (Iin != Isrc)
    ImageFree(&Iin) ;

  if (Iout != Idst)
  {
    /* should worry about using fmin and fmax to rescale byte images here... */
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int        
ImageBuildExponentialFilter(IMAGE *gradImage, int wsize, float k,
                           IMAGE *offsetImage, 
                           IMAGE *filterSequence)
{
  int    rows, cols, x, y, whalf, xc, yc, x0, y0,
         dx, dy, frame ;
  float  fpix, *g, norm, val, *filterPix ;
  static float *gaussian = NULL ;
  static int   w = 0 ;

  rows = gradImage->rows ;
  cols = gradImage->cols ;

  whalf = (wsize-1)/2 ;

  if (wsize != w)
  {
    w = wsize ;
    free(gaussian) ;
    gaussian = NULL ;
  }

  if (!gaussian)  /* allocate a gaussian bump */
  {
    float den ;

    gaussian = (float *)calloc(wsize*wsize, sizeof(float)) ;
    den = wsize*wsize + wsize + 1 ;
    norm = 0.0f ;
    for (g = gaussian, y = 0 ; y < wsize ; y++)
    {
      yc = y - whalf ;
      for (x = 0 ; x < wsize ; x++, g++) 
      {
        xc = x - whalf ;
        *g = (float)exp(-36.0 * sqrt((double)(xc*xc+yc*yc)) / (double)den) ;
        norm += *g ;
      }
    }

    /* normalize gaussian */
    for (g = gaussian, y = 0 ; y < wsize ; y++)
    {
      for (x = 0 ; x < wsize ; x++, g++) 
        *g /= norm ;
    }
  }

/*
  x and y are in window coordinates, while xc and yc are in image
  coordinates.
*/
  for (frame = y0 = 0 ; y0 < rows ; y0++)
  {
    for (x0 = 0 ; x0 < cols ; x0++, frame++)
    {
      if (offsetImage)
      {
        dx = nint(*IMAGEFpix(offsetImage, x0, y0)) ;
        dy = nint(*IMAGEIseq_pix(offsetImage, x0, y0, 1)) ;
      }
      else
        dx = dy = 0 ;

      norm = 0.0f ;
      filterPix = IMAGEFseq_pix(filterSequence, 0, 0, frame) ;

      if (x0 == 5 && y0 == 10)
        x0 = 70 ;

      for (g = gaussian, y = -whalf ; y <= whalf ; y++)
      {
        /* reflect across the boundary */
        yc = y + y0 + dy ;
        if (yc < 0)
          yc = -yc ;
        else if (yc >= rows)
          yc = rows - (yc - rows + 1) ;
        
        for (x = -whalf ; x <= whalf ; x++)
        {
          xc = x0 + x + dx ;
          if (xc < 0)
            xc = -xc ;
          else if (xc >= cols)
            xc = cols - (xc - cols + 1) ;
          
          fpix = *IMAGEFpix(gradImage, xc, yc) ;
          val = (float)exp((double)(-fpix*fpix / k))/*  * *g++ */ ;
          norm += val ;
          *filterPix++ = val ;
        }
      }

      if (FZERO(norm))
        continue ;

      /* normalize kernel weights to sum to 1 */
      filterPix = IMAGEFseq_pix(filterSequence, 0, 0, frame) ;
      for (y = 0 ; y < wsize ; y++)
      {
        for (x = 0 ; x < wsize ; x++)
          *filterPix++ /= norm ;
      }
    }
  }

/*  ImageWrite(filterSequence, "filter.hipl") ;*/
  return(0) ;
}
#if 0
/*----------------------------------------------------------------------
            Parameters:
           Description:
----------------------------------------------------------------------*/
int ImageCalculateMomentOffset(IMAGE *gradImage, int wsize, float c,
                          IMAGE *offsetImage)
{
  return(0) ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageSpaceVariantFilter(IMAGE *inImage, IMAGE *filterSequence, 
                                  IMAGE *outImage)
{
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int        
ImageExponentialFilter(IMAGE *inImage, IMAGE *gradImage, 
                      int wsize, float k,
                           IMAGE *offsetImage, 
                           IMAGE *outImage)
{
  int    rows, cols, x, y, whalf, xc, yc, x0, y0, dx, dy ;
  float  fpix, *g, norm, val, *filterPix, *filter, *outPix ;
  static float w, *gaussian ;

  xc = yc = 0 ;   /* eliminate compiler warning */
  filter = (float *)calloc(wsize*wsize, sizeof(float)) ;
  if (!filter)
    ErrorReturn(ERROR_NO_MEMORY, 
                (ERROR_NO_MEMORY,
                 "ImageExponentialFilter: could not allocate filter")) ;


  rows = gradImage->rows ;
  cols = gradImage->cols ;

  whalf = (wsize-1)/2 ;

  if (wsize != w)
  {
    w = wsize ;
    free(gaussian) ;
    gaussian = NULL ;
  }

  if (!gaussian)  /* allocate a gaussian bump */
  {
    float den ;

    gaussian = (float *)calloc(wsize*wsize, sizeof(float)) ;
    den = wsize*wsize + wsize + 1 ;
    norm = 0.0f ;
    for (g = gaussian, y = 0 ; y < wsize ; y++)
    {
      yc = y - whalf ;
      for (x = 0 ; x < wsize ; x++, g++) 
      {
        xc = x - whalf ;
        *g = (float)exp(-36.0 * sqrt((double)(xc*xc+yc*yc)) / (double)den) ;
        norm += *g ;
      }
    }

    /* normalize gaussian */
    for (g = gaussian, y = 0 ; y < wsize ; y++)
    {
      for (x = 0 ; x < wsize ; x++, g++) 
        *g /= norm ;
    }
  }

/*
  x and y are in window coordinates, while xc and yc are in image
  coordinates.
*/
  outPix = IMAGEFpix(outImage, 0, 0) ;
  for (y0 = 0 ; y0 < rows ; y0++)
  {
    for (x0 = 0 ; x0 < cols ; x0++)
    {
      if (offsetImage)
      {
        dx = nint(*IMAGEFpix(offsetImage, x0, y0)) ;
        dy = nint(*IMAGEFseq_pix(offsetImage, x0, y0, 1)) ;
      }
      else
        dx = dy = 0 ;

      norm = 0.0f ;
      filterPix = filter ;

      for (g = gaussian, y = -whalf ; y <= whalf ; y++)
      {
        /* reflect across the boundary */
        yc = y + y0 + dy ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= rows)
          yc = rows - 1 ;
        for (x = -whalf ; x <= whalf ; x++)
        {
          xc = x0 + x + dx ;
          if (xc < 0)
            xc = 0 ;
          else if (xc >= cols)
            xc = cols - 1 ;
          
          fpix = *IMAGEFpix(gradImage, xc, yc) ;
          val = (float)exp((double)(-fpix*fpix / k))/*  * *g++ */ ;
          norm += val ;
          *filterPix++ = val ;
        }
      }

      if (FZERO(norm)) /* neigborhood is all zeros */
      {
        *outPix++ = 0.0f ;
        continue ;
      }

      /* normalize kernel weights to sum to 1 */
      filterPix = filter ;
      for (y = 0 ; y < wsize ; y++)
      {
        for (x = 0 ; x < wsize ; x++)
          *filterPix++ /= norm ;
      }

/* 
      now apply filter to this point in the image, taking possible 
      offset into accound 
*/
      filterPix = filter ;
      val = 0.0f ;
      for (y = -whalf ; y <= whalf ; y++)
      {
        /* reflect across the boundary */
        yc = y + y0 + dy ;
        if (yc < 0)
          yc = 0 ;
        else if (yc >= rows)
          yc = rows - 1 ;

        for (x = -whalf ; x <= whalf ; x++)
        {
          xc = x0 + x + dx ;
          if (xc < 0)
            xc = 0 ;
          else if (xc >= cols)
            xc = cols - 1 ;
          
          fpix = *IMAGEFpix(inImage, xc, yc) ;
          val += fpix * *filterPix++ ;
        }
      }

      if (isnan(val))
      {
        fprintf(stderr, "(%d, %d) = NaN!\n", xc, yc) ;
      }
      *outPix++ = val ;
    }
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int compare_sort_array(const void *pf1, const void *pf2) ;
int        
ImageMedianFilter(IMAGE *inImage, int wsize, 
                            IMAGE *offsetImage, IMAGE *outImage)
{
  static float *sort_array = NULL ;
  static int   sort_size = 0 ;
  int    ecode, x0, y0, rows, cols, x, y, whalf, yc, dx, dy, frame,wsq,median_index;
  float  *sptr, *outPix, min_val, max_val, *inPix ;
  byte   *in_image, *out_image ;
  IMAGE  *Iout, *Iin ;
  register int xc ;

  rows = inImage->rows ;
  cols = inImage->cols ;
  if (!offsetImage)
  {
    /* h_median only takes byte formatted input and output images */
    if (inImage->pixel_format != PFBYTE)
    {
      Iin = ImageAlloc(rows, cols,PFBYTE, inImage->num_frame);
      ImageValRange(inImage, &min_val, &max_val) ;
      ImageScale(inImage, inImage, 0.0f, 255.0f) ;
      ImageCopy(inImage, Iin) ;
      ImageScale(inImage, inImage, min_val, max_val) ;  /* restore old image */
    }
    else
      Iin = inImage ;

    if (inImage->pixel_format != PFBYTE)
      Iout = ImageAlloc(rows, cols, PFBYTE, outImage->num_frame);
    else
      Iout = outImage ;

    out_image = Iout->image ;
    in_image = Iin->image ;
    for (frame = 0 ; frame < inImage->num_frame ; frame++)
    {
      ecode = h_median(Iin, Iout, wsize) ;
      if (ecode != HIPS_OK)
        ErrorReturn(-1, (ecode, "ImageMedian: h_median failed (%d)", ecode)) ;

      Iout->firstpix += Iout->sizeimage ;
      Iout->image += Iout->sizeimage ;

      Iin->firstpix += Iin->sizeimage ;
      Iin->image += Iin->sizeimage ;
    }

    Iout->firstpix = Iout->image = out_image ;
    Iin->firstpix = Iin->image = in_image ;
    if (Iin != inImage)
      ImageFree(&Iin) ;
    if (outImage != Iout)
    {
      ImageCopy(Iout, outImage) ;
      ImageScale(outImage, outImage, min_val, max_val) ;
      ImageFree(&Iout) ;
    }
    return(NO_ERROR) ;
  }

  median_index = wsize*wsize/2 ;
  wsq = wsize*wsize ;
  whalf = (wsize-1)/2 ;

  /* create a static array for sorting pixels in */
  if (wsize > sort_size)
  {
    sort_size = wsize ;
    if (sort_array)
      sort_array = NULL ;
  }

  if (!sort_array)
    sort_array = (float *)calloc(wsq, sizeof(float)) ;

  for (frame = 0 ; frame < inImage->num_frame ; frame++)
  {
    outPix = IMAGEFseq_pix(outImage, 0, 0, frame) ;
    for (y0 = 0 ; y0 < rows ; y0++)
    {
      for (x0 = 0 ; x0 < cols ; x0++)
      {
/*
         x and y are in window coordinates, while xc and yc are in image
         coordinates.
 */
        if (offsetImage)
        {
          dx = nint(*IMAGEFpix(offsetImage, x0, y0)) ;
          dy = nint(*IMAGEFseq_pix(offsetImage, x0, y0, 1)) ;
        }
        else
          dx = dy = 0 ;
        
        for (sptr = sort_array, y = -whalf ; y <= whalf ; y++)
        {
          /* reflect across the boundary */
          yc = y + y0 + dy ;
#if 0
          if (yc < 0)
            yc = -yc ;
          else if (yc >= rows)
            yc = rows - (yc - rows + 1) ;
#else
          if (yc < 0)
            yc = 0 ;
          else if (yc >= rows)
            yc = rows - 1 ;
#endif

          inPix = IMAGEFseq_pix(inImage, 0, yc, frame) ;
          for (x = -whalf ; x <= whalf ; x++)
          {
            xc = x0 + x + dx ;
#if 0
            if (xc < 0)
              xc = -xc ;
            else if (xc >= cols)
              xc = cols - (xc - cols + 1) ;
#else
            if (xc < 0)
              xc = 0 ;
            else if (xc >= cols)
              xc = cols - 1 ;
#endif

#if 0            
            *sptr++ = *IMAGEFseq_pix(inImage, xc, yc, frame) ;
#else
            *sptr++ = *(inPix + xc) ;
#endif
          }
        }
        qsort(sort_array, wsq, sizeof(float), compare_sort_array) ;
        *outPix++ = sort_array[median_index] ;
      }
    }
  }
  return(0) ;
}
static int
compare_sort_array(const void *pf1, const void *pf2)
{
  register float f1, f2 ;

  f1 = *(float *)pf1 ;
  f2 = *(float *)pf2 ;

/*  return(f1 > f2 ? 1 : f1 == f2 ? 0 : -1) ;*/
  if (f1 > f2)
    return(1) ;
  else if (f1 < f2)
    return(-1) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageConvolveGaussian(IMAGE *Isrc,IMAGE *gImage, IMAGE *outImage,
                     int dst_frameno)
{
  static IMAGE     *tmpImage = NULL ;
  int              ksize ;
  float            *kernel, *buf ;

  if (!ImageCheckSize(Isrc, tmpImage, 0, 0, 0))
  {
    if (tmpImage)
      ImageFree(&tmpImage) ;
    tmpImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
  }
  if (!outImage)
    outImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;

  ImageSetSize(tmpImage, Isrc->rows, Isrc->cols) ;

  kernel = IMAGEFpix(gImage, 0, 0) ;
  ksize = gImage->cols ;
  ImageConvolve1d(Isrc, tmpImage, kernel, ksize, IMAGE_VERTICAL) ;

  buf = IMAGEFpix(outImage, 0, 0) ;
  outImage->image = (byte *)IMAGEFseq_pix(outImage, 0, 0, dst_frameno) ;
  ImageConvolve1d(tmpImage, outImage, kernel, ksize, IMAGE_HORIZONTAL) ;

  outImage->image = (byte *)buf ;
  return(outImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
ImageConvolve1d(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int           x, y, width, height, halflen ;
  register int  xi,yi, i ;
  float         total, *ki, *outPix, *inPix, *inBase ;

  width = I->cols ;
  height = I->rows ;

  halflen = len/2 ;
  outPix = IMAGEFpix(J, 0, 0) ;
  if (axis == IMAGE_HORIZONTAL)
  {
    for (y = 0 ; y < height ; y++)
    {
      inBase = IMAGEFpix(I, 0, y) ;
      for (x = 0 ; x < width ; x++)
      {
        total = 0.0f ;

        for (ki = k, i = 0 ; i < len ; i++)
        {
          xi = x + i - halflen ;
          if (xi < 0)
            xi = 0 ;
          else if (xi >= width)
            xi = width - 1 ;

          inPix = inBase + xi ;
          total += *ki++ * *inPix ;
        }

        *outPix++ = total ;
      }
    }
  }
  else
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        inBase = IMAGEFpix(I, x, 0) ;
        total = 0.0f ;

        for (ki = k, i = 0 ; i < len ; i++)
        {
          yi = y + i - halflen ;
          if (yi < 0)
            yi = 0 ;
          else if (yi >= height)
            yi = height - 1 ;

          inPix = inBase + yi*width ;
          total += *ki++ * *inPix ;
        }
        *outPix++ = total ;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              use the kernel described by Burt and Adelson to reduce
              the pyramid 1 level.
----------------------------------------------------------------------*/
#define KERNEL_SIZE 5
#define K_A         0.4f
static float kernel[KERNEL_SIZE]  = 
{
  0.25f - K_A/2.0f, .25f, K_A, 0.25f, 0.25f-K_A/2.0f
} ;

IMAGE *
ImageReduce(IMAGE *Isrc, IMAGE *outImage)
{
  int  rows, cols ;
  static IMAGE *tmpImage = NULL ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;
  if (!ImageCheckSize(Isrc, tmpImage, rows, cols, 0) && tmpImage)
  {
    ImageFree(&tmpImage) ;
    tmpImage = NULL ;
  }

  if (!tmpImage)
  {
    tmpImage = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    if (!tmpImage)
      return(NULL) ;
  }
  else
  {
    ImageSetSize(tmpImage,rows,cols) ;
    ImageClearArea(tmpImage, 0, 0, -1, -1, 0.0f) ;
  }

  rows /= 2 ;
  cols /= 2 ;

  if (!outImage)
    outImage = ImageAlloc(rows, cols, PFFLOAT, 1) ;
  else
  {
    if (!ImageCheckSize(Isrc, outImage, rows, cols, 0))
      ErrorReturn(outImage, (ERROR_NO_MEMORY,
                             "ImageReduce: output image is too small\n")) ;
  }

  /* blur vertically */
  ImageConvolve1d(Isrc, tmpImage, kernel, KERNEL_SIZE, IMAGE_VERTICAL) ;
  ImageReduce1d(tmpImage, outImage, kernel, KERNEL_SIZE, IMAGE_HORIZONTAL) ;
#if 0
{
  char str[100] ;
  sprintf(str, "tmp%d.hipl", rows*2) ;
  ImageWrite(tmpImage, str) ;
  sprintf(str, "out%d.hipl", rows*2) ;
  ImageWrite(outImage, str) ;
}
#endif
  return(outImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             use k[] to scale the image down by 2.
----------------------------------------------------------------------*/
void
ImageReduce1d(IMAGE *I, IMAGE *J, float k[], int len, int axis)
{
  int    x, y, i, Jwidth, Jheight, xi,yi, halflen, Iwidth, Iheight ;
  float  total ;

  Jwidth = J->cols ;
  Jheight = J->rows ;
  Iwidth = I->cols ;
  Iheight = I->rows ;

  halflen = (len-1)/2 ;
  if (axis == IMAGE_HORIZONTAL)
  {
    for (y = 0 ; y < Jheight ; y++)
    {
      yi = 2*y ;
      for (x = 0 ; x < Jwidth ; x++)
      {
        total = 0.0f ;

        if (x >= 3 && x <= 6)
          i = 0 ;

        for (i = 0 ; i < len ; i++)
        {
          /* Neumann boundary conditions */
          xi = 2*x + i - halflen ;
#if 0
          if (xi < 0)
            xi = -xi ;
          else if (xi >= Iwidth)
            xi = Iwidth - (xi - Iwidth+1) ;
#else
          if (xi < 0)
            xi = 0 ;
          else if (xi >= Iwidth)
            xi = Iwidth - 1 ;
#endif
          
          total = total + k[i] * *IMAGEFpix(I, xi, yi) ;
        }
        *IMAGEFpix(J, x, y) = total ;
      }
    }
  }
  else
  {
    for (y = 0 ; y < Jheight ; y++)
    {
      for (x = 0 ; x < Jwidth ; x++)
      {
        total = 0.0f ;
        xi = 2*x ;

#if 0
        if (((x == 9) && (y == 127)) ||
            ((y == 9) && (x == 127)))
          i = 0 ;
#endif

        for (i = 0 ; i < len ; i++)
        {
          /* Neumann boundary conditions */
          yi = 2*y + i - halflen ;
#if 0
          if (yi < 0)
            yi = -yi ;
          else if (yi >= Iheight)
            yi = Iheight - (yi - Iheight+1) ;
#else
          if (yi < 0)
            yi = 0 ;
          else if (yi >= Iheight)
            yi = Iheight - 1 ;
#endif
          
          total = total + k[i] * *IMAGEFpix(I, xi, yi) ;
        }
        *IMAGEFpix(J, x, y) = total ;
      }
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a gaussian bump which tails to 0. Returns
             an image which is ((4*xsigma)+1) x ((4*ysigma)+1)
----------------------------------------------------------------------*/
IMAGE *
ImageGaussian(float xsigma, float ysigma)
{
  IMAGE *image ;
  float     norm, ytwo_sigma, xtwo_sigma, fx, fy, k ;
  int       x, y, xlen, ylen, xhalf, yhalf ;

  /* build the kernel in k */
  xlen = (int)nint(6.0f * xsigma)+1 ;
  ylen = (int)nint(6.0f * ysigma)+1 ;
  if (EVEN(xlen))
    xlen++ ;
  if (EVEN(ylen))
    ylen++ ;
  xhalf = xlen/2 ;
  yhalf = ylen/2 ;
  image = ImageAlloc(xlen, ylen, PFFLOAT, 1) ;

  norm = 0.0f ;
  xtwo_sigma = 2.0f * xsigma ;
  ytwo_sigma = 2.0f * ysigma ;

  for (x = 0 ; x < xlen ; x++)
  {
    fx = (float)(x-xhalf) ;
    if (fabs(fx) <= xtwo_sigma)
      k = (float)exp((double)(-fx*fx/(xtwo_sigma*xsigma))) ; 
    else if (xtwo_sigma < (float)fabs(fx) && (float)fabs(fx) <= 4.0f*xsigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * 
        (float)pow(4.0 - fabs(fx)/(double)xsigma, 4.0) ;
    else
      k = 0 ;

    for (y = 0 ; y < ylen ; y++)
      *IMAGEFpix(image, x, y) = k ;
  }

  for (y = 0 ; y < ylen ; y++)
  {
    fy = (float)(y-yhalf) ;
    if (fabs(fy) <= ytwo_sigma)
      k = (float)exp((double)(-fy*fy/(ytwo_sigma*ysigma))) ;
    else if (ytwo_sigma < fabs(fy) && fabs(fy) <= 4*ysigma)
      k = 1.0f / (16.0f * (float)(M_E * M_E)) * (float)pow(4.0 - fabs(fy)/(double)ysigma, 4.0) ;
    else
      k = 0 ;

    for (x = 0 ; x < xlen ; x++)
    {
      *IMAGEFpix(image, x, y) *= k ;
      norm += *IMAGEFpix(image, x, y) ;
    }
  }

  /* normalize kernel to sum to 1 */
  for (x = 0 ; x < xlen ; x++)
  {
    for (y = 0 ; y < ylen ; y++)
      *IMAGEFpix(image, x, y) /= norm ;
  }

  return(image) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             construct a splined gaussian bump which tails to 0. 
             Returns an image which is (8*sigma)+1 
             (Nitzberg and Shiota, 1993)
----------------------------------------------------------------------*/
IMAGE *
ImageGaussian1d(float sigma, int max_len)
{
  IMAGE *image ;
  float     norm, two_sigma, fx, k ;
  int       x, half, len ;

  /* build the kernel in k */
  len = (int)nint(8.0f * sigma)+1 ;
  if (max_len && max_len > len)
    len = max_len ;
  half = len/2 ;
  image = ImageAlloc(1, len, PFFLOAT, 1) ;

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

    *IMAGEFpix(image, x, 0) = k ;
    norm += k ;
  }

  /* normalize kernel to sum to 1 */
  for (x = 0 ; x < len ; x++)
    *IMAGEFpix(image, x, 0) /= norm ;

  return(image) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int        
ImageSobel(IMAGE *Isrc, IMAGE *gradImage, 
                     IMAGE *dxImage, IMAGE *dyImage)
{
  static IMAGE *xImage = NULL, *yImage = NULL ;
  IMAGE        *Iin ;
  int          x, y, rows, cols ;
  float        *xpix, *ypix, *gradpix = NULL, xval, yval, gval ;

  if (Isrc->pixel_format != PFFLOAT)
  {
    Iin = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
    ImageCopy(Isrc, Iin) ;
  }
  else
    Iin = Isrc ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (!dxImage)
  {
    if (!ImageCheckSize(Isrc, xImage, 0, 0, 0))
    {
      if (xImage)
        ImageFree(&xImage) ;
      xImage = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    }
    else
      ImageSetSize(xImage, rows, cols) ;

    dxImage = xImage ;
  }

  if (!dyImage)
  {
    if (!ImageCheckSize(Isrc, yImage, 0, 0, 0))
    {
      if (yImage)
        ImageFree(&yImage) ;
      yImage = ImageAlloc(rows, cols, PFFLOAT, 1) ;
    }
    else
      ImageSetSize(yImage, rows, cols) ;

    dyImage = yImage ;
  }

  
  ImageSetSize(dxImage, rows, cols) ;
  ImageSetSize(dyImage, rows, cols) ;
#if 0
  ImageConvolve3x3(Isrc, sx, dxImage) ;
  ImageConvolve3x3(Isrc, sy, dyImage) ;
#else
  ImageSobelX(Isrc, dxImage) ;
  ImageSobelY(Isrc, dyImage) ;
#endif
  if (gradImage)
  {
    ImageSetSize(gradImage, rows, cols) ;
    gradpix = IMAGEFpix(gradImage, 0, 0) ;

    xpix = IMAGEFpix(dxImage, 0, 0) ;
    ypix = IMAGEFpix(dyImage, 0, 0) ;
    for (y = 0 ; y < rows ; y++)
    {
      for (x = 0 ; x < cols ; x++)
      {
        xval = *xpix++ ;
        yval = *ypix++ ;
        gval = (float)sqrt((double)(xval * xval + yval * yval)) ;
        *gradpix++ = gval ;
      }
    }
  }

  if (Iin != Isrc)
    ImageFree(&Iin) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
/*
   -0.25    0    0.25
   -0.50    0    0.50
   -0.25    0    0.25
*/
int
ImageSobelX(IMAGE *Isrc, IMAGE *xImage)
{
  register float *tl_pix, *ml_pix, *bl_pix, *tr_pix, *mr_pix, *br_pix, *outPtr;
  int     rows, cols, row, col ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;
  outPtr = IMAGEFpix(xImage, 1, 1) ;
  tl_pix = IMAGEFpix(Isrc, 0, 0) ;
  ml_pix = IMAGEFpix(Isrc, 0, 1) ;
  bl_pix = IMAGEFpix(Isrc, 0, 2) ;
  tr_pix = IMAGEFpix(Isrc, 2, 0) ;
  mr_pix = IMAGEFpix(Isrc, 2, 1) ;
  br_pix = IMAGEFpix(Isrc, 2, 2) ;

  /* don't apply sobel to outer ring to pixels to avoid border effects */
  rows-- ;
  cols-- ;
  for (row = 1 ; row < rows ; row++)
  {
    for (col = 1 ; col < cols ; col++)
    {
      *outPtr++ =
        -.25f * *tl_pix++ - .5f * *ml_pix++ - .25f * *bl_pix++ +
         .25f * *tr_pix++ + .5f * *mr_pix++ + .25f * *br_pix++ ;
    }
    outPtr += 2 ;
    tl_pix += 2 ;
    ml_pix += 2 ;
    bl_pix += 2 ;
    tr_pix += 2 ;
    mr_pix += 2 ;
    br_pix += 2 ;
  }

  return(0) ;
}
/*
   -0.25   -.50   -0.25
    0       0      0
    0.25    .50    0.25
*/
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageSobelY(IMAGE *Isrc, IMAGE *yImage)
{
  register float *tl_pix, *tm_pix,*tr_pix, *bl_pix, *bm_pix, *br_pix, *outPtr;
  int     rows, cols, row, col ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;
  outPtr = IMAGEFpix(yImage, 1, 1) ;
  tl_pix = IMAGEFpix(Isrc, 0, 0) ;
  tm_pix = IMAGEFpix(Isrc, 1, 0) ;
  tr_pix = IMAGEFpix(Isrc, 2, 0) ;
  bl_pix = IMAGEFpix(Isrc, 0, 2) ;
  bm_pix = IMAGEFpix(Isrc, 1, 2) ;
  br_pix = IMAGEFpix(Isrc, 2, 2) ;

  /* don't apply sobel to outer ring to pixels to avoid border effects */
  rows-- ;
  cols-- ;
  for (row = 1 ; row < rows ; row++)
  {
    for (col = 1 ; col < cols ; col++)
    {
      *outPtr++ =
        -.25f * *tl_pix++ - .5f * *tm_pix++ - .25f * *tr_pix++ +
         .25f * *bl_pix++ + .5f * *bm_pix++ + .25f * *br_pix++ ;
    }
    outPtr += 2 ;
    tl_pix += 2 ;
    tm_pix += 2 ;
    tr_pix += 2 ;
    bl_pix += 2 ;
    bm_pix += 2 ;
    br_pix += 2 ;
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageXDerivative(IMAGE *Isrc, IMAGE *xImage)
{
  if (!xImage)
    xImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;

  ImageConvolve3x3(Isrc, sx, xImage) ;

  return(xImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageYDerivative(IMAGE *Isrc, IMAGE *yImage)
{
  if (!yImage)
    yImage = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;

  ImageConvolve3x3(Isrc, sy, yImage) ;

  return(yImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageConvolve3x3(IMAGE *Isrc, float kernel[], IMAGE *outImage)
{
  int     rows, cols, x, y, xk, yk, k, xi, yi, frame ;
  float   *fkpix, sum, *fopix, *fipix ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;

  for (frame = 0 ; frame < Isrc->num_frame ; frame++)
  {
    switch (Isrc->pixel_format)
    {
    case PFFLOAT:
      fopix = (float *)IMAGEFseq_pix(outImage, 0, 0, frame) ;
      for (y = 0 ; y < rows ; y++)
      {
        for (x = 0 ; x < cols ; x++, fopix++)
        {
          fkpix = kernel ;
          for (sum = 0.0f, k = 0, yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y + yk ;    /* image coordinate */
            if (yi < 0)
              yi = 0 ;
            else if (yi >= rows)
              yi = rows-1 ;
            
            for (xk = -1 ; xk <= 1 ; xk++, k++, fkpix++)
            {
              xi = x + xk ;   /* image coordinate */
              if (xi < 0)
                xi = 0 ;
              else if (xi >= cols)
                xi = cols-1 ;
              fipix = IMAGEFseq_pix(Isrc, xi, yi, frame) ;
              sum = sum + *fipix * *fkpix ;
            }
          }
          *fopix = sum ;
        }
      }
      break ;
    default:
      fprintf(stderr, "ImageConvolve3x3: unsupported pixel format %d\n",
              Isrc->pixel_format) ;
      exit(-1);
      break ;
    }
  }

  return(0) ;
}
