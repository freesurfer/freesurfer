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
