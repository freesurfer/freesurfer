/*
 *       FILE NAME:   image.c
 *
 *       DESCRIPTION: image processing utilities
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
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
#include <fcntl.h>
#include <unistd.h>

#include <hipl_format.h>

#include "image.h"
#include "error.h"
#include "matrix.h"
#include "matfile.h"
#include "utils.h"
#include "macros.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
static int alloc_image(struct header *hd) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE   *
ImageAlloc(int rows, int cols, int format, int nframes)
{
  IMAGE *I ;

  I = (IMAGE *)calloc(1, sizeof(IMAGE)) ;
  if (!I)
    ErrorExit(ERROR_NO_MEMORY, "ImageAlloc: could not allocate header\n") ;

  init_header(I, "orig", "seq", nframes, "today", rows,cols,format,1, "temp");
  alloc_image(I) ;
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           stolen from hips2 code and modified to allocate multiple frames.
------------------------------------------------------*/
static int
alloc_image(struct header *hd)
{
  int fcb,cb;

  if (hd->sizeimage == 0) 
  {
    hd->imdealloc = FALSE;
    return(HIPS_OK);
  }
  if ((hd->image = hmalloc(hd->sizeimage*hd->num_frame)) == (byte *)HIPS_ERROR)
    return(HIPS_ERROR);
  if (hd->pixel_format == PFMSBF || hd->pixel_format == PFLSBF) 
  {
    fcb = hd->fcol/8;
    cb = (hd->ocols + 7)/8;
    hd->firstpix = hd->image + ((cb * hd->frow) + fcb);
  }
  else
    hd->firstpix = hd->image +
      ((hd->ocols * hd->frow) + hd->fcol) * hd->sizepix;
  hd->imdealloc = TRUE;
  memset(hd->image, 0, hd->sizeimage*hd->num_frame) ;
  return(HIPS_OK);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageFree(IMAGE **pI)
{
  IMAGE *I = *pI ;

  if (!I)
    ErrorExit(ERROR_BADPARM, "ImageFree: null pointer") ;

  free_header(I) ;
  *pI = NULL ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageWrite(IMAGE *I, char *fname)
{
  FILE  *fp ;
  int   ecode ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    ErrorExit(ERROR_NO_FILE, "ImageWrite(%s) failed\n", fname) ;

  ecode = ImageFWrite(I, fp, fname) ;
  fclose(fp) ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageFWrite(IMAGE *I, FILE *fp, char *fname)
{
  int ecode, frame ;
  char *image ;

  if (!fname)
    fname = "ImageFWrite" ;

  ecode = fwrite_header(fp,I,"fwrite") ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageFWrite: fwrite_header failed (%d)\n",ecode);

  image = I->image ;
  for (frame = 0 ; frame < I->num_frame ; frame++)
  {
    ecode = fwrite_image(fp, I, frame, "fwrite") ;
    if (ecode != HIPS_OK)
      ErrorExit(ERROR_NO_FILE, 
              "ImageFWrite: fwrite_image frame %d failed (%d)\n",ecode,frame);
    I->image += I->sizeimage ;  /* next frame */
  }
  I->image = image ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageFRead(FILE *fp, char *fname, int frame)
{
  int    ecode, end_frame ;
  IMAGE  *I ;
  char   *startpix ;

  if (!fname)
    fname = "ImageFRead" ;

  I = (IMAGE *)calloc(1, sizeof(IMAGE)) ;
  if (!I)
    ErrorExit(ERROR_NO_MEMORY,"ImageFRead: could not allocate header\n") ;

  ecode = fread_header(fp, I, fname) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_header failed (%d)\n",ecode);

  if (frame < 0)    /* read all frames */
  {
    end_frame = I->num_frame-1 ;
    frame = 0 ;
  }
  else              /* read only specified frame */
  {
    if (fseek(fp, I->sizeimage*frame, SEEK_CUR) < 0)
    {
      ImageFree(&I) ;
      ErrorReturn(NULL, 
                  (ERROR_BADFILE, 
                   "ImageFRead(%s, %d) - could not seek to specified frame",
                   fname, frame)) ;
    }
    end_frame = frame ;
  }
  if (end_frame >= I->num_frame)
    ErrorReturn(NULL, 
                (ERROR_BADFILE,
                 "ImageFRead(%s, %d) - frame out of bounds", fname, frame)) ;
  I->num_frame = end_frame - frame + 1 ;
  alloc_image(I) ;

  startpix = I->image ;
  for ( ; frame <= end_frame ; frame++)
  {
    ecode = fread_image(fp, I, frame, fname) ;
    if (ecode != HIPS_OK)
      ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_image failed (%d)\n", ecode);
    I->image += I->sizepix ;
  }
  I->image = startpix ;
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageRead(char *fname)
{
  IMAGE   *I = NULL ;
  MATRIX  *mat ;
  FILE    *fp ;
  int     type, frame ;
  char    buf[100] ;

  strcpy(buf, fname) ;   /* don't destroy callers string */
  fname = buf ;
  ImageUnpackFileName(fname, &frame, &type, fname) ;

  switch (type)
  {
  case MATLAB_IMAGE:
    mat = MatlabRead(fname) ;
    if (!mat)
      ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s) failed\n", fname)) ;
    I = ImageFromMatrix(mat, NULL) ;
    h_invert(I,I) ;
    MatrixFree(&mat) ;
    break ;
  case HIPS_IMAGE:
    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s, %d) failed\n", 
                         fname, frame)) ;
    I = ImageFRead(fp, fname, frame) ;
    fclose(fp) ;
    break ;
  default:
    break ;
  }
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static unsigned char thicken_se[9] = { 255, 1, 255, 1, 1, 1, 255, 1, 255 } ;
IMAGE *
ImageDilate(IMAGE *Isrc, IMAGE *Idst, int which)
{
  int    ecode, center_row, center_col, gray ;
  IMAGE  Ise, *Iin, *Itmp, *Iout ;

  init_header(&Ise, "orig", "seq", 1, "today", 3, 3, PFBYTE, 1, "temp");

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;

  if (Isrc->pixel_format != PFBYTE)
  {
    Itmp = ImageScale(Isrc, NULL, 0, 255) ;
    Iin = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;
    ImageCopy(Itmp, Iin) ;
    ImageFree(&Itmp) ;
  }
  else
    Iin = Isrc ;

  if (Idst->pixel_format != PFBYTE)
    Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;
  else
    Iout = Idst ;

  ImageReplace(Iin, Iin, 255.0f, 1.0f) ;
  ImageReplace(Iin, Iin, 0.0f, 255.0f) ;

  switch (which)
  {
  case MORPH_THICKEN:
    Ise.image = Ise.firstpix = thicken_se ;
    center_row = Ise.rows / 2 ;
    center_col = Ise.cols / 2 ;
    ecode = h_morphdil(Iin, &Ise, Iout, center_row, center_col, 128) ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageMorph: unsupported morphological operation %d\n", which) ;
    break ;
  }

  if (Iin != Isrc)   /* source image was not in proper format */
    ImageFree(&Iin) ;
  else    /* translate back to original format */
  {
    ImageReplace(Isrc, Isrc, 255.0f, 0.0f) ;
    ImageReplace(Isrc, Isrc, 1.0f, 255.0f) ;
  }

  if (Idst != Iout)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }

  ImageReplace(Idst, Idst, 255.0f, 0.0f) ;
  ImageReplace(Idst, Idst, 1.0f, 255.0f) ;
  
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE   *
ImageErode(IMAGE *Isrc, IMAGE *Idst, int which)
{
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageThreshold(IMAGE *Isrc, IMAGE *Idst, float threshold)
{
  int   ecode ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols,Isrc->pixel_format,
                      Isrc->num_frame);

  ecode = h_softthresh(Isrc, Idst, threshold) ;
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageDFT(IMAGE *Isrc, IMAGE *Idst)
{
  float    loglen ;
  int      ecode ;
  IMAGE    *Itmp ;
  Pixelval p ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFCOMPLEX, 1) ;

  if (Isrc->pixel_format == PFBYTE)  /* must convert to float */
  {
    Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
    ImageCopy(Isrc, Itmp) ;
    Isrc = Itmp ;
  }
  else
    Itmp = NULL ;

  hips_rtocplx = CPLX_RVI0 ;
  ecode = h_toc(Isrc, Idst) ;
  if  (ecode != HIPS_OK)
  {
    ImageFree(&Idst) ;
    ErrorReturn(NULL, (ecode, "ImageDFT: h_fourtr failed (%d)\n", ecode)) ;
  }

  loglen = log2(Isrc->rows) ;
  if ((Isrc->rows == Isrc->cols) && (floor(loglen) == loglen)) /* FFT */
  {
  }
  else   /* not a power of 2 - use DFT */
  {
  }

  ecode = h_fourtr(Idst) ;
  if  (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageDFT: h_fourtr failed (%d)\n", ecode) ;

#if 0
  /* h_divscale can't handle complex quantities */
  p.v_complex[REAL_PIX] = 1.0f/(float)(Idst->numpix) ;
  p.v_complex[IMAG_PIX] = 0.0f ;
  ImageMulScale(Idst, Idst, &p) ;
#endif

  if (Itmp)
    ImageFree(&Itmp) ;

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           perform the inverse fourier transform of the input image,
           scaling the output by 1/n.
------------------------------------------------------*/
IMAGE *
ImageInverseDFT(IMAGE *Isrc, IMAGE *Idst)
{
  float    loglen ;
  IMAGE    *Itmp ;
  int      ecode ;
  Pixelval p ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
  Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFCOMPLEX, 1) ;

  ImageCopy(Isrc, Itmp) ;

  loglen = log2(Isrc->rows) ;
  if ((Isrc->rows == Isrc->cols) && (floor(loglen) == loglen)) /* FFT */
  {
  }
  else   /* not a power of 2 - use DFT */
  {
  }

  ecode = h_invfourtr(Itmp) ;
  if  (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageInverseDFT: h_invfourtr failed (%d)\n", ecode) ;

  hips_cplxtor = CPLX_REAL ;
  h_tof(Itmp, Idst) ;
#if 1
  p.v_float = (float)Idst->numpix ;
  ecode = h_divscale(Idst, Idst, &p) ;
  if  (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageInverseDFT: h_divscale failed (%d)\n", ecode) ;
#endif
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageMul(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst)
{
  int   ecode ;

  if (!Idst)
    Idst = ImageAlloc(Isrc1->rows, Isrc1->cols, 
                      Isrc1->pixel_format, Isrc1->num_frame) ;

  ecode = h_mul(Isrc1, Isrc2, Idst) ;
  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageMul: h_mul failed (%d)\n", ecode) ;

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageResize(IMAGE *Isrc, IMAGE *Idst, int drows, int dcols)
{
  float x_scale, y_scale ;
  int   ecode ;

  if (!Idst)
    Idst = ImageAlloc(drows, dcols, Isrc->pixel_format, Isrc->num_frame) ;

  x_scale = (float)dcols / (float)Isrc->cols ;
  y_scale = (float)drows / (float)Isrc->rows ;

  if (FEQUAL(x_scale, 1.0f))
    ImageCopy(Isrc, Idst) ;
  else
    switch (Isrc->pixel_format)
    {
    case PFBYTE:
#if 0
      ecode = h_affine(Isrc, Idst, x_scale, 0.0f, 0.0f, 0.0f, y_scale, 0.0f) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode,
                  "ImageResize: h_affine(%2.3f, %2.3f) returned %d\n",ecode);
#else
      if (x_scale > 1.0f)
        ecode = h_enlarge(Isrc, Idst, x_scale, y_scale) ;
      else 
        ecode = h_reduce(Isrc, Idst, nint(1.0f/x_scale), nint(1.0f/y_scale)) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode,
                  "ImageResize: h_%s(%2.3f, %2.3f) returned %d\n",
                  x_scale > 1.0f ? "enlarge" : "reduce", ecode);
#endif
      break ;
    default:
      if (x_scale > 1.0f)
        ecode = h_enlarge(Isrc, Idst, nint(x_scale), nint(y_scale)) ;
      else 
        ecode = h_reduce(Isrc, Idst, nint(1.0f/x_scale), nint(1.0f/y_scale)) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode,
                  "ImageResize: h_%s(%2.3f, %2.3f) returned %d\n",
                  x_scale > 1.0f ? "enlarge" : "reduce", ecode);
      break ;
    }
  
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageCopy(IMAGE *Isrc, IMAGE *Idst)
{
  int  old, ecode ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  if (Idst->pixel_format == Isrc->pixel_format)
  {
    ecode = h_copy(Isrc, Idst) ;
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageCopy: h_copy failed (%d)\n", ecode) ;
  }
  else
  {
    switch (Idst->pixel_format)
    {
    case PFFLOAT:
      old = hips_cplxtor ;
      hips_cplxtor = CPLX_REAL ;
      ecode = h_tof(Isrc, Idst) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode, "ImageCopy: h_tof failed (%d)\n", ecode) ;
      hips_cplxtor = old ;
      break ;
    case PFCOMPLEX:
      old = hips_rtocplx ;
      hips_rtocplx = CPLX_RVI0 ;
      ecode = h_toc(Isrc, Idst) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode, "ImageCopy: h_toc failed (%d)\n", ecode) ;
      hips_rtocplx = old ;
      break ;
    case PFDBLCOM:
      old = hips_rtocplx ;
      hips_rtocplx = CPLX_RVI0 ;
      ecode = h_todc(Isrc, Idst) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode, "ImageCopy: h_todc failed (%d)\n", ecode) ;
      hips_rtocplx = old ;
      break ;
    case PFINT:
      ecode = h_toi(Isrc, Idst) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode, "ImageCopy: h_toi failed (%d)\n", ecode) ;
      break ;
    case PFBYTE:
      ecode = h_tob(Isrc, Idst) ;
      if (ecode != HIPS_OK)
        ErrorExit(ecode, "ImageCopy: h_tob failed (%d)\n", ecode) ;
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED,
                "ImageCopy %d-->%d, unsupported conversion\n",
                Isrc->pixel_format, Idst->pixel_format) ;
      break ;
    }
  }

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#define DISCEDGE_VARCRIT  0.0f
#define DISCEDGE_SIZE     7

IMAGE   *
ImageEdgeDetect(IMAGE *Isrc, IMAGE *Idst, float sigma, int wsize, float lthresh, float uthresh,
                int dothin)
{
  int    ecode ;
  IMAGE  *Itmp, *Iout ;
  char   *command_str = "" ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;

  ImageScale(Isrc, Isrc, 0.0f, 255.0f) ;
  if (Idst->pixel_format != PFBYTE)
    Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;
  else
    Iout = Idst ;

  ecode = h_canny(Isrc, Iout, sigma, wsize, lthresh, uthresh, dothin) ;
  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE   *
ImageCorrelate(IMAGE *Itemplate, IMAGE *Isrc, int zeropad, IMAGE *Icorr)
{
  IMAGE *Iconj, *Ifcorr, *Ifsrc ;
  int   ecode ;

  /* assumes the template as already been FTed */
  Iconj = ImageConjugate(Itemplate, NULL) ;
  Ifsrc = ImageDFT(Isrc, NULL) ;
  Ifcorr = ImageMul(Iconj, Ifsrc, NULL) ;
  Icorr = ImageInverseDFT(Ifcorr, Icorr) ;

  ecode = h_flipquad(Icorr, Icorr) ;
  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode) ;

#if 0
ImageWrite(Itemplate, "Itemplate.hipl") ;
ImageWrite(Isrc, "Isrc.hipl") ;
ImageWrite(Ifsrc, "Ifsrc.hipl") ;
ImageWrite(Iconj, "Iconj.hipl") ;
ImageWrite(Ifcorr, "Ifcorr.hipl") ;
ImageWrite(Icorr, "Iflip.hipl") ;
#endif

  ImageFree(&Iconj) ;
  ImageFree(&Ifcorr) ;
  ImageFree(&Ifsrc) ;

  return(Icorr) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE   *
ImageCopyArea(IMAGE *Isrc, IMAGE *Idst, int srow, int scol,
            int drow, int dcol, int rows, int cols)
{
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageClearArea(IMAGE *I, int r0, int c0, int rows, int cols, float val)
{
  float   *fptr ;
  int     row, col ;
  char    *cptr, cval ;

  rows = MIN(I->rows, r0+rows) ;
  cols = MIN(I->cols, c0+cols) ;
  
  for (row = r0 ; row < rows ; row++)
  {
    switch (I->pixel_format)
    {
    case PFFLOAT:
      fptr = IMAGEFpix(I, c0, row) ;
      for (col = c0 ; col < cols ; col++)
        *fptr++ = val ;
      break ;
    case PFBYTE:
      cptr = IMAGEpix(I, c0, row) ;
      cval = (char)val ;
      for (col = c0 ; col < cols ; col++)
        *cptr++ = cval ;
      break ;
    default:
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "ImageClearArea: only handles PFFLOAT")) ;
      break ;
    }
  }
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float
ImageFindPeak(IMAGE *I, int *prow, int *pcol, float *pval)
{
  float  max_val, *fpix, val ;
  int    max_row, max_col, row, col, rows, cols ;

  if (I->pixel_format != PFFLOAT)
    ErrorReturn(0.0f, (ERROR_UNSUPPORTED, "ImageFindPeak: only supports PFFLOAT")) ;

  rows = I->rows ;
  cols = I->cols ;

  fpix = IMAGEFpix(I, 0, 0) ;
  max_val = -1000000.0f ;
  for (row = 0 ; row < rows ; row++)
  {
    for (col = 0 ; col < cols ; col++)
    {
      val = *fpix++ ;
      if (val >= max_val)
      {
        max_val = val ;
        max_row = row ;
        max_col = col ;
      }
    }
  }

  *prow = max_row ;
  *pcol = max_col ;
  return(max_val) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImagePowerSpectrum(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE  *Idft, *Iconj ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, Isrc->num_frame) ;

  Idft = ImageAlloc(Isrc->rows, Isrc->cols, PFCOMPLEX, Isrc->num_frame) ;

  if (Isrc->pixel_format != PFCOMPLEX)  /* not FFT'd yet */
    ImageDFT(Isrc, Idft) ;
  else
    ImageCopy(Isrc, Idft) ;

  Iconj = ImageConjugate(Idft, NULL) ;

  ImageMul(Idft, Iconj, Iconj) ;
  ImageCopy(Iconj, Idst) ;       /* change it to floating point */

  ImageFree(&Iconj) ;

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
imageLargeEnough(IMAGE *Isrc, IMAGE *Idst)
{
  if (Isrc->num_frame > Idst->num_frame)
    return(0) ;
  if (Isrc->numpix > Idst->numpix)
    return(0) ;

  return(1) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageNormalizePix(IMAGE *Isrc, IMAGE *Idst)
{
  float  scale, fmin, fmax ;
  int    ecode ;
  byte   bmin, bmax ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  switch (Isrc->pixel_format)
  {
  case PFBYTE:
    ecode = h_minmax(Isrc, &bmin, &bmax, 0) ;
    fmin = (float)bmin ;
    fmax = (float)bmax ;
    break ;
  case PFFLOAT:
    ecode = h_minmax(Isrc, &fmin, &fmax, 0) ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageNormalize: unsupported pixel format %d\n",
              Isrc->pixel_format) ;
    break ;
  }

  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageNormalize: h_minmax failed (%d)\n", ecode) ;

  scale = 1.0f / (fmax - fmin) ;
  fmin = -fmin * scale ;
  ecode = h_linscale(Isrc, Idst, scale, fmin) ;
  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageNormalize: h_linscale failed (%d)\n", ecode) ;

  ecode = h_minmax(Isrc, &fmin, &fmax, 0) ;
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageConjugate(IMAGE *Isrc, IMAGE *Idst)
{
  CPIX       *spix, *dpix ;
  int        npix, i ;

  npix = Isrc->orows * Isrc->ocols * Isrc->num_frame ;
  switch (Isrc->pixel_format)
  {
  case PFCOMPLEX:
    if (!Idst)
      Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFCOMPLEX, Isrc->num_frame) ;
    spix = (CPIX *)IMAGECpix(Isrc, 0, 0) ;
    dpix = (CPIX *)IMAGECpix(Idst, 0, 0) ;
    for (i = 0 ; i < npix ; i++, spix++, dpix++)
    {
      dpix->real = spix->real ;
      dpix->imag = -spix->imag ;
    }
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageConjugate: unsupported pixel format %d\n", 
              Isrc->pixel_format) ;
    break ;
  }

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
Pixelval
ImageAccum(IMAGE *Isrc)
{
  IMAGE    *Itiny ;
  Pixelval retval ;
  int      ecode, row, col, rows, cols, endrow, endcol ;
  float    ftotal, *fpix, real, imag ;
  CPIX     *cpix ;

#if 0
  Itiny = ImageAlloc(1,1,Isrc->pixel_format, 1) ;
  ecode = h_reduce(Isrc, Itiny, (float)Isrc->cols, (float)Isrc->rows) ;
#endif

  endrow = Isrc->frow + Isrc->rows ;
  endcol = Isrc->fcol + Isrc->cols ;
  switch (Isrc->pixel_format)
  {
  case PFCOMPLEX:
    real = imag = 0.0f ;
    for (row = Isrc->frow ; row < endrow ; row++)
    {
      cpix = IMAGECpix(Isrc, row, Isrc->fcol) ;
      for (col = Isrc->fcol ; col < endcol ; col++, cpix++)
      {
        real += cpix->real ;
        imag += cpix->imag ;
      }
    }
    retval.v_complex[REAL_PIX] = real ;
    retval.v_complex[IMAG_PIX] = imag ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, "ImageAccum: unsupported pixel format %d\n",
              Isrc->pixel_format) ;
  }

  return(retval) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MATRIX *
ImageToMatrix(IMAGE *I)
{
  MATRIX *mat ;
  int    format, bytes ;

  switch (I->pixel_format)
  {
  case PFCOMPLEX:
    format = MATRIX_COMPLEX ;
    break ;
  case PFFLOAT:
    format = MATRIX_REAL ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, "ImageToMatrix: unsupported image type %d",
              I->pixel_format) ;
    break ;
  }

  mat = MatrixAlloc(I->rows, I->cols, format) ;
  bytes = mat->rows * mat->cols * sizeof(float) ;
  if (mat->type == MATRIX_COMPLEX)
    bytes *= 2 ;

  memcpy(mat->data, I->image, bytes) ; 
  return(mat) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageFromMatrix(MATRIX *matrix, IMAGE *I)
{
  int    format, bytes ;

  format = (matrix->type == MATRIX_COMPLEX) ? PFCOMPLEX : PFFLOAT ;

  if (!I)
    I = ImageAlloc(matrix->rows, matrix->cols, format, 1) ;
  else
    if (I->rows != matrix->rows || I->cols != matrix->cols)
      ErrorExit(ERROR_BADPARM, "ImageFromMatrix: size mismatch") ;

  bytes = matrix->rows * matrix->cols * sizeof(float) ;
  if (matrix->type == MATRIX_COMPLEX)
    bytes *= 2 ;

  memcpy(I->image, matrix->data, bytes) ; 
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageInverse(IMAGE *Isrc, IMAGE *Idst)
{
  MATRIX *mat, *mat_inverse ;

  mat = ImageToMatrix(Isrc) ;
  mat_inverse = MatrixInverse(mat, NULL) ;
  Idst = ImageFromMatrix(mat_inverse, Idst) ;
  MatrixFree(&mat) ;
  MatrixFree(&mat_inverse) ;
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageMatrixMul(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst)
{
  MATRIX *mat1, *mat2, *mat_dst ;

  if (Isrc2->rows != Isrc1->cols)
    ErrorExit(ERROR_BADPARM, 
              "ImageMatrixMul: inner dimensions must agree (%d, %d)",
              Isrc1->cols, Isrc2->rows) ;

  mat1 = ImageToMatrix(Isrc1) ;
  mat2 = ImageToMatrix(Isrc2) ;
  mat_dst = MatrixMultiply(mat1, mat2, NULL) ;
  Idst = ImageFromMatrix(mat_dst, Idst) ;
  MatrixFree(&mat1) ;
  MatrixFree(&mat2) ;
  MatrixFree(&mat_dst) ;

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageType(char *fname)
{
  char *dot, buf[200] ;

  strcpy(buf, fname) ;
  dot = strrchr(buf, '.') ;

  if (dot)
  {
    dot = StrUpper(dot+1) ;
    if (!strcmp(dot, "MAT"))
      return(MATLAB_IMAGE) ;
  }

  return(HIPS_IMAGE) ;
} 
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageFrame(char *fname)
{
  char *colon, buf[200] ;
  int   frame ;

  strcpy(buf, fname) ;
  colon = strrchr(buf, ':') ;

  if (colon)
  {
    sscanf(colon+1, "%d", &frame) ;
    *colon = 0 ;
  }
  else
    frame = 0 ;

  return(frame) ;
} 
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageScale(IMAGE *Isrc, IMAGE *Idst, float new_min, float new_max)
{
  float  scale, old_min, old_max ;
  int    ecode, imin, imax ;
  byte   bmin, bmax ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  switch (Isrc->pixel_format)
  {
  case PFBYTE:
    ecode = h_minmax(Isrc, &bmin, &bmax, 0) ;
    old_min = (float)bmin ;
    old_max = (float)bmax ;
    break ;
  case PFINT:
    ecode = h_minmax(Isrc, &imin, &imax, 0) ;
    old_min = (float)imin ;
    old_max = (float)imax ;
    break ;
  case PFFLOAT:
    ecode = h_minmax(Isrc, &old_min, &old_max, 0) ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageScale: unsupported pixel format %d\n",
              Isrc->pixel_format) ;
    break ;
  }

  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageScale: h_minmax failed (%d)\n", ecode) ;

  scale = (new_max - new_min) / (old_max - old_min) ;
  ecode = h_linscale(Isrc, Idst, scale, new_min - old_min * scale) ;
  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageScale: h_linscale failed (%d)\n", ecode) ;

  ecode = h_minmax(Isrc, &old_min, &old_max, 0) ;
  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ImageCheckSize(IMAGE *inImage,IMAGE *outImage, int rows, int cols, int nframes)
{
  int  inPix, outPix ;

  if (!outImage)
    return(0) ;

  if (!rows)
    rows = inImage->rows ;
  if (!cols)
    cols = inImage->cols ;
  if (!nframes)
    nframes = inImage->num_frame ;

  inPix = rows * cols * nframes * inImage->sizepix ;
  outPix = outImage->numpix * outImage->sizepix * outImage->num_frame ;

  return(outPix >= inPix) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageCopyFrames(IMAGE *inImage, IMAGE *outImage,int start, int nframes, 
               int dst_frame)
{
  UCHAR *cIn, *cOut ;
  UINT  *iIn, *iOut ;
  float *fsrc, *fdst ;
  int   size, frameno, pix_per_frame, needed, end ;

  if (!ImageCheckSize(inImage, outImage, 0, 0, dst_frame + nframes))
  {
    fprintf(stderr, "ImageCopyFrames: outImage not large enough\n") ;
    return(-1) ;
  }

  end = start + nframes - 1 ;
  pix_per_frame = inImage->rows * inImage->cols ;
  for (frameno = start ; frameno <= end ; frameno++)
  {
    size = inImage->rows * inImage->cols ;
    switch (inImage->pixel_format)
    {
    case PFFLOAT:
      if (outImage->pixel_format != PFFLOAT)
        ErrorExit(ERROR_UNSUPPORTED, 
                  "ImageCopyFrames: unsupported image pixel format %d -> %d\n",
                  inImage->pixel_format, outImage->pixel_format) ;

      fsrc = IMAGEFseq_pix(inImage, 0, 0, frameno) ;
      fdst = IMAGEFseq_pix(outImage, 0, 0, dst_frame+frameno-start) ;
      memcpy((char *)fdst, (char *)fsrc, pix_per_frame*sizeof(float)) ;
      break ;
    case PFBYTE:
      if (outImage->pixel_format == PFBYTE)
        memcpy(IMAGEseq_pix(inImage,0,0,frameno),
               IMAGEseq_pix(outImage, 0, 0,frameno),
               pix_per_frame * sizeof(char)) ;
      else
      {
        size = inImage->rows * inImage->cols ;
        cIn = IMAGEseq_pix(inImage, 0, 0, frameno) ;
        switch (outImage->pixel_format)
        {
        case PFFLOAT:
          fdst = IMAGEFseq_pix(outImage, 0, 0, frameno) ;
          while (size--)
            *fdst++ = (float)*cIn++ ;
          break ;
        case PFBYTE:
          cOut = IMAGEseq_pix(outImage, 0, 0, frameno) ;
          while (size--)
            *cOut++ = *cIn++ ;
          break ;
        case PFINT:
          iOut = IMAGEIseq_pix(outImage, 0, 0, frameno) ;
          while (size--)
            *iOut++ = (UINT)*cIn++ ;
          break ;
      default:
          ErrorExit(ERROR_UNSUPPORTED, 
                    "ImageCopyFrames: unsupported image pixel format %d -> %d\n"
                    , inImage->pixel_format, outImage->pixel_format) ;
          break ;
        }
        break ;
      }
    
    case PFINT:
      iIn = IMAGEIpix(inImage, 0, 0) + pix_per_frame * frameno ;
      switch (outImage->pixel_format)
      {
      case PFINT:
        iOut = IMAGEIpix(outImage, 0, 0) + pix_per_frame * frameno ;
        memcpy((char *)iOut, (char *)iIn, pix_per_frame*sizeof(int)) ;
        break ;
      case PFBYTE:
        cOut = IMAGEpix(outImage, 0, 0) + pix_per_frame * frameno ;
        while (size--)
          *cOut++ = (UCHAR)*iIn++ ;
        break ;
      default:
        ErrorExit(ERROR_UNSUPPORTED, 
                  "ImageCopyFrames: unsupported image pixel format %d -> %d\n",
                  inImage->pixel_format, outImage->pixel_format) ;
        break ;
      }
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED, 
                "ImageCopyFrames: unsupported image pixel format %d -> %d\n",
                inImage->pixel_format, outImage->pixel_format) ;
      break ;
    }
  }

  outImage->rows = inImage->rows ;
  outImage->cols = inImage->cols ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int  
ImageScaleRange(IMAGE *image, float fmin, float fmax, int low, int high)
{
  int    min_val, max_val, size, val ;
  UCHAR *csrc, cmin_val, cmax_val, cval ;
  int   *isrc, imin_val, imax_val, ival ;
  float *fsrc, fmin_val, fmax_val, fval, norm ;

  size = image->cols * image->rows ;
  switch (image->pixel_format)
  {
  case PFBYTE:
    cmax_val = (UCHAR)fmax ;
    cmin_val = (UCHAR)fmin ;
    size = image->cols * image->rows ;
    csrc = IMAGEpix(image, 0, 0) ;
    norm = ((float)high - (float)low) / 
      ((float)cmax_val - (float)cmin_val) ;
    while (size--)
    {
      cval = *csrc ;
      cval = (UCHAR)((float)((float)cval - (float)cmin_val) * norm) + low ;
      *csrc++ = cval ;
    }
    break ;

  case PFINT:
    imax_val = (int)fmax ;
    imin_val = (int)fmin ;
    size = image->cols * image->rows ;
    isrc = (int *)IMAGEIpix(image, 0, 0) ;
    norm = ((float)high - (float)low) / ((float)imax_val-(float)imin_val);
    while (size--)
    {
      ival = *isrc ;
      ival = (int)((float)((float)ival - (float)imin_val) * norm) + low ;
      *isrc++ = ival ;
    }
    break ;

  case PFFLOAT:
    size = image->cols * image->rows ;
    fsrc = (float *)IMAGEFpix(image, 0, 0) ;
    norm = ((float)high - (float)low) / (fmax-fmin);
    while (size--)
    {
      fval = *fsrc ;
      fval = ((fval - fmin)*norm) + (float)low ;
      *fsrc++ = fval ;
    }
    break ;

  default:
    fprintf(stderr, "ImageScale: unsupported format %d\n", image->pixel_format);
    exit(1) ;
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:
              fname - the name of the file to read from

           Description:
              read a hips image from a file, and allocate an image
              header and data space for it.  Returns the newly
              allocated image.
----------------------------------------------------------------------*/
int
ImageReadInto(char *fname, IMAGE *I, int image_no)
{
  FILE   *fp ;
  int    ecode ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadInto(%s) failed\n", fname)) ;
  
  ecode = fread_header(fp, I, fname) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageReadInto(%s): fread_header failed (%d)\n", 
              fname, ecode) ;
  ecode = fread_image(fp, I, image_no, fname) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageReadInto(%s): fread_image failed (%d)\n", 
              fname, ecode) ;
  
  fclose(fp) ;
  
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:
              image - the image to write
              fname - the name of the file to write to.

           Description:
              write a hips image to file 'fname'
----------------------------------------------------------------------*/
int
ImageWriteFrames(IMAGE *image, char *fname, int start, int nframes)
{
  IMAGE  *tmp_image ;
  
  tmp_image = ImageAlloc(image->cols, image->rows,image->pixel_format,nframes);
  ImageCopyFrames(image, tmp_image, start, nframes, 0) ;
  ImageWrite(tmp_image, fname) ;
  ImageFree(&tmp_image) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageRescale(IMAGE *inImage, IMAGE *outImage, float scale)
{
  if (!outImage)
    outImage = ImageAlloc(inImage->rows, inImage->cols, inImage->pixel_format, 
                          inImage->num_frame) ;

  if (scale == 1)
    ImageCopy(inImage, outImage) ;
  else if (scale > 1)
    ImageScaleUp(inImage, outImage, scale) ;
  else
    ImageScaleDown(inImage, outImage, scale) ;

  return(outImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageScaleDown(IMAGE *inImage, IMAGE *outImage, float scale)
{
  int    inRow, inCol, outRow, outCol, inCols, inRows, endCol, endRow,
         outRows, outCols ;
  UCHAR  *inPix, *outPix ;
  UINT   *inIPix, *outIPix ;
  float  *finPix, *foutPix ;

  if (!ImageCheckSize(inImage, outImage, inImage->rows*scale,
                          inImage->cols*scale, 0))
  {
    fprintf(stderr, "ImageScaleDown: output image not big enough\n") ;
    return(-1) ;
  }

  outImage->cols = inImage->cols * scale ;
  outImage->rows = inImage->rows * scale ;
  inRows = inImage->rows ;
  inCols = inImage->cols ;
  outRows = outImage->rows ;
  outCols = outImage->cols ;

  switch (inImage->pixel_format)
  {
  case PFBYTE:
    switch (outImage->pixel_format)
    {
    case PFFLOAT:  /* byte --> float */
      foutPix = IMAGEFpix(outImage, 0, 0) ;
      for (outRow = 0 ; outRow < outRows ; outRow++)
        for (outCol = 0 ; outCol < outCols ; outCol++, foutPix++)
        {
          /* map center point to this output point */
          inRow = outRow / scale ;
          inCol = outCol / scale ;
          *foutPix = (float)*IMAGEpix(inImage, inCol, inRow) ;
        }
      break ;
    case PFBYTE:  /* byte --> byte */
      outPix = IMAGEpix(outImage, 0, 0) ;
      for (outRow = 0 ; outRow < outRows ; outRow++)
        for (outCol = 0 ; outCol < outCols ; outCol++, outPix++)
        {
          /* map center point to this output point */
          inRow = outRow / scale ;
          inCol = outCol / scale ;
          if (inRow >= inRows || outRow >= outRows ||
              inCol >= inCols || outCol >= outCols)
          {
            fprintf(stderr, "in: %d, %d --> out: %d, %d!\n",
                    inRow, inCol, outRow, outCol) ;
            exit(2) ;
          }
          *outPix = *IMAGEpix(inImage, inCol, inRow) ;
        }
      break ;
    default:
      fprintf(stderr, "ImageScaleDown: unsupported output pixel format %d\n", 
              outImage->pixel_format) ;
      return(-1) ;
      break ;
    }
    break ;
  case PFFLOAT:   /* float --> byte */
    switch (outImage->pixel_format)
    {
    case PFBYTE:
      outPix = IMAGEpix(outImage, 0, 0) ;
      for (outRow = 0 ; outRow < outRows ; outRow++)
        for (outCol = 0 ; outCol < outCols ; outCol++, outPix++)
        {
          /* map center point to this output point */
          inRow = outRow / scale ;
          inCol = outCol / scale ;
          *outPix = (UCHAR)*IMAGEFpix(inImage, inCol, inRow) ;
        }
      break ;
    case PFFLOAT:
      foutPix = IMAGEFpix(outImage, 0, 0) ;
      for (outRow = 0 ; outRow < outRows ; outRow++)
        for (outCol = 0 ; outCol < outCols ; outCol++, foutPix++)
        {
          /* map center point to this output point */
          inRow = outRow / scale ;
          inCol = outCol / scale ;
          *foutPix = (float)*IMAGEFpix(inImage, inCol, inRow) ;
        }
      break ;
    default:
      fprintf(stderr, "ImageScaleDown: unsupported output pixel format %d\n", 
              outImage->pixel_format) ;
      return(-1) ;
      break ;
    }
    break ;
  case PFINT:

  default:
    fprintf(stderr, "ImageScaleDown: unsupported pixel format %d\n", 
            inImage->pixel_format) ;
    return(-2) ;
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageScaleUp(IMAGE *inImage, IMAGE *outImage, float scale)
{
  int    needed, inRow, inCol, outRow, outCol, inCols, inRows, endCol, endRow,
         outRows, outCols ;
  UCHAR  *inPix, *outPix ;
  UINT   *inIPix, *outIPix ;
  float  *finPix, *foutPix ;

  if (!ImageCheckSize(inImage, outImage, nint(inImage->rows*scale),
                          nint(inImage->cols*scale), 0))
  {
    fprintf(stderr, 
          "ImageScaleUp: output image not large enough %d x %d -> %d x %d\n",
            inImage->rows, inImage->cols, outImage->rows, outImage->cols);
    return(-1) ;
  }

  outCols = outImage->cols = nint((float)inImage->cols * scale) ;
  outRows = outImage->rows = nint((float)inImage->rows * scale) ;

  inRows = inImage->rows ;
  inCols = inImage->cols ;

  switch (inImage->pixel_format)
  {
  case PFBYTE:
    switch (outImage->pixel_format)
    {
    case PFBYTE:
      outPix = IMAGEpix(outImage, 0, 0) ;
      for (outRow = 0 ; outRow < outRows ; outRow++)
      {
        for (outCol = 0 ; outCol < outCols ; outCol++)
        {
          inCol = (int)((float)outCol / scale) ;
          inRow = (int)((float)outRow / scale) ;
          inPix = IMAGEpix(inImage, inCol, inRow) ;
          *outPix++ = *inPix ;
        }
      }
      break ;
    case PFFLOAT:
      inPix = IMAGEpix(inImage, 0, 0) ;
      for (inRow = 0 ; inRow < inRows ; inRow++)
        for (inCol = 0 ; inCol < inCols ; inCol++, inPix++)
        {
          /* fill in a scale x scale area in the output image */
          endRow = inRow * scale + scale ;
          endCol = inCol * scale + scale ;
          for (outRow = inRow * scale ; outRow < endRow ; outRow++)
          {
            foutPix = IMAGEFpix(outImage, nint((float)inCol * scale), outRow) ;
            
            for (outCol = inCol * scale ; outCol < endCol ; outCol++,foutPix++)
              *foutPix = (float)(*inPix) ;
          }
        }
      break ;
    default:
      fprintf(stderr, "ImageScaleUp: unsupported output pixel format %d\n", 
              outImage->pixel_format) ;
      return(-1) ;
      break ;
    }
    break ;
  case PFFLOAT:
    switch (outImage->pixel_format)
    {
    case PFBYTE:
      finPix = IMAGEFpix(inImage, 0, 0) ;
      for (inRow = 0 ; inRow < inRows ; inRow++)
        for (inCol = 0 ; inCol < inCols ; inCol++, finPix++)
        {
          /* fill in a scale x scale area in the output image */
          endRow = inRow * scale + scale ;
          endCol = inCol * scale + scale ;
          for (outRow = inRow * scale ; outRow < endRow ; outRow++)
          {
            outPix = IMAGEpix(outImage, nint((float)inCol * scale), outRow) ;
            
            for (outCol = inCol * scale ; outCol < endCol ; outCol++, outPix++)
              *outPix = (UCHAR)(*finPix) ;
          }
        }
      break ;
    case PFFLOAT:
      finPix = IMAGEFpix(inImage, 0, 0) ;
      for (inRow = 0 ; inRow < inRows ; inRow++)
        for (inCol = 0 ; inCol < inCols ; inCol++, finPix++)
        {
          /* fill in a scale x scale area in the output image */
          endRow = inRow * scale + scale ;
          endCol = inCol * scale + scale ;
          for (outRow = inRow * scale ; outRow < endRow ; outRow++)
          {
            foutPix = IMAGEFpix(outImage, nint((float)inCol * scale), outRow) ;
            
            for (outCol = inCol * scale ; outCol < endCol ; outCol++,foutPix++)
              *foutPix = *finPix ;
          }
        }
      break ;
    default:
      fprintf(stderr, "ImageScaleUp: unsupported output pixel format %d\n", 
              outImage->pixel_format) ;
      return(-1) ;
      break ;
    }
    break ;
  case PFINT:
    inIPix = IMAGEIpix(inImage, 0, 0) ;
    inRows = inImage->rows ;
    inCols = inImage->cols ;
    for (inRow = 0 ; inRow < inRows ; inRow++)
      for (inCol = 0 ; inCol < inCols ; inCol++, inIPix++)
      {
        /* fill in a scale x scale area in the output image */
        endRow = inRow * scale + scale ;
        endCol = inCol * scale + scale ;
        for (outRow = inRow * scale ; outRow < endRow ; outRow++)
        {
          outIPix = IMAGEIpix(outImage, nint((float)inCol * scale), outRow) ;
          
          for (outCol = inCol * scale ; outCol < endCol ; outCol++, outIPix++)
            *outIPix = *inIPix ;
        }
      }
    break ;
  default:
    fprintf(stderr, "ImageScaleUp: unsupported pixel format %d\n", 
            inImage->pixel_format) ;
    return(-2) ;
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageReflect(IMAGE *inImage, IMAGE *outImage, int how)
{
  int           x, y, ymax ;
  unsigned char *src, *dst, *tmp ;
  unsigned int  *isrc, *idst, *itmp ;

  if (!ImageCheckSize(inImage, outImage, 0, 0, 0))
  {
    fprintf(stderr, "ImageReflect: output image not large enough\n") ;
    return(-1) ;
  }
  outImage->rows = inImage->rows ;
  outImage->cols = inImage->cols ;
  switch(inImage->pixel_format)
  {
  case PFBYTE:
    switch (how)
    {
    case IMAGE_REFLECT_AROUND_X_AXIS:
      ymax = inImage->rows - 1 ;
      src = inImage->image ;
      tmp = outImage->image ;
      for (y = 0 ; y < inImage->rows ; y++)
      {
        for (x = 0 ; x < inImage->cols ; x++)
        {
          dst = IMAGEpix(outImage, x, ymax - y) ;
          *dst = *src++ ;
        }
      }
      break ;
    case IMAGE_REFLECT_AROUND_Y_AXIS:
      break ;
    default:
      fprintf(stderr, "ImageReflect: unknown how parm (%d)\n", how) ;
      exit(1) ;
    }
    break ;
  case PFINT:
    switch (how)
    {
    case IMAGE_REFLECT_AROUND_X_AXIS:
      ymax = inImage->rows - 1 ;
      isrc = (unsigned int *)inImage->image ;
      itmp = (unsigned int *)outImage->image ;
      for (y = 0 ; y < inImage->rows ; y++)
      {
        for (x = 0 ; x < inImage->cols ; x++)
        {
          idst = IMAGEIpix(outImage, x, ymax - y) ;
          *idst = *isrc++ ;
        }
      }
      break ;
    case IMAGE_REFLECT_AROUND_Y_AXIS:
      break ;
    default:
      fprintf(stderr, "ImageReflect: unknown how parm (%d)\n", how) ;
      exit(1) ;
    }
    break ;

  default:
    fprintf(stderr, "ImageReflect: unsupported image format %d\n", 
            inImage->pixel_format) ;
    break ;
  }

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              add multiplicative "speckle" noise to an image.
----------------------------------------------------------------------*/
int
ImageAddSpeckleNoise(IMAGE *inImage,IMAGE *outImage, float amp)
{
  long   idum = 1.0l ;
  int    npix ;
  float  *inPix, *outPix, gnoise ;

  if (inImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported input format %d\n",
            inImage->pixel_format) ;
    return(-1) ;
  }
  if (outImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported output format %d\n",
            outImage->pixel_format) ;
    return(-1) ;
  }

  npix = inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = randomNumber(1.0f-amp, 1.0f+amp) ;
    *outPix++ += *inPix++ * gnoise ;
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              corrupt an image with salt & pepper noise: randomly 
              generated 0s and 1s.
----------------------------------------------------------------------*/
int
ImageAddSaltNoise(IMAGE *inImage,IMAGE *outImage, float density)
{
  long   idum = 1.0l ;
  int    npix ;
  float  *inPix, *outPix, gnoise, in ;

  if (inImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported input format %d\n",
            inImage->pixel_format) ;
    return(-1) ;
  }
  if (outImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported output format %d\n",
            outImage->pixel_format) ;
    return(-1) ;
  }

  npix = inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = randomNumber(0.0f, 1.0f) ;
    in = *inPix++ ;
    if (gnoise < density)
    {
      if (gnoise < density/2.0f)
        in = 0.0f ;
      else
        in = 1.0f ;
    }
    *outPix++ = in ;
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             corrupt an image with additive zero mean gaussian noise.
----------------------------------------------------------------------*/
int
ImageAddNoise(IMAGE *inImage, IMAGE *outImage, float amp)
{
  long   idum = 1.0l ;
  int    npix ;
  float  *inPix, *outPix, gnoise ;

  if (inImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported input format %d\n",
            inImage->pixel_format) ;
    return(-1) ;
  }
  if (outImage->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageAddNoise: unsupported output format %d\n",
            outImage->pixel_format) ;
    return(-1) ;
  }

  npix = inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = randomNumber(-amp, amp) ;
    *outPix++ = *inPix++ + gnoise ;
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int  
ImageValRange(IMAGE *image, void *pvmin, void *pvmax)
{
  float  *pfmax, *pfmin, fmax, fmin, *fpix ;
  int    size ;

  size = image->rows * image->cols * image->num_frame ;
  switch (image->pixel_format)
  {
  case PFFLOAT:
    pfmax = (float *)pvmax ;
    pfmin = (float *)pvmin ;
    fpix = IMAGEFpix(image, 0, 0) ;
    fmax = fmin = *fpix ;
    while (size--)
    {
      if (*fpix > fmax)
        fmax = *fpix ;
      if (*fpix < fmin)
        fmin = *fpix ;
      fpix++ ;
    }
    *pfmax = fmax ;
    *pfmin = fmin ;
    break ;
  default:
    fprintf(stderr, "ImageValRange: unsupported pixel format %d\n",
            image->pixel_format) ;
    return(-1);
  }

    
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageCatSeq(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst)
{
  IMAGE *Itmp ;
  int   num_frame, frameno ;

  if ((Isrc1->rows != Isrc2->rows) || (Isrc1->cols != Isrc2->cols))
    return(NULL) ;

  if (!Idst)
    Idst = ImageAlloc(Isrc1->rows, Isrc1->cols, Isrc1->pixel_format,
                      Isrc1->num_frame + Isrc2->num_frame) ;

  num_frame = Isrc2->num_frame ;
  if (Isrc1)
    num_frame += Isrc1->num_frame ;
  Itmp = ImageAlloc(Isrc2->rows, Isrc2->cols, Isrc2->pixel_format, num_frame) ;
  if (Isrc1)
  {
    ImageCopyFrames(Isrc1, Itmp, 0, Isrc1->num_frame, 0) ;
    frameno = Isrc1->num_frame ;
  }
  else
    frameno = 0 ;
  ImageCopyFrames(Isrc2, Itmp, 0, Isrc2->num_frame, frameno) ;

#if 0
  if (frameno > 0)
    ImageWrite(Itmp, "Itmp.hipl") ;
  if (frameno > 0)
    ImageWrite(Isrc2, "Isrc2.hipl") ;
{
  IMAGE *Itmp3 ;

  static int n = 0 ;
  char fname[100] ;
  sprintf(fname, "new%d.hipl", n) ;
  ImageWrite(Isrc2, fname) ;
  sprintf(fname, "add%d.hipl", n++) ;
  ImageWrite(Itmp, fname) ;

  Itmp3 = ImageAlloc(Itmp->rows, Itmp->cols, Itmp->pixel_format, 1) ;
  ImageCopyFrames(Itmp, Itmp3, Itmp->num_frame-1, 1, 0) ;
  ImageWrite(Itmp3, "Ilast.hipl") ;
  ImageFree(&Itmp3) ;
}
#endif

  if (Isrc1)
    ImageFree(&Isrc1) ;
  return(Itmp) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageMulScale(IMAGE *Isrc, IMAGE *Idst, Pixelval *p)
{
  int    ecode, size ;
  float  real, imag, sreal, simag ;
  CPIX   *csrc, *cdst ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 1);

  switch (Isrc->pixel_format)
  {
  case PFCOMPLEX:
    csrc = IMAGECpix(Isrc, 0, 0) ;
    cdst = IMAGECpix(Idst, 0, 0) ;
    real = p->v_complex[REAL_PIX] ;
    imag = p->v_complex[IMAG_PIX] ;
    size = Isrc->numpix ;
    while (size--)
    {
      simag = csrc->imag ;
      sreal = csrc->real ;
      cdst->real = real * sreal - imag * simag ;
      cdst->imag = real * simag + sreal * imag ;
      csrc++ ;
      cdst++ ;
    }
    break ;
  default:
    ecode = h_mulscale(Isrc, Idst, p) ;
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageMulScale: h_mulscale failed (%d)", ecode) ;
    break ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageAddScalar(IMAGE *Isrc, IMAGE *Idst, float scalar)
{
  int    size ;
  float  *fpix ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame);

  switch (Isrc->pixel_format)
  {
  case PFFLOAT:
    size = Isrc->numpix * Isrc->num_frame ;
    fpix = IMAGEFpix(Isrc, 0, 0) ;
    while (size--)
      *fpix++ += scalar ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, "ImageAddScalar: unsupported pixel type %d",
              Isrc->pixel_format) ;
    break ;
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE    *
ImageReplace(IMAGE *Isrc, IMAGE *Idst, float inpix, float outpix)
{
  float  *fin, *fout ;
  char   *cin, *cout, cinpix, coutpix ;
  int    npix ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, Isrc->num_frame) ;

  if (Idst->pixel_format != Isrc->pixel_format)
    ErrorReturn(NULL, (ERROR_BADPARM, "ImageReplace: src and dst formats must match")) ;

  npix = Isrc->numpix * Isrc->num_frame ;
  switch (Isrc->pixel_format)
  {
  case PFFLOAT:
    fin = IMAGEFpix(Isrc, 0, 0) ;
    fout = IMAGEFpix(Idst, 0, 0) ;
    while (npix--)
    {
      if (*fin == inpix)
      {
        *fout++ = outpix ;
        fin++;
      }
      else
        *fout++ = *fin++ ;
    }
    break ;
  case PFBYTE:
    cinpix = (char)inpix ;
    coutpix = (char)outpix ;
    cin = IMAGEpix(Isrc, 0, 0) ;
    cout = IMAGEpix(Idst, 0, 0) ;
    while (npix--)
    {
      if (*cin == cinpix)
      {
        *cout++ = coutpix ;
        cin++;
      }
      else
        *cout++ = *cin++ ;
    }
    break ;
  default:
    ErrorReturn(NULL, 
                (ERROR_UNSUPPORTED, "ImageReplace: unsupported pixel format %d", Isrc->pixel_format));
    break ;
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              decompose a file name, extracting the type and the frame #.
----------------------------------------------------------------------*/
int
ImageUnpackFileName(char *inFname, int *pframe, int *ptype, char *outFname)
{
  char *colon, *dot, buf[100] ;

  strcpy(outFname, inFname) ;
  colon = strrchr(outFname, ':') ;
  dot = strchr(outFname, '.') ;

  if (colon)   /* : in filename indicates frame # */
  {
    if (sscanf(colon+1, "%d", pframe) < 1)
      *pframe = 0 ;
    *colon = 0 ;
  }
  else
    *pframe = 0 ;

  if (dot)
  {
    dot = StrUpper(strcpy(buf, dot+1)) ;
    if (!strcmp(dot, "MAT"))
      *ptype = MATLAB_IMAGE ;
    else
      *ptype = HIPS_IMAGE ;
  }
  else
    *ptype = HIPS_IMAGE ;

  return(NO_ERROR) ;
}

