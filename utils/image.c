/*
 *       FILE NAME:   image.c
 *
 *       DESCRIPTION: image processing utilities
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
 * replaced alloc_image with alloc_image_buffer (dng, 3/12/96).
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
 #include <unistd.h>  /* for SEEK_ constants */

#include "hmem.h"
#include <hipl_format.h>

#include "hips.h"

#include "image.h"
#include "error.h"
#include "matrix.h"
#include "matfile.h"
#include "utils.h"
#include "macros.h"
#include "machine.h"
#include "proto.h"
#include "diag.h"
#include "canny.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/
static int alloc_image_buffer(struct header *hd) ;
static void break_here(void) ;
static void break_here(void){ fprintf(stderr, "break point encountered!\n") ;}

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
  int   ecode ;
  
  I = (IMAGE *)calloc(1, sizeof(IMAGE)) ;
  if (!I)
    ErrorExit(ERROR_NO_MEMORY, "ImageAlloc: could not allocate header\n") ;

  init_header(I, "orig", "seq", nframes, "today", rows,cols,format,1, "temp");
  ecode = alloc_image_buffer(I) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_MEMORY, "ImageAlloc: could not allocate %dx%d buffer\n",
              rows,cols) ;
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           stolen from hips2 code and modified to allocate multiple frames.
------------------------------------------------------*/
static int
alloc_image_buffer(struct header *hd)
{
  int fcb,cb;
  long  npix ;

  if (hd->sizeimage == (hsize_t) 0)  /*dng*/
  {
    hd->imdealloc = (h_boolean) FALSE; /*dng*/
    return(HIPS_OK);
  }
  npix = (long)hd->sizeimage*(long)hd->num_frame ;
  if((hd->image = hcalloc(npix, sizeof(byte))) == (byte *)NULL)    
    return(HIPS_ERROR);
  if (hd->pixel_format == PFMSBF || hd->pixel_format == PFLSBF) 
  {
    fcb = hd->fcol/8;
    cb = (hd->ocols + 7)/8;
    hd->firstpix = hd->image + ((cb * hd->frow) + fcb);
  }
  else
    hd->firstpix = hd->image +
      (((long)hd->ocols * (long)hd->frow) + (long)hd->fcol) * hd->sizepix;
  hd->imdealloc = TRUE;
  hmemset(hd->image, 0, hd->sizeimage*hd->num_frame) ;
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
    ErrorReturn(-1, (ERROR_NO_FILE, "ImageWrite(%s) failed\n", fname)) ;

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
  byte *image ;

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
int
ImageUpdateHeader(IMAGE *I, char *fname)
{
  FILE  *fp ;
  int   ecode ;

  fp = fopen(fname, "r+b") ;
  if (!fp)
    ErrorReturn(-1, (ERROR_NO_FILE, "ImageUpdateHeader(%s) failed\n", fname)) ;

  ecode = fwrite_header(fp, I, fname) ;
  if (ecode != HIPS_OK)
    ErrorReturn(-1, (ERROR_NO_FILE, 
                  "ImageUpdateHeader: fwrite_header failed (%d)\n",ecode)) ;

  fclose(fp) ;
  return(0) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageFRead(FILE *fp, char *fname, int start, int nframes)
{
  int    ecode, end_frame, frame ;
  IMAGE  *I ;
  byte   *startpix ;

  if (!fname)
    fname = "ImageFRead" ;

  I = (IMAGE *)calloc(1, sizeof(IMAGE)) ;
  if (!I)
    ErrorExit(ERROR_NO_MEMORY,"ImageFRead: could not allocate header\n") ;

  ecode = fread_header(fp, I, fname) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_header failed (%d)\n",ecode);

  if (start < 0)    /* read all frames */
  {
    start = 0 ;
    nframes = I->num_frame ;
  }
  else              /* read only specified frames */
  {
    if (fseek(fp, (long)I->sizeimage*(long)start, SEEK_CUR) < 0)
    {
      ImageFree(&I) ;
      ErrorReturn(NULL, 
                  (ERROR_BADFILE, 
                   "ImageFRead(%s, %d) - could not seek to specified frame",
                   fname, start)) ;
    }
  }

  if (nframes < 0)
    nframes = I->num_frame - start + 1 ;

  end_frame = start + nframes - 1 ;
  if (end_frame >= I->num_frame)
    ErrorReturn(NULL, 
                (ERROR_BADFILE,
                 "ImageFRead(%s, %d) - frame out of bounds", fname,end_frame));
  I->num_frame = nframes ;
  alloc_image_buffer(I) ;

  startpix = I->image ;
  for (frame = start ; frame <= end_frame ; frame++)
  {
    ecode = fread_image(fp, I, frame, fname) ;
    if (ecode != HIPS_OK)
      ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_image failed (%d)\n", ecode);
    I->image += I->sizeimage ;
  }
  I->image = startpix ;

  if (!ImageValid(I))
  {
    DCPIX  *dcpix, dcval ;
    CPIX   *cpix, cval ;
    double *dpix, dval ;
    float  *fpix, fval ;
    long   npix ;

    npix = (long)I->numpix * (long)I->num_frame ;

    switch (I->pixel_format)
    {
    case PFDBLCOM:
      dcpix = IMAGEDCpix(I, 0, 0) ;
      while (npix--)
      {
        dcval = *dcpix ;
        dcpix->real = swapDouble(dcval.real) ;
        dcpix->imag = swapDouble(dcval.imag) ;
        dcpix++ ;
      }
      break ;
    case PFCOMPLEX:
      cpix = IMAGECpix(I, 0, 0) ;
      while (npix--)
      {
        cval = *cpix ;
        cpix->real = swapDouble(cval.real) ;
        cpix->imag = swapDouble(cval.imag) ;
        cpix++ ;
      }
      break ;
    case PFDOUBLE:
      dpix = IMAGEDpix(I, 0, 0) ;
      while (npix--)
      {
        dval = *dpix ;
        *dpix++ = swapDouble(dval) ;
      }
      break ;
    case PFFLOAT:
      fpix = IMAGEFpix(I, 0, 0) ;
      while (npix--)
      {
        fval = *fpix ;
        *fpix++ = swapFloat(fval) ;
      }
      break ;
    default:
      ErrorReturn(NULL, 
                  (ERROR_UNSUPPORTED, 
                   "ImageFRead: unsupported type %d\n", I->pixel_format)) ;
    }
  }
  
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageReadFrames(char *fname, int start, int nframes)
{
  IMAGE *I ;
  FILE  *fp ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadFrames(%s) fopen failed\n", 
                       fname)) ;

  I = ImageFRead(fp, fname, start, nframes) ;
  fclose(fp) ;
  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageReadHeader(char *fname)
{
  IMAGE   *I = NULL ;
  FILE    *fp ;
  int     type, frame ;
  char    buf[100] ;

  strcpy(buf, fname) ;   /* don't destroy callers string */
  fname = buf ;
  ImageUnpackFileName(fname, &frame, &type, fname) ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadHeader(%s, %d) failed\n", 
                       fname, frame)) ;

  I = ImageFReadHeader(fp, fname) ;
  fclose(fp) ;

  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageFReadHeader(FILE *fp, char *fname)
{
  IMAGE   *I = NULL ;
  int     ecode ;
  int     type, frame ;
  char    buf[100] ;

  strcpy(buf, fname) ;   /* don't destroy callers string */
  fname = buf ;
  ImageUnpackFileName(fname, &frame, &type, fname) ;

  I = (IMAGE *)calloc(1, sizeof(IMAGE)) ;
  if (!I)
    ErrorExit(ERROR_NO_MEMORY, "ImageReadHeader: could not allocate header\n");

  switch (type)
  {
  case MATLAB_IMAGE:
  {
    MATFILE mf ;
    
    MatReadHeader(fp, &mf) ;
    init_header(I, "matlab", "seq", 1, "today", (int)mf.mrows, (int)mf.ncols,
               mf.imagf ? PFCOMPLEX : PFFLOAT, 1, "temp") ;
  }
    break ;
  case HIPS_IMAGE:
  default:
    ecode = fread_header(fp, I, fname) ;
    if (ecode != HIPS_OK)
    {
      fclose(fp) ;
      ErrorReturn(NULL, (ERROR_NO_FILE, 
                         "ImageReadHeader(%s): fread_header failed (%d)\n", 
                         fname, ecode)) ;
    }
    break ;
  }

  return(I) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          read an image from file and convert it to the specified
          format.
------------------------------------------------------*/
IMAGE *
ImageReadType(char *fname, int pixel_format)
{
  IMAGE *Itmp, *I ;

  Itmp = ImageRead(fname) ;          
  if (!Itmp)
    ErrorReturn(NULL, (ERROR_NO_FILE, 
      "ImageReadType(%s, %d): could not read image",
      fname, pixel_format)) ;
  if (Itmp->pixel_format != pixel_format)
  {
    I = ImageAlloc(Itmp->rows, Itmp->cols, pixel_format, Itmp->num_frame) ;
    ImageCopy(Itmp, I) ;
    ImageFree(&Itmp) ;
  }
  else
    I = Itmp ;

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
    DiagPrintf(DIAG_WRITE, 
               "ImageRead: fname=%s, frame=%d, type=%d (M=%d,H=%d)\n",
               fname, frame, type , MATLAB_IMAGE, HIPS_IMAGE);
    mat = MatlabRead(fname) ;
    if (!mat)
      ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s) failed\n", fname)) ;
    I = ImageFromMatrix(mat, NULL) ;
    ImageInvert(I, I) ;
    MatrixFree(&mat) ;
    break ;
  case HIPS_IMAGE:
    fp = fopen(fname, "rb") ;
    if (!fp)
      ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s, %d) failed\n", 
                         fname, frame)) ;
    I = ImageFRead(fp, fname, frame, 1) ;
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
IMAGE *
ImageThreshold(IMAGE *Isrc, IMAGE *Idst, float threshold)
{
  int      ecode ;
  Pixelval p ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols,Isrc->pixel_format,
                      Isrc->num_frame);

  p.v_float = threshold ;
  ecode = h_softthresh(Isrc, Idst, &p) ;
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
  /*float    loglen ;*/
  int      ecode ;
  IMAGE    *Itmp ;
  /*Pixelval p ;*/

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

  if (Isrc->pixel_format != PFCOMPLEX)
  {
    hips_rtocplx = CPLX_RVI0 ;
    ecode = h_toc(Isrc, Idst) ;
    if  (ecode != HIPS_OK)
    {
      ImageFree(&Idst) ;
      ErrorReturn(NULL, (ecode, "ImageDFT: h_toc failed (%d)\n", ecode)) ;
    }
  }
  else
    ImageCopy(Isrc, Idst) ;

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
  /*float    loglen ;*/
  IMAGE    *Itmp ;
  int      ecode ;
  Pixelval p ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
  Itmp = ImageAlloc(Isrc->rows, Isrc->cols, PFCOMPLEX, 1) ;

  ImageCopy(Isrc, Itmp) ;

#if 0
  loglen = log2(Isrc->rows) ;
  if ((Isrc->rows == Isrc->cols) && (floor(loglen) == loglen)) /* FFT */
  {
  }
  else   /* not a power of 2 - use DFT */
  {
  }
#endif
  ecode = h_invfourtr(Itmp) ;
  if  (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageInverseDFT: h_invfourtr failed (%d)\n", ecode) ;

  if (Idst->pixel_format != PFCOMPLEX)
  {
    hips_cplxtor = CPLX_REAL ;
    h_tof(Itmp, Idst) ;
  }
  else
    ImageCopy(Itmp, Idst) ;

#if 1
  p.v_float = (float)Idst->numpix ;
  ecode = h_divscale(Idst, Idst, &p) ;
  if  (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageInverseDFT: h_divscale failed (%d)\n", ecode) ;
#endif
  ImageFree(&Itmp) ;
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
  IMAGE *Itmp ;

  if (!Idst)
    Idst = ImageAlloc(drows, dcols, Isrc->pixel_format, Isrc->num_frame) ;

  x_scale = (float)dcols / (float)Isrc->cols ;
  y_scale = (float)drows / (float)Isrc->rows ;

  if (drows == dcols && Isrc->rows == Isrc->cols && (!ISINT(x_scale) ||
      !ISINT(y_scale)))
    return(ImageRescale(Isrc, Idst, (float)drows / (float)Isrc->rows)) ;

  if (!Idst)
    Idst = ImageAlloc(drows, dcols, Isrc->pixel_format, Isrc->num_frame) ;

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
        ecode = h_enlarge(Isrc, Idst, nint(x_scale), nint(y_scale) ) ;
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
      {
        float scale ;
        
        scale = 1.0f / x_scale ;
        if ((x_scale == y_scale) && ISPOW2(scale))
        {
          int   reductions, i ;
      
          reductions = (int)log2(scale) ;
          fprintf(stderr, "reducing %d times\n", reductions) ;
          for (i = 0 ; i < reductions ; i++)
          {
            Itmp = Isrc ;
            Isrc = ImageReduce(Itmp, NULL) ;
            if (i)   /* first one is real source image, don't free it */
              ImageFree(&Itmp) ;
          }
          ImageCopy(Isrc, Idst) ;
        }
        else
        {
          ecode = h_reduce(Isrc, Idst, nint(1.0f/x_scale), nint(1.0f/y_scale));
          if (ecode != HIPS_OK)
            ErrorExit(ecode,
                      "ImageResize: h_%s(%2.3f, %2.3f) returned %d\n",
                      x_scale > 1.0f ? "enlarge" : "reduce", ecode);
        }
      }
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
  int  old, ecode, frame, nframes ;
  byte *src_image, *dst_image ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  if (Idst->numpix != Isrc->numpix || Idst->num_frame < Isrc->num_frame)
    ErrorReturn(NULL, (ERROR_BADPARM, "ImageCopy: dst not big enough")) ;

  src_image = Isrc->image ;
  dst_image = Idst->image ;
  nframes = Isrc->num_frame ;
  Isrc->num_frame = Idst->num_frame = 1 ;
  
  for (frame = 0 ; frame < nframes ; frame++)
  {
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
      case PFDOUBLE:
        old = hips_cplxtor ;
        hips_cplxtor = CPLX_REAL ;
        ecode = h_tod(Isrc, Idst) ;
        if (ecode != HIPS_OK)
          ErrorExit(ecode, "ImageCopy: h_tod failed (%d)\n", ecode) ;
        hips_cplxtor = old ;
        break ;
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
    Isrc->firstpix += Isrc->sizeimage ;
    Isrc->image += Isrc->sizeimage ;
    if (Idst != Isrc)
    {
      Idst->firstpix += Idst->sizeimage ;
      Idst->image += Idst->sizeimage ;
    }
  }

  Isrc->firstpix = Isrc->image = src_image ;
  Idst->firstpix = Idst->image = dst_image ;
  Isrc->num_frame =  Idst->num_frame = nframes ;

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
ImageEdgeDetect(IMAGE *Isrc, IMAGE *Idst, float sigma, int wsize, 
                float lthresh, float uthresh, int dothin)
{
  int    ecode ;
  IMAGE  *Iout ;
  float  fmin, fmax ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;

  ImageValRange(Isrc, &fmin, &fmax) ;
  ImageScale(Isrc, Isrc, 0.0f, 255.0f) ;
  if (Idst->pixel_format != PFBYTE)
    Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, Isrc->num_frame) ;
  else
    Iout = Idst ;

  ecode = h_canny(Isrc, Iout, sigma, wsize, lthresh, uthresh, dothin) ;
  ImageScale(Isrc, Isrc, fmin, fmax) ;
  if (ecode != NO_ERROR)
    ErrorPrintf(ecode, "h_canny returned error code %d", ecode) ;

  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    if (Iout->rows == 0)
      break_here() ;
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
  IMAGE *Iconj, *Ifcorr, *Ifsrc, *Ireal, *Iimag ;
  int   ecode ;

#if 0
  if (zeropad)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, 
                       "ImageCorrelate: zero padding unsupported")) ;
#else
  if (zeropad)
    Isrc = ImageZeroPad(Isrc, NULL) ;
#endif

  /* assumes the template as already been FTed */
  Iconj = ImageConjugate(Itemplate, NULL) ;
  Ifsrc = ImageDFT(Isrc, NULL) ;
  Ifcorr = ImageMul(Iconj, Ifsrc, NULL) ;
  Icorr = ImageInverseDFT(Ifcorr, Icorr) ;

  if (Icorr->pixel_format == PFCOMPLEX)
  {
    /* flipquad can't handle complex images, do it separately */
    Ireal = ImageAlloc(Icorr->rows, Icorr->cols, PFFLOAT, 1) ;
    Iimag = ImageAlloc(Icorr->rows, Icorr->cols, PFFLOAT, 1) ;
    ImageSplit(Icorr, Ireal, Iimag) ;
    ecode = h_flipquad(Ireal, Ireal) ;
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode) ;
    ecode = h_flipquad(Iimag, Iimag) ;
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode) ;
    ImageCombine(Ireal, Iimag, Icorr) ;
    ImageFree(&Ireal) ;
    ImageFree(&Iimag) ;
  }
  else
  {
    ecode = h_flipquad(Icorr, Icorr) ;
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageCorrelate: h_flipquad failed (%d)\n", ecode) ;
  }

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
           perform a spatial correlation of Ikernel with Isrc around the
           region of the point x0,y0. The size of the correlation region is 
           given by wsize 
------------------------------------------------------*/
IMAGE *
ImageCorrelateRegion(IMAGE *Isrc, IMAGE *Ikernel, IMAGE *Idst, int row0, 
                     int col0, int wsize)
{
  int    col_offset, row_offset, row, col, rows, cols, whalf, krow0, kcol0,
         rowstart, rowend, colstart, colend, krow, kcol, drow,dcol ;
  CPIX   *src, *kernel ;
  float  sreal, simag, kreal, kimag, val, *dst, total ;
  
  rows = Isrc->rows ;
  cols = Isrc->cols ;

  if (!Idst)
    Idst = ImageAlloc(rows, cols, PFFLOAT, 1) ;

  kcol0 = col0 - Isrc->cols/2 ;
  krow0 = row0 - Isrc->rows/2 ;

  whalf = (wsize-1)/2 ;
  rowstart = MAX(row0 - whalf, 0) ;
  colstart = MAX(col0 - whalf, 0) ;
  rowend = MIN(rows-1, row0+wsize) ;
  colend = MIN(cols-1, col0+wsize) ;

  for (row_offset = -whalf ; row_offset <= whalf ; row_offset++)
  {
    drow = row0 + row_offset ;
    if (drow < 0 || drow >= rows)
      continue ;

    dst = IMAGEFpix(Idst, col0-whalf, drow) ;
    for (col_offset = -whalf ; col_offset <= whalf ; col_offset++, dst++)
    {
      dcol = col0 + col_offset ;
      if (dcol < 0 || dcol >= cols)
        continue ;
      total = 0.0f ;
      for (row = 0 ; row < rows ; row++)
      {
        krow = row + row_offset + krow0 ;
        if (krow < 0 || krow >= rows)
          continue ;

        src = IMAGECpix(Isrc, 0, row) ;
        for (col = 0 ; col < cols ; col++, src++)
        {
          kcol = col + col_offset + kcol0 ;
          if (kcol < 0 || kcol >= cols)
            continue ;
          kernel = IMAGECpix(Ikernel, krow, kcol) ;
          kreal = kernel->real ;
          kimag = kernel->imag ;
          sreal = src->real ;
          simag = src->imag ;
          val = kreal*sreal + kimag*simag ;  /* real part */
          total += val ;
        }
      }
      *dst = total ;
    }
  }

  return(Idst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
#if 0
IMAGE   *
ImageCopyArea(IMAGE *Isrc, IMAGE *Idst, int srow, int scol,
            int drow, int dcol, int rows, int cols)
{
  return(Idst) ;
}
#endif
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
  byte    *cptr, cval ;

  if (r0 < 0) r0 = 0 ;
  if (c0 < 0) c0 = 0 ;

  if (rows < 0) rows = I->rows ;
  if (cols < 0) cols = I->cols ;

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
      ErrorReturn(ERROR_UNSUPPORTED, 
                  (ERROR_UNSUPPORTED, 
                   "ImageClearArea: unsupported image format %d",
                   I->pixel_format)) ;
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
  int    max_row = -1, max_col = -1, row, col, rows, cols ;

  if (I->pixel_format != PFFLOAT)
    ErrorReturn(0.0f, 
                (ERROR_UNSUPPORTED, "ImageFindPeak: only supports PFFLOAT")) ;

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
  *pval = max_val ;
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
#if 0
static int
imageLargeEnough(IMAGE *Isrc, IMAGE *Idst)
{
  if (Isrc->num_frame > Idst->num_frame)
    return(0) ;
  if (Isrc->numpix > Idst->numpix)
    return(0) ;

  return(1) ;
}
#endif
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *
ImageNormalizePix(IMAGE *Isrc, IMAGE *Idst)
{
  float  scale, fmin = 0.0f, fmax = 0.0f ;
  int    ecode ;
  Pixelval  pmin, pmax ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  ecode = h_minmax(Isrc, &pmin, &pmax, 0) ;
  switch (Isrc->pixel_format)
  {
  case PFBYTE:
    fmin = (float)pmin.v_byte ;
    fmax = (float)pmax.v_byte ;
    break ;
  case PFFLOAT:
    fmin = pmin.v_float ;
    fmax = pmax.v_float ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageNormalize: unsupported pixel format %d\n",
              Isrc->pixel_format) ;
    break ;
  }

  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode,"ImageNormalize: h_minmax failed (%d)\n", ecode));

  if (FEQUAL(fmax, fmin))
#if 1
    {
      ImageCopy(Isrc, Idst) ;
      return(Idst) ;
    }
#else
    ErrorReturn(NULL, (ERROR_BADPARM, "ImageNormalize: constant image")) ;
#endif

  scale = 1.0f / (fmax - fmin) ;
  fmin = -fmin * scale ;

  ecode = h_linscale(Isrc, Idst, scale, fmin) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, 
                (ecode, "ImageNormalize: h_linscale failed (%d)\n", ecode)) ;

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
  long        npix,i;

  npix = (long)Isrc->orows * Isrc->ocols * Isrc->num_frame ;
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
  Pixelval retval ;
  int      row, col, endrow, endcol ;
  float    real, imag ;
  CPIX     *cpix ;

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
  int    format = 0;
  long   bytes ;

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
  bytes = (long) mat->rows * mat->cols * sizeof(float) ;
  if (mat->type == MATRIX_COMPLEX)
    bytes *= 2 ;

  hmemcpy(mat->data, I->image, bytes) ; 
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
  int    format;
  long   bytes ;

  format = (matrix->type == MATRIX_COMPLEX) ? PFCOMPLEX : PFFLOAT ;

  if (!I)
    I = ImageAlloc(matrix->rows, matrix->cols, format, 1) ;
  else
    if (I->rows != matrix->rows || I->cols != matrix->cols)
      ErrorExit(ERROR_BADPARM, "ImageFromMatrix: size mismatch") ;

  bytes = (long)matrix->rows * matrix->cols * sizeof(float) ;
  if (matrix->type == MATRIX_COMPLEX)
    bytes *= 2 ;

  hmemcpy(I->image, matrix->data, bytes) ; 
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
  char *number, buf[200] ;
  int   frame ;

  strcpy(buf, fname) ;
  number = strrchr(buf, '#') ;

  if (number)
  {
    sscanf(number+1, "%d", &frame) ;
    *number = 0 ;
  }
  else
    frame = 0 ;

  return(frame) ;
} 
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           since h_linscale only outputs FLOAT images (and doesn't
           complain otherwise!), we must allocate a temp. image for
           output purposes unless the supplied one is already in
           float format.
------------------------------------------------------*/
IMAGE *
ImageScale(IMAGE *Isrc, IMAGE *Idst, float new_min, float new_max)
{
  float  scale, old_min, old_max ;
  int    ecode, nframes, frame ;
  Pixelval pmin, pmax ;
  IMAGE  *Iout ;
  byte   *src_image, *out_image ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame) ;

  if (FEQUAL(new_max, new_min))
#if 1
    {
      ImageCopy(Isrc, Idst) ;
      return(Idst) ;
    }
#else
    ErrorReturn(NULL, 
            (ERROR_BADPARM, "ImageScale: specified min and max are equal"));
#endif

  old_min = old_max = 0.0f ;  /* remove warning */

  if (Idst->pixel_format != PFFLOAT)  /* h_linscale outputs floats */
    Iout = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, Isrc->num_frame) ;
  else
    Iout = Idst ;


  src_image = Isrc->image ;
  out_image = Iout->image ;
  nframes = Isrc->num_frame ;
  Isrc->num_frame =  Iout->num_frame = 1 ;
  for (frame = 0 ; frame < nframes ; frame++)
  {
    ecode = h_minmax(Isrc, &pmin, &pmax, 0) ;
    switch (Isrc->pixel_format)
    {
    case PFBYTE:
      old_min = (float)pmin.v_byte ;
      old_max = (float)pmax.v_byte ;
      break ;
    case PFINT:
      old_min = (float)pmin.v_int ;
      old_max = (float)pmax.v_int ;
      break ;
    case PFFLOAT:
      old_min = pmin.v_float ;
      old_max = pmax.v_float ;
      break ;
    case PFDOUBLE:
      old_min = pmin.v_double ;
      old_max = pmax.v_double ;
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED, 
                "ImageScale: unsupported pixel format %d\n",
                Isrc->pixel_format) ;
      break ;
    }
    
    if (ecode != HIPS_OK)
      ErrorExit(ecode, "ImageScale: h_minmax failed (%d)\n", ecode) ;

    if (FEQUAL(old_max, old_min))
#if 1
      {
        ImageCopy(Isrc, Idst) ;
        return(Idst) ;
      }
#else
      ErrorReturn(NULL, (ERROR_BADPARM, "ImageScale: constant image")) ;
#endif
    
    scale = (new_max - new_min) / (old_max - old_min) ;


    ecode = h_linscale(Isrc, Iout, scale, new_min - old_min * scale) ;
    if (ecode != HIPS_OK)
      ErrorReturn(NULL, (ecode,"ImageScale: h_linscale failed (%d)\n",ecode));

    Isrc->firstpix += Isrc->sizeimage ;
    Isrc->image += Isrc->sizeimage ;
    if (Isrc != Iout)
    {
      Iout->firstpix += Iout->sizeimage ;
      Iout->image += Iout->sizeimage ;
    }
  }

  Isrc->firstpix = Isrc->image = src_image ;
  Iout->firstpix = Iout->image = out_image ;
  Isrc->num_frame =  Iout->num_frame = nframes ;

  if (Iout != Idst)  /* copy it to desired format */
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
int
ImageCheckSize(IMAGE *inImage,IMAGE *outImage, int rows, int cols, int nframes)
{
  long  inPix, outPix ;

  if (!outImage)
    return(0) ;

  if (!rows)
    rows = inImage->rows ;
  if (!cols)
    cols = inImage->cols ;
  if (!nframes)
    nframes = inImage->num_frame ;

  inPix = (long)rows * (long)cols * (long)nframes * (long)inImage->sizepix ;
  outPix = (long)outImage->numpix * (long)outImage->sizepix * 
           (long)outImage->num_frame ;

  return(outPix >= inPix) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
               change the size of an image
----------------------------------------------------------------------*/
int    
ImageSetSize(IMAGE *I, int rows, int cols)
{
  if (!ImageCheckSize(I, I, rows, cols, 0))
    return(0) ;

  I->frow = I->fcol = 0 ;
  I->rows = I->orows = rows ;
  I->cols = I->ocols = cols ;
  I->numpix = rows*cols ;
  I->sizeimage = rows*cols*I->sizepix ;
  return(1) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageCopyFrames(IMAGE *inImage, IMAGE *outImage,int start, int nframes, 
               int dst_frame)
{
  byte  *cIn, *cOut ;
  unsigned int  *iIn, *iOut ;
  float *fsrc, *fdst ;
  double *dsrc, *ddst ;
  int   size, frameno, pix_per_frame, end ;

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
    case PFDOUBLE:
      if (outImage->pixel_format != PFDOUBLE)
        ErrorExit(ERROR_UNSUPPORTED, 
                  "ImageCopyFrames: unsupported image pixel format %d -> %d\n",
                  inImage->pixel_format, outImage->pixel_format) ;

      dsrc = IMAGEDseq_pix(inImage, 0, 0, frameno) ;
      ddst = IMAGEDseq_pix(outImage, 0, 0, dst_frame+frameno-start) ;
      hmemcpy((char *)ddst, (char *)dsrc, pix_per_frame*sizeof(double)) ;
      break ;
    case PFFLOAT:
      if (outImage->pixel_format != PFFLOAT)
        ErrorExit(ERROR_UNSUPPORTED, 
                  "ImageCopyFrames: unsupported image pixel format %d -> %d\n",
                  inImage->pixel_format, outImage->pixel_format) ;

      fsrc = IMAGEFseq_pix(inImage, 0, 0, frameno) ;
      fdst = IMAGEFseq_pix(outImage, 0, 0, dst_frame+frameno-start) ;
      hmemcpy((char *)fdst, (char *)fsrc, pix_per_frame*sizeof(float)) ;
      break ;
    case PFBYTE:
      if (outImage->pixel_format == PFBYTE)
        hmemcpy(IMAGEseq_pix(inImage,0,0,frameno),
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
          cOut = (byte *)IMAGEseq_pix(outImage, 0, 0, frameno) ;
          while (size--)
            *cOut++ = *cIn++ ;
          break ;
        case PFINT:
          iOut = (unsigned int *)IMAGEIseq_pix(outImage, 0, 0, frameno) ;
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
        hmemcpy((char *)iOut, (char *)iIn, pix_per_frame*sizeof(int)) ;
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

  ImageSetSize(outImage, inImage->rows, inImage->cols) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int  
ImageScaleRange(IMAGE *image, float fmin, float fmax, int low, int high)
{
  int    size;
  byte   *csrc, cmin_val, cmax_val, cval ;
  int    *isrc, imin_val, imax_val, ival ;
  float  *fsrc,  fval, norm ;
  double *dsrc,  dval, dmin, dmax, dnorm, dlow ;

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
      fval = (float)(cval - cmin_val) * norm ;
      cval = (byte)((byte)fval +(byte)low) ;
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
  case PFDOUBLE:
    dmin = (double)fmin ;
    dmax = (double)fmax ;
    size = image->cols * image->rows ;
    dsrc = IMAGEDpix(image, 0, 0) ;
    dnorm = ((double)high - (double)low) / (dmax-dmin);
    dlow = (double)low ;
    while (size--)
    {
      dval = *dsrc ;
      dval = ((dval - dmin)*dnorm) + (double)low ;
      *dsrc++ = dval ;
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
    fprintf(stderr, "ImageScale: unsupported format %d\n",image->pixel_format);
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
    ErrorPrintf(ERROR_NO_FILE, "ImageReadInto(%s) failed\n", fname) ;
  
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
  
  tmp_image = ImageAlloc(image->rows, image->cols,image->pixel_format,nframes);
  ImageCopyFrames(image, tmp_image, start, nframes, 0) ;
  ImageWrite(tmp_image, fname) ;
  ImageFree(&tmp_image) ;
  return(NO_ERROR);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageRescale(IMAGE *inImage, IMAGE *outImage, float scale)
{
  int rows, cols ;

  cols = nint((float)inImage->cols * scale) ;
  rows = nint((float)inImage->rows * scale) ;
  if (!outImage)
   outImage = ImageAlloc(rows, cols, inImage->pixel_format,inImage->num_frame);

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
  int    inRow, inCol, outRow, outCol, inCols, inRows, 
         outRows, outCols, frame ;
  UCHAR  *outPix ;
  byte   *in_image, *out_image ;
  float  *foutPix ;

  if (!ImageCheckSize(inImage, outImage, nint((float)inImage->rows*scale),
                        nint((float)inImage->cols*scale), inImage->num_frame))
    ErrorReturn(-1, (ERROR_NO_MEMORY,
                     "ImageScaleDown: output image not big enough\n")) ;

  ImageSetSize(outImage, nint((float)inImage->rows * scale), 
               nint((float)inImage->cols * scale)) ;

  inRows = inImage->rows ;
  inCols = inImage->cols ;
  outRows = outImage->rows ;
  outCols = outImage->cols ;

  in_image = inImage->image ;
  out_image = outImage->image ;
  for (frame = 0 ; frame < inImage->num_frame ; frame++)
  {
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
            inRow = nint((float)outRow / scale) ;
            inCol = nint((float)outCol / scale) ;
            *foutPix = (float)*IMAGEpix(inImage, inCol, inRow) ;
          }
        break ;
      case PFBYTE:  /* byte --> byte */
        outPix = IMAGEpix(outImage, 0, 0) ;
        for (outRow = 0 ; outRow < outRows ; outRow++)
          for (outCol = 0 ; outCol < outCols ; outCol++, outPix++)
          {
            /* map center point to this output point */
            inRow = nint((float)outRow / scale) ;
            inCol = nint((float)outCol / scale) ;
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
        ErrorReturn(-1, (ERROR_UNSUPPORTED,
                       "ImageScaleDown: unsupported output pixel format %d\n", 
                         outImage->pixel_format)) ;
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
            inRow = nint((float)outRow / scale) ;
            inCol = nint((float)outCol / scale) ;
            *outPix = (UCHAR)*IMAGEFpix(inImage, inCol, inRow) ;
          }
        break ;
      case PFFLOAT:
        /* if scale is a power of 2, use reduce */
        /* if scale is a power of 2, use reduce */

        if (ISPOW2(1.0f/scale))
        {
          int   reductions, i ;
          IMAGE *Itmp ;
      
          reductions = (int)log2(1.0f/scale) ;
          for (i = 0 ; i < reductions ; i++)
          {
            Itmp = inImage ;
            inImage = ImageReduce(Itmp, NULL) ;
            if (i)  /* first one is real source image, don't free it */
              ImageFree(&Itmp) ;
          }
          ImageCopy(inImage, outImage) ;
        }
        else
        {
          foutPix = IMAGEFpix(outImage, 0, 0) ;
          for (outRow = 0 ; outRow < outRows ; outRow++)
            for (outCol = 0 ; outCol < outCols ; outCol++, foutPix++)
            {
              /* map center point to this output point */
              inRow = nint((float)outRow / scale) ;
              inCol = nint((float)outCol / scale) ;
              *foutPix = (float)*IMAGEFpix(inImage, inCol, inRow) ;
            }
        }
        break ;
      default:
        ErrorReturn(-1, (ERROR_UNSUPPORTED,
                       "ImageScaleDown: unsupported output pixel format %d\n", 
                         outImage->pixel_format)) ;
        break ;
      }
      break ;
    case PFINT:
      
    default:
      ErrorReturn(-2, (ERROR_UNSUPPORTED,
                       "ImageScaleDown: unsupported pixel format %d\n", 
                       inImage->pixel_format)) ;
    }
    inImage->image += inImage->sizeimage ;
    inImage->firstpix += inImage->sizeimage ;
    if (inImage != outImage)
    {
      outImage->image += outImage->sizeimage ;
      outImage->firstpix += outImage->sizeimage ;
    }
  }
  inImage->image = inImage->firstpix = in_image ;
  outImage->image = outImage->firstpix = out_image ;
  

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
ImageScaleUp(IMAGE *inImage, IMAGE *outImage, float scale)
{
  int    inRow, inCol, outRow, outCol, inCols, inRows, endCol, endRow,
         outRows, outCols, frame ;
  UCHAR  *inPix, *outPix ;
  UINT   *inIPix, *outIPix ;
  float  *finPix, *foutPix ;
  byte   *in_image, *out_image ;

  if (!ImageCheckSize(inImage, outImage, nint(inImage->rows*scale),
                          nint(inImage->cols*scale), inImage->num_frame))
    ErrorReturn(-1, (ERROR_NO_MEMORY,
          "ImageScaleUp: output image not large enough %d x %d -> %d x %d\n",
            inImage->rows, inImage->cols, outImage->rows, outImage->cols));

  outCols = outImage->cols = nint((float)inImage->cols * scale) ;
  outRows = outImage->rows = nint((float)inImage->rows * scale) ;

  inRows = inImage->rows ;
  inCols = inImage->cols ;

  in_image = inImage->image ;
  out_image = outImage->image ;
  for (frame = 0 ; frame < inImage->num_frame ; frame++)
  {
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
            endRow = nint( (float)inRow * scale + scale) ;
            endCol = nint( (float)inCol * scale + scale) ;
            for (outRow = nint( (float)inRow *scale);outRow < endRow ;outRow++)
            {
              foutPix = IMAGEFpix(outImage, nint((float)inCol * scale),outRow);
              
              for (outCol = nint( (float)inCol * scale);outCol < endCol ; 
                   outCol++,foutPix++)
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
            endRow = nint( (float)inRow * scale + scale) ;
            endCol = nint( (float)inCol * scale + scale) ;
            for (outRow = nint( (float)inRow*scale);outRow < endRow ;outRow++)
            {
              outPix = IMAGEpix(outImage, nint((float)inCol * scale), outRow) ;
              
              for (outCol = nint( (float)inCol * scale) ; outCol < endCol ; 
                   outCol++, outPix++)
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
            endRow = nint( (float)inRow * scale + scale) ;
            endCol = nint( (float)inCol * scale + scale) ;
            for (outRow = nint( (float)inRow*scale);outRow < endRow ;outRow++)
            {
              foutPix = IMAGEFpix(outImage, nint((float)inCol * scale),outRow);
              
              for (outCol = nint( (float)inCol * scale) ; outCol < endCol ; 
                   outCol++,foutPix++)
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
          endRow = nint( (float)inRow * scale + scale) ;
          endCol = nint( (float)inCol * scale + scale) ;
          for (outRow = nint( (float)inRow * scale);outRow < endRow ; outRow++)
          {
            outIPix = IMAGEIpix(outImage, nint((float)inCol * scale), outRow) ;
            
            for (outCol = nint( (float)inCol * scale); outCol < endCol ; 
                 outCol++, outIPix++)
              *outIPix = *inIPix ;
          }
        }
      break ;
    default:
      fprintf(stderr, "ImageScaleUp: unsupported pixel format %d\n", 
              inImage->pixel_format) ;
      return(-2) ;
    }
    inImage->image += inImage->sizeimage ;
    inImage->firstpix += inImage->sizeimage ;
    if (inImage != outImage)
    {
      outImage->image += outImage->sizeimage ;
      outImage->firstpix += outImage->sizeimage ;
    }
  }

  inImage->image = inImage->firstpix = in_image ;
  outImage->image = outImage->firstpix = out_image ;
  
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
  ImageSetSize(outImage, inImage->rows, inImage->cols) ;

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
  long    npix ;
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

  npix = (long)inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = (float)randomNumber(1.0-(double)amp, 1.0+(double)amp) ;
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
  long     npix ;
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

  npix = (long)inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = (float)randomNumber(0.0, 1.0) ;
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
  long    npix ;
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

  npix = (long)inImage->rows * inImage->cols * inImage->num_frame ;
  inPix = IMAGEFpix(inImage, 0, 0) ;
  outPix = IMAGEFpix(outImage, 0, 0) ;
  while (npix--)
  {
    gnoise = (float)randomNumber(-(double)amp, (double)amp) ;
    *outPix++ = *inPix++ + gnoise ;
  }
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int  
ImageValRange(IMAGE *image, float *pfmin, float *pfmax)
{
  float  fmax, fmin, *fpix ;
  double dmax, dmin, *dpix ;
  unsigned int    size, imax, imin, *ipix ; /* "unsiged" added dng */
  byte   bmin, bmax, *bpix ;

  size = image->rows * image->cols * image->num_frame ;
  switch (image->pixel_format)
  {
  case PFCOMPLEX:
  case PFDBLCOM:
    *pfmin = 0.0f ;
    *pfmax = 1.0f ;
    break ;
  case PFDOUBLE:
    dpix = IMAGEDpix(image, 0, 0) ;
    if (!isnan(*dpix))
      dmax = dmin = *dpix ;
    else
      dmax = dmin = 0.0f ;
    while (size--)
    {
      if (isnan(*dpix))
        continue ;
      if (*dpix > dmax)
        dmax = *dpix ;
      if (*dpix < dmin)
        dmin = *dpix ;
      dpix++ ;
    }
    *pfmax = (double)dmax ;
    *pfmin = (double)dmin ;
    break ;
  case PFFLOAT:
    fpix = IMAGEFpix(image, 0, 0) ;
    if (!isnan(*fpix))
      fmax = fmin = *fpix ;
    else
      fmax = fmin = 0.0f ;
    while (size--)
    {
      if (isnan((double)*fpix))
        continue ;
      if (*fpix > fmax)
        fmax = *fpix ;
      if (*fpix < fmin)
        fmin = *fpix ;
      fpix++ ;
    }
    *pfmax = fmax ;
    *pfmin = fmin ;
    break ;
  case PFINT:
    ipix = IMAGEIpix(image, 0, 0) ;
    imax = imin = *ipix ;
    while (size--)
    {
      if (*ipix > imax)
        imax = *ipix ;
      if (*ipix < imin)
        imin = *ipix ;
      ipix++ ;
    }
    *pfmax = (float)imax ;
    *pfmin = (float)imin ;
    break ;
  case PFBYTE:
    bpix = IMAGEpix(image, 0, 0) ;
    bmax = bmin = *bpix ;
    while (size--)
    {
      if (*bpix > bmax)
        bmax = *bpix ;
      if (*bpix < bmin)
        bmin = *bpix ;
      bpix++ ;
    }
    *pfmax = (float)bmax ;
    *pfmin = (float)bmin ;
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
  int    ecode;
  hsize_t size ;
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
  hsize_t    size ;
  float  *fpix ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 
                      Isrc->num_frame);

  switch (Isrc->pixel_format)
  {
  case PFFLOAT:
    size = Isrc->numpix * (hsize_t)Isrc->num_frame ;
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
              replace pixels of value 'inpix' with the value 'outpix'
----------------------------------------------------------------------*/
IMAGE    *
ImageReplace(IMAGE *Isrc, IMAGE *Idst, float inpix, float outpix)
{
  float  *fin, *fout ;
  byte   *cin, *cout, cinpix, coutpix ;
  hsize_t    npix ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows,Isrc->cols,Isrc->pixel_format,Isrc->num_frame);

  if (Idst->pixel_format != Isrc->pixel_format)
    ErrorReturn(NULL, (ERROR_BADPARM, "ImageReplace: src and dst formats must match")) ;

  npix = Isrc->numpix * (hsize_t)Isrc->num_frame ;
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
    cinpix = (byte)inpix ;
    coutpix = (byte)outpix ;
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
                (ERROR_UNSUPPORTED, 
                 "ImageReplace: unsupported pixel format %d", Isrc->pixel_format));
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
  char *number, *dot, buf[100] ;

  strcpy(outFname, inFname) ;
  number = strrchr(outFname, '#') ;
  dot = strrchr(outFname, '.') ;

  if (number)   /* : in filename indicates frame # */
  {
    if (sscanf(number+1, "%d", pframe) < 1)
      *pframe = -1 ;
    *number = 0 ;
  }
  else
    *pframe = -1 ;

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
/*----------------------------------------------------------------------
            Parameters:

           Description:
              compare 2 images. Return values are:
              0  - images are the same
              1  - images are linearly independent
              -1 - images are linearly dependent
----------------------------------------------------------------------*/
int
ImageCmp(IMAGE *Isrc, IMAGE *Idst)
{
  int      ret, ecode ;
  IMAGE    *Idiv ;
  Pixelval pmin, pmax ;
  float    fmin = 0.0f, fmax = 0.0f ;

  Idiv = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 1) ;

  ecode = h_div(Isrc, Idst, Idiv) ;
  if (ecode != HIPS_OK)
    ErrorReturn(-1, (ecode, "ImageCmp: h_div returned %d", ecode)) ;

  ecode = h_minmax(Idiv, &pmin, &pmax, 0) ;
  if (ecode != HIPS_OK)
    ErrorExit(ecode, "ImageCmp: h_minmax failed (%d)\n", ecode) ;

  switch (Isrc->pixel_format)
  {
  case PFBYTE:
    fmin = (float)pmin.v_byte ;
    fmax = (float)pmax.v_byte ;
    break ;
  case PFFLOAT:
    fmin = pmin.v_float ;
    fmax = pmax.v_float ;
    break ;
  default:
    ErrorExit(ERROR_UNSUPPORTED, 
              "ImageCmp: unsupported pixel format %d\n", Isrc->pixel_format);
    break ;
  }

  if (fmin != fmax)  /* if Idiv is constant - they are linearly dependent */
    ret = 1 ;
  else
  {
    if (FEQUAL(fmin, 1.0f))
      ret = 0 ;
    else
      ret = -1 ;
  }

  return(ret);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              return the number of frames stored in the file 'fname'
----------------------------------------------------------------------*/
int
ImageNumFrames(char *fname)
{
  IMAGE  I ;
  FILE   *fp ;
  int    frame, type, ecode, nframes ;
  char   buf[100] ;

  ImageUnpackFileName(fname, &frame, &type, buf) ;
  fname = buf ;
  if ((frame >= 0) || (type != HIPS_IMAGE))
    return(1) ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, 
                (ERROR_NO_FILE, "ImageNumFrame(%s) could not open file\n", 
                 fname)) ;

  ecode = fread_header(fp, &I, fname) ;
  if (ecode != HIPS_OK)
    ErrorReturn(ERROR_NO_FILE,
                (ERROR_NO_FILE, 
                 "ImageNumFrame: fread_header failed (%d)\n",ecode));

  nframes = I.num_frame ;
  fclose(fp) ;
  free_hdrcon(&I) ;
  return(nframes) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              append an image to the end of a hips sequence file, incrementing
              the number of frames recorded in the header.
----------------------------------------------------------------------*/
int
ImageAppend(IMAGE *I, char *fname)
{
  FILE   *fp ;
  int    ecode, frame = 0, nframes ;
  IMAGE  Iheader, *Iframe ;
  char   tmpname[200] ;

  fp = fopen(fname, "r+b") ;
#if 0
  if (!fp)
    ErrorReturn(-1, (ERROR_NO_FILE, "ImageAppend(%s) failed\n", fname)) ;
#endif

  if (!fp)
    return(ImageWrite(I, fname)) ;

  ecode = fread_header(fp, &Iheader, fname) ;
  if (ecode != HIPS_OK)
    ErrorExit(ERROR_NO_FILE, "ImageAppend: fread_header failed (%d)\n",ecode);

  /* increment # of frames, and update file header */
  Iheader.num_frame++ ;
  nframes = Iheader.num_frame ;
  if (nframes == 10 || nframes == 100 || nframes == 1000)
  {
    /* header size will grow by 1 byte, must copy whole file (ughhh) */
    fclose(fp) ;
    strcpy(tmpname,FileTmpName(NULL)) ;
    FileRename(fname, tmpname) ;

    /* write out new header */
    fp = fopen(fname, "wb") ;
    if (!fp)
      ErrorReturn(-1, (ERROR_NO_FILE, "ImageAppend(%s) failed\n", fname)) ;
    
    ecode = fwrite_header(fp, &Iheader, fname) ;
    if (ecode != HIPS_OK)
     ErrorExit(ERROR_NO_FILE,"ImageAppend: fwrite_header failed (%d)\n",ecode);

    nframes = Iheader.num_frame - 1 ;
    for (frame = 0 ; frame < nframes ; frame++)
    {
      Iframe = ImageReadFrames(tmpname, frame, 1) ;
      if (!Iframe)
        ErrorReturn(-3, (ERROR_BADFILE, 
                         "ImageAppend: could not read %dth frame", frame)) ;
      ecode = fwrite_image(fp, Iframe, frame, fname) ;
      if (ecode != HIPS_OK)
        ErrorReturn(-4, (ERROR_BADFILE, 
                          "ImageAppend: fwrite_image frame %d failed (%d)\n",
                          ecode,frame));
    }
    unlink(tmpname) ;
  }
  else    /* seek back to start and increment # of frames */
  {
    if (fseek(fp, 0L, SEEK_SET) < 0)
      ErrorReturn(-2,(ERROR_BADFILE,"ImageAppend(%s): could not seek to end"));
    ecode = fwrite_header(fp, &Iheader, fname) ;
    if (ecode != HIPS_OK)
     ErrorExit(ERROR_NO_FILE,"ImageAppend: fwrite_header failed (%d)\n",ecode);
  }

  if (fseek(fp, 0L, SEEK_END) < 0)
    ErrorReturn(-2, (ERROR_BADFILE, "ImageAppend(%s): could not seek to end"));

  ecode = fwrite_image(fp, I, frame, "fwrite") ;
  if (ecode != HIPS_OK)
    ErrorReturn(-1, (ERROR_BADFILE, 
              "ImageAppend: fwrite_image frame %d failed (%d)\n",ecode,frame));

  free_hdrcon(&Iheader) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageAbs(IMAGE *inImage, IMAGE *outImage)
{
  UCHAR *cIn, *cOut ;
  UINT  *iIn, *iOut ;
  float *fIn, *fOut ;
  int   size, nframes, frameno, pix_per_frame, rows, cols ;

  rows = inImage->rows ;
  cols = inImage->cols ;
  if (!outImage)
  {
    outImage = ImageAlloc(rows, cols, 1, inImage->pixel_format) ;
  }

  nframes = inImage->num_frame ;
  pix_per_frame = inImage->rows * inImage->cols ;
  for (frameno = 0 ; frameno < nframes ; frameno++)
  {
    size = inImage->rows * inImage->cols ;
    switch (inImage->pixel_format)
    {
    case PFFLOAT:
      fIn = IMAGEFpix(inImage, 0, 0) + pix_per_frame * frameno ;
      fOut = IMAGEFpix(outImage, 0, 0) + pix_per_frame * frameno  ;
      while (size--)
        *fOut++ = (float)fabs(*fIn++) ;
      break ;
    case PFBYTE:
      cIn = IMAGEpix(inImage, 0, 0) + pix_per_frame * frameno ;
      switch (outImage->pixel_format)
      {
      case PFBYTE:
        cOut = IMAGEpix(outImage, 0, 0) + pix_per_frame * frameno ;
        while (size--)
          *cOut++ = (byte)abs((int)(*cIn++)) ;
        break ;
      case PFINT:
        iOut = IMAGEIpix(outImage, 0, 0) + pix_per_frame * frameno ;
        while (size--)
          *iOut++ = (UINT)abs(*cIn++) ;
        break ;
      default:
        ErrorExit(ERROR_BADPARM,
                  "ImageAbs: unsupported output image pixel format (%d)\n",
                  outImage->pixel_format) ;
        return(NULL) ;
        break ;
      }
      break ;
    case PFINT:
      iIn = IMAGEIpix(inImage, 0, 0) + pix_per_frame * frameno ;
      switch (outImage->pixel_format)
      {
      case PFINT:
        iOut = IMAGEIpix(outImage, 0, 0) + pix_per_frame * frameno ;
        while (size--)
          *iOut++ = abs(*iIn++) ;
        break ;
        break ;
      default:
        ErrorExit(ERROR_BADPARM,
                  "ImageAbs: unsupported output image pixel format (%d)\n",
                  outImage->pixel_format) ;
        return(NULL) ;
        break ;
      }
      break ;
    default:
      ErrorExit(ERROR_BADPARM,
                "ImageAbs: unsupported input image pixel format (%d)\n",
                inImage->pixel_format) ;
      return(NULL) ;
      break ;
    }
  }

  ImageSetSize(outImage, inImage->rows, inImage->cols) ;
  
  return(outImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageAdd(IMAGE *Is1, IMAGE *Is2, IMAGE *Idst)
{
  int ecode ;

  if (!Idst)
    Idst = ImageAlloc(Is1->rows, Is1->cols,Is1->pixel_format, Is1->num_frame);

  ecode = h_add(Is1, Is2, Idst) ;
  if (ecode != HIPS_OK)
    ErrorPrintf(ecode, "ImageAdd: h_add failed (%d)", ecode) ;

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Idst = Is1 - Is2
----------------------------------------------------------------------*/
IMAGE *
ImageSubtract(IMAGE *Is1, IMAGE *Is2, IMAGE *Idst)
{
  int ecode ;

  if (!Idst)
    Idst = ImageAlloc(Is1->rows, Is1->cols,Is1->pixel_format, Is1->num_frame);

  ecode = h_diff(Is1, Is2, Idst) ;
  if (ecode != HIPS_OK)
    ErrorPrintf(ecode, "ImageSubtract: h_diff failed (%d)", ecode) ;

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageExtractInto(IMAGE *Isrc, IMAGE *Idst, int x0, int y0, 
                     int dx, int dy, int xdst, int ydst)
{
  CPIX     *cpsrc, *cpdst ;
  UCHAR    *csrc, *cdst ;
  float    *fsrc, *fdst ;
  int      xin, yin, yout, x1, y1, yend, xend ;

  if (!Idst)
    Idst = ImageAlloc(dy, dx, Isrc->pixel_format, Isrc->num_frame) ;
  else
    if (Isrc->pixel_format != Idst->pixel_format)
      ErrorReturn(NULL, (ERROR_BADPARM, 
                         "ImageExtractInto: out format must match"
                         "input format\n")) ;

  x1 = x0 + dx ;
  y1 = y0 + dy ;
  xend = Isrc->cols-1 ;
  yend = Isrc->rows-1 ;

  switch (Isrc->pixel_format)
  {
  case PFCOMPLEX:
    yout = ydst ;
    for (yin = y0 ; yin < y1 ; yin++, yout++)
    {
      cpsrc = IMAGECpix(Isrc, x0, yin) ;
      cpdst = IMAGECpix(Idst, xdst, yout) ;
      for (xin = x0 ; xin < x1 ; xin++, cpdst++, cpsrc++)
      {
        if (xin < 0 || xin > xend || yin < 0 || yin > yend)
          cpdst->real = cpdst->imag = 0.0f ;
        else
          *cpdst = *cpsrc ;
      }
    }
    break ;
  case PFBYTE:
    yout = ydst ;
    for (yin = y0 ; yin < y1 ; yin++, yout++)
    {
      csrc = IMAGEpix(Isrc, x0, yin) ;
      cdst = IMAGEpix(Idst, xdst, yout) ;
      for (xin = x0 ; xin < x1 ; xin++, cdst++, csrc++)
      {
        if (xin < 0 || xin > xend || yin < 0 || yin > yend)
          *cdst = 0 ;
        else
          *cdst = *csrc ;
      }
    }
    break ;
  case PFFLOAT:
    yout = ydst ;
    for (yin = y0 ; yin < y1 ; yin++, yout++)
    {
      fsrc = IMAGEFpix(Isrc, x0, yin) ;
      fdst = IMAGEFpix(Idst, xdst, yout) ;
      for (xin = x0 ; xin < x1 ; xin++, fdst++, fsrc++)
      {
        if (xin < 0 || xin > xend || yin < 0 || yin > yend)
          *fdst = 0.0f ;
        else
          *fdst = *fsrc ;
      }
    }
    break ;
  default:
    ErrorReturn(Idst, (ERROR_UNSUPPORTED, 
                       "ImageExtractInto: unsupported image format %d\n", 
                       Isrc->pixel_format)) ;
    break ;
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageExtract(IMAGE *Isrc, IMAGE *Idst, int x0, int y0,int dx, int dy)
{
  UCHAR    *csrc, *cdst ;
  float    *fsrc, *fdst ;
  int      xin, yin, xout, yout, x1, y1, yend, xend ;

  if (!Idst)
    Idst = ImageAlloc(dy, dy, Isrc->pixel_format, Isrc->num_frame) ;

  x1 = x0 + dx ;
  y1 = y0 + dy ;
  xend = Isrc->cols-1 ;
  yend = Isrc->rows-1 ;

  switch (Isrc->pixel_format)
  {
  case PFBYTE:
    yout = xout = 0 ;
    cdst = IMAGEpix(Idst, 0, 0) ;
    for (yin = y0 ; yin < y1 ; yin++, yout++)
    {
      csrc = IMAGEpix(Isrc, x0, yin) ;
      cdst = IMAGEpix(Idst, 0, yout) ;
      for (xout = 0, xin = x0 ; xin < x1 ; xin++, xout++, cdst++, csrc++)
      {
        if (xin < 0 || xin > xend || yin < 0 || yin > yend)
          *cdst = 0 ;
        else
          *cdst = *csrc ;
      }
    }
    break ;
  case PFFLOAT:
    yout = xout = 0 ;
    fdst = IMAGEFpix(Idst, 0, 0) ;
    for (yin = y0 ; yin < y1 ; yin++, yout++)
    {
      fsrc = IMAGEFpix(Isrc, x0, yin) ;
      fdst = IMAGEFpix(Idst, 0, yout) ;
      for (xout = 0, xin = x0 ; xin < x1 ; xin++, xout++, fdst++, fsrc++)
      {
        if (xin < 0 || xin > xend || yin < 0 || yin > yend)
          *fdst = 0.0f ;
        else
          *fdst = *fsrc ;
      }
    }
    break ;
  default:
    ErrorReturn(Idst, (ERROR_UNSUPPORTED, 
                       "ImageExtract: unsupported image format %d\n", 
                       Isrc->pixel_format)) ;
    break ;
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
ImageZeroMean(IMAGE *Isrc, IMAGE *Idst)
{
  int    frameno, rows, cols, row, col, pix_per_frame, nframes ;
  float  ftotal, fmean, *fSrcPtr, *fDstPtr, *fSrcBase, *fDstBase ;

  if (Isrc->pixel_format != PFFLOAT || 
      (Idst && Idst->pixel_format != PFFLOAT))
  {
    fprintf(stderr, "ImageZeroMean: unsupported pixel format %d\n",
            Isrc->pixel_format) ;
    return(NULL) ;
  }

  if (!Idst)
    Idst = ImageClone(Isrc) ;

  nframes = Isrc->num_frame ;
  rows = Isrc->rows ;
  cols = Isrc->cols ;
  pix_per_frame = rows * cols ;
  fSrcBase = IMAGEFpix(Isrc, 0, 0) ;
  fDstBase = IMAGEFpix(Idst, 0, 0) ;

  /* mean of each pixel across all frames */
  for (row = 0 ; row < rows ; row++)
  {
    for (col = 0 ; col < cols ; col++, fSrcBase++, fDstBase++)
    {
      ftotal = 0.0f ;
      fSrcPtr = fSrcBase ;
      for (frameno = 0 ; frameno < nframes ; frameno++)
      {
        ftotal += *fSrcPtr ;
        fSrcPtr += pix_per_frame ;
      }
      fmean = ftotal / (float)nframes ;

      /* subtract mean from this pixel in every frame */
      fDstPtr = fDstBase ;
      fSrcPtr = fSrcBase ;
      for (frameno = 0 ; frameno < nframes ; frameno++)
      {
        *fDstPtr = *fSrcPtr - fmean ;
        fDstPtr += pix_per_frame ;
        fSrcPtr += pix_per_frame ;
      }
    }
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              form the covariance matrix treating each frame of
              Isrc as an observation.
----------------------------------------------------------------------*/
IMAGE *
ImageCovarMatrix(IMAGE *image, float **pmeans)
{
  static IMAGE *zimage = NULL ;    /* zero-mean version of image */
  IMAGE  *cimage ;
  int        rows, cols, row, col, crow, ccol, crows, ccols, pix_per_frame,
             nframes, frameno, i ;
  float      *flPtr, *frPtr, *fDstPtr, *flBase, *frBase, ftotal,
             *means, *meanPtr, *fSrcPtr, *fSrcBase ;

  if (image->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImageCovarMatrix: input image must be FLOAT\n") ;
    return(NULL) ;
  }

  if (!ImageCheckSize(image, zimage, 0, 0, 1))
  {
    if (zimage)
      ImageFree(&zimage) ;
    zimage = ImageClone(image) ;
  }

  rows = image->rows ;
  cols = image->cols ;
  nframes = image->num_frame ;
  pix_per_frame = rows * cols ;

  /* first calculate mean across observations (frames) */
  means = (float *)calloc((unsigned int)pix_per_frame, sizeof(float)) ;
  if (!means)
  {
    fprintf(stderr, "ImageCovarMatrix: could not allocate mean vector\n") ;
    return(NULL) ;
  }

  /* mean of each pixel across all frames */
  fSrcBase = IMAGEFpix(image, 0, 0) ;
  for (i = row = 0 ; row < rows ; row++)
  {
    for (col = 0 ; col < cols ; col++, fSrcBase++, i++)
    {
      ftotal = 0.0f ;
      fSrcPtr = fSrcBase ;
      for (frameno = 0 ; frameno < nframes ; frameno++)
      {
        ftotal += *fSrcPtr ;
        fSrcPtr += pix_per_frame ;
      }
      means[i] = ftotal / (float)nframes ;
    }
  }

  /* zimage will hold the centered (0 mean) version of image */
  fDstPtr = IMAGEFpix(zimage, 0, 0) ;
  fSrcPtr = IMAGEFpix(image, 0, 0) ;
  for (frameno = 0 ; frameno < nframes ; frameno++)
  {
    meanPtr = means ;
    for (i = row = 0 ; row < rows ; row++)
    {
      for (col = 0 ; col < cols ; col++, i++)
      {
        *fDstPtr++ = *fSrcPtr++ - *meanPtr++ ;
      }
    }
  }

/* 
  now form covariance matrix, treating each frame of the input sequence
  as a row in a num_frame x (rows*cols) sized matrix and forming the 
  outer product of that matrix with itself. The covariance matrix will then
  have the dimension (rows*cols) x (rows*cols).
*/
  ccols = crows = rows*cols ;
  cimage = ImageAlloc(crows, ccols, PFFLOAT, 1) ;
  if (!cimage)
  {
    fprintf(stderr, 
            "ImageCovarMatrix: could not allocate %d x %d covariance matrix\n",
            crows, ccols) ;
    free(means) ;
    return(NULL) ;
  }
  fDstPtr = IMAGEFpix(cimage, 0, 0) ;
  for (crow = 0 ; crow < crows ; crow++)
  {
    for (ccol = 0 ; ccol < ccols ; ccol++)
    {
/* 
      Calculate value of this entry in covariance matrix by multiplying
      the crow'th position in each image by the ccol'th position in each image.
*/
      ftotal = 0.0f ;
      flPtr = flBase = (float *)zimage->image + crow ;
      frPtr = frBase = (float *)zimage->image + ccol ;
      for (i = 0 ; i < nframes ; i++)
      {
        ftotal += *flPtr * *frPtr ;
        flPtr += pix_per_frame ;
        frPtr += pix_per_frame ;
      }

      *fDstPtr++ = ftotal / (nframes-1) ;  /* unbiased covariance */
    }
  }

  if (pmeans)
    *pmeans = means ;
  else
    free(means) ;
  return(cimage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
#include "eigen.h"
#include "matfile.h"

static int compare_evalues(const void *l1, const void *l2)  ;


typedef struct
{
  int   eno ;
  float evalue ;
} EIGEN_VALUE, EVALUE ;

static int
compare_evalues(const void *l1, const void *l2)
{
  EVALUE *e1, *e2 ;

  e1 = (EVALUE *)l1 ;
  e2 = (EVALUE *)l2 ;
  return(fabs(e1->evalue) < fabs(e2->evalue) ? 1 : -1) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
              perform the Hotelling or Karhunen-Loeve transform on a
              sequence of images. The output image will contain the
              sorted eigenvectors of the covariance matrix (the principal
              components) in decreasing order of eigenvalue size. The 
              num_frame-2st frame contains the mean vectors, while the 
              last frame contains the eigenvalues.

              Given an input image x where each frame is an observation of
              the cols*rows variables, we form the matrix A whose rows are
              orthonormal eigenvectors of the covariance matrix of x, then
              construct the principal components y via:

              y = A (x - mx)

              where mx is the mean vector of the x variables.
----------------------------------------------------------------------*/
IMAGE *
ImagePrincipalComponents(IMAGE *image, int nterms, IMAGE **pcoefImage)
{
  IMAGE *cimage, *pcImage, *zImage, *coefImage ;
  float     *evalues, *evectors ;
  int       frameno, nframes, row, col, rows, cols, pix_per_frame, i, nevalues,
             frame ;
  float     *fSrcPtr, *fDstPtr, *mean_vector, *fPcPtr ;
  EVALUE    *eigen_values ;

  if (image->pixel_format != PFFLOAT)
  {
    fprintf(stderr, "ImagePrincipalComponents: input image must be FLOAT\n") ;
    return(NULL) ;
  }

  cimage = ImageCovarMatrix(image, &mean_vector) ;

  nevalues = nframes = cimage->rows ;
  if (!nterms)
    nterms = nevalues ;

  evalues = (float *)calloc((unsigned int)nevalues, sizeof(float)) ;
  pix_per_frame = image->rows * image->cols ;
  eigen_values = (EVALUE *)calloc((UINT)nevalues, sizeof(EIGEN_VALUE));
  evectors = (float *)calloc((UINT)(nevalues*pix_per_frame), sizeof(float)) ;

  /* 2 extra frames - 1 for mean vector, and one for eigenvalues */
  pcImage = ImageAlloc(image->rows, image->cols, PFFLOAT, nterms+2) ;
  rows = pcImage->rows ;
  cols = pcImage->cols ;
  
  EigenSystem((float *)cimage->image, cimage->rows, evalues, evectors) ;
#if 0
  MatFileWrite("evectors.mat", evectors, cimage->rows, cimage->rows,
               "evectors") ;
#endif

/* 
  sort eigenvalues in order of decreasing absolute value. The
  EIGEN_VALUE structure is needed to record final order of eigenvalue so 
  we can also sort the eigenvectors.
*/
  for (i = 0 ; i < nevalues ; i++)
  {
    eigen_values[i].eno = i ;
    eigen_values[i].evalue = evalues[i] ;
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues) ;
  for (i = 0 ; i < nevalues ; i++)
    evalues[i] = eigen_values[i].evalue ;

  /* columns of evectors are eigenvectors. */

  /* copy eigenvectors into each frame in eigenvalue order */
  for (frameno = 0 ; frameno < nterms ; frameno++)
  {
    fDstPtr = IMAGEFseq_pix(pcImage,0,0, frameno) ;
    fSrcPtr = evectors + eigen_values[frameno].eno ;
    for (row = 0 ; row < rows ; row++)
    {
      for (col = 0 ; col < cols ; col++)
      {
        *fDstPtr++ = *fSrcPtr ;
        fSrcPtr += pix_per_frame ;
      }
    }
  }

  /* copy mean vector into next to last frame */
  fDstPtr = IMAGEFseq_pix(pcImage, 0, 0, pcImage->num_frame-2) ;
  for (i = 0 ; i < pix_per_frame ; i++)
    *fDstPtr++ = mean_vector[i] ;

  /* copy eigenvalues into last frame */
  fDstPtr = IMAGEFseq_pix(pcImage, 0, 0, pcImage->num_frame-1) ;
  for (i = 0 ; i < nevalues ; i++)
    *fDstPtr++ = evalues[i] ;

#if 0
ImageWrite(pcImage, "pc.hipl") ;
MatFileWrite("pc.mat", (float *)pcImage->image, pcImage->num_frame, 
                        pcImage->cols, "pcs") ;
#endif

  if (pcoefImage)  /* also return coefficients for reconstruction */
  {
/*
    the coefficient image will contain the same number of frames as the
    # of principal components to be used in reconstruction. Each frame 
    contains one coefficient per frame in the input image.
*/
    *pcoefImage = coefImage = ImageAlloc(1, image->num_frame, PFFLOAT,nterms);
    if (!coefImage)
    {
      fprintf(stderr, 
              "ImagePrincipalComponents: could not allocated coef image\n");
      exit(3) ;
    }
    zImage = ImageZeroMean(image, NULL) ;
/*
    the coefficients for reconstruction of the observation vectors from
    the means and eigenvectors are given by:
    
    y = A . (x - mx)
*/
    fPcPtr = IMAGEFpix(pcImage, 0, 0) ;
    fSrcPtr = IMAGEFpix(zImage, 0, 0) ;
    fDstPtr = IMAGEFpix(coefImage, 0, 0) ;
    for (frame = 0 ; frame < nterms ; frame++)
    {
/* 
      for first set of coefficients. frame is also the number of the
      eigenvector used in the construction.
*/
      fDstPtr = IMAGEFseq_pix(coefImage, 0, 0, frame) ;

      for (col = 0 ; col < image->num_frame ; col++)
      {
/* 
        for each column in the transposed centered image matrix, multiply
        the col'th frame by the frame'th row.
*/
        fPcPtr = IMAGEFseq_pix(pcImage, 0, 0, frame) ;
        fSrcPtr = IMAGEFseq_pix(zImage, 0, 0, col) ;
        for (i = 0 ; i < nevalues ; i++)
          *fDstPtr += *fPcPtr++ * *fSrcPtr++ ;

        fDstPtr++ ;
      }
    }

    ImageFree(&zImage) ;
  }

  free(mean_vector) ;
  free(evalues) ;
  free(evectors) ;
  free(eigen_values) ;
  ImageFree(&cimage) ;
  return(pcImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              reconstruct an image sequence from a set of principal
              components.

              note that the coefficient image may represent the 
              reconstruction coefficients of more than one input
              image. In that case, this function may be called
              repeatedly, each call using an 'nframes' subset of
              coefImage starting at 'start'.
----------------------------------------------------------------------*/
IMAGE *
ImageReconstruct(IMAGE *pcImage, IMAGE *coefImage, IMAGE *xrImage,
                int start, int nframes)
{
  int     rows, cols, frame, nterms, term, i, pix_per_frame,pix_per_coef_frame;
  float   *fPcPtr, *fXPtr, *fCoefPtr, *means, *Mx, ftotal ;
  float   *fBasePcPtr, *fBaseCoefPtr ;

  rows = pcImage->rows ;
  cols = pcImage->cols ;
  pix_per_frame = rows * cols ;
  pix_per_coef_frame = coefImage->cols * coefImage->rows ;

  
  /* one coefficient for each frame to be reconstructed */
  if (!nframes)
    nframes = coefImage->cols ; 

  nterms = coefImage->num_frame ; /* # of terms in expansion */

  if (!xrImage)
    xrImage = ImageAlloc(rows, cols, PFFLOAT, nframes) ;
  if (!xrImage)
  {
    fprintf(stderr, "ImageReconstruct: could not allocate image!\n") ;
    exit(-2) ;
  }

/*
  reconstruct x image based on a (possibly) reduced set of eigenvectors 
  (in pcImage) and coefficients (in coefImage).

  xr = (Ak' * y)' + mx

  where ' denotes transpose

  We have one coefficient for each frame to be reconstructed, and one frame of
  coefficients for each principal component to use in the reconstruction.
*/
  means = IMAGEFseq_pix(pcImage, 0, 0, pcImage->num_frame-2) ;

  fBaseCoefPtr = IMAGEFpix(coefImage, start, 0) ;
  for (frame = 0 ; frame < nframes ; frame++, fBaseCoefPtr++)
  {
    /* reconstruct the frame #'th observation */
    fXPtr = IMAGEFseq_pix(xrImage, 0, 0, frame) ;
    Mx = means ;
    fBasePcPtr = IMAGEFpix(pcImage, 0, 0) ;

    for (i = 0 ; i < pix_per_frame ; i++)
    {
      /* use i'th component of every term */
      fPcPtr = fBasePcPtr++ ;
      fCoefPtr = fBaseCoefPtr ;
#if 0
      fPcPtr = IMAGEFpix(pcImage, i, 0) ;  
      fCoefPtr = IMAGEFpix(coefImage, frame+start, 0) ;
#endif

      for (ftotal = 0.0f, term = 0 ; term < nterms ; term++)
      {
        /* use the term #'th column of the Ak matrix */
        ftotal += *fPcPtr * *fCoefPtr ;

        fCoefPtr += pix_per_coef_frame ;
        fPcPtr += pix_per_frame ;
      }

      *fXPtr++ = ftotal + *Mx++ ;
    }
  }

  return(xrImage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             normalize the amplitude of a complex image.
----------------------------------------------------------------------*/
IMAGE *
ImageNormalizeComplex(IMAGE *Isrc, IMAGE *Idst, float thresh)
{
  int    rows, cols ;
  long   npix ;
  CPIX   *src, *dst ;
  float  real, imag, mag ;

  if (Isrc->pixel_format != PFCOMPLEX)
    ErrorReturn(NULL,(ERROR_BADPARM,"ImageNormalizeComplex: Isrc not complex"));
  rows = Isrc->rows ;
  cols = Isrc->cols ;
  if (!Idst)
    Idst = ImageAlloc(rows, cols, PFCOMPLEX, 1) ;

  src = IMAGECpix(Isrc, 0, 0) ;
  dst = IMAGECpix(Idst, 0, 0) ;

  npix = Isrc->numpix ;
  if (FZERO(thresh))
    thresh = 0.00001f ;

  while (npix--)
  {
    real = src->real ;
    imag = src->imag ;
    src++ ;
    mag = (float)sqrt((double)real*real + imag*imag) ;
    if (mag < thresh)
      dst->real = dst->imag = 0.0f ;
    else
    {
      dst->real = real / mag ;
      dst->imag = imag / mag ;
    }
    dst++ ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              treat each frame in the sequence as a rows x cols dimensional
              vector and normalize it so its length is 1.
----------------------------------------------------------------------*/
int
ImageNormalizeFrames(IMAGE *Isrc, IMAGE *Idst)
{
  float  flen, *fsrcPtr, *fdstPtr, fval ;
  int    frameno, rows, cols;
  long   npix, pix_per_frame ;

  rows = Isrc->rows ;
  cols = Isrc->cols ;
  pix_per_frame = (long)rows * cols ;
  
  switch (Isrc->pixel_format)
  {
  case PFFLOAT:
    for (frameno = 0 ; frameno < Isrc->num_frame ; frameno++)
    {
      fsrcPtr = IMAGEFpix(Isrc, 0, 0) + frameno * pix_per_frame ;
      npix = pix_per_frame ;
      flen = 0.0f ;
      while (npix--)
      {
        fval = *fsrcPtr++ ;
        flen += fval * fval ;
      }
      flen = (float)sqrt(flen) ;

      if (FZERO(flen))
        return(-2) ;

      fsrcPtr = IMAGEFpix(Isrc, 0, 0) + frameno * pix_per_frame ;
      fdstPtr = IMAGEFpix(Idst, 0, 0) + frameno * pix_per_frame ;
      npix = pix_per_frame ;

      while (npix--)
        *fdstPtr++ = *fsrcPtr++ / flen ;
    }
    break ;
  default:
    fprintf(stderr, "ImageNormalizeFrames: unsupported pixel format %d\n",
            Isrc->pixel_format) ;
    return(-1);
  }
  
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
               generate a complex image from a real and an imaginary one
----------------------------------------------------------------------*/
IMAGE *
ImageCombine(IMAGE *Ireal, IMAGE *Iimag, IMAGE *Idst)
{
  int    x, y, rows, cols ;
  float  *real, *imag ;
  CPIX   *dst ;

  rows = Ireal->rows ;
  cols = Ireal->cols ;

  if (!Idst)
    Idst = ImageAlloc(rows, cols, PFCOMPLEX, 1) ;

  if (Idst->pixel_format != PFCOMPLEX)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, 
                       "ImageCombine: destination must be complex")) ;

  dst = IMAGECpix(Idst, 0, 0) ;
  real = IMAGEFpix(Ireal, 0, 0) ;
  imag = IMAGEFpix(Iimag, 0, 0) ;

  for (y = 0 ; y < rows ; y++)
  {
    for (x = 0 ; x < cols ; x++)
    {
      dst->imag = *imag++ ;
      dst->real = *real++ ;
      dst++ ;
    }
  }

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              flip an image about its horizontal axis
----------------------------------------------------------------------*/
IMAGE *
ImageInvert(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE *Ireal, *Iimag ;
  int   ecode ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, Isrc->pixel_format, 1) ;

  switch (Isrc->pixel_format)
  {
  case PFCOMPLEX:
    Iimag = ImageAlloc(Isrc->rows, Isrc->cols, PFFLOAT, 1) ;
    Ireal = ImageSplit(Isrc, NULL, Iimag) ;
    ecode = h_invert(Ireal, Ireal) ;
    if (ecode != HIPS_OK)
      ErrorReturn(NULL, (ecode, "ImageInvert: h_invert failed %d", ecode)) ;
    ecode = h_invert(Iimag, Iimag) ;
    if (ecode != HIPS_OK)
      ErrorReturn(NULL, (ecode, "ImageInvert: h_invert failed %d", ecode)) ;
    ImageCombine(Ireal, Iimag, Idst) ;
    break ;
  default:
    ecode = h_invert(Isrc, Idst) ;
    if (ecode != HIPS_OK)
      ErrorReturn(NULL, (ecode, "ImageInvert: h_invert failed %d", ecode)) ;
    break ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
               generate a complex image from a real and an imaginary one
----------------------------------------------------------------------*/
IMAGE *
ImageSplit(IMAGE *Icomp, IMAGE *Ireal, IMAGE *Iimag)
{
  int    x, y, rows, cols ;
  float  *real, *imag = NULL ;
  double *dreal, *dimag = NULL ;
  CPIX   *cpix ;
  DCPIX  *dcpix ;

  rows = Icomp->rows ;
  cols = Icomp->cols ;

  if (!Ireal)
    Ireal = ImageAlloc(rows, cols, PFFLOAT, 1) ;

  if (!COMPLEX_IMAGE(Icomp))
    ErrorReturn(ImageCopy(Icomp, Ireal), (ERROR_UNSUPPORTED, 
                       "ImageSplit: source must be complex")) ;

  real = IMAGEFpix(Ireal, 0, 0) ;
  if (Iimag)
    imag = IMAGEFpix(Iimag, 0, 0) ;

  switch (Icomp->pixel_format)
    {
    case PFCOMPLEX:
      cpix = IMAGECpix(Icomp, 0, 0) ;
      
      for (y = 0 ; y < rows ; y++)
        {
          for (x = 0 ; x < cols ; x++)
            {
              if (Iimag)
                *imag++ = cpix->imag ;
              *real++ = cpix->real ;
              cpix++ ;
            }
        }
      break ;
    case PFDBLCOM:
      dcpix = IMAGEDCpix(Icomp, 0, 0) ;
      
      switch (Ireal->pixel_format)
      {
      case PFFLOAT:
        for (y = 0 ; y < rows ; y++)
        {
          for (x = 0 ; x < cols ; x++)
          {
            if (Iimag)
              *imag++ = (float)dcpix->imag ;
            *real++ = (float)dcpix->real ;
            dcpix++ ;
          }
        }
        break ;
      case PFDOUBLE:
        dreal = IMAGEDpix(Ireal, 0, 0) ;
        if (Iimag)
          dimag = IMAGEDpix(Iimag, 0, 0) ;

        for (y = 0 ; y < rows ; y++)
        {
          for (x = 0 ; x < cols ; x++)
          {
            if (Iimag)
              *dimag++ = dcpix->imag ;
            *dreal++ = dcpix->real ;
            dcpix++ ;
          }
        }
        break ;
      }
      break ;
    default:
      break ;
    }

  return(Ireal) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
               shrink an image using Gaussian blurred sampling.
----------------------------------------------------------------------*/
IMAGE *
ImageShrink(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE   *Iin, *Iout, *Igaussian ;
  int     srows, scols, drows, dcols, x, y, xc, yc, xhalf, yhalf,xk,yk,ys ;
  float   smax, smin, xscale, yscale, *dpix ;
  register int xs ;
  register float *spix, *gpix, total ;
  

  srows = Isrc->rows ;
  scols = Isrc->cols ;
  drows = Idst->rows ;
  dcols = Idst->cols ;

  if (Isrc->pixel_format != PFFLOAT)
  {
    Iin = ImageAlloc(srows, scols, PFFLOAT, 1) ;
    ImageCopy(Isrc, Iin) ;
    ImageValRange(Isrc, &smin, &smax) ;
  }
  else
    Iin = Isrc ;

  if (Idst->pixel_format != PFFLOAT)
    Iout = ImageAlloc(drows, dcols, PFFLOAT, 1) ;
  else
    Iout = Idst ;

/* 
  create a Gaussian sampling function whose sigma depends on the 
  amount of scaling.
*/
  xscale = (float)scols / (float)dcols ;
  yscale = (float)srows / (float)drows ;
#if 0
  Igaussian = ImageGaussian(xscale/4.0f, yscale/4.0f) ;
#else
  Igaussian = ImageGaussian(xscale/10.0f, yscale/10.0f) ;
/*  fprintf(stderr, "grows,gcols = %d,%d\n", Igaussian->rows, Igaussian->cols);*/
#endif

  xhalf = (Igaussian->cols-1) / 2 ;
  yhalf = (Igaussian->rows-1) / 2 ;

#if 0
fprintf(stderr, 
        "ImageShrink: xscale %2.3f, yscale %2.3f, grows,gcols = %d,%d\n",
        xscale, yscale, Igaussian->rows, Igaussian->cols) ;
#endif

  /* for each dst pixel, Gaussian sample region around it */
  dpix = IMAGEFpix(Iout, 0, 0) ;
  for (y = 0 ; y < drows ; y++)
  {
    for (x = 0 ; x < dcols ; x++)
    {
      total = 0.0f ;

      xc = nint((float)x * xscale) ;  /* center of Gaussian in source coords */
      yc = nint((float)y * yscale) ;

      /* apply kernel */
      gpix = IMAGEFpix(Igaussian, 0, 0) ;
      for (yk = -yhalf ; yk <= yhalf ; yk++)
      {
        ys = yc + yk ;
        if (ys < 0)
          ys = 0 ;
        else if (ys >= srows)
          ys = srows - 1 ;

        for (xk = -xhalf ; xk <= xhalf ; xk++)
        {
          xs = xc + xk ;
          if (xs < 0)
            xs = 0 ;
          else if (xs >= scols)
            xs = scols - 1 ;

          spix = IMAGEFpix(Iin, xs, ys) ;
          total += *spix * *gpix++ ;
        }
      }

      *dpix++ = total ;
    }
  }
  
  if (Iin != Isrc)
    ImageFree(&Iin) ;
  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
  }
  ImageFree(&Igaussian) ;
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             perform histogram equalization on an image
----------------------------------------------------------------------*/
IMAGE *
ImageHistoEqualize(IMAGE *Isrc, IMAGE *Idst)
{
  IMAGE     *Iin, *Iout ;
  struct hips_histo histogram ;
  int       ecode, count ;
  float     fmin, fmax ;
  Pixelval  crap ;
  byte      map[256] ;

  if (Isrc->pixel_format != PFBYTE)
    Iin = ImageConvertToByte(Isrc, NULL) ;
  else
    Iin = Isrc ;

  ImageValRange(Isrc, &fmin, &fmax) ;   /* scale output image properly*/

  if (Idst->pixel_format != PFBYTE)
    Iout = ImageAlloc(Idst->rows, Idst->cols, PFBYTE, 1) ;
  else
    Iout = Idst ;

  alloc_histo(&histogram, &crap, &crap, 256, PFBYTE) ;
  ecode = h_clearhisto(&histogram) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode, "ImageHistoEqualize: h_clearhisto failed")) ;
  ecode = h_histo(Iin, &histogram, 0, &count) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode, "ImageHistoEqualize: h_histo failed")) ;
  ecode = h_histoeq(&histogram, count, map) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode, "ImageHistoEqualize: h_histoeq failed")) ;
  ecode = h_pixmap(Iin, Iout, map) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode, "ImageHistoEqualize: h_pixmap failed")) ;

  free(histogram.histo) ;

  if (Iin != Isrc)
    ImageFree(&Iin) ;
  if (Iout != Idst)
  {
    ImageCopy(Iout, Idst) ;
    ImageFree(&Iout) ;
    ImageScale(Idst, Idst, fmin, fmax) ;
  }
  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              convert an image to byte format, scaling its intensity values to
              0-255
----------------------------------------------------------------------*/
#include "rescale.h"
IMAGE *
ImageConvertToByte(IMAGE *Isrc, IMAGE *Idst)
{
  int MinPoint[2], MaxPoint[2], ecode ;

  if (!Idst)
    Idst = ImageAlloc(Isrc->rows, Isrc->cols, PFBYTE, 1) ;

  ecode = h_rescale(Isrc, 0.0f, 255.0f, MinPoint, MaxPoint, Idst) ;
  if (ecode != HIPS_OK)
    ErrorReturn(NULL, (ecode, "ImageConvert: h_rescale failed\n")) ;

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              calculate the mean and variance of the intensity in an image.
----------------------------------------------------------------------*/
static void break_now(void) ;
static void break_now(void) {}

int
ImageStatistics(IMAGE *Isrc, float *pmean, float *pvar)
{
  long   npix ;
  float  *pix, total, dif, mean ;
  IMAGE  *I ;

  if (Isrc->pixel_format != PFFLOAT)
  {
    I = ImageAlloc(I->rows, I->cols, PFFLOAT, 1) ;
    ImageCopy(Isrc, I) ;
  }
  else
    I = Isrc ;

  npix = I->numpix ;
  pix = IMAGEFpix(I, 0, 0) ;

  /* compute mean */
  total = 0.0f ;
  while (npix--)
    total += *pix++ ;

  npix = I->numpix ;
  mean = *pmean = total / (float)npix ;
  pix = IMAGEFpix(I, 0, 0) ;

  /* compute variance */
  total = 0.0f ;
  while (npix--)
  {
    dif = *pix++ - mean ;
    total += dif * dif ;
  }

  *pvar = total / (float)I->numpix ;

  if (isnan(*pvar) || *pvar < 0.0f)
  {
    break_now() ;
    ImageStatistics(I, pmean, pvar) ;
  }

  if (I != Isrc)
    ImageFree(&I) ;

  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              pad an image with zeros out to the next power of 2.
----------------------------------------------------------------------*/
IMAGE *
ImageZeroPad(IMAGE *Isrc, IMAGE *Idst)
{
  int  drows, dcols, logrows, logcols, scols, srows, dcol, drow ;

  scols = Isrc->cols ;
  srows = Isrc->rows ;
  logcols = (int)log2((double)scols) ;
  logrows = (int)log2((double)srows) ;

  dcols = (int)exp2((double)(logcols+1)) ;
  drows = (int)exp2((double)(logrows+1)) ;

  if (!Idst)
    Idst = ImageAlloc(drows, dcols, Isrc->pixel_format, Isrc->num_frame) ;

  drow = (drows - srows) / 2 ;
  dcol = (dcols - scols) / 2 ;
  ImageExtractInto(Isrc, Idst, 0, 0, scols, srows, dcol, drow) ;

  return(Idst) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             extract the interior region of a zero-padded image
----------------------------------------------------------------------*/
IMAGE *
ImageUnpad(IMAGE *Isrc, IMAGE *Idst, int rows, int cols)
{
  int   row0, col0 ;

  if (!Idst)
    Idst = ImageAlloc(rows, cols, Isrc->pixel_format, Isrc->num_frame) ;

  if (!rows)
    rows = Idst->rows ;
  if (!cols)
    cols = Idst->cols ;

  row0 = (Isrc->rows - rows) / 2 ;
  col0 = (Isrc->cols - cols) / 2 ;
  ImageExtractInto(Isrc, Idst, col0, row0, cols, rows, 0, 0) ;
  return(Idst) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
             determine whether the values in an image are valid or not
----------------------------------------------------------------------*/
int
ImageValid(IMAGE *I)
{
  int    size ;
  float  *fpix ;
  double *dpix, exponent, val ;

  size = I->rows * I->cols * I->num_frame ;

  switch (I->pixel_format)
  {
  case PFFLOAT:
    fpix = IMAGEFpix(I, 0, 0); 
    while (size--)
    {
      val = *fpix++ ;
      if (val == 0.0)
        continue ;
      exponent = log(fabs(val)) ;
      if (fabs(exponent > 20))
        return(0) ;
    }
    break ;
  case PFDOUBLE:
    dpix = IMAGEDpix(I, 0, 0); 
    while (size--)
    {
      val = *dpix++ ;
      if (val == 0.0)
        continue ;
      exponent = log(fabs(val)) ;
      if (fabs(exponent > 20))
        return(0) ;
    }
    break ;
  case PFDBLCOM:
    break ;
  }

  return(1) ;
}

