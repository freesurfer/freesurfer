/**
 * @brief image processing prototypes
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

#ifndef OPEN_GL_CODE

#ifndef IMAGE_H
#define IMAGE_H

#include "matrix.h"

#include "hips.h"

typedef struct header IMAGE ;

/* include this after definition of IMAGE */
#include "kernel.h"

/* allocation */
IMAGE   *ImageAlloc(int rows, int cols, int format, int nframes) ;
IMAGE   *ImageAllocHeader(int rows, int cols, int format, int nframes) ;
int     ImageAllocBuffer(IMAGE *image) ;
int     ImageFree(IMAGE **pI) ;

/* resizing */
int      ImageScaleDown(IMAGE *inImage, IMAGE *outImage, float scale) ;
IMAGE    *ImageDifferentialScale(IMAGE *inImage, IMAGE *outImage,
                                 int outRows, int outCols) ;
int      ImageDifferentialScaleUp(IMAGE *inImage, IMAGE *outImage,
                                  int outRows, int outCols) ;
int      ImageDifferentialScaleDown(IMAGE *inImage, IMAGE *outImage,
                                    int outRows, int outCols) ;
int      ImageScaleUp(IMAGE *inImage, IMAGE *outImage, float scale) ;
IMAGE   *ImageResize(IMAGE *Isrc, IMAGE *Idst, int drows, int dcols) ;


/* file I/O */

/* reading */
IMAGE   *ImageRead(const char *fname) ;
IMAGE   *ImageReadType(const char *fname, int pixel_format) ;
IMAGE    *ImageReadFrames(const char *fname, int start, int nframes) ;
int      ImageReadInto(const char *fname, IMAGE *image, int image_no) ;
IMAGE   *ImageFRead(FILE *fp,const  char *fname, int start, int nframes) ;
IMAGE   *ImageReadHeader(const char *fname) ;
IMAGE   *ImageFReadHeader(FILE *fp,const  char *fname) ;
IMAGE   *ImageInvert(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageCombine(IMAGE *Ireal, IMAGE *Iimag, IMAGE *Idst) ;
IMAGE   *ImageSplit(IMAGE *Icomp, IMAGE *Ireal, IMAGE *Iimag) ;
IMAGE   *ImageShrink(IMAGE *Isrc, IMAGE *Idst) ;

/* writing */
int     ImageWrite(IMAGE *image,const char *fname) ;
int     ImageFWrite(IMAGE *image, FILE *fp,const  char *fname) ;
int      ImageWriteFrames(IMAGE *image,const  char *fname, int start, int nframes) ;
int      ImageAppend(IMAGE *image,const  char *fname) ;
int      ImageUpdateHeader(IMAGE *image,const  char *fname) ;

/* morphology */
IMAGE   *ImageDilate(IMAGE *Isrc, IMAGE *Idst, int which) ;
IMAGE   *ImageErode(IMAGE *Isrc, IMAGE *Idst, int which) ;
IMAGE   *ImageMorph(IMAGE *Isrc, IMAGE *Idst, int which) ;
IMAGE    *ImageClose(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE    *ImageOpen(IMAGE *Isrc, IMAGE *Idst) ;

IMAGE   *ImageThreshold(IMAGE *Isrc, IMAGE *Idst, float threshold) ;
IMAGE   *ImageDFT(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageInverseDFT(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageMul(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst) ;
IMAGE   *ImageMulScale(IMAGE *Isrc, IMAGE *Idst, Pixelval *p) ;
IMAGE   *ImageCopy(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageEdgeDetect(IMAGE *Isrc, IMAGE *Idst, float sigma, int wsize,
                         float lthresh, float uthresh, int dothin);
IMAGE   *ImageCorrelate(IMAGE *Itemplate, IMAGE *Isrc, int zpad,IMAGE *Icorr) ;
IMAGE   *ImageCopyArea(IMAGE *Isrc, IMAGE *Idst, int srow, int scol,
                       int drow, int dcol, int rows, int cols) ;
int     ImageClearArea(IMAGE *image, int row, int col, int rows, int cols,
                       float val, int frame) ;
float   ImageFindPeak(IMAGE *image, int *prow, int *pcol, float *pval) ;
IMAGE   *ImagePowerSpectrum(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageNormalizePix(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE   *ImageConjugate(IMAGE *Isrc, IMAGE *Idst) ;
Pixelval ImageAccum(IMAGE *Isrc) ;
MATRIX  *ImageToMatrix(IMAGE *image) ;
IMAGE   *ImageFromMatrix(MATRIX *matrix, IMAGE *Idst) ;
int     ImageType(const char *fname) ;
#define ImageSize(image)   ((image)->orows * (image)->ocols)
#define ImageBytes(image)  ((image)->sizeimage)
IMAGE   *ImageScale(IMAGE *Isrc, IMAGE *Idst, float fmin, float fmax) ;
int     ImageCheckSize(IMAGE *inImage,IMAGE *outImage, int rows,
                       int cols, int nframes) ;
int     ImageFrame(const char *fname) ;
int     ImageUnpackFileName(const char *inFname, int *pframe, int *ptype,
                            char *outFname) ;
int      ImageCopyFrames(IMAGE *inImage, IMAGE *outImage,int start,
                         int nframes, int dst_frame);
int      ImageScaleRange(IMAGE *image, float fmin, float fmax,
                         int low, int high) ;
IMAGE    *ImageRescale(IMAGE *inImage, IMAGE *outImage, float scale) ;
int      ImageReflect(IMAGE *inImage, IMAGE *outImage, int how) ;
int      ImageAddNoise(IMAGE *inImage, IMAGE *outImage, float amp) ;
int      ImageAddSaltNoise(IMAGE *inImage,IMAGE *outImage, float density) ;
int      ImageAddSpeckleNoise(IMAGE *inImage,IMAGE *outImage, float amp) ;
int      ImageValRange(IMAGE *image, float *pfmin, float *pfmax) ;
IMAGE    *ImageCatSeq(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst) ;
IMAGE    *ImageInverse(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE    *ImageMatrixMul(IMAGE *Isrc1, IMAGE *Isrc2, IMAGE *Idst) ;
IMAGE    *ImageAddScalar(IMAGE *Isrc, IMAGE *Idst, float scalar) ;
IMAGE    *ImageAdd(IMAGE *Is1, IMAGE *Is2, IMAGE *Idst) ;
IMAGE    *ImageSubtract(IMAGE *Is1, IMAGE *Is2, IMAGE *Idst) ;
IMAGE    *ImageReplace(IMAGE *Isrc, IMAGE *Idst, float inpix, float outpix) ;
int      ImageCmp(IMAGE *Isrc, IMAGE *Idst) ;
int      ImageNumFrames(const char *fname) ;
float    ImageRMSDifference(IMAGE *I1, IMAGE *I2) ;

/* offset stuff */
IMAGE    *ImageOffsetMagnitude(IMAGE *Isrc, IMAGE *Idst, int maxsteps) ;
int      ImageAdjustOffsetLen(IMAGE *inImage, IMAGE *outImage,
                              int wsize) ;
IMAGE    *ImageCalculateNitShiOffset(IMAGE *Ix, IMAGE *Iy, int wsize,
                                     float mu, float c, IMAGE *offsetImage) ;
IMAGE    *ImageCalculateOffsetDirection(IMAGE *Ix, IMAGE *Iy, int wsize,
                                        IMAGE *offsetImage) ;
IMAGE    *ImageOffsetDirectionMap(IMAGE *Ix, IMAGE *Iy, int wsize,
                                  IMAGE *Iorient, IMAGE *Idir, IMAGE *Ioffset);
IMAGE    *ImageOffsetScale(IMAGE *Isrc, IMAGE *Idst) ;
int      ImageCalculateMomentOffset(IMAGE *gradImage, int wsize, float c,
                                    IMAGE *offsetImage);
IMAGE    *ImageCalculateOffset(IMAGE *Ix, IMAGE *Iy, int wsize,
                               IMAGE *offsetImage);
IMAGE    *ImageOffsetOrientation(IMAGE *Ix,IMAGE *Iy,int wsize,IMAGE *Iorient);
IMAGE    *ImageOffsetDirection(IMAGE *Ix,IMAGE *Iy,int wsize,IMAGE *Iorient,
                               IMAGE *Ioffset);
IMAGE    *ImageNitshiOffsetDirection(IMAGE *Ix,IMAGE *Iy,int wsize,
                                     IMAGE *Iorient, IMAGE *Ioffset);
IMAGE    *ImageOffsetDirectionMagnitude(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy,
                                        int wsize, IMAGE *Idst,int maxsteps);

IMAGE    *ImageFilterMinMax(IMAGE *Imin, IMAGE *Imax, IMAGE *Ioffset,
                            IMAGE *Idir, IMAGE *Idst) ;


/* various filters */
int      ImageMedianFilter(IMAGE *inImage, int wsize,
                           IMAGE *offsetImage, IMAGE *outImage) ;
IMAGE    *ImageSigmaFilter(IMAGE *Isrc, int wsize, float nsigma,
                           IMAGE *Ioffset, IMAGE *Idst) ;
int      ImageBuildExponentialFilter(IMAGE *inImage, int wsize, float k,
                                     IMAGE *offsetImage,
                                     IMAGE *filterSequnce) ;
int      ImageExponentialFilter(IMAGE *inImage, IMAGE *gradImage,
                                int wsize, float k, IMAGE *offsetImage,
                                IMAGE *outImage) ;
int      ImageSpaceVariantFilter(IMAGE *inImage, IMAGE *filterSequence,
                                 IMAGE *outImage) ;
int      ImageDiffuse(IMAGE *inImage, IMAGE *outImage, double k,
                      int niter, int which, double slope, KIMAGE *kimage) ;
int      ImageDiffuseCurvature(IMAGE *inImage,IMAGE *outImage,
                               double A,int niter, double slope,KIMAGE *kimage) ;
int      ImageDiffusePerona(IMAGE *inImage, IMAGE *outImage,
                            double k, int niter,double slope, KIMAGE *kimage);
int      ImageDiffuseHV(IMAGE *inImage, IMAGE *outImage, double k,
                        int niter, double slope, KIMAGE *kimage) ;

int      ImageCurvature(IMAGE *inImage, float A, IMAGE *gradImage) ;
IMAGE    *ImageAbs(IMAGE *inImage, IMAGE *outImage) ;
int      ImageSobel(IMAGE *inImage, IMAGE *gradImage,
                    IMAGE *xImage, IMAGE *yImage) ;
int      ImageSobelX(IMAGE *inImage, IMAGE *xImage) ;
int      ImageSobelY(IMAGE *inImage, IMAGE *yImage) ;
int      ImageConvolve3x3(IMAGE *inImage, float kernel[], IMAGE *outImage) ;
IMAGE    *ImageXDerivative(IMAGE *inImage, IMAGE *xImage) ;
IMAGE    *ImageYDerivative(IMAGE *inImage, IMAGE *yImage) ;
void    ImageCircularConvolve1d(IMAGE *imageI, IMAGE *J, float k[], int len,
                                int axis) ;
void    ImageConvolve1d(IMAGE *imageI, IMAGE *J,
                        float k[], int len, int axis) ;
void    ImageConvolve1dByte(IMAGE *imageI, IMAGE *J,
                            float k[], int len, int axis) ;
IMAGE   *ImageGaussian(float xsigma, float ysigma) ;
IMAGE   *ImageGaussian1d(float sigma, int max_len) ;
IMAGE   *ImageCircularConvolveGaussian(IMAGE *Isrc,IMAGE *gImage, IMAGE *Iout,
                                       int dst_frameno) ;
void    ImageCircularConvolve1d(IMAGE *imageI, IMAGE *J,
                                float k[], int len, int axis) ;
IMAGE   *ImageConvolveGaussianByte(IMAGE *inImage, IMAGE *gImage,
                                   IMAGE *outImage, int dst_frameno) ;
IMAGE   *ImageConvolveGaussian(IMAGE *inImage, IMAGE *gImage,
                               IMAGE *outImage, int dst_frameno) ;
IMAGE   *ImageCircularConvolveGaussian(IMAGE *inImage, IMAGE *gImage,
                                       IMAGE *outImage, int dst_frameno) ;
IMAGE   *ImageConvolveGaussianFrames(IMAGE *inImage, IMAGE *gImage,
                                     IMAGE *outImage) ;
IMAGE   *ImageLaplacian(IMAGE *inImage, IMAGE *outImage) ;
IMAGE   *ImageMeanFilter(IMAGE *Isrc, int wsize, IMAGE *Idst) ;

IMAGE   *ImageExtractInto(IMAGE *inImage, IMAGE *outImage, int x0,
                          int y0, int dx, int dy, int xdst, int ydst) ;
IMAGE   *ImageExtract(IMAGE *inImage, IMAGE *outImage, int x0,
                      int y0,int dx,int dy) ;
IMAGE   *ImageReduce(IMAGE *inImage, IMAGE *outImage) ;
void   ImageReduce1d(IMAGE *imageI, IMAGE *J, float k[], int len, int axis) ;
IMAGE   *ImageCovarMatrix(IMAGE *image, float **pmeans) ;
IMAGE   *ImagePrincipalComponents(IMAGE *image, int nterms,
                                  IMAGE **pcoefImage) ;
IMAGE   *ImageReconstruct(IMAGE *pcImage, IMAGE *coefImage,
                          IMAGE *xrImage, int start, int nframes) ;
IMAGE  *ImageZeroMean(IMAGE *inImage, IMAGE *outImage) ;
IMAGE  *ImageNormalizeComplex(IMAGE *Isrc, IMAGE *Idst, float thresh) ;
int    ImageNormalizeFrames(IMAGE *inImage, IMAGE *outImage) ;
int    ImageSetSize(IMAGE *image, int rows, int cols) ;
IMAGE  *ImageNormalizeOffsetDistances(IMAGE *Isrc, IMAGE *Idst, int maxstep) ;
IMAGE  *ImageSmoothOffsets(IMAGE *Isrc, IMAGE *Idst, int wsize) ;
IMAGE  *ImageApplyOffset(IMAGE *Isrc, IMAGE *offsetImage, IMAGE *Idst) ;
IMAGE  *ImageOffsetMedialAxis(IMAGE *offsetImage, IMAGE *Iedge) ;
IMAGE  *ImageHistoEqualize(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *ImageConvertToByte(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *ImageCorrelateRegion(IMAGE *Isrc,IMAGE *Ikernel,IMAGE *Idst,
                             int row0, int col0, int wsize);
int    ImageStatistics(IMAGE *image, float *pmean, float *pvar) ;
IMAGE  *ImageZeroPad(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *ImageUnpad(IMAGE *Isrc, IMAGE *Idst, int rows, int cols) ;

int    ImageValid(IMAGE *image) ;

double ImageEntropy(IMAGE *image, int pairflag) ;

IMAGE  *ImageDownsample2Horizontal(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *ImageDownsample2(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE  *ImageUpsample2(IMAGE *Isrc, IMAGE *Idst) ;

// hips functions
int init_header(IMAGE *I,const char *onm,const char *snm,int nfr,const char *odt,int rw,int cl,int pfmt,int nc,const char *desc);
int h_copy(IMAGE *Isrc, IMAGE *Idst);
int free_header(IMAGE *I);

/*
  constants which specify the morphological operation carried out
  by ImageMorph routine.
*/
#define MORPH_BRIDGE       0
#define MORPH_CLEAN        1
#define MORPH_DIAG         2
#define MORPH_DILATE       3
#define MORPH_ERODE        4
#define MORPH_FATTEN       5
#define MORPH_FILL         6
#define MORPH_H_BREAK      7
#define MORPH_MAJORITY     8
#define MORPH_REMOVE       9
#define MORPH_SHRINK1      10
#define MORPH_SHRINK2      11
#define MORPH_SKELETONIZE1 12
#define MORPH_SKELETONIZE2 13
#define MORPH_THIN1        14
#define MORPH_THIN2        15
#define MORPH_THICKEN1     16
#define MORPH_THICKEN2     17
#define MORPH_SPUR         18
#define MORPH_PERIM4       19
#define MORPH_PERIM8       20

#define MORPH_THICKEN      21


/* image formats */
#define IMAGE_BYTE     PFBYTE
#define IMAGE_SHORT    PFSHORT
#define IMAGE_FLOAT    PFFLOAT
#define IMAGE_DOUBLE   PFDOUBLE
#define IMAGE_C_FLOAT  PFCOMPLEX
#define IMAGE_C_DOUBLE PFDBLCOM


#define REAL_PIX      0
#define IMAG_PIX      1

#define ImageClone(I)  ImageAlloc(I->rows,I->cols,I->pixel_format,I->num_frame)
#define ImageClear(I)  ImageClearArea(I, 0, 0, I->rows, I->cols, 0.0f,-1)

/* IMAGEpix is for byte size images */
#define IMAGEpix(im, x, y)           ((im->image) + (((long)y) * im->ocols) \
                                        + (x))

#define IMAGEIpix(im, x, y)          (((unsigned int *)im->image) \
                                     + (((long)y) * im->ocols) + (x))
#define IMAGEFpix(im, x, y)          (((float *)im->image) \
                                     + (((long)y) * im->ocols) + (x))
#define IMAGESpix(im, x, y)          (((short *)im->image) \
                                     + (((long)y) * im->ocols) + (x))

#define IMAGEDpix(im, x, y)          (((double *)im->image) \
                                     + (((long)y) * im->ocols) + (x))
/* for complex images */
#define IMAGECpix(im, x, y)          (((CPIX *)im->image) \
                                     + ((((long)y) * im->ocols)) + (x))

/* for double complex images */
#define IMAGEDCpix(im, x, y)          (((DCPIX *)im->image) \
                                     + ((((long)y) * im->ocols)) + (x))

/* IMAGEpix is for byte size images */
#define IMAGEseq(im, n)           ((im->image) + (((long)n) *im->ocols \
                                                  *im->orows))

/* for float images */
#define IMAGEFseq(im, n)          ((float *)((im)->image) + \
                                  (((long)n) * (im)->ocols * (im)->orows))

/* for double images */
#define IMAGEDseq(im, n)          ((double *)((im)->image) + \
                                  (((long)n) * (im)->ocols * (im)->orows))

/* for integer images */
#define IMAGEIseq(im, n)          ((int *)((im)->image) + \
                                  ((n) * (im)->ocols * (im)->orows))
/* for complex images */
#define IMAGECseq(im, n)          ((float *)(im->image) + \
                                  (((long)n) * 2 * im->ocols * im->orows))

// RGB images are 3 bytes per pixel
#define IMAGERGBpix(im, x, y)           ((im->image) + (((int) y) * im->ocols * 3) + (x*3))


#define IMAGEDseq_pix(im,x,y,n)   ((IMAGEDseq(im,n)) \
                                     + (((long)y) * im->ocols) + (x))
#define IMAGEFseq_pix(im,x,y,n)   ((IMAGEFseq(im,n)) \
                                     + (((long)y) * im->ocols) + (x))
#define IMAGEseq_pix(im,x,y,n)   ((IMAGEseq(im,n)) \
                                     + (((long)y) * im->ocols) + (x))
#define IMAGEIseq_pix(im,x,y,n)   ((IMAGEIseq(im,n)) \
                                     + (((long)y) * im->ocols) + (x))


#define COMPLEX_IMAGE(i)  (((i)->pixel_format == PFCOMPLEX) || \
                                           ((i)->pixel_format == PFDBLCOM))
typedef struct
{
  float  real ;
  float  imag ;
}
CPIX ;

typedef struct
{
  double  real ;
  double  imag ;
}
DCPIX ;


#define HIPS_IMAGE     0
#define MATLAB_IMAGE   1
#define TIFF_IMAGE     2
#define TIF_IMAGE      TIFF_IMAGE
#define JPEG_IMAGE     3
#define PGM_IMAGE      4
#define PPM_IMAGE      5
#define PBM_IMAGE      6 /* Write not yet implemented */
#define RGBI_IMAGE     7

#define END_UNDEF 0
#define END_BIG   1
#define END_SMALL 2

#define IMAGE_REFLECT_AROUND_X_AXIS   0
#define IMAGE_REFLECT_AROUND_Y_AXIS   1

/* for use in ImageConvolve1d */
#define IMAGE_VERTICAL         0
#define IMAGE_HORIZONTAL       1

extern float sx[9], sy[9] ;  /* sobel coefficients */

#define  ImageInitHeader(I, cols, rows, format, nframes) \
  init_header(I, "orig", "seq", nframes, "today", rows,cols,format,1, "temp")

#include "filter.h"

#define DIFFUSE_PERONA         FILTER_DIFFUSE_GRAD
#define DIFFUSE_CURVATURE      FILTER_DIFFUSE_CURV
#define DIFFUSE_HV             FILTER_DIFFUSE_HV

#endif

#endif
