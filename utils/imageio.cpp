/**
 * @brief image processing utilities
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
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

/*
 * replaced alloc_image with alloc_image_buffer (dng, 3/12/96).
 */

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* for SEEK_ constants */

#include "hmem.h"
#include "hips.h"

#include "diag.h"
#include "error.h"
#include "image.h"
#include "machine.h"
#include "macros.h"
#include "matfile.h"
#include "matrix.h"
#include "proto.h"
#include "rgb.h"
#include "tiffio.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static IMAGE *TiffReadImage(const char *fname, int frame);
static IMAGE *TiffReadHeader(const char *fname, IMAGE *I);
static int TiffWriteImage(IMAGE *I, const char *fname, int frame);
static IMAGE *JPEGReadImage(const char *fname);
static IMAGE *JPEGReadHeader(FILE *fp, IMAGE *);
static int JPEGWriteImage(IMAGE *I, const char *fname, int frame);
static int RGBwrite(IMAGE *I, char *fname, int frame);
static IMAGE *RGBReadImage(char *fname);
static IMAGE *RGBReadHeader(char *fname, IMAGE *I);
static IMAGE *PGMReadImage(const char *fname);
static IMAGE *PGMReadHeader(FILE *fp, IMAGE *);
static int PGMWriteImage(IMAGE *I, const char *fname, int frame);
static int PPMWriteImage(IMAGE *I, const char *fname, int frame);
static IMAGE *PPMReadImage(const char *fname);
static IMAGE *PPMReadHeader(FILE *fp, IMAGE *);
static IMAGE *PBMReadImage(const char *fname);
static IMAGE *PBMReadHeader(FILE *fp, IMAGE *);
static byte FindMachineEndian(void);
static void ImageSwapEndian(IMAGE *I);

static byte endian = END_UNDEF;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ImageWrite(IMAGE *I, const char *fname)
{
  FILE *fp;
  // int ecode;

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(-1, (ERROR_NO_FILE, "ImageWrite(%s) failed\n", fname));

  // ecode =
  ImageFWrite(I, fp, fname);
  fclose(fp);
  return (0);
}

static byte FindMachineEndian(void)
{
  short int word = 0x0001;
  char *bite = (char *)&word;
  return (bite[0] ? END_SMALL : END_BIG);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ImageFWrite(IMAGE *I, FILE *fp, const char *fname)
{
  byte *image;
  int ecode, type, frame;
  char *buf;

  if (!fname) fname = "ImageFWrite";
  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname); /* don't destroy callers string */

  ImageUnpackFileName(buf, &frame, &type, buf);

  switch (type) {
    case RGBI_IMAGE:
      RGBwrite(I, buf, frame);
      break;
    case TIFF_IMAGE:
      TiffWriteImage(I, buf, frame);
      break;
    default:
    case JPEG_IMAGE:
      JPEGWriteImage(I, buf, frame);
      break;
    case PGM_IMAGE:
      PGMWriteImage(I, buf, frame);
      break;
    case PPM_IMAGE:
      PPMWriteImage(I, buf, frame);
      break;
    case HIPS_IMAGE:
      if (endian == END_UNDEF) endian = FindMachineEndian();

      if (findparam(I, "endian")) clearparam(I, "endian"); /* Clear any existing endian parameter */
      setparam(I, "endian", PFBYTE, 1, endian);

      ecode = fwrite_header(fp, I, "fwrite");
      if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageFWrite: fwrite_header failed (%d)\n", ecode);

      image = I->image;
      for (frame = 0; frame < I->num_frame; frame++) {
        ecode = fwrite_image(fp, I, frame, "fwrite");
        if (ecode != HIPS_OK)
          ErrorExit(ERROR_NO_FILE, "ImageFWrite: fwrite_image frame %d failed (%d)\n", ecode, frame);
        I->image += I->sizeimage; /* next frame */
      }
      I->image = image;
      break;
  }

  free(buf);
  return (0);
}

static void ImageSwapEndian(IMAGE *I)
{
  DCPIX *dcpix, dcval;
  CPIX *cpix, cval;
  double *dpix, dval;
  float *fpix, fval;
  long npix;

  npix = (long)I->numpix * (long)I->num_frame;

  switch (I->pixel_format) {
    case PFDBLCOM:
      dcpix = IMAGEDCpix(I, 0, 0);
      while (npix--) {
        dcval = *dcpix;
        dcpix->real = swapDouble(dcval.real);
        dcpix->imag = swapDouble(dcval.imag);
        dcpix++;
      }
      break;
    case PFCOMPLEX:
      cpix = IMAGECpix(I, 0, 0);
      while (npix--) {
        cval = *cpix;
        cpix->real = swapFloat(cval.real);
        cpix->imag = swapFloat(cval.imag);
        cpix++;
      }
      break;
    case PFDOUBLE:
      dpix = IMAGEDpix(I, 0, 0);
      while (npix--) {
        dval = *dpix;
        *dpix++ = swapDouble(dval);
      }
      break;
    case PFFLOAT:
      fpix = IMAGEFpix(I, 0, 0);
      while (npix--) {
        fval = *fpix;
        *fpix++ = swapFloat(fval);
      }
      break;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "ImageFRead: unsupported type %d\n", I->pixel_format);
  }
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageFRead(FILE *fp, const char *fname, int start, int nframes)
{
  int ecode, end_frame, frame, count = 1;
  IMAGE *I;
  byte *startpix, end = END_UNDEF;

  if (!fname) fname = "ImageFRead";

  I = (IMAGE *)calloc(1, sizeof(IMAGE));
  if (!I) ErrorExit(ERROR_NO_MEMORY, "ImageFRead: could not allocate header\n");

  ecode = fread_header(fp, I, fname);
  if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_header failed (%d)\n", ecode);

  if (endian == END_UNDEF) endian = FindMachineEndian();

  if (findparam(I, "endian")) getparam(I, "endian", PFBYTE, &count, &end);

  if (start < 0) /* read all frames */
  {
    start = 0;
    nframes = I->num_frame;
  }
  else /* read only specified frames */
  {
    if (fseek(fp, (long)I->sizeimage * (long)start, SEEK_CUR) < 0) {
      ImageFree(&I);
      ErrorReturn(NULL, (ERROR_BADFILE, "ImageFRead(%s, %d) - could not seek to specified frame", fname, start));
    }
  }

  if (nframes < 0) nframes = I->num_frame - start + 1;

  end_frame = start + nframes - 1;
  if (end_frame >= I->num_frame)
    ErrorReturn(NULL, (ERROR_BADFILE, "ImageFRead(%s, %d) - frame out of bounds", fname, end_frame));
  I->num_frame = nframes;
  if (ImageAllocBuffer(I) != NO_ERROR) ErrorExit(Gerror, "ImageAllocBuffer failed");

  startpix = I->image;
  for (frame = start; frame <= end_frame; frame++) {
    ecode = fread_image(fp, I, frame, fname);
    if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageFRead: fread_image failed (%d)\n", ecode);
    I->image += I->sizeimage;
  }
  I->image = startpix;

  /* We only swap endians if there wasn't an endian parameter and the image is
     invalid (ie.  values in the image seem to be extreme) OR the endian of
     the machine does not match the endian of the image */

  switch (end) {
    case END_UNDEF:
      if (!ImageValid(I)) ImageSwapEndian(I);
      break;
    case END_BIG:
    case END_SMALL:
      if (end != endian) ImageSwapEndian(I);
      break;
  }

  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageReadFrames(const char *fname, int start, int nframes)
{
  IMAGE *I;
  FILE *fp;

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadFrames(%s) fopen failed\n", fname));

  I = ImageFRead(fp, fname, start, nframes);
  fclose(fp);
  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageReadHeader(const char *fname)
{
  IMAGE *I = NULL;
  FILE *fp;
  int type, frame;
  char *buf;

  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname); /* don't destroy callers string */

  ImageUnpackFileName(buf, &frame, &type, buf);

  fp = fopen(buf, "rb");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadHeader(%s, %d) failed\n", buf, frame));

  I = ImageFReadHeader(fp, buf);
  fclose(fp);

  free(buf);
  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageFReadHeader(FILE *fp, const char *fname)
{
  IMAGE *I = NULL;
  int ecode;
  int type, frame;
  char* buf;

  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname); /* don't destroy callers string */

  ImageUnpackFileName(buf, &frame, &type, buf);

  I = (IMAGE *)calloc(1, sizeof(IMAGE));
  if (!I) ErrorExit(ERROR_NO_MEMORY, "ImageReadHeader: could not allocate header\n");

  switch (type) {
    case TIFF_IMAGE:
      if (TiffReadHeader(buf, I) == NULL) return (NULL);
      break;
    case RGBI_IMAGE:
      RGBReadHeader(buf, I);
      break;
    case MATLAB_IMAGE: {
      MATFILE mf;

      MatReadHeader0(fp, &mf);
      init_header(
          I, "matlab", "seq", 1, "today", (int)mf.mrows, (int)mf.ncols, mf.imagf ? PFCOMPLEX : PFFLOAT, 1, "temp");
    } break;
    case JPEG_IMAGE:
      JPEGReadHeader(fp, I);
      break;
    case PGM_IMAGE:
      PGMReadHeader(fp, I);
      break;
    case PPM_IMAGE:
      PPMReadHeader(fp, I);
      break;
    case PBM_IMAGE:
      PBMReadHeader(fp, I);
      break;
    case HIPS_IMAGE:
    default:
      ecode = fread_header(fp, I, buf);
      if (ecode != HIPS_OK) {
        fclose(fp);
        ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadHeader(%s): fread_header failed (%d)\n", buf, ecode));
      }
      break;
  }

  free(buf)
  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          read an image from file and convert it to the specified
          format.
------------------------------------------------------*/
IMAGE *ImageReadType(const char *fname, int pixel_format)
{
  IMAGE *Itmp, *I;

  Itmp = ImageRead(fname);
  if (!Itmp) ErrorReturn(NULL, (ERROR_NO_FILE, "ImageReadType(%s, %d): could not read image", fname, pixel_format));
  if (Itmp->pixel_format != pixel_format) {
    I = ImageAlloc(Itmp->rows, Itmp->cols, pixel_format, Itmp->num_frame);
    ImageCopy(Itmp, I);
    ImageFree(&Itmp);
  }
  else
    I = Itmp;

  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
IMAGE *ImageRead(const char *fname)
{
  IMAGE *I = NULL;
  MATRIX *mat;
  FILE *fp;
  int type, frame;
  char* buf;

  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname); /* don't destroy callers string */

  ImageUnpackFileName(buf, &frame, &type, buf);

  switch (type) {
    case TIFF_IMAGE:
      I = TiffReadImage(buf, frame);
      if (I == NULL) return (NULL);
      break;
    case MATLAB_IMAGE:
      DiagPrintf(
          DIAG_WRITE, "ImageRead: buf=%s, frame=%d, type=%d (M=%d,H=%d)\n", buf, frame, type, MATLAB_IMAGE, HIPS_IMAGE);
      mat = MatlabRead(buf);
      if (!mat) ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s) failed\n", buf));
      I = ImageFromMatrix(mat, NULL);
      ImageInvert(I, I);
      MatrixFree(&mat);
      break;
    case HIPS_IMAGE:
      fp = fopen(buf, "rb");
      if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "ImageRead(%s, %d) failed\n", buf, frame));
      I = ImageFRead(fp, buf, frame, 1);
      fclose(fp);
      break;
    case JPEG_IMAGE:
      I = JPEGReadImage(buf);
      break;
    case PGM_IMAGE:
      I = PGMReadImage(buf);
      break;
    case PPM_IMAGE:
      I = PPMReadImage(buf);
      break;
    case PBM_IMAGE:
      I = PBMReadImage(buf);
      break;
    case RGBI_IMAGE:
      I = RGBReadImage(buf);
    default:
      break;
  }

  free(buf);
  return (I);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ImageType(const char *fname)
{
  char *dot, *buf;

  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname);
  dot = strrchr(buf, '.');

  if (dot) {
    dot = StrUpper(dot + 1);
    if (!strcmp(dot, "MAT")) return (MATLAB_IMAGE);
  }

  free(buf);
  return (HIPS_IMAGE);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ImageFrame(const char *fname)
{
  char *number, *buf;
  int frame;

  buf = (char*)malloc(strlen(fname) + 1);
  strcpy(buf, fname);
  number = strrchr(buf, '#');

  if (number) {
    sscanf(number + 1, "%d", &frame);
    *number = 0;
  }
  else
    frame = 0;

  free(buf);
  return (frame);
}
/*----------------------------------------------------------------------
            Parameters:
              fname - the name of the file to read from

           Description:
              read a hips image from a file, and allocate an image
              header and data space for it.  Returns the newly
              allocated image.
----------------------------------------------------------------------*/
int ImageReadInto(const char *fname, IMAGE *I, int image_no)
{
  FILE *fp;
  int ecode;

  fp = fopen(fname, "rb");
  if (!fp) ErrorPrintf(ERROR_NO_FILE, "ImageReadInto(%s) failed\n", fname);

  ecode = fread_header(fp, I, fname);
  if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageReadInto(%s): fread_header failed (%d)\n", fname, ecode);
  ecode = fread_image(fp, I, image_no, fname);
  if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageReadInto(%s): fread_image failed (%d)\n", fname, ecode);

  fclose(fp);

  return (0);
}
/*----------------------------------------------------------------------
            Parameters:
              image - the image to write
              fname - the name of the file to write to.

           Description:
              write a hips image to file 'fname'
----------------------------------------------------------------------*/
int ImageWriteFrames(IMAGE *image, const char *fname, int start, int nframes)
{
  IMAGE *tmp_image;

  tmp_image = ImageAlloc(image->rows, image->cols, image->pixel_format, nframes);
  ImageCopyFrames(image, tmp_image, start, nframes, 0);
  ImageWrite(tmp_image, fname);
  ImageFree(&tmp_image);
  return (NO_ERROR);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              decompose a file name, extracting the type and the frame #.
----------------------------------------------------------------------*/
int ImageUnpackFileName(const char *inFname, int *pframe, int *ptype, char *outFname)
{
  char *number, *dot, *buf;

  buf = (char*)malloc(strlen(inFname) + strlen(outFname) + 1);
  if (inFname != outFname) strcpy(outFname, inFname);
  number = strrchr(outFname, '#');
  dot = strrchr(outFname, '.');

  if (number) /* : in filename indicates frame # */
  {
    if (sscanf(number + 1, "%d", pframe) < 1) *pframe = -1;
    *number = 0;
  }
  else
    *pframe = -1;

  if (dot) {
    dot = StrUpper(strcpy(buf, dot + 1));
    if (!strcmp(dot, "MAT"))
      *ptype = MATLAB_IMAGE;
    else if (!strcmp(dot, "TIF") || !strcmp(dot, "TIFF"))
      *ptype = TIFF_IMAGE;
    else if (!strcmp(dot, "JPG") || !strcmp(dot, "JPEG"))
      *ptype = JPEG_IMAGE;
    else if (!strcmp(dot, "PGM"))
      *ptype = PGM_IMAGE;
    else if (!strcmp(dot, "PPM"))
      *ptype = PPM_IMAGE;
    else if (!strcmp(dot, "PBM"))
      *ptype = PBM_IMAGE;
    else if (!strcmp(dot, "RGB"))
      *ptype = RGBI_IMAGE;
    else
      *ptype = HIPS_IMAGE;
  }
  else
    *ptype = HIPS_IMAGE;

  free(buf);
  return (NO_ERROR);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              return the number of frames stored in the file 'fname'
----------------------------------------------------------------------*/
int ImageNumFrames(const char *fname)
{
  IMAGE I;
  FILE *fp;
  int frame, type, ecode, nframes;
  char *buf;

  buf = (char*)malloc(strlen(fname) + 1);

  ImageUnpackFileName(fname, &frame, &type, buf);
  fname = buf;
  if ((frame >= 0) || (type != HIPS_IMAGE)) return (1);

  fp = fopen(fname, "rb");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "ImageNumFrame(%s) could not open file\n", fname));

  ecode = fread_header(fp, &I, fname);
  if (ecode != HIPS_OK) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "ImageNumFrame: fread_header failed (%d)\n", ecode));

  nframes = I.num_frame;
  fclose(fp);
  free_hdrcon(&I);
  free(buf);
  return (nframes);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              append an image to the end of a hips sequence file, incrementing
              the number of frames recorded in the header.
----------------------------------------------------------------------*/
int ImageAppend(IMAGE *I, const char *fname)
{
  FILE *fp;
  int ecode, frame = 0, nframes;
  IMAGE Iheader, *Iframe;
  // PW 2020/12/18 not ideal, but won't cause an overflow as long string size
  // here matches with `char *FileTmpName(const char *basename)` @ utils.cpp:629
  char tmpname[STR_LEN];

  fp = fopen(fname, "r+b");
#if 0
  if (!fp)
    ErrorReturn(-1, (ERROR_NO_FILE, "ImageAppend(%s) failed\n", fname)) ;
#endif

  if (!fp) return (ImageWrite(I, fname));

  ecode = fread_header(fp, &Iheader, fname);
  if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageAppend: fread_header failed (%d)\n", ecode);

  /* increment # of frames, and update file header */
  Iheader.num_frame++;
  nframes = Iheader.num_frame;
  if (nframes == 10 || nframes == 100 || nframes == 1000) {
    /* header size will grow by 1 byte, must copy whole file (ughhh) */
    fclose(fp);
    strcpy(tmpname, FileTmpName(NULL));
    FileRename(fname, tmpname);

    /* write out new header */
    fp = fopen(fname, "wb");
    if (!fp) ErrorReturn(-1, (ERROR_NO_FILE, "ImageAppend(%s) failed\n", fname));

    ecode = fwrite_header(fp, &Iheader, fname);
    if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageAppend: fwrite_header failed (%d)\n", ecode);

    nframes = Iheader.num_frame - 1;
    for (frame = 0; frame < nframes; frame++) {
      Iframe = ImageReadFrames(tmpname, frame, 1);
      if (!Iframe) ErrorReturn(-3, (ERROR_BADFILE, "ImageAppend: could not read %dth frame", frame));
      ecode = fwrite_image(fp, Iframe, frame, fname);
      if (ecode != HIPS_OK)
        ErrorReturn(-4, (ERROR_BADFILE, "ImageAppend: fwrite_image frame %d failed (%d)\n", ecode, frame));
    }
    unlink(tmpname);
  }
  else /* seek back to start and increment # of frames */
  {
    if (fseek(fp, 0L, SEEK_SET) < 0) ErrorReturn(-2, (ERROR_BADFILE, "ImageAppend(%s): could not seek to end"));
    ecode = fwrite_header(fp, &Iheader, fname);
    if (ecode != HIPS_OK) ErrorExit(ERROR_NO_FILE, "ImageAppend: fwrite_header failed (%d)\n", ecode);
  }

  if (fseek(fp, 0L, SEEK_END) < 0) ErrorReturn(-2, (ERROR_BADFILE, "ImageAppend(%s): could not seek to end"));

  ecode = fwrite_image(fp, I, frame, "fwrite");
  if (ecode != HIPS_OK)
    ErrorReturn(-1, (ERROR_BADFILE, "ImageAppend: fwrite_image frame %d failed (%d)\n", ecode, frame));

  free_hdrcon(&Iheader);
  fclose(fp);
  return (NO_ERROR);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
             Read a TIFF image from a file.
----------------------------------------------------------------------*/

static IMAGE *TiffReadImage(const char *fname, int frame0)
{
  IMAGE *I;
  TIFF *tif = TIFFOpen(fname, "r");
  int type = PFBYTE;  // just make compiler happy
  int width, height, ret, row;
  short nsamples, bits_per_sample;
  int nframe, frame;
  byte *iptr;
  tdata_t *buf;
  short photometric;
  int photometricInt;
  short fillorder;
  short compression;
  int compressionInt;
  short orientation;
  short resunit;
#if 0  // we used to translate RGB image into grey scale
//unsigned char     *buffer;
  int      skip;
  int      i;
  float    r, g, b, y;
  float    *pf;
#endif
  unsigned int scanlinesize;  //, extra_samples;
  int index = 0;
  float xres, yres, res;

  if (!tif) return (NULL);

  /* Find out how many frames we have */
  nframe = 1;  // note that TIFFOpen reads the 1st directory
  while (TIFFReadDirectory(tif)) nframe++;

  // some tif image just cannot be handled
  if (nframe == 0) nframe = 1;

  /* Go back to the beginning */
  TIFFSetDirectory(tif, 0);

  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH, &width);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH, &height);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_PHOTOMETRIC, &photometric);
  // fill order is LSB or MSB TIFFReadScanLine() handles it automatically
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_FILLORDER, &fillorder);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &compression);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_RESOLUTIONUNIT, &resunit);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_XRESOLUTION, &xres);
  if (ret == 0) xres = resunit == 1 ? .1 : (10 * 2.54);
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_YRESOLUTION, &yres);
  if (ret == 0) yres = resunit == 1 ? .1 : (10 * 2.54);

  // orientation
  // #define TIFFTAG_ORIENTATION             274     /* +image orientation */
  //    ORIENTATION_TOPLEFT         1       /* row 0 top, col 0 lhs */
  //    ORIENTATION_TOPRIGHT        2       /* row 0 top, col 0 rhs */
  //    ORIENTATION_BOTRIGHT        3       /* row 0 bottom, col 0 rhs */
  //    ORIENTATION_BOTLEFT         4       /* row 0 bottom, col 0 lhs */
  //    ORIENTATION_LEFTTOP         5       /* row 0 lhs, col 0 top */
  //    ORIENTATION_RIGHTTOP        6       /* row 0 rhs, col 0 top */
  //    ORIENTATION_RIGHTBOT        7       /* row 0 rhs, col 0 bottom */
  //    ORIENTATION_LEFTBOT         8       /* row 0 lhs, col 0 bottom */
  ret = TIFFGetFieldDefaulted(tif, TIFFTAG_ORIENTATION, &orientation);
  if (DIAG_VERBOSE_ON) {
    fprintf(stderr, "\ntiff info\n");
    fprintf(stderr, "         size: (%d, %d)\n", width, height);
    fprintf(stderr, "samples/pixel: %d\n", nsamples);
    fprintf(stderr, "  bits/sample: %d\n", bits_per_sample);
    fprintf(stderr, "  orientation: %d\n", orientation);
    photometricInt = photometric;  // used in 'case' statement to avoid
    // gcc warning 'case value out of range'
    switch (photometricInt) {
      case PHOTOMETRIC_MINISWHITE:
        fprintf(stderr, "  photometric: min value is white.\n");
        break;
      case PHOTOMETRIC_MINISBLACK:
        fprintf(stderr, "  photometric: min value is black.\n");
        break;
      case PHOTOMETRIC_RGB:
        fprintf(stderr, "  photometric: RGB color model.\n");
        break;
      case PHOTOMETRIC_PALETTE:
        fprintf(stderr, "  photometric: use palette.\n");
        break;
      case PHOTOMETRIC_MASK:
        fprintf(stderr, "  photometric: $holdout mask.\n");
        break;
      case PHOTOMETRIC_SEPARATED:
        fprintf(stderr, "  photometric: color separations.\n");
        break;
      case PHOTOMETRIC_YCBCR:
        fprintf(stderr, "  photometric: YCbCr6 CCIR 601.\n");
        break;
      case PHOTOMETRIC_CIELAB:
        fprintf(stderr, "  photometric: 1976 CIE L*a*b* \n");
        break;
      case PHOTOMETRIC_ITULAB:
        fprintf(stderr, "  photometric: ITU L*a*b* \n");
        break;
      case PHOTOMETRIC_LOGL:
        fprintf(stderr, "  photometric: CIE Log2(L) \n");
        break;
      case PHOTOMETRIC_LOGLUV:
        fprintf(stderr, "  photometric: CIE Log2(L) (u',v') \n");
        break;
      default:
        fprintf(stderr, "  photometric: unknown type\n");
        break;
    }
    compressionInt = compression;  // used in 'case' statement to avoid
    // gcc warning 'case value out of range'
    // we are not supporting compression at this time
    switch (compressionInt) {
      case COMPRESSION_NONE:
        fprintf(stderr, "  compression: no compression\n");
        break;
      case COMPRESSION_LZW:
        fprintf(stderr, "  compression: Lempel-Ziv & Welch\n");
        break;
      case COMPRESSION_JPEG:
        fprintf(stderr, "  compression: JPEG DCT compression\n");
        break;
      case COMPRESSION_PACKBITS:
        fprintf(stderr, "  compression: Macintosh RLE\n");
        break;
      default:
        fprintf(stderr, "  compression: %d see /usr/include/tiff.h for meaning\n", compression);
        break;
    }
  }
  // extra_samples = 0;
  switch (nsamples) {
    case 1:
      switch (bits_per_sample) /* not valid - I don't know why */
      {
        default:
        case 8:
          type = PFBYTE;
          break;
        case 16:
          type = PFSHORT;
          break;
        case 32:
          type = PFFLOAT;
          break;
        case 64:
          type = PFDOUBLE;
          break;
      }
      break;
    case 4:
      // The lack of break is deliberate and the 'fall through' comment is recognised by GCC
      // extra_samples = 1;
      nsamples = 3;
      // fall through
    case 3:
      switch (bits_per_sample) {
        default:
        case 8:
          type = PFRGB;
          break;
      }
      break;
    default:
      ErrorExit(ERROR_BADPARM,
                "IMAGE: nsamples=%d, bits_per_sample=%d.  only grey scale or RGB image is supported\n",
                nsamples,
                bits_per_sample);
  }
  // nsamples not handled here
  if (nsamples != 1 && nsamples != 3)
    ErrorExit(ERROR_BADPARM, "IMAGE: nsamples = %d.  only grey scale or RGB image is supported\n", nsamples);

  // type can be grey scale or RGB
  if (frame0 < 0)
    I = ImageAlloc(height, width, type, nframe);
  else
    I = ImageAlloc(height, width, type, 1);

  res = (xres + yres) / 2;
  switch (resunit) {
    case RESUNIT_CENTIMETER:  // 3 - cm
    case RESUNIT_NONE:        // 1 - no units
      I->sizepix = 100 / res;
      I->xsize = 100.0 / xres;  // mm
      I->ysize = 100.0 / yres;  // mm
      break;
    default:
    case RESUNIT_INCH:  // 2 - inches
      I->sizepix = 2.54 / res;
      I->xsize = 10.0 * 2.54 / xres;  // mm
      I->ysize = 10.0 * 2.54 / yres;  // mm
      break;
  }

  iptr = I->image;

  for (frame = 0; frame < nframe; frame++) {
    int planar_config, fillorder;
    TIFFSetDirectory(tif, frame);

    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_FILLORDER, &fillorder);
    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &planar_config);
    if (planar_config == PLANARCONFIG_SEPARATE)
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "TiffReadImage:  PLANARCONFIG_SEPARATE unsupported"));
    // else planar_config ==  PLANARCONFIG_CONTIG
    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH, &width);
    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH, &height);
    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
    ret = TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    scanlinesize = TIFFScanlineSize(tif);
    for (row = 0; row < height; row++) {
      // get the pointer at the first column of a row
      // note that the orientation is column, row
      switch (orientation) {
        case ORIENTATION_TOPLEFT:
          index = height - row - 1;
          break;
        case ORIENTATION_BOTLEFT:
          index = row;
          break;
        default:
          ErrorExit(ERROR_BADPARM, "IMAGE: orientation = %d. we support only topleft or botleft\n", orientation);
      }

      if (nsamples == 1) {
        switch (bits_per_sample) {
          default:
          case 8:
            buf = (tdata_t *)IMAGEpix(I, 0, index);
            break;
          case 16:
            buf = (tdata_t *)IMAGESpix(I, 0, index);
            break;
          case 32:
            buf = (tdata_t *)IMAGEFpix(I, 0, index);
            break;
          case 64:
            buf = (tdata_t *)IMAGEDpix(I, 0, index);
            break;
        }
        if (TIFFReadScanline(tif, buf, row, 0) < 0)  // row must be sequentially read for compressed data
          ErrorReturn(NULL, (ERROR_BADFILE, "TiffReadImage:  TIFFReadScanline returned error"));
        if (bits_per_sample == 1)  // unpack bitmap
        {
          unsigned char *bitmap, bitmask;
          unsigned int byte_, col, bit, b;
          bitmap = (unsigned char *)calloc(scanlinesize, sizeof(unsigned char));

          memmove(bitmap, buf, scanlinesize);
          for (col = b = 0; b < scanlinesize; b++, col += 8) {
            byte_ = bitmap[b];
            if (byte_ > 0) DiagBreak();
            if (fillorder == FILLORDER_LSB2MSB) {
              for (bitmask = 0x01, bit = 0; bit < 8; bit++) {
                if (col + bit == (unsigned)Gx && index == Gy) DiagBreak();
                *IMAGEpix(I, col + bit, index) = ((byte_ & bitmask) > 0);
                bitmask = bitmask << 1;
              }
            }
            else  // fillorder == FILLORDER_MSB2LSB
            {
              for (bitmask = 0x01 << 7, bit = 0; bit < 8; bit++) {
                if (col + bit == (unsigned)Gx && index == Gy) DiagBreak();
                *IMAGEpix(I, col + bit, index) = ((byte_ & bitmask) > 0);
                bitmask = bitmask >> 1;
              }
            }
          }

          free(bitmap);
        }
      }
      else if (nsamples == 4)  // RGB model + alpha
      {
        int s;
        unsigned char *ipix;
        ipix = (unsigned char *)calloc(scanlinesize, sizeof(unsigned char));
        buf = (tdata_t *)ipix;
        switch (bits_per_sample) {
          default:
          case 8:
            //          buf = (tdata_t*) IMAGERGBpix(I, 0, index);
            if (TIFFReadScanline(tif, buf, row, 0) < 0)  // row must be sequentially read for compressed data
              ErrorReturn(NULL, (ERROR_BADFILE, "TiffReadImage:  TIFFReadScanline returned error"));
        }
        for (s = 0; s < width; s++) {
          unsigned char *opix;
          opix = IMAGERGBpix(I, s, index);
          *opix++ = *ipix;
          *opix++ = *(ipix + 1);
          *opix++ = *(ipix + 2);
          ipix += nsamples;
        }
        free(buf);
      }
      else if (nsamples == 3)  // RGB model
      {
        switch (bits_per_sample) {
          default:
          case 8:
            buf = (tdata_t *)IMAGERGBpix(I, 0, index);
            if (TIFFReadScanline(tif, buf, row, 0) < 0)  // row must be sequentially read for compressed data
              ErrorReturn(NULL, (ERROR_BADFILE, "TiffReadImage:  TIFFReadScanline returned error"));
        }

#if 0
        ////////////////////////////////////////////////////////////
        // we used to translate into the grey value
        // now translate it into Y value
        // RGB range 0 to 1.0
        // then YIQ is
        //     Y   =  0.299  0.587   0.114  R
        //     I      0.596 -0.275  -0.321  G
        //     Q      0.212 -0.523   0.311  B
        // andd use Y for grey scale (this is color tv signal into bw tv
        switch (bits_per_sample)
        {
        default:
        case 8:
          skip = 3; //
          for (i = 0; i < width; ++i)
          {
            r = (float) buf[i*skip];
            g = (float) buf[i*skip+1];
            b = (float) buf[i*skip+2];
            y = (0.299*r + 0.587*g + 0.114*b);
            *IMAGEpix(I, i, row) = (unsigned char) y;
          }
          break; // 3 bytes at a time
        case 32:
          skip = 12; // 3x4 bytes at a time
          for (i=0; i < width ; ++i)
          {
            pf = (float *) &buf[i*skip];
            r = *pf;
            pf = (float *) &buf[i*skip+4];
            g = *pf;
            pf = (float *) &buf[i*skip+8];
            b = *pf;
            y = (0.299*r + 0.587*g + 0.114*b);
            *IMAGEFpix(I, i, row) = y;
          }
        case 64:
          ErrorExit(ERROR_BADPARM, "At this time we don't support RGB double valued tiff.\n");
        }
        free(buffer);
#endif
      }
    }
    if (frame0 < 0)
      I->image += I->sizeimage;
    else if (frame == frame0) /* only interested in one frame */
      break;
  }
  I->image = iptr;

  TIFFClose(tif);

  return (I);
}

#ifndef Darwin
#ifndef Windows_NT
/* unresolved in libtiff for some reason... */
void __eprintf(void);

void __eprintf(void) {}
#endif
#endif

/*----------------------------------------------------------------------
            Parameters:

           Description:
             Read the header info from a tiff image.
----------------------------------------------------------------------*/
static IMAGE *TiffReadHeader(const char *fname, IMAGE *I)
{
  TIFF *tif = TIFFOpen(fname, "r");
  int width, height, bits_per_sample;  //, ret, extra_samples;
  short nsamples;
  int type = PFBYTE;  // just make compiler happy
  if (!tif) return (NULL);

  TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetFieldDefaulted(tif, TIFFTAG_IMAGELENGTH, &height);
  TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &nsamples);
  // ret =
  TIFFGetFieldDefaulted(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);

  // extra_samples = 0;
  switch (nsamples) {
    case 1:
      switch (bits_per_sample) {
        case 32:
          type = PFFLOAT;
          break;
        case 16:
          type = PFSHORT;
          break;
        case 8:
        default:
          type = PFBYTE;
          break;
      }
      break;
    case 4:
    // extra_samples = 1;
    // no break
    case 3:
      switch (bits_per_sample) {
        default:
        case 8:
          type = PFRGB;
          break;
      }
      break;
    default:
      ErrorExit(ERROR_BADPARM, "IMAGE: nsamples = %d.  only grey scale or RGB image is supported\n", nsamples);
  }

  if (!I)
    I = ImageAlloc(height, width, type, 1);
  else
    init_header(I, "orig", "seq", 1, "today", height, width, type, 1, "temp");

  return (I);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Write an image to disk in TIFF format.
----------------------------------------------------------------------*/
static int TiffWriteImage(IMAGE *I, const char *fname, int frame)
{
  TIFF *out;
  short bits_per_sample, samples_per_pixel, sample_format;
  int row, frames;
  byte *timage;
  tdata_t *buf;

  out = TIFFOpen(fname, "w");
  if (out == NULL) return (ERROR_NOFILE);

  switch (I->pixel_format) {
    case PFBYTE:
      sample_format = SAMPLEFORMAT_INT;
      samples_per_pixel = 1;
      bits_per_sample = sizeof(byte) * 8;
      break;
    case PFSHORT:
      sample_format = SAMPLEFORMAT_INT;
      samples_per_pixel = 1;
      bits_per_sample = sizeof(short) * 8;
      break;
    case PFINT:
      sample_format = SAMPLEFORMAT_INT;
      samples_per_pixel = 1;
      bits_per_sample = sizeof(int) * 8;
      break;
    case PFFLOAT:
      sample_format = SAMPLEFORMAT_IEEEFP;
      samples_per_pixel = 1;
      bits_per_sample = sizeof(float) * 8;
      break;
    case PFDOUBLE:
      sample_format = SAMPLEFORMAT_IEEEFP;
      samples_per_pixel = 1;
      bits_per_sample = sizeof(double) * 8;
      break;
    default:
      ErrorReturn(ERROR_UNSUPPORTED,
                  (ERROR_UNSUPPORTED, "TiffWrite: pixel format %d not supported currently supported", I->pixel_format));
      samples_per_pixel = 3;
      sample_format = SAMPLEFORMAT_UINT;
      bits_per_sample = 8;
      break;
  }

  timage = I->image;

  for (frames = 0; frames < I->num_frame; frames++) {
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, (uint32)sample_format);
    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32)I->cols);
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32)I->rows);
    // orientation is bot-left
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
    if (!FZERO(I->xsize)) TIFFSetField(out, TIFFTAG_XRESOLUTION, (double)nint(2.54 * 10.0 / I->xsize));
    if (!FZERO(I->ysize)) TIFFSetField(out, TIFFTAG_YRESOLUTION, (double)nint(2.54 * 10.0 / I->ysize));

    TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    /* write out the data, line by line */
    for (row = 0; row < I->rows; row++) {
      switch (I->pixel_format) {
        case PFINT:
          buf = (tdata_t *)IMAGEIpix(I, 0, row);
          break;
        case PFBYTE:
          buf = (tdata_t *)IMAGEpix(I, 0, row);
          break;
        case PFSHORT:
          buf = (tdata_t *)IMAGESpix(I, 0, row);
          break;
        case PFFLOAT:
          buf = (tdata_t *)IMAGEFpix(I, 0, row);
          break;
        case PFDOUBLE:
          buf = (tdata_t *)IMAGEDpix(I, 0, row);
          break;
        default:
          ErrorReturn(ERROR_UNSUPPORTED, (ERROR_BADFILE, "TiffWrite: unsupported pixel format %d", I->pixel_format));
      }
      if (TIFFWriteScanline(out, buf, row, 0) < 0)
        ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "TiffWrite: TIFFWriteScanline returned error"));
    }
    TIFFWriteDirectory(out);
    I->image += I->sizeimage;
  }
  //  printf("closing TIFF file\n") ;
  TIFFClose(out);
  I->image = timage;

  return (NO_ERROR);
}

#ifndef HAVE_NETPBM

static IMAGE *PGMReadImage(const char *fname) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "pgm not supported")); }
static IMAGE *PGMReadHeader(FILE *fp, IMAGE *I) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "pgm not supported")); }
static int PGMWriteImage(IMAGE *I, const char *fname, int frame) { ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "pgm not supported")); }
static IMAGE *PPMReadImage(const char *fname) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ppm not supported")); }
static IMAGE *PPMReadHeader(FILE *fp, IMAGE *I) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "ppm not supported")); }
static IMAGE *PBMReadImage(const char *fname) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "pbm not supported")); }
static IMAGE *PBMReadHeader(FILE *fp, IMAGE *I) { ErrorReturn(NULL, (ERROR_UNSUPPORTED, "pbm not supported")); }
static int PPMWriteImage(IMAGE *I, const char *fname, int frame) { ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "ppm not supported")); }

#else  // HAVE_NETPBM

#include <netpbm/pbm.h>
#include <netpbm/pgm.h>
#include <netpbm/ppm.h>

static IMAGE *PGMReadHeader(FILE *fp, IMAGE *I)
{
  int rows, cols, format;
  gray maxval;

  pgm_readpgminit(fp, &cols, &rows, &maxval, &format);

  if (!I)
    I = ImageAlloc(rows, cols, PFBYTE, 1);
  else
    init_header(I, "orig", "seq", 1, "today", rows, cols, PFBYTE, 1, "temp");

  return I;
}

static IMAGE *PGMReadImage(const char *fname)
{
  FILE *infile;
  IMAGE *I;
  int rows, cols, format, i;
  gray maxval;

  if ((infile = pm_openr(fname)) == NULL) ErrorExit(ERROR_NO_FILE, "PGMReadImage:  Input file does not exist\n");

  pgm_readpgminit(infile, &cols, &rows, &maxval, &format);

  I = ImageAlloc(rows, cols, PFBYTE, 1);
  for (i = rows - 1; i >= 0; i--) pgm_readpgmrow(infile, IMAGEpix(I, 0, i), cols, maxval, format);

  pm_close(infile);

  return I;
}

static int PPMWriteImage(IMAGE *I, const char *fname, int frame)
{
  FILE *outf;
  int i, j;
  byte pval;
  pixel *cpix, *pp;
  if (I->pixel_format != PFBYTE)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "PPMWrite: only PFBYTE currently supported"));

  if ((outf = fopen(fname, "wb")) == NULL)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "PPMWrite: only PFBYTE currently supported"));

  cpix = (pixel *)malloc(sizeof(pixel) * I->ocols);
  if (!cpix) ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "Could not allocate color pixel buffer"));

  ppm_writeppminit(outf, I->ocols, I->orows, 255, 0);
  for (i = I->orows - 1; i >= 0; i--) {
    for (pp = cpix, j = 0; j < I->ocols; j++, pp++) {
      pval = *IMAGEpix(I, j, i);
      PPM_ASSIGN(*pp, pval, pval, pval);
    }
    ppm_writeppmrow(outf, cpix, I->ocols, 255, 0);
  }

  fclose(outf);

  free(cpix);

  return NO_ERROR;
}

static int PGMWriteImage(IMAGE *I, const char *fname, int frame)
{
  FILE *outf;
  int i;

  if (I->pixel_format != PFBYTE)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "PGMWrite: only PFBYTE currently supported"));

  if ((outf = fopen(fname, "wb")) == NULL)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "PGMWrite: could not open file"));

  pgm_writepgminit(outf, I->ocols, I->orows, 255, 0);
  for (i = I->orows - 1; i >= 0; i--) pgm_writepgmrow(outf, IMAGEpix(I, 0, i), I->ocols, 255, 0);

  fclose(outf);

  return NO_ERROR;
}

static IMAGE *PBMReadImage(const char *fname)
{
  FILE *infile;
  IMAGE *I;
  int rows, cols, i, j;
  bit **inbits;
  byte *ptr;

  if ((infile = pm_openr(fname)) == NULL) ErrorExit(ERROR_NO_FILE, "PGMReadImage:  Input file does not exist\n");

  inbits = pbm_readpbm(infile, &cols, &rows);

  pm_close(infile);

  I = ImageAlloc(rows, cols, PFBYTE, 1);
  ptr = I->image;

  for (j = rows - 1; j >= 0; j--)
    for (i = 0; i < cols; i++) *ptr++ = (byte)((inbits[j][i] == PBM_WHITE) ? 255 : 0);

  return I;
}

static IMAGE *PPMReadImage(const char *fname)
{
  FILE *infile;
  IMAGE *I;
  int rows, cols, format, i, j;
  pixval maxval;
  pixel *pixelrow, *pptr;
  byte *ptr;

  if ((infile = pm_openr(fname)) == NULL) ErrorExit(ERROR_NO_FILE, "PGMReadImage:  Input file does not exist\n");

  ppm_readppminit(infile, &cols, &rows, &maxval, &format);

  pixelrow = ppm_allocrow(cols);

  I = ImageAlloc(rows, cols, PFBYTE, 1);

  for (i = rows - 1; i >= 0; i--) {
    ppm_readppmrow(infile, pixelrow, cols, maxval, format);
    pptr = pixelrow;
    ptr = IMAGEpix(I, 0, i);
    for (j = 0; j < cols; j++, pptr++, ptr++) *ptr = (byte)(PPM_LUMIN(*pptr) + 0.5);
  }

  pm_close(infile);

  return I;
}

static IMAGE *PBMReadHeader(FILE *fp, IMAGE *I)
{
  int rows, cols, format;

  pbm_readpbminit(fp, &cols, &rows, &format);

  if (!I)
    I = ImageAlloc(rows, cols, PFBYTE, 1);
  else
    init_header(I, "orig", "seq", 1, "today", rows, cols, PFBYTE, 1, "temp");

  return I;
}

static IMAGE *PPMReadHeader(FILE *fp, IMAGE *I)
{
  int rows, cols, format;
  gray maxval;

  ppm_readppminit(fp, &cols, &rows, &maxval, &format);

  if (!I)
    I = ImageAlloc(rows, cols, PFBYTE, 1);
  else
    init_header(I, "orig", "seq", 1, "today", rows, cols, PFBYTE, 1, "temp");

  return I;
}

#endif  // HAVE_NETPBM

#define JPEG_INTERNALS

extern "C" {
#include "jinclude.h"
#include "jpeglib.h"
#include "jmorecfg.h"
}

static IMAGE *JPEGReadImage(const char *fname)
{
  FILE *infile;
  IMAGE *I;
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW  ptr;
  int rowctr ;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  if ((infile = fopen(fname, "rb")) == NULL) ErrorExit(ERROR_NO_FILE, "JPEGReadImage:  Input file does not exist\n");

  jpeg_stdio_src(&cinfo, infile);
  jpeg_read_header(&cinfo, TRUE);

  cinfo.out_color_space = JCS_GRAYSCALE;
  cinfo.out_color_space = JCS_RGB ;

  jpeg_start_decompress(&cinfo);
  I = ImageAlloc(cinfo.output_height, cinfo.output_width, PFRGB, 1);

  rowctr = I->orows - 1;
  while (cinfo.output_scanline < cinfo.output_height) {
    
    ptr = IMAGERGBpix(I, 0, rowctr--);
    jpeg_read_scanlines(&cinfo, &ptr, 1);
  }

  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  return I;
}

static IMAGE *JPEGReadHeader(FILE *fp, IMAGE *I)
{
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);

  if (!I)
    I = ImageAlloc(cinfo.image_height, cinfo.image_width, PFBYTE, 1);
  else
    init_header(I, "orig", "seq", 1, "today", cinfo.image_height, cinfo.image_width, PFBYTE, 1, "temp");

  jpeg_destroy_decompress(&cinfo);

  return (I);
}

static int JPEGWriteImage(IMAGE *I, const char *fname, int frame)
{
  FILE *outf;
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPROW ptr;

  if (I->pixel_format != PFBYTE)
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "JPEGWrite: only PFBYTE currently supported"));

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  if ((outf = fopen(fname, "wb")) == NULL) ErrorExit(ERROR_BAD_FILE, "JPEGWriteImage:  Could not open file");

  jpeg_stdio_dest(&cinfo, outf);

  cinfo.image_width = I->ocols;
  cinfo.image_height = I->orows;
  cinfo.input_components = 1;
  cinfo.in_color_space = JCS_GRAYSCALE;

  jpeg_set_defaults(&cinfo);

  jpeg_start_compress(&cinfo, TRUE);

  while (cinfo.next_scanline < cinfo.image_height) {
    ptr = IMAGEpix(I, 0, (I->rows - cinfo.next_scanline - 1));
    jpeg_write_scanlines(&cinfo, &ptr, 1);
  }

  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  fclose(outf);

  return NO_ERROR;
}

static IMAGE *RGBReadHeader(char *fname, IMAGE *I) {
  RGB_IMAGE *rgb;

  rgb = iopen(fname, "r", 0, 0, 0, 0, 0);

  if (!I)
    I = ImageAlloc(rgb->ysize, rgb->xsize, PFRGB, 1) ;
  else
    init_header(I, "orig", "seq", 1, "today", rgb->ysize,
                rgb->xsize, PFRGB, 1, "temp");

  iclose(rgb);

  return(I) ;
}

static IMAGE *RGBReadImage(char *fname) {
  IMAGE *I;
  RGB_IMAGE *rgb;
  unsigned short rows,cols,*r,*g,*b,i,j,*tr,*tg,*tb;
  byte *iptr;

  rgb = iopen(fname, "r", 0, 0, 0, 0, 0);
  rows = rgb->ysize;
  cols = rgb->xsize;

  if (rgb->zsize>3)
    ErrorReturn(NULL, (ERROR_BAD_PARM,
                       "Too many color planes in RGBReadImage (%s)\n",fname));

  I = ImageAlloc(rows, cols, PFRGB, 1);

  if ((r = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  if ((g = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  if ((b = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  iptr = I->image;

  for (i=0;i<rows;i++) {
    getrow(rgb,r,i,0); /* Red */
    getrow(rgb,g,i,1); /* Green */
    getrow(rgb,b,i,2); /* Blue */

    /* Translate color planes to RGB format */
    tr = r;
    tg = g;
    tb = b;
    for (j=0;j<cols;j++) {
      *iptr++ = *tr++;
      *iptr++ = *tg++;
      *iptr++ = *tb++;
    }
  }

  free(r);
  free(g);
  free(b);
  iclose(rgb);

  return I;
}

static int RGBwrite(IMAGE *I, char *fname, int frame) {
  RGB_IMAGE  *image ;
  int    x, y ;
  unsigned short *r ;

#ifndef Linux
  image = iopen(fname,"w",RLE(1), 2, I->cols, I->rows, 1);
#else
  image = iopen(fname,"w",UNCOMPRESSED(1), 2, I->cols, I->rows, 1);
#endif
  r = (unsigned short *)calloc(I->cols, sizeof(unsigned short)) ;
  for (y = 0 ; y < I->rows; y++) {
    for (x = 0 ; x < I->cols ; x++)
      r[x] = (unsigned short)(*IMAGEpix(I, x, y)) ;

    /* fill rbuf, gbuf, and bbuf with pixel values */
    putrow(image, r, y, 0);    /* red row */
    putrow(image, r, y, 1);    /* green row */
    putrow(image, r, y, 2);    /* blue row */
  }
  iclose(image);
  free(r) ;
  return(NO_ERROR) ;
}
