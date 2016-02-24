/**
 * @file  hips_header.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/24 16:23:27 $
 *    $Revision: 1.5 $
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


#ifndef HIPS_HEADER_H
#define HIPS_HEADER_H

/*
 * Copyright (c) 1991 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 */

/*
 * hips_header.h - definitions related to the HIPS header
 *
 * Michael Landy - 12/28/90
 */

/* The HIPS header as stored in memory */

struct header
{
  char *orig_name; /* The originator of this sequence */
  h_boolean ondealloc; /* If nonzero, free orig_name when requested */
  char *seq_name; /* The name of this sequence */
  h_boolean sndealloc; /* If nonzero, free seq_name when requested */
  int num_frame; /* The number of frames in this sequence */
  char *orig_date; /* The date the sequence was originated */
  h_boolean oddealloc; /* If nonzero, free orig_date when requested */
  int orows;  /* The number of rows in each stored image */
  int ocols;  /* The number of columns in each stored image */
  int rows;  /* The number of rows in each image (ROI) */
  int cols;  /* The number of columns in each image (ROI) */
  int frow;  /* The first ROI row */
  int fcol;  /* The first ROI col */
  int pixel_format; /* The format of each pixel */
  int numcolor; /* The number of color frames per image */
  int numpix;  /* The number of pixels per stored frame */
  fs_hsize_t sizepix; /* The number of bytes per pixel */
  fs_hsize_t sizeimage; /* The number of bytes per stored frame */
  byte *image;  /* The image itself */
  h_boolean imdealloc; /* if nonzero, free image when requested */
  byte *firstpix; /* Pointer to first pixel (for ROI) */
  int sizehist; /* Number of bytes in history (excluding
             null, including <newline>) */
  char *seq_history; /* The sequence's history of transformations */
  h_boolean histdealloc; /* If nonzero, free history when requested */
  int sizedesc; /* Number of bytes in description (excluding
             null, including <newline>) */
  char *seq_desc; /* Descriptive information */
  h_boolean seqddealloc; /* If nonzero, free desc when requested */
  int numparam; /* Count of additional parameters */
  h_boolean paramdealloc; /* If nonzero, free param structures and/or
             param values when requested */
  struct extpar *params; /* Additional parameters */
  float xsize ;
  float ysize ;
};

struct hips_roi
{
  int rows;  /* The number of rows in the ROI */
  int cols;  /* The number of columns in the ROI */
  int frow;  /* The first ROI row */
  int fcol;  /* The first ROI col */
};

/*
 * Pixel Format Codes
 */

#define PFBYTE  0 /* Bytes interpreted as unsigned integers */
#define PFSHORT  1 /* Short integers (2 bytes) */
#define PFINT  2 /* Integers (4 bytes) */
#define PFFLOAT  3 /* Float's (4 bytes)*/
#define PFCOMPLEX  4 /* 2 Float's interpreted as (real,imaginary) */
#define PFASCII  5 /* ASCII rep, with linefeeds after each row */
#define PFDOUBLE  6 /* Double's (8 byte floats) */
#define PFDBLCOM  7 /* Double complex's (2 Double's) */
#define PFQUAD  10 /* quad-tree encoding (Mimaging) */
#define PFQUAD1  11 /* quad-tree encoding */
#define PFHIST  12 /* histogram of an image (using ints) */
#define PFSPAN  13 /* spanning tree format */
#define PLOT3D  24 /* plot-3d format */
#define PFMSBF  30 /* packed, most-significant-bit first */
#define PFLSBF  31 /* packed, least-significant-bit first */
#define PFSBYTE  32 /* signed bytes */
#define PFUSHORT 33 /* unsigned shorts */
#define PFUINT  34 /* unsigned ints */
#define PFRGB  35 /* RGB RGB RGB bytes */
#define PFRGBZ  36 /* RGB0 RGB0 RGB0 bytes */
#define PFZRGB  37 /* 0RGB 0RGB 0RGB bytes */
#define PFMIXED  40 /* multiple frames in different pixel formats */
#define PFBGR  41 /* BGR BGR BGR bytes */
#define PFBGRZ  42 /* BGR0 BGR0 BGR0 bytes */
#define PFZBGR  43 /* 0BGR 0BGR 0BGR bytes */
#define PFINTPYR 50 /* integer pyramid */
#define PFFLOATPYR 51 /* float pyramid */
#define PFPOLYLINE 100 /* 2D points */
#define PFCOLVEC 101 /* Set of RGB triplets defining colours */
#define PFUKOOA  102 /* Data in standard UKOOA format */
#define PFTRAINING 104 /* Set of colour vector training examples */
#define PFTOSPACE 105 /* TOspace world model data structure */
#define PFSTEREO 106 /* Stereo sequence (l, r, l, r, ...) */
#define PFRGPLINE 107 /* 2D points with regions */
#define PFRGISPLINE 108 /* 2D points with regions and interfaces */
#define PFCHAIN  200 /* Chain code encoding (Mimaging) */
#define PFLUT  300 /* LUT format (uses Ints) (Mimaging) */
#define PFAHC  400 /* adaptive hierarchical encoding */
#define PFOCT  401 /* oct-tree encoding */
#define PFBT  402 /* binary tree encoding */
#define PFAHC3  403 /* 3-d adaptive hierarchical encoding */
#define PFBQ  404 /* binquad encoding */
#define PFRLED  500 /* run-length encoding */
#define PFRLEB  501 /* run-length encoding, line begins black */
#define PFRLEW  502 /* run-length encoding, line begins white */
#define PFPOLAR  600 /* rho-theta format (Mimaging) */
#define PFGRLE  601 /* gray scale run-length encoding */
#define PFSRLE  602 /* monochrome run-scale encoding */
#define PFVFFT3D 701 /* float complex 3D virtual-very fast FT */
#define PFVFFT2D 702 /* float complex 2D virtual-very fast FT */
#define PFDVFFT3D 703 /* double complex 3D VFFT */
#define PFDVFFT2D 704 /* double complex 2D VFFT */
#define PFVVFFT3D 705 /* float 3D VFFT in separated planes */
#define PFDVVFFT3D 706 /* double 3D VVFFT in separated planes */


/* The Extended Parameters Structure */

struct extpar
{
  char *name;  /* name of this variable */
  int format;  /* format of values (PFBYTE, PFINT, etc.) */
  int count;  /* number of values */
  union {
    byte v_b; /* PFBYTE/PFASCII, count = 1 */
    int v_i; /* PFINT, count = 1 */
    short v_s; /* PFSHORT, count = 1 */
    float v_f; /* PFFLOAT, count = 1 */
    byte *v_pb; /* PFBYT/PFASCIIE, count > 1 */
    int *v_pi; /* PFINT, count > 1 */
    short *v_ps; /* PFSHORT, count > 1 */
    float *v_pf; /* PFFLOAT, count > 1 */
  } val;
  h_boolean dealloc; /* if nonzero, free memory for val */
  struct extpar *nextp; /* next parameter in list */
};

#define FBUFLIMIT 30000  /* increase this if you use large PLOT3D
files */
#define LINELENGTH 200  /* max characters per line in header vars */
#define NULLPAR ((struct extpar *) 0)

#endif
