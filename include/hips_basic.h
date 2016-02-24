/**
 * @file  hips_basic.h
 * @brief basic definitions for HIPS
 *
 */
/*
 * Original Author: Michael Landy - 12/28/90
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2016/02/24 16:23:27 $
 *    $Revision: 1.6 $
 *
 * Copyright (c) 1991 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 *
 */


#ifndef HIPS_BASIC_H
#define HIPS_BASIC_H

/*
 * Machine-dependent portion
 *
 * The next lines are the only portion of the file which should be tailored
 * to an individual installation.
 */

typedef unsigned long fs_hsize_t; /* variable which can hold the size of an
        image in bytes */
#ifdef AIX
#define HPUXORAIX
#endif

#ifdef HPUX
#define HPUXORAIX
#endif

#ifdef HPUXORAIX

#define H__RANDOM lrand48 /* routine to call for random numbers */
#define H__RANDTYPE long /* type of H__RANDOM() */
#define H__SRANDOM srand48 /* routine to call to set the seed */
#define H__MAXRAND (0x7fffffff) /* maximum random number */
#define H__RANDBITS (31) /* number of random bits returned */

#else

#define H__RANDOM random /* routine to call for random numbers */
#define H__RANDTYPE long /* type of H__RANDOM() */
#define H__SRANDOM srandom /* routine to call to set the seed */
#define H__MAXRAND (0x7fffffff) /* maximum random number */
#define H__RANDBITS (31) /* number of random bits returned */

#endif

/* *******************END OF MACHINE-DEPENDENT PORTION*********************/

typedef unsigned char byte;
typedef char  sbyte;
typedef unsigned short h_ushort;
typedef unsigned int h_uint;
typedef float  h_complex[2];
typedef double  h_dblcom[2];
typedef const char *  Filename;

union pixelval {
  byte v_byte;
  sbyte v_sbyte;
  short v_short;
  h_ushort v_ushort;
  int v_int;
  h_uint v_uint;
  float v_float;
  double v_double;
  h_complex v_complex;
  h_dblcom v_dblcom;
};

typedef union pixelval Pixelval;

/*
 * For general readability
 */

#ifndef TRUE
# define TRUE 1
# define FALSE 0
#endif

typedef int h_boolean;

/*
 * Histogram structure
 *
 * The zero-th bin is underflows.  The last bin is overflows. So there are
 * nbins+2 slots.  The n'th bin counts pixels such that:
 *
 * min + ((n-1)*binwidth) <= value < min + n*binwidth
 *
 * For complex and double complex images, the complex magnitude is
 * histogrammed, and min/binwidth are either floats (for complex images) or
 * doubles (for double complex images).
 */

struct hips_histo
{
  int nbins;
  int *histo;
  fs_hsize_t sizehist;
  h_boolean histodealloc;
  int pixel_format;
  Pixelval minbin;
  Pixelval binwidth;
};

/*
 * Statistics structure
 *
 * The variable nelem counts the number of pixels that contributed to  these
 * image statistics (which might be less than the number of pixels in the
 * region-of-interest if zero-valued pixels aren't included.
 */

struct hips_stats
{
  int nelem;
  int pixel_format;
  Pixelval statmin;
  Pixelval statmax;
  double sum,ssq,mean,var,stdev;
};

/*
 * Convolution mask set structure
 */

struct hips_mask
{
  char *name;  /* name of this mask set */
  int nmasks;  /* number of masks */
  int func_num;  /* function applied to mask outputs */
  int pixel_format; /* format of mask elements */
  union {
    float **f_values; /* float mask pointers */
    int **i_values;  /* int mask pointers */
  } vals;
  int *mask_rows;  /* number of rows in each mask */
  int *mask_cols;  /* number of columns in each mask */
  int *row_offset; /* row number of mask value overlying image
             pixel */
  int *col_offset; /* column number of mask value overlying image
             pixel */
};

/*
 * Mask function numbers
 */

#define MASKFUN_MAXABS          1
#define MASKFUN_MEANSQ          2
#define MASKFUN_SUMABS          3
#define MASKFUN_MAX             4
#define MASKFUN_MAXFLR          5
#define MASKFUN_MXASFLR         6
#define MASKFUN_MUL             7
#define MASKFUN_NORM            8
#define MASKFUN_DIFF            9
#define MASKFUN_ORIENT          10
#define MASKFUN_IDENT           11
#define MASKFUN_MAXMASKS        11

/*
 * Filter types and structure
 *
 * A bandpass filter is a concatenation (i.e., product) of a lowpass (using
 * highcut/highorder) and a highpass (using lowcut/loworder) filter.  A band
 * reject filter is one minus the corresponding bandpass filter.
 */

#define FILTMETHOD_IDEAL 1
#define FILTMETHOD_BUTTERWORTH 2
#define FILTMETHOD_EXPONENTIAL 3

#define FILTDIST_ROW  1
#define FILTDIST_COL  2
#define FILTDIST_BOTH  3

#define FILTTYPE_LOWPASS 1
#define FILTTYPE_HIGHPASS 2
#define FILTTYPE_BANDPASS 3
#define FILTTYPE_BANDREJ 4

struct hips_filter
{
  int method;  /* Ideal/Butterworth/Exponential */
  int disttype;  /* scale by number of rows, columns or both */
  int ftype;  /* lowpass/highpass/bandpass/bandreject */
  double dmetric;  /* Minkowdki metric */
  double lowcut;
  int loworder;
  double highcut;
  int highorder;
};

char *strsave(),*memalloc(),*formatheader(),*formatheadera();
char *formatheaderc(),*hformatname(),*hformatname_f(),*hformatname_t();
byte *halloc(),*hmalloc();
fs_hsize_t hsizepix();
struct extpar *findparam(),*grepparam();
FILE *hfopenr(),*ffopen(),*ffreopen();
/* h_boolean getline(),swallownl(),hfgets(); */

/*
 * avoid hassles of including string.h or strings.h
 */

#if 0
extern char *strcat(),*strncat(),*strcpy(),*strncpy(),*index(),*rindex();
extern char *strchr(),*strdup(),*strstr();
extern int strcmp(),strncmp();
#else
#include <string.h>
#endif

/* omit strlen so that it defaults to int, but can also be size_t as in gcc */

/*
 * image and pyramid type declarations for the pyramid routines.
 *
 * The pyramid utilities are derived from code originally written by
 * Raj Hingorani at SRI/David Sarnoff Research Institute.  The original
 * Gaussian and Laplacian pyramid algorithms were designed by Peter Burt (also
 * currently at SRI/DSRC).  See:  Computer Graphics and Image Processing,
 * Volume 16, pp. 20-51, 1981, and IEEE Transactions on Communications,
 * Volume COM-31, pp. 532-540, 1983.
 */

#define MAXLEV 12

typedef struct
{
  float **ptr;
  int nr;
  int nc;
}
FIMAGE;

typedef struct
{
  int **ptr;
  int nr;
  int nc;
}
IIMAGE;

typedef FIMAGE FPYR[MAXLEV];
typedef IIMAGE IPYR[MAXLEV];

typedef struct
{
  float *k;
  int taps2;  /* the number of taps from the center rightward,
            total number is 2*taps2-1 */
}
FILTER;

/* function definitions */

float  **_read_fimgstr();
int  **_read_iimgstr();
float  **_alloc_fimage();
int  **_alloc_iimage();

/* externals */

extern int Image_border;

/* image macros */

#ifndef MAX
# define MAX(A,B)  ((A) > (B) ? (A) : (B))
#endif
#ifndef MIN
# define MIN(A,B)  ((A) < (B) ? (A) : (B))
#endif
#ifndef ABS
# define ABS(A)    ((A) > 0 ? (A) : (-(A)))
#endif
#ifndef BETWEEN
# define BETWEEN(A,B,C) (((A) < (B)) ? (B) : (((A) > (C)) ? (C) : (A)))
#endif
#ifndef SIGN
# define SIGN(A,B) (((B) > 0) ? (A) : (-(A)))
#endif
#ifndef TOascii
# define TOascii(c) ((c) & 0x7f)
#endif

/* compatibilities, type lists, etc. */

#define LASTTYPE -1 /* the last type in a type list */

#define CM_ROWS  01 /* check compatibility: ROI rows */
#define CM_COLS  02 /* check compatibility: ROI cols */
#define CM_FRAMES 04 /* check compatibility: frames & depths */
#define CM_FORMAT 010 /* check compatibility: pixel format */
#define CM_NUMCOLOR 020 /* check compatibility: numcolor */
#define CM_NUMLEV 040 /* check compatibility: pyramid levels */
#define CM_OROWS 0100 /* check compatibility: total # rows */
#define CM_OCOLS 0200 /* check compatibility: total # cols */
#define CM_FRAMESC 0400 /* check compatibility: check frames & depths
if numcolor != 1 or numdepth != 1 */
#define CM_NUMCOLOR3 01000 /* check compatibility: numcolor (treat
RGB, etc. as if numcolor=3) */
#define CM_DEPTH 02000 /* check compatibility: numdepth */

  /* converting to packed formats */

#define MSBF_PACKING 1 /* use most significant bit first packing */
#define LSBF_PACKING 2 /* use least significant bit first packing */

  /* converting complex numbers to single-valued numbers */

#define CPLX_MAG 1 /* complex magnitude */
#define CPLX_PHASE 2 /* complex phase */
#define CPLX_REAL 3 /* complex - real part only */
#define CPLX_IMAG 4 /* complex - imaginary part only */

  /* converting single-valued numbers to complex numbers */

#define CPLX_RVI0 1 /* real part = value, imaginary = 0 */
#define CPLX_R0IV 2 /* real part = 0, imaginary = value */
#define CPLX_RVIV 3 /* real part = value, imaginary = same value */

  /*
   * type conversion methods
   *
   * Note: because find_method returns a method identifier, or METH_IDENT,
   * or the negative of a method identifier (for conversion via PFINT), it is
   * essential that none of these possible values be identical to HIPS_ERROR so
   * that it also can give a normal hips error return.
   */

#define METH_IDENT 2
#define METH_BYTE 3
#define METH_COMPLEX 4
#define METH_DOUBLE 5
#define METH_DBLCOM 6
#define METH_FLOAT 7
#define METH_INT 8
#define METH_LSBF 9
#define METH_MSBF 10
#define METH_SHORT 11
#define METH_SBYTE 12
#define METH_UINT 13
#define METH_USHORT 14
#define METH_RGB 15
#define METH_RGBZ 16
#define METH_ZRGB 17
#define METH_BGR 18
#define METH_BGRZ 19
#define METH_ZBGR 20

  /* conversion-related variables */

  extern int hips_packing;
    extern byte hips_lchar;
      extern byte hips_hchar;
        extern int hips_cplxtor;
          extern int hips_rtocplx;
            extern int hips_convback;
              extern int hips_lclip,hips_hclip;
                extern int hips_zdiv;
                  extern h_boolean hips_oldhdr;

                    /* header handling */

                    extern int hips_fullhist;
                      extern int hips_fulldesc;
                        extern int hips_fullxpar;

                          struct h_types
                          { /* the type names structure defined in htypes.c */
                            int h_pfmt;  /* pixel format */
char *h_fmtname; /* sprintf string for error code */
};

struct h_convs
{ /* the conversion names structure defined in htype.c */
  int h_cnvtype;  /* pixel format */
char *h_cnvname; /* sprintf string for error code */
};

/* useful constants and functions */

#define RND01           ((double)H__RANDOM()/(double)H__MAXRAND)
#define RNDPM           (RND01 > 0.5 ? 1.0 : -1.0)
#ifndef SQR
#define SQR(x)          ((x)*(x))
#endif
#define H_E  2.7182818284590452354
#define H_LOG2E  1.4426950408889634074
#define H_LOG10E 0.43429448190325182765
#define H_LN2  0.69314718055994530942
#define H_LN10  2.30258509299404568402
#define H_PI  3.14159265358979323846
#define H_2PI  6.28318530717958647692
#define H_PI_2  1.57079632679489661923
#define H_PI_4  0.78539816339744830962
#define H_3PI_2         4.71238898038468985769
#define H_ONE_PI 0.31830988618379067154
#define H_TWO_PI 0.63661977236758134308
#define H_TWO_SQRTPI 1.12837916709551257390
#define H_180_PI 57.2957795130823208768 /* degrees to radians */
#define H_SQRT2  1.41421356237309504880
#define H_SQRT1_2 0.70710678118654752440
#define H_SQRT3OVER2 0.86602540378443864676

#endif
