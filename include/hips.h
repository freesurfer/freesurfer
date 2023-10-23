#ifndef HIPS_H
#define HIPS_H

// ------ typedefs ------

typedef unsigned char ubyte;
typedef char sbyte;
typedef unsigned short h_ushort;
typedef unsigned int h_uint;
typedef float h_complex[2];
typedef double h_dblcom[2];
typedef unsigned long fs_hsize_t;
typedef const char *  Filename;
typedef int h_boolean;

// ------ macros ------

#define MSBF_PACKING 1 /* use most significant bit first packing */
#define LSBF_PACKING 2 /* use least significant bit first packing */

#define CPLX_MAG 1 /* complex magnitude */
#define CPLX_PHASE 2 /* complex phase */
#define CPLX_REAL 3 /* complex - real part only */
#define CPLX_IMAG 4 /* complex - imaginary part only */

#define CPLX_RVI0 1 /* real part = value, imaginary = 0 */
#define CPLX_R0IV 2 /* real part = 0, imaginary = value */

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

#define HIPS_ERROR -1 /* error-return from hips routines */
#define HIPS_OK  0 /* one possible nonerror-return value */

#define FBUFLIMIT 30000  /* increase this if you use large PLOT3D files */
#define LINELENGTH 200  /* max characters per line in header vars */
#define NULLPAR ((struct extpar *) 0)

// ------ math macros ------

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

// ------ structs ------

union pixelval {
  ubyte v_byte;
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

struct header {
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
  ubyte *image;  /* The image itself */
  h_boolean imdealloc; /* if nonzero, free image when requested */
  ubyte *firstpix; /* Pointer to first pixel (for ROI) */
  int sizehist; /* Number of bytes in history (excluding null, including <newline>) */
  char *seq_history; /* The sequence's history of transformations */
  h_boolean histdealloc; /* If nonzero, free history when requested */
  int sizedesc; /* Number of bytes in description (excluding null, including <newline>) */
  char *seq_desc; /* Descriptive information */
  h_boolean seqddealloc; /* If nonzero, free desc when requested */
  int numparam; /* Count of additional parameters */
  h_boolean paramdealloc; /* If nonzero, free param structures and/or param values when requested */
  struct extpar *params; /* Additional parameters */
  float xsize ;
  float ysize ;
};

struct hips_roi {
  int rows;  /* The number of rows in the ROI */
  int cols;  /* The number of columns in the ROI */
  int frow;  /* The first ROI row */
  int fcol;  /* The first ROI col */
};

struct extpar {
  char *name;  /* name of this variable */
  int format;  /* format of values (PFBYTE, PFINT, etc.) */
  int count;  /* number of values */
  union {
    ubyte v_b; /* PFBYTE/PFASCII, count = 1 */
    int v_i; /* PFINT, count = 1 */
    short v_s; /* PFSHORT, count = 1 */
    float v_f; /* PFFLOAT, count = 1 */
    ubyte *v_pb; /* PFBYT/PFASCIIE, count > 1 */
    int *v_pi; /* PFINT, count > 1 */
    short *v_ps; /* PFSHORT, count > 1 */
    float *v_pf; /* PFFLOAT, count > 1 */
  } val;
  h_boolean dealloc; /* if nonzero, free memory for val */
  struct extpar *nextp; /* next parameter in list */
};

struct hips_histo {
  int nbins;
  int *histo;
  fs_hsize_t sizehist;
  h_boolean histodealloc;
  int pixel_format;
  Pixelval minbin;
  Pixelval binwidth;
};

// ------ functions ------

double h_entropy(int *table,int count,int pairflag);
int alloc_histo(struct hips_histo *histo,union pixelval *min,union pixelval *max,int nbins,int format);
int alloc_histobins(struct hips_histo *histo);
int clearparam(struct header *hd, const char *name);
int fread_header(FILE *fp,struct header *hd,const char *fname);
int fread_image(FILE *fp,struct header *hd,int fr,const char *fname);
int free_hdrcon(struct header *hd);
// int free_header(struct header *hd);
int fwrite_header(FILE *fp,struct header *hd,const char *fname);
int fwrite_image(FILE *fp,struct header *hd,int fr,const char *fname);
int getparam(struct header *hda,...);
int h_add(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_clearhisto(struct hips_histo *histogram);
// int h_copy(struct header *hdi,struct header *hdo);
int h_diff(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_div(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_divscale(struct header *hdi,struct header *hdo,union pixelval *b);
int h_enlarge(struct header *hdi,struct header *hdo,int xf,int yf);
int h_entropycnt(struct header *hd,int *table,int pairflag);
int h_flipquad(struct header *hdi,struct header *hdo);
int h_fourtr(struct header *hd);
int h_histo(struct header *hd,struct hips_histo *histogram,int nzflag,int *count);
int h_histoeq(struct hips_histo *histogram,int count,unsigned char *map);
int h_invert(struct header *hdi,struct header *hdo);
int h_invfourtr(struct header *hd);
int h_linscale(struct header *hdi,struct header *hdo,float b,float c);
int h_median(struct header *hdi,struct header *hdo,int size);
int h_minmax(struct header *hd,union pixelval *minval,union pixelval *maxval,int nzflag);
int h_morphdil(struct header *hdi,struct header *hde,struct header *hdo,int centerr,int centerc,int gray);
int h_mul(struct header *hdi1,struct header *hdi2,struct header *hdo);
int h_mulscale(struct header *hdi,struct header *hdo,union pixelval *b);
int h_pixmap(struct header *hdi,struct header *hdo,unsigned char *map);
int h_reduce(struct header *hdi,struct header *hdo,int xf,int yf);
int h_softthresh(struct header *hdi,struct header *hdo,union pixelval *thresh);
int h_tob(struct header *hdi,struct header *hdo);
int h_toc(struct header *hdi,struct header *hdo);
int h_tod(struct header *hdi,struct header *hdo);
int h_todc(struct header *hdi,struct header *hdo);
int h_tof(struct header *hdi,struct header *hdo);
int h_toi(struct header *hdi,struct header *hdo);
int setparam(struct header *hda,...);
int update_header(struct header *hd,int argc,char **argv );
struct extpar *findparam(struct header *hd, const char *name);

// ------ externs ------

extern int hips_rtocplx;
extern int hips_cplxtor;

// ------ canny -----

void canny(int *magmax, int *hthresh, int *lthresh, int *image, int *xsize, int *ysize, short *shortim, int *windowsize, double *sigma, int *bordermode, double *hfrac, double *lfrac, int *pflag,
           short *gx, short *gy, short *mag, int *hist, int *histsize, unsigned char *nms, unsigned char *edgemap, float *gm, float *gmp,short *temp) ;
void cleanup(unsigned char *map, int xsize, int ysize) ;
void find_edges(unsigned char *map, short *mag, int xsize, int ysize, int maxmag, float hpixel_fraction, float lpixel_fraction, int *hgram, int hsize, int *actual_hthresh, int *actual_lthresh);
void follow_edges(unsigned char *edgemapptr, short *edgemagptr) ;
void clear_borders(unsigned char *charimage, int xsize, int ysize) ;
void gauss_filter(short *inimage, int inx, int iny, int direction, int boundary, int masksize, float sigma, short *grad, int *outx, int *outy, float *gmask, float *gprimemask, short *tempimage) ;
void copyimage(int *charimage, int ncols, int nrows, short *shortimage) ;
void thin(unsigned char *edges, int height, int width) ;
int h_canny(struct header *Isrc, struct header *Idst, double sigma, int mask_size,double lfrac,double hfrac,int dothin);

#endif
