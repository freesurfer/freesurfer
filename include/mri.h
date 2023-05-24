/**
 * @brief prototypes and structures for working with MRI volumes.
 *
 *  prototypes and structures for working with MRI volumes.
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

#ifndef MRI_H
#define MRI_H

#include <vector>
#include <array>
#include <string>

#include "faster_variants.h"

#include "minc.h"
#include "const.h"
#include "matrix.h"
#include "dmatrix.h"
#include "machine.h"
#include "colortab.h"
#include "affine.h"
#include "fnvhash.h"
#include "itkImage.h"


#define BUFTYPE  unsigned char

#define SAMPLE_NEAREST       0
#define SAMPLE_TRILINEAR     1
#define SAMPLE_SINC          2
#define SAMPLE_CUBIC         3 /*E*/
#define SAMPLE_WEIGHTED      4
#define SAMPLE_CUBIC_BSPLINE 5
#define SAMPLE_VOTE          6   // for downsampling aseg volumes

#define MRI_UCHAR   0
#define MRI_INT     1
#define MRI_LONG    2
#define MRI_FLOAT   3
#define MRI_SHORT   4
#define MRI_BITMAP  5
#define MRI_TENSOR  6
#define MRI_FLOAT_COMPLEX  7
#define MRI_DOUBLE_COMPLEX  8
#define MRI_RGB             9
#define MRI_USHRT          10

#define NEAREST_NEIGHBOR_FACE   1
#define NEAREST_NEIGHBOR_EDGE   2
#define NEAREST_NEIGHBOR_CORNER 3

#define MAX_CMDS 1000

#define SEQUENCE_MPRAGE      1
#define SEQUENCE_EPI         2

#define FRAME_TYPE_ORIGINAL             0
#define FRAME_TYPE_DIFFUSION_AUGMENTED  1

#define MRInvox(mri)  ((mri)->width * (mri)->height * (mri)->depth * (mri)->nframes)

#define BVEC_SPACE_UNKNOWN 0
#define BVEC_SPACE_SCANNER 1
#define BVEC_SPACE_VOXEL   2

#define MB_RADIAL 0
#define MB_TANGENTIAL 1

#define FS_COORDS_UNKNOWN      0
#define FS_COORDS_TKREG_RAS    1 // same as surfaceRAS
#define FS_COORDS_SCANNER_RAS  2 // same as surface "RealRAS"
#define FS_COORDS_VOXEL        3

// standard itk image type
typedef itk::Image<float, 3> ITKImageType;

typedef struct
{
  int     type ;           // code for what is stored in this frame
  float   TE ;             // echo time
  float   TR ;             // recovery time
  float   flip ;           // flip angle
  float   TI ;             // time-to-inversion
  float   TD ;             // delay time
  int     sequence_type ;  // see SEQUENCE* constants
  float   echo_spacing ;
  float   echo_train_len ; // length of the echo train
  float   read_dir[3] ;    // read-out direction in RAS coords
  float   pe_dir[3] ;      // phase-encode direction in RAS coords
  float   slice_dir[3] ;   // slice direction in RAS coords
  int     label ;          // index into CLUT
  char    name[STRLEN] ;   // human-readable description of frame contents
  int     dof ;            // for stat maps (e.g. # of subjects)
  MATRIX  *m_ras2vox ;
  float   thresh ;
  int     units ;          // e.g. UNITS_PPM,  UNITS_RAD_PER_SEC, ...

  // for Herr Dr. Prof. Dr. Dr. Witzel
  // directions: maybe best in both reference frames
  // or just 3 coordinates and a switch which frame it is ?
  double DX ;
  double DY ;
  double DZ ;

  double DR ;
  double DP ;
  double DS ;

// B-value
  double bvalue ;

// Mixing time
  double TM ;

// What kind of diffusion scan is this (can be an enum)
// stejskal-tanner,trse,steam etc....

  long diffusion_type ;

// Gradient values
  long D1_ramp ;
  long D1_flat ;
  double D1_amp ;

  long D2_ramp ;
  long D2_flat ;
  double D2_amp ;

  long D3_ramp ;
  long D3_flat ;
  double D3_amp ;
  long D4_ramp ;
  long D4_flat ;
  double D4_amp ;

} MRI_FRAME ;


typedef struct
{
  int  x ;
  int  y ;
  int  z ;
  int  dx ;
  int  dy ;
  int  dz ;
}
MRI_REGION ;


struct VOL_GEOM;
MATRIX *MRIxfmCRS2XYZ( const VOL_GEOM *mri, int base ); /* Native Vox2RAS Matrix (scanner and xfm too) */
MATRIX *MRIxfmCRS2XYZtkreg( const VOL_GEOM *mri );      // TkReg  Vox2RAS Matrix
MATRIX *VGras2tkreg(VOL_GEOM *vg, MATRIX *ras2tkreg);
MATRIX *VGtkreg2RAS(VOL_GEOM *vg, MATRIX *tkreg2ras);

struct VOL_GEOM
{
  // NOTE: VOL_GEOM is subclassed by MRI, so take that into account when changing
  int           valid;   /* whether this is a valid info or not (1 valid, 0 not valid) */
  int           width ;
  int           height ;
  int           depth ;
  float         xsize = 1; // Added =1 because that is what MRI does
  float         ysize = 1;
  float         zsize = 1;

  int ras_good_flag = 0;  // indicates whether the RAS coordinates are accurate
  float x_r = -1, x_a = 0, x_s =  0;
  float y_r =  0, y_a = 0, y_s = -1;
  float z_r =  0, z_a = 1, z_s =  0;
  float c_r =  0, c_a = 0, c_s =  0;
  char          fname[STRLEN];  // volume filename

  // i_to_r__, r_to_i__, register_mat were moved here from MRI in this commit:
  // https://github.com/freesurfer/freesurfer/commit/d7b80d733e96380543093e1bb0c11f0123f97132
  AffineMatrix *i_to_r__ = nullptr;  // cached i->r transform
  MATRIX *r_to_i__ = nullptr;        // cached r->i transform
  MATRIX *register_mat = nullptr;
  
  // The functions below compute these matrices on-the-fly
  // RAS = scanner RAS (sometimes known as "real" RAS in surface contexts)
  // TkregRAS = RAS used by tkregister; surface coords are by default in this TkregRAS space
  MATRIX *get_Vox2RAS(void){       return(MRIxfmCRS2XYZ(this,0));}
  MATRIX *get_RAS2Vox(void){       return(MatrixInverse(get_Vox2RAS(),NULL));}
  MATRIX *get_Vox2TkregRAS(void){  return(MRIxfmCRS2XYZtkreg(this));}
  MATRIX *get_TkregRAS2Vox(void){  return(MatrixInverse(get_Vox2TkregRAS(),NULL));}
  MATRIX *get_RAS2TkregRAS(void){  return(VGras2tkreg(this, NULL));}
  MATRIX *get_TkregRAS2RAS(void){  return(VGtkreg2RAS(this, NULL));}
};

typedef VOL_GEOM VG;

class MRI : public VOL_GEOM
{
public:

  class Shape
  {
  public:
    Shape() {};
    Shape(const std::vector<int>& shape);
    Shape(const std::vector<ssize_t>& shape) : Shape(std::vector<int>(shape.begin(), shape.end())) {};
    Shape(const std::initializer_list<int>& shape) : Shape(std::vector<int>(shape)) {}
    ssize_t width, height, depth, nframes, size;

    operator std::vector<ssize_t>() const { return {width, height, depth, nframes}; }

    friend bool operator == (const Shape &l, const Shape &r) { return (std::vector<ssize_t>(l) == std::vector<ssize_t>(r)); }
    friend bool operator != (const Shape &l, const Shape &r) { return !(l == r); }
  };

  MRI(const Shape volshape, int dtype, bool alloc = true);
  MRI(const std::string& filename);
  ~MRI();

  void initIndices();
  void initSlices();
  void write(const std::string& filename);
  FnvHash hash();

  // ITK image conversions
  ITKImageType::Pointer toITKImage(int frame = 0);
  void loadITKImage(ITKImageType::Pointer image, int frame = 0);

  // ---- image geometry ----
  //int width;        // number of columns // Now inherited from VOL_GEOM
  //int height;       // number of rows    // Now inherited from VOL_GEOM
  //int depth;        // number of slices  // Now inherited from VOL_GEOM
  int nframes;      // number of frames
  Shape shape;      // volume shape
  int imnr0;        // starting image number
  int imnr1;        // ending image number
  float xstart;     // starting x (in xsize units)
  float ystart;     // starting y (in ysize units)
  float zstart;     // starting z (in zsize units)
  float xend;       // ending x (in xsize units)
  float yend;       // ending y (in ysize units)
  float zend;       // ending z (in zsize units)
  //float xsize = 1;  // size of a voxel in the x direction // Now inherited from VOL_GEOM
  //float ysize = 1;  // size of a voxel in the y direction // Now inherited from VOL_GEOM
  //float zsize = 1;  // size of a voxel in the z direction // Now inherited from VOL_GEOM
  float thick = 1;
  int scale = 1;
  float ps = 1;
  float fov;

  // ---- indices to handle boundary conditions ----
  int *xi = nullptr;
  int *yi = nullptr;
  int *zi = nullptr;

  // ---- RAS distances ---- // Now inherited from VOL_GEOM
  //float x_r = -1, x_a = 0, x_s =  0;
  //float y_r =  0, y_a = 0, y_s = -1;
  //float z_r =  0, z_a = 1, z_s =  0;
  //float c_r =  0, c_a = 0, c_s =  0;
  //int ras_good_flag = 0;  // indicates whether the RAS coordinates are accurate
  
  // ---- transforms ----
  char transform_fname[STRLEN];
  General_transform transform;
  Transform *linear_transform = nullptr;
  Transform *inverse_linear_transform = nullptr;
  int free_transform = 0;
  MATRIX *AutoAlign = nullptr;       // for Andre

  // ---- volume metadata ----
  double outside_val = 0;
  double mean;
  int brightness = 1;
  int yinvert = 1;              // for converting between MNC and coronal slices
  int dof = 1;
  MRI_FRAME *frames = nullptr;
  COLOR_TABLE *ct = nullptr;
  MRI_REGION roi;

  // ---- scan parameters ----
  float tr = 0;                 // time to recovery
  float te = 0;                 // time to echo
  float ti = 0;                 // time to inversion
  double flip_angle = 0;        // flip angle in radians
  float FieldStrength = 0;      // field strength
  char *pedir = nullptr;        // phase enc direction: ROW, COL, etc
  MATRIX *origRas2Vox = nullptr ;         // to get to original voxel grid from ras
  float location = 0;           // NOT USED

  // ---- DTI ----
  int bvec_space = 0;           // 0: unknown, 1: scanner, 2: voxel
  MATRIX *bvals = nullptr;
  MATRIX *bvecs = nullptr;

  // ---- file metadata ----
  //char fname[STRLEN];           // filename // Now inherited from VOL_GEOM
  char fnamePostFixes[STRLEN];    // used in MRIwrite(), append to output file name
  int  len_fnamePostFixes;
  char fname_format[STRLEN];    // file extension
  char subject_name[STRLEN];    // fs subject name
  char path_to_t1[STRLEN];      // NOT USED
  char gdf_image_stem[STRLEN];
  char *cmdlines[MAX_CMDS];     // command line provenance
  int ncmds = 0;                // number of commands run previously
  void *tag_data = nullptr;     // saved tag data
  int tag_data_size = 0;        // size of tag data

  // ---- image buffer ----
  int type;                     // image data type
  int ptype = 2;                // NOT USED
  size_t bytes_per_vox = 0;     // bytes per voxel
  size_t bytes_total = 0;       // total bytes in buffer
  size_t vox_per_row = 0;       // number of voxels per volume row
  size_t vox_per_slice = 0;     // number of voxels per volume slice
  size_t vox_per_vol = 0;       // number of voxels per volume frame
  size_t vox_total = 0;         // total number of voxels in the volume
  int ischunked;                // indicates whether the buffer is chunked (contiguous)
  bool owndata = true;          // indicates ownership of the chunked buffer data
  BUFTYPE ***slices = nullptr;  // fallback non-contiguous storage for 3D-indexed image data
  void *chunk = nullptr;        // default contiguous storage for image data
};


typedef struct
{
  int type; // MB_RADIAL or MB_TANGENTIAL
  double offset,slope;
  int c0,r0; // center of motion in full volume space
  int cR,rR; // col and row of first voxel of region in full volume space
  double DeltaD; // sample spacing along radius
  double cutoff; // number of stddevs to cut off kernel
  int Interp; // SAMPLE_NEAREST or SAMPLE_TRILINEAR
  MRI *d0;    // Image of distance to voxel from center
  MRI *theta; // Image of angle of voxel
  MRI *fwhm;  // Image of FWHM at voxel
  MRI *dmin;  // Image of start distance (cut off)
  MRI *nd;    // Image of number of samples
} MOTIONBLUR2D, MB2D;


MRI *MRImotionBlur2D(MRI *src, MB2D *mb, MRI *out);
int MB2Dfree(MB2D **pmb);
MB2D *MB2Dcopy(MB2D *src, int CopyMRI, MB2D *copy);
MRI *MB2Dgrid(MRI *mbtemplate, int skip, MRI *outvol);

MATRIX *MRIcopyFramesToMatrixRows(MRI *mri, MATRIX *m_dst, int start_frame, int nframes, int dst_row) ;

MATRIX *vg_i_to_r(const VOL_GEOM *vg);
MATRIX *vg_r_to_i(const VOL_GEOM *vg);
#define vg_getRasToVoxelXform vg_r_to_i
#define vg_getVoxelToRasXform vg_i_to_r
MATRIX *TkrVox2RASfromVolGeom(const VOL_GEOM *vg);
MATRIX *TkrRAS2VoxfromVolGeom(const VOL_GEOM *vg);



MATRIX *MRIxfmCRS2XYZfsl(VOL_GEOM *mri);        // FSL/FLIRT  Vox2RAS Matrix

int MRIsetVox2RASFromMatrix(VOL_GEOM *mri, MATRIX *m_vox2ras);
int MRIsetVox2RASFromMatrixUnitTest(MRI *mri);
MATRIX *MRItkRegMtxFromVox2Vox(VOL_GEOM *ref, VOL_GEOM *mov, MATRIX *vox2vox);//ras2ras from vox2vox
MATRIX *MRItkReg2Native(VOL_GEOM *ref, VOL_GEOM *mov, MATRIX *R); /* tkreg2native (scanner and xfm too) */
MATRIX *MRItkRegMtx(VOL_GEOM *ref, VOL_GEOM *mov, MATRIX *D); /* native2tkreg (scanner and xfm too) */
MATRIX *MRIvoxToVoxFromTkRegMtx(VOL_GEOM *mov, VOL_GEOM *targ, MATRIX *tkR);
MATRIX *MRIfsl2TkReg(VOL_GEOM *ref, VOL_GEOM  *mov, MATRIX *FSLRegMat);
MATRIX *MRItkreg2FSL(VOL_GEOM *ref, VOL_GEOM *mov, MATRIX *tkRegMat);
MATRIX *MtxCRS1toCRS0(MATRIX *Q);
int MRIp0ToCRAS(VOL_GEOM *mri, double r0, double a0, double s0);
MATRIX *MRIfixTkReg(VOL_GEOM *mov, MATRIX *R);
MATRIX *MRImatrixOfDirectionCosines(VOL_GEOM *mri, MATRIX *Mdc);
MATRIX *MRImatrixOfVoxelSizes(VOL_GEOM *mri, MATRIX *D);
MATRIX *MRImatrixOfTranslations(VOL_GEOM *mri, MATRIX *P0);

int MRIhfs2Sphinx(VOL_GEOM *mri);

float  MRIgetVoxDx(MRI *mri, int c, int r, int s, int f);
float  MRIgetVoxDy(MRI *mri, int c, int r, int s, int f);
float  MRIgetVoxDz(MRI *mri, int c, int r, int s, int f);

float  MRIgetVoxVal( const MRI *mri, int c, int r, int s, int f);
float  MRIgetVoxVal2( const MRI *mri, int c, int r, int s, int f);
int    MRIsetVoxVal(MRI *mri, int c, int r, int s, int f, float voxval);
int    MRIsetVoxVal2(MRI *mri, int c, int r, int s, int f, float voxval);
void   MRIdbl2ptr(double v, void *pmric, int mritype);
double MRIptr2dbl(void *pmric, int mritype);

size_t MRIsizeof(int mritype);

const char * MRIprecisionString(int PrecisionCode);
int MRIprecisionCode(const char *PrecisionString);

MRI *MRImakeMosaic(MRI **mri, int nimages, int rectify) ;

int MRIareNonzeroInNbhd(MRI *mri, int wsize, int x, int y, int z) ;
float  MRIfindNearestNonzero(MRI *mri,
                             int wsize,
                             double x0, double y0, double z0,
                             float max_dist) ;
float  MRIfindNearestNonzeroLocation(MRI *mri, int wsize,
                                     double xr, double yr, double zr,
                                     int *pxv, int *pyv, int *pzv) ;
/* single pixel filtering */
float MRIvoxelMedian(MRI *mri, int x0, int y0, int z0, int wsize) ;
float MRIvoxelMean( const MRI *mri, int x, int y, int z, int wsize, int frame) ;
float MRIvoxelMin(MRI *mri, int x0, int y0, int z0, int wsize) ;
float MRIvoxelMax(MRI *mri, int x0, int y0, int z0, int wsize) ;
float MRIvoxelStd(MRI *mri, int x, int y, int z, float mean, int wsize) ;
float MRIvoxelZscore(MRI *mri, int x, int y, int z, int wsize) ;
float MRIvoxelDx(MRI *mri, int x, int y, int z) ;
float MRIvoxelDy(MRI *mri, int x, int y, int z) ;
float MRIvoxelDz(MRI *mri, int x, int y, int z) ;
float MRIvoxelGradient(MRI *mri, int x, int y, int z, float *pdx, float *pdy,
                       float *pdz) ;
float MRIvoxelDirection(MRI *mri, int x, int y, int z, int wsize) ;
MRI *MRI2ndDirectionalDerivative(MRI *mri_src, MRI *mri_deriv, float nx, float ny, float nz) ;
float MRIvoxelGradientDir2ndDerivative(MRI *mri, int x0, int y0, int z0,
                                       int wsize) ;
MRI  * MRIgradientDir2ndDerivative(MRI *mri_src, MRI *mri_dst, int wsize) ;


double MRImeanFrameThresh(MRI *mri, int frame, float thresh);



/* use these constants for MRIreorder */
#define XDIM  1
#define YDIM  2
#define ZDIM  3
/* ch ov */
/*
  MRI  *MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim);
*/

/* I/O functions */
/* ch ov */
/*
  int    MRIwrite(MRI *mri, char *fpref) ;
*/
int    MRIappend(MRI *mri,const char *fpref) ;
int    MRIwriteInfo(MRI *mri,const char *fpref) ;
/* ch ov */
/*
  MRI   *MRIread(char *fpref) ;
  MRI   *MRIreadInfo(char *fpref) ;
*/

/* memory allocation routines */
int   MRIfree(MRI **pmri) ;
MRI   *MRIalloc(int width, int height, int depth, int type) ;
MRI   *MRIallocSequence(int width, int height,int depth,int type,int nframes);
MRI   *MRIallocHeader(int width, int height, int depth, int type, int nframes) ;
int   MRIsetResolution(MRI *mri, float xres, float yres, float zres) ;
int   MRIsetTransform(MRI *mri,   General_transform *transform) ;


/* correlation routines */
MRI   *MRIxcorr(MRI *mri_ref, MRI *mri_in, MRI *mri_dst) ;
MRI   *MRIxcorrWindow(MRI *mri_ref, MRI *mri_in,MRI *mri_dst,int window_size) ;
MRI   *MRInxcorr(MRI *mri_ref, MRI *mri_in, MRI *mri_dst) ;
MRI   *MRInxcorrWindow(MRI *mri_ref,MRI *mri_in,MRI *mri_dst,int window_size) ;
long  MRIcorrelate(MRI *mri_ref, MRI *mri_in, int xoff, int yoff, int zoff) ;


int   MRIpeak(MRI *mri, int *px, int *py, int *pz) ;
int   MRIcompareHeaders(MRI const *mri1, MRI const *mri2) ;
MRI   *MRIcopyHeader( const MRI *mri_src, MRI *mri_dst) ;
int   MRIcopyPulseParameters(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcopy(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIreslice(MRI *mri_src, MRI *mri_dst, int slice_direction) ;
int   MRIboundingBox(MRI *mri, int thresh, MRI_REGION *region) ;
int   MRIlabelBoundingBox(MRI *mri, int label, MRI_REGION *region) ;
int   MRIfindApproximateSkullBoundingBox(MRI *mri, int thresh,
    MRI_REGION *region) ;
int   MRIboundingBoxNbhd(MRI *mri, int thresh, int wsize,MRI_REGION *region) ;
MRI *MRIsetBoundingBox(MRI *mri_template,
                       MRI_REGION *region,
                       double InVal,
                       double OutVal);

/* coordinate transforms */
MRI   *MRItranslate(MRI *mri_src, MRI *mri_dst,
                    double dx, double dy, double dz) ;
MRI   *MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;
MRI   *MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz) ;
MRI   *MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB) ;
MRI   *MRIinverseLinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA) ;
MRI   *MRIlinearTransformInterp(MRI *mri_src, MRI *mri_dst, MATRIX *mA,
                                int InterpMethod);
MRI   *MRIlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA) ;
MRI   *MRIapplyRASlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA) ;
MRI   *MRIapplyRASinverseLinearTransform(MRI *mri_src, MRI *mri_dst,
    MATRIX *mA) ;
MRI   *MRIapplyRASlinearTransformInterp(MRI *mri_src, MRI *mri_dst,
                                        MATRIX *mA, int interpMethod) ;
MRI   *MRIapplyRASinverseLinearTransformInterp(MRI *mri_src, MRI *mri_dst,
    MATRIX *mA, int interpMethod) ;

int MRIinterpCode(const char *InterpString);
const char * MRIinterpString(int InterpCode);
MRI   *MRIinterpolate(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIconfThresh(MRI *mri_src, MRI *mri_probs, MRI *mri_classes,
                     MRI *mri_dst,float thresh, int min_target,int max_target);

/* debugging */
int   MRIdump(MRI *mri, FILE *fp) ;
int   MRIdumpBuffer(MRI *mri, FILE *fp) ;

/* arithmetic operations */
double MRIrmsDifferenceNonzero(MRI *mri1, MRI *mri2) ;
MRI   *MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabsdiff(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIaddScalar(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI   *MRIaverage(MRI *mri_src, int dof, MRI *mri_dst) ;
MRI   *MRIaddToFrame(MRI *mri_src, MRI *mri_to_add, MRI *mri_dst, int src_frame_no, int dst_frame_no) ;
MRI   *MRIdivide(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
  MRI   *MRIdivideFrames(MRI *mri1, MRI *mri2, int frame1, int frame2, MRI *mri_dst) ;
MRI   *MRImultiply(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIscaleAndMultiply(MRI *mri1, float scale, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabs(MRI *mri, MRI *mri_dst) ;
MRI   *MRIneg(MRI *mri_src, MRI *mri_dst);
MRI   *MRIpos(MRI *mri_src, MRI *mri_dst);
MRI   *MRIlinearScale(MRI *mri_src,
                      MRI *mri_dst,
                      float scale,
                      float offset,
                      int only_nonzer) ;
MRI   *MRIscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;
MRI   *MRIscalarMulFrame(MRI *mri_src, MRI *mri_dst, float scalar, int frame) ;
void  MRIrms(MRI *in, MRI *out);

/* filtering */
int   MRIcpolvAllQuadrantsFilled(MRI *mri,
                                 int x, int y, int z,int vertex,
                                 int wsize) ;
MRI   *MRIremoveIslands(MRI *mri_src, MRI*mri_dst, int wsize, int thresh) ;
MRI   *MRIresegmentThinWMStrands(MRI *mri_src, MRI *mri_dst, int thickness);
MRI   *MRIthickenThinWMStrands(MRI *mri_T1,
                               MRI *mri_src, MRI *mri_dst,
                               int thickness, int nsegments, float wm_hi) ;
MRI   *MRIfindThinWMStrands(MRI *mri_src, MRI *mri_dst, int wsize);
MRI   *MRIcentralPlaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst,
    int wsize);
MRI   *MRIplaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize) ;
int   MRIcpolvMaxWhiteAtVoxel(MRI *mri, int x, int y, int z, int wsize) ;
MRI   *MRIpolvZscore(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize) ;
MRI   *MRIpolvNormalCurvature(MRI *mri_src, MRI *mri_dst, MRI *mri_polv,
                              int wsize) ;
MRI   *MRIpolvMean(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize) ;
MRI   *MRIpolvMedian(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize) ;
MRI   *MRIpolvOrder(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize,
                    int thresh) ;
MRI   *MRIpolvCount(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize,
                    int low_lim, int hi_lim) ;
MRI   *MRIorderThreshold(MRI *mri_src, MRI *mri_dst, MRI *mri_order, int num) ;

MRI   *MRIpolvMeanRegion(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize,
                         MRI_REGION *region);
MRI   *MRIpolvMedianRegion(MRI *mri_src, MRI *mri_dst,MRI *mri_polv,int wsize,
                           MRI_REGION *region);

MRI   *MRIdivergence(MRI *mri_src, MRI *mri_divergence) ;
MRI   *MRIlaplacian(MRI *mri_src, MRI *mri_laplacian);
MRI   *MRIsobelFrame(MRI *mri_src, MRI *mri_grad, MRI *mri_mag, int frame) ;
MRI   *MRIsobel(MRI *mri_src, MRI *mri_grad, MRI *mri_mag);
MRI   *MRIxSobel(MRI *mri_src, MRI *mri_x, int frame) ;
MRI   *MRIxSobelForAllTypes(MRI *mri_src, MRI *mri_x, int frame) ;
MRI   *MRIySobel(MRI *mri_src, MRI *mri_y, int frame) ;
MRI   *MRIySobelForAllTypes(MRI *mri_src, MRI *mri_y, int frame) ;
MRI   *MRIzSobel(MRI *mri_src, MRI *mri_z, int frame) ;
MRI   *MRIzSobelForAllTypes(MRI *mri_src, MRI *mri_z, int frame) ;
MRI   *MRIsobelRegion(MRI *mri_src, MRI *mri_grad, int domag,
                      MRI_REGION *region);
MRI   *MRIxSobelRegion(MRI *mri_src, MRI *mri_x, int frame,MRI_REGION *region);
MRI   *MRIySobelRegion(MRI *mri_src, MRI *mri_y, int frame,MRI_REGION *region);
MRI   *MRIzSobelRegion(MRI *mri_src, MRI *mri_z, int frame,MRI_REGION *region);

MRI   *MRIreduce(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIreduce2D(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIreduceSlice(MRI *mri_src, MRI *mri_dst,
                      float *k, int len, int axis) ;
MRI   *MRIreduceByte(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIconvolve1dFloat(MRI *mri_src, MRI *mri_dst,
                        float *k, int len, int axis,
                        int src_frame, int dst_frame);
MRI   *MRIconvolve1dShort(MRI *mri_src, MRI *mri_dst, float *k, int len,
                          int axis, int src_frame, int dst_frame) ;
MRI   *MRIconvolve1dInt(MRI *mri_src, MRI *mri_dst, float *k, int len,
                          int axis, int src_frame, int dst_frame) ;
MRI   *MRIconvolve1dByte(MRI *mri_src, MRI *mri_dst, float *k, int len,
                         int axis, int src_frame, int dst_frame) ;
MRI   *MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *kernel,
                     int len, int axis, int src_frame, int dst_frame) ;
MRI   *MRIreduce1d(MRI *mri_src, MRI *mri_dst,float *kernel,int len,int axis);
MRI   *MRIreduce1dByte(MRI *mri_src, MRI *mri_dst,float *kernel,int len,
                       int axis);
double MRIrmsDiff(MRI *mri1, MRI *mri2) ;
MRI   *MRIdiffuse(MRI *mri_src, MRI *mri_dst, double k,
                  int niter, int which, double slope) ;
MRI   *MRIdiffuseCurvature(MRI *mri_src, MRI *mri_dst,
                           double A,int niter, double slope) ;
MRI   *MRIdiffusePerona(MRI *mri_src, MRI *mri_dst,
                        double k, int niter,double slope);
MRI   *MRIdirectionMap(MRI *mri_grad, MRI *mri_direction, int wsize);
MRI   *MRIdirectionMapUchar(MRI *mri_grad, MRI *mri_direction, int wsize);
void  MRIcalcCRASforSampledVolume(MRI *src, MRI *sampled,
                                  double *pr, double *pa, double *ps);
void  MRIcalcCRASforExtractedVolume(MRI *src, MRI *dst,
                                    int x0, int y0, int z0,
                                    int x1, int y1, int z1,
                                    double *pr, double *pa, double *ps);
// 0 is the src extract position start
// 1 is the dst extracted region start
MRI   *MRIsrcTransformedCentered(MRI *src, MRI *dst,
                                 MATRIX *stod_voxtovox, int interp_method);
MRI   *MRITransformedCenteredMatrix(MRI *src, MRI *orig_dst, MATRIX *m_L) ;

/* offset stuff */
MRI   *MRIoffsetDirection(MRI *mri_grad, int wsize, MRI *mri_direction,
                          MRI *mri_dir);
MRI   *MRIoffsetMagnitude(MRI *mri_src, MRI *mri_dst, int maxsteps) ;
MRI   *MRIapplyOffset(MRI *mri_src, MRI *mri_dst, MRI *mri_offset) ;


 /* it just copies the header info, not image data */
MRI   *MRIclone( const MRI *mri_src, MRI *mri_dst );
MRI   *MRIcloneDifferentType(MRI *mri_src, int type) ;
MRI   *MRIcloneRoi(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcloneBySpace(MRI *mri_src, int type, int nframes);
MRI   *MRIthreshold(MRI *mri_src, MRI *mri_dst, float threshold) ;
MRI   *MRIthresholdAllFrames(MRI *mri_src, MRI *mri_dst, float threshold) ;
MRI   *MRIupperthresholdAllFrames(MRI *mri_src, MRI *mri_dst, float threshold) ;
MRI   *MRIthresholdFrame(MRI *mri_src, MRI *mri_dst, float threshold, int frame) ;
MRI   *MRIupperthresholdFrame(MRI *mri_src, MRI *mri_dst, float threshold, int frame) ;
MRI   *MRIinvert(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIinvertContrast(MRI *mri_src, MRI *mri_dst, float threshold) ;
MRI   *MRIbinarizeNoThreshold(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIbinarize(MRI *mri_src, MRI *mri_dst, float threshold,
                   float low_val, float hi_val) ;
MRI   *MRIthresholdRangeInto(MRI *mri_src, MRI *mri_dst,
                             BUFTYPE low_val, BUFTYPE hi_val) ;
int   MRIprincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues,
                             double *means, BUFTYPE theshold) ;
int   MRIcenterOfMass(MRI *mri,double *means, BUFTYPE threshold) ;
int   MRIbinaryPrincipleComponents(MRI *mri, MATRIX *mEvectors,
                                   float *evalues,
                                   double *means, BUFTYPE theshold) ;
int   MRIprincipleComponentsRange(MRI *mri,
                                  MATRIX *mEvectors,
                                  float *evalues,
                                  double *means,
                                  float low_thresh,
                                  float hi_thresh) ;
int   MRIclear(MRI *mri_src) ;

/* these routines use trilinear interpolation */
MRI   *MRIrotateX_I(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY_I(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ_I(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate_I(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;

/* extraction routines */
MRI   *MRIextract(MRI *mri_src, MRI *mri_dst,
                  int x0, int y0, int z0,
                  int dx, int dy, int dz) ;
MRI   *MRIextractInto(MRI *mri_src, MRI *mri_dst,
                      int x0, int y0, int z0,
                      int dx, int dy, int dz,
                      int x1, int y1, int z1) ;
MRI   *MRIextractIntoRegion(MRI *mri_src, MRI *mri_dst,
                            int x0, int y0, int z0,
                            MRI_REGION *region) ;

MRI *MRIextractRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *region) ;
MRI *MRIinsertRegion(MRI *regionvol, MRI_REGION *region, MRI *temp, MRI *out);

MRI   *MRIextractPolvPlane(MRI *mri_src, MRI *mri_dst, MRI *mri_polv,
                           int wsize, int x, int y, int z);
MRI   *MRIextractCpolv(MRI *mri_src, MRI *mri_dst, MRI *mri_polv,
                       int wsize, int x, int y, int z);
MRI   *MRIextractCpolvCoords(MRI *mri_src, int *px, int *py, int *pz,
                             MRI *mri_polv, int x, int y,int z,int wsize);
MRI   *MRIextractValues(MRI *mri_src, MRI *mri_dst, float min_val,
                        float max_val) ;
MRI   *MRIwmfilter(MRI *mri_src, MRI *mri_cpolv, MRI *mri_dst,
                   float nslope, float pslope) ;
MRI   *MRIorder(MRI *mri_src, MRI *mri_dst, int wsize, float pct) ;
#if 1
MRI   *MRIremoveHoles(MRI *mri_src, MRI*mri_dst, int wsize, float pct,
                      int use_all) ;
#else
MRI   *MRIremoveHoles(MRI *mri_src, MRI*mri_dst, int wsize, float pct) ;
#endif

MRI   *MRInormalizeFrameVectorLength(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcomputeFrameVectorLength(MRI *mri_src, MRI *mri_dst);
MRI   *MRIcomputeFrameVectorL1Length(MRI *mri_src, MRI *mri_dst);

/* morphology */
MRI   *MRImorph(MRI *mri_src, MRI *mri_dst, int which) ;
MRI   *MRIerode(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIerodeNN(MRI *in, MRI *out, int NNDef);
MRI   *MRIerodeLabels(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIerodeThresh(MRI *mri_src, MRI *mri_intensity, double thresh,
                      MRI *mri_dst) ;
MRI * MRIdilate6Thresh(MRI *mri_src, MRI *mri_intensity, double thresh,
                       MRI *mri_dst) ;
MRI   *MRIdilateThresh(MRI *mri_src, MRI *mri_intensity, double thresh,
                      MRI *mri_dst) ;
MRI   *MRIerodeZero(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIerode2D(MRI *mri_src, MRI *mri_dst);
MRI   *MRIerodeRegion(MRI *mri_src, MRI *mri_dst,int wsize,MRI_REGION *region);
MRI   *MRIerodeSegmentation(MRI *seg, MRI *out, int nErodes, int nDiffThresh);
MRI *MRIdilateSegmentation(MRI *seg, MRI *out, int nDils, MRI *mask,
			   int maskframe, double maskthresh, int *pnchanges);
MRI   *MRIdilate(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIdilateUchar(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIopen(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIopenN(MRI *mri_src, MRI *mri_dst, int order) ;
MRI   *MRIclose(MRI *mri_src, MRI *mri_dst) ;
/* the following use 4 (or 6 in 3-D) connectivity */
MRI   *MRIerode6(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIdilate6(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIopen6(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIclose6(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIunion(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIintersect(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIcomplement(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIxor(MRI *mri1, MRI *mri2, MRI *mri_dst, int t1, int t2) ;
MRI   *MRIand(MRI *mri1, MRI *mri2, MRI *mri_dst, int thresh) ;
MRI   *MRIandVal(MRI *mri1, MRI *mri2, MRI *mri_dst, int thresh) ;
MRI   *MRIor(MRI *mri1, MRI *mri2, MRI *mri_dst, int thresh) ;
MRI   *MRIorVal(MRI *mri1, MRI *mri2, MRI *mri_dst, int thresh) ;
MRI   *MRInot(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcomputeResidual(MRI *mri1, MRI *mri2, MRI *mri_dst, int t1, int t2) ;

/* filtering operations */
MRI   *MRImodeFilter(MRI *mri_src, MRI *mri_dst, int niter) ;
MRI   *MRImodeFilterWithControlPoints(MRI *mri_src,
                                      MRI *mri_ctrl,
                                      MRI *mri_dst,
                                      int niter) ;
MRI *MRIcombineDistanceTransforms(MRI *mri_src1, MRI *mri_src2, MRI *mri_dst);
MRI   *MRIthreshModeFilter(MRI *mri_src, MRI *mri_dst, int niter,float thresh);
MRI   *MRImin(MRI *mri1, MRI *mri2, MRI *mri_min);
MRI   *MRIminmax(MRI *mri_src, MRI *mri_dst, MRI *mri_dir, int wsize) ;
MRI   *MRIgaussian1d(float sigma, int max_len) ;
MRI   *MRIconvolveGaussian(MRI *mri_src, MRI *mri_dst, MRI *mri_gaussian) ;
MRI   *MRIgaussianSmooth(MRI *src, double std, int norm, MRI *targ);
MRI   *MRImaskedGaussianSmooth(MRI *src, MRI *binmask, double std, MRI *targ);
MRI   *MRIconvolveGaussianMeanAndStdByte(MRI *mri_src, MRI *mri_dst,
    MRI *mri_gaussian) ;
MRI *MRIgaussianSmoothNI(MRI *src, double cstd, double rstd, double sstd,
			 MRI *targ);

/* frequency filtering*/
MRI* MRI_fft(MRI *mri_src, MRI* dst);
MRI *MRI_ifft(MRI *src, MRI *dst, int w, int h, int d);
MRI *MRI_fft_gaussian(MRI *src, MRI *dst, float std, int norm);
MRI *MRI_fft_lowpass(MRI *src, MRI *dst, int percent);
MRI *MRI_fft_highpass(MRI *src, MRI *dst, int percent);

MRI  *MRIcomputeMeanAndStandardDeviation(MRI *mri_src, MRI *mri_dst, int dof) ;
MRI *MRIscaleMeanIntensities(MRI *mri_src, MRI *mri_ref, MRI *mri_dst) ;
MRI   *MRIscaleIntensities(MRI *mri_src, MRI *mri_dst, float scale, float offset) ;
MRI   *MRImedian(MRI *mri_src, MRI *mri_dst, int wsize, MRI_REGION *box) ;
MRI   *MRImean(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRIminAbs(MRI *mri_src1, MRI *mri_src2, MRI *mri_dst) ;
MRI   *MRImeanInMask(MRI *mri_src, MRI *mri_dst, MRI *mri_mask, int wsize) ;
MRI   *MRIstdInMask(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_mask, int wsize) ;
double MRImeanInLabel(MRI *mri_src, MRI *mri_labeled, int label) ;
MATRIX *MRIcovarianceInLabelMultispectral(MRI *mri_src, MRI *mri_labeled, VECTOR *v_means, int label) ;
VECTOR *MRImeanInLabelMultispectral(MRI *mri_src, MRI *mri_labeled, int label) ;
double MRIstdInLabel(MRI *mri_src, MRI *mri_labeled, MRI *mri_mean, int label) ;
double MRImeanInLabelInRegion(MRI *mri_src, MRI *mri_labeled,
                              int label, int x0, int y0, int z0, int whalf, double *sigma);
MRI   *MRImeanByte(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRIstd(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize) ;
MRI   *MRIstdNonzero(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize) ;
MRI   *MRIzScore(MRI *mri_src, MRI *mri_dst, MRI *mri_mean, MRI *mri_std) ;

MRI   *MRIdirectionMapRegion(MRI *mri_grad, MRI *mri_direction, int wsize,
                             MRI_REGION *region);
MRI   *MRImeanRegion(MRI *mri_src, MRI *mri_dst, int wsize,MRI_REGION *region);
MRI   *MRIstdRegion(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize,
                    MRI_REGION *region) ;
MRI   *MRIzScoreRegion(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, MRI *mri_std,
                       MRI_REGION *region) ;

int   MRIcheckSize(MRI *mri_src, MRI *mri_check, int width, int height,
                   int depth) ;
/* ch ov */
/*
  MRI   *MRIreadRaw(FILE *fp, int width, int height, int depth, int type) ;
*/
int   MRIreInitCache(MRI *mri); /* when header is modified,
                                   you must call this function
                                   to update cached info */
int   MRIvoxelToWorld(MRI *mri, double xv, double yv, double zv,
                      double *xw, double *yw, double *zw) ;
int   MRIworldToVoxel(MRI *mri, double xw, double yw, double zw,
                      double *pxv, double *pyv, double *pzv) ;
int   MRIworldToVoxelIndex(MRI *mri, double xw, double yw, double zw,
                           int *pxv, int *pyv, int *pzv) ;
MRI *MRItoTalairach(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfromTalairach(MRI *mri_src, MRI *mri_dst) ;
int   MRIworldToTalairachVoxel(MRI *mri, double xw, double yw, double zw,
                               double *pxv, double *pyv, double *pzv) ;
int   MRIvoxelToTalairachVoxel(MRI *mri, double xv, double yv, double zv,
                               double *pxt, double *pyt, double *pzt) ;
int   MRIvoxelToVoxel(MRI *mri_src, MRI *mri_dst,
                      double xv, double yv, double zv,
                      double *pxt, double *pyt, double *pzt) ;
int   MRItalairachVoxelToVoxel(MRI *mri, double xt, double yt, double zt,
                               double *pxv, double *pyv, double *pzv) ;
int   MRItalairachVoxelToWorld(MRI *mri, double xt, double yt, double zt,
                               double *pxw, double *pyw, double *pzw) ;
int   MRIvoxelToTalairach(MRI *mri, double xv, double yv, double zv,
                          double *pxt, double *pyt, double *pzt) ;
int   MRItalairachToVoxel(MRI *mri, double xt, double yt, double zt,
                          double *pxv, double *pyv, double *pzv) ;
int MRIworldToTalairach(MRI *mri, double xw, double yw, double zw, double *pxt, double *pyt, double *pzt);

int   MRItransformRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *src_region,
                         MRI_REGION *dst_region) ;
MRI  *MRIextractArbitraryPlane(MRI *mri_src, MRI *mri_dst,
                               double e1_x, double e1_y, double e1_z,
                               double e2_x, double e2_y, double e2_z,
                               int x, int y, int z, int wsize);
MRI   *MRIextractTalairachPlane(MRI *mri_src, MRI *mri_dst, int orientation,
                                int x, int y, int z, int size) ;
int   MRIeraseTalairachPlane(MRI *mri, MRI *mri_mask, int orientation,
                             int x, int y, int z,int size,int fill_val);
int   MRIeraseTalairachPlaneNew(MRI *mri, MRI *mri_mask, int orientation,
                                int x, int y, int z,int size,int fill_val);

MRI   *MRIextractPlane(MRI *mri_src, MRI *mri_dst, int orientation,int where);
MRI   *MRIfillPlane(MRI *mri_mask, MRI *mri_dst,
                    int orientation, int where, int fillval);
int   MRIerasePlane(MRI *mri, float x0, float y0, float z0,
                    float dx, float dy, float dz, int fill_val);

int   MRIeraseBorders(MRI *mri, int width) ;
int   MRIindexNotInVolume( const MRI *mri,
			   const double col, const double row, const double slice );
int   MRIsampleVolume( const MRI *mri,double x, double y, double z, double *pval );
DMATRIX *MRIgradTrilinInterp(const MRI *mri, double x, double y, double z, DMATRIX *grad);

double *MRItrilinKernel(MRI *mri,
                        double c,
                        double r,
                        double s,
                        double *kernel);
int   MRIinterpolateIntoVolume(MRI *mri,
                               double x, double y, double z,
                               double val) ;
int   MRIinterpolateIntoVolumeFrame(MRI *mri,
                               double x, double y, double z,
                                    int frame, double val) ;
int   MRIsampleVolumeSlice(MRI *mri,
                           double x, double y, double z,
                           double *pval,
                           int slice_direction) ;
int   MRIsampleSeqVolume( const MRI *mri,
                          double x, double y, double z,
                          float *valvect,
                          int firstframe, int lastframe);
int   MRIsampleSeqVolumeType(MRI *mri,
			     double x, double y, double z,
			     float *valvect,
			     int firstframe, int lastframe, int type);
int   MRIsampleVolumeType( const MRI *mri,
                           double x, double y, double z,
                           double *pval,
                           int type );
int   MRIsampleLabeledVolume(MRI *mri,
                             double x, double y, double z,
                             double *pval,
                             unsigned char ucharLabel);
int   MRIsampleVolumeFrame( const MRI *mri,
			    double x, double y, double z,
			    const int frame,
			    double *pval);
#ifdef FASTER_MRI_EM_REGISTER
int   MRIsampleVolumeFrame_xyzInt_nRange_floats(const MRI *mri,
                            int x, int y, int z, 
			    const int frameBegin,
			    const int frameEnd,		// [frameBegin] .. [frameEnd-1] done
			    float *valForEachFrame);	// vals loaded into [0] .. [frameEnd-1 - frameBegin]
#endif
int   MRIsampleVolumeFrameType( const MRI *mri,
				const double x, const double y, const double z,
				const int frame,
				int interp_type,
				double *pval );
#ifdef FASTER_MRI_EM_REGISTER
int   MRIsampleVolumeFrameType_xyzInt_nRange_SAMPLE_NEAREST_floats(const MRI *mri,
                            int x, int y, int z, 
			    const int frameBegin,
			    const int frameEnd,		// [frameBegin] .. [frameEnd-1] done
			    float *valForEachFrame);	// vals loaded into [0] .. [frameEnd-1 - frameBegin]
#endif
int   MRIsampleVolumeGradient(MRI *mri, double x, double y, double z,
                              double *pdx, double *pdy, double *pdz) ;
int   MRIsampleVolumeGradientFrame( const MRI *mri,
				    double x, double y, double z,
				    double *pdx, double *pdy, double *pdz,
				    int frame ) ;
int   MRIsampleVolumeDerivative(MRI *mri,
                                double x, double y, double z,
                                double dx, double dy, double dz,
                                double *pmag) ;
int   MRIsampleVolumeDerivativeScale(MRI *mri,
                                     double x, double y, double z,
                                     double dx, double dy, double dz,
                                     double *pmag,
                                     double sigma) ;
int   MRIsampleVolumeDirectionScale(MRI *mri,
                                    double x, double y, double z,
                                    double dx, double dy, double dz,
                                    double *pmag,
                                    double sigma) ;
float MRIsampleCardinalDerivative(MRI *mri, int x, int y, int z,
                                  int xk, int yk, int zk) ;
float MRIsampleXDerivative(MRI *mri, int x, int y, int z, int dir) ;
float MRIsampleYDerivative(MRI *mri, int x, int y, int z, int dir) ;
float MRIsampleZDerivative(MRI *mri, int x, int y, int z, int dir) ;
MRI   *MRIxDerivative(MRI *mri_src, MRI *mri_dx) ;
MRI   *MRIyDerivative(MRI *mri_src, MRI *mri_dy) ;
MRI   *MRIzDerivative(MRI *mri_src, MRI *mri_dz) ;

/* resampling routines */
MRI   *MRIupsample2(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIupsampleN(MRI *mri_src, MRI *mri_dst, int N) ;
MRI   *MRIupsampleNConserve(MRI *mri_src, MRI *mri_dst, int N);
MRI   *MRIdownsampleN(MRI *src, MRI *dst, int Fc, int Fr, int Fs, int KeepType);
MRI   *MRIdownsample2(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIdownsample2LabeledVolume(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIresize(MRI *mri, double xsize, double ysize, double zsize, int nframes);

/* surfaceRAS and voxel routines */
MATRIX *surfaceRASFromVoxel_(MRI *mri);
MATRIX *voxelFromSurfaceRAS_(MRI *mri);
MATRIX *surfaceRASFromRAS_(MRI const *mri);
MATRIX *RASFromSurfaceRAS_(MRI const *mri, MATRIX *RASFromSRAS);

int MRIscannerRASToVoxel(MRI *mri, double xr, double yr, double zr, double *xv, double *yv, double *zv);
int MRIvoxelToSurfaceRAS(MRI *mri, double xv, double yv, double zv,
                         double *xs, double *ys, double *zs);
int MRIsurfaceRASToVoxel(MRI *mri, double xr, double yr, double zr,
                         double *xv, double *yv, double *zv);
int MRIsurfaceRASToVoxelCached(MRI *mri, double xr, double yr, double zr,
                               double *xv, double *yv, double *zv);
int MRIRASToSurfaceRAS(MRI *mri, double xr, double yr, double zr,
                       double *xsr, double *ysr, double *zsr);
int MRIsurfaceRASToRAS(MRI *mri, double xsr, double ysr, double zsr,
                       double *xr, double *yr, double *zr);

/* bitmap image access macros */
#define MRIset_bit(mri,x,y,z)    MRIvox(mri,(x)/8,y,z) |= (0x001 << ((x)%8))
#define MRItest_bit(mri,x,y,z)   (MRIvox(mri,(x)/8,(y),(z))&(0x001 << ((x)%8)))
#define MRIclear_bit(mri,x,y,z)  MRIvox(mri,(x)/8,y,z) &= ~(0x001 << ((x)%8))

#define MRISvox(mri,x,y,z)  (((short *)mri->slices[z][y])[x])
#define MRIUSvox(mri,x,y,z) (((unsigned short *)mri->slices[z][y])[x])
#define MRIFvox(mri,x,y,z)  (((float *)(mri->slices[z][y]))[x])
#define MRIvox(mri,x,y,z)   (((BUFTYPE *)mri->slices[z][y])[x])
#define MRISCvox(mri,x,y,z) (((signed char *)mri->slices[z][y])[x])
#define MRIIvox(mri,x,y,z)  (((int *)mri->slices[z][y])[x])
#define MRILvox(mri,x,y,z)  (((long32 *)mri->slices[z][y])[x])

#define MRISseq_vox(mri,x,y,z,n)  (((short*)\
mri->slices[z+(n)*mri->depth][y])[x])
#define MRIUSseq_vox(mri,x,y,z,n) (((unsigned short*)\
mri->slices[z+(n)*mri->depth][y])[x])
#define MRISCseq_vox(mri,x,y,z,n) (((signed char*)\
mri->slices[z+(n)*mri->depth][y])[x])
#define MRIFseq_vox(mri,x,y,z,n)  (((float*)\
(mri->slices[z+((n)*mri->depth)][y]))[x])
#define MRIseq_vox(mri,x,y,z,n)   (((BUFTYPE *)\
mri->slices[z+(n)*mri->depth][y])[x])
#define MRIIseq_vox(mri,x,y,z,n)  (((int *)\
mri->slices[z+(n)*mri->depth][y])[x])
#define MRILseq_vox(mri,x,y,z,n)  (((long32 *)\
mri->slices[z+(n)*mri->depth][y])[x])

#define MRI_HEIGHT      0
#define MRI_WIDTH       1
#define MRI_DEPTH       2

#define MRI_UNDEFINED   -1
#define MRI_CORONAL     0
#define MRI_SAGITTAL    1
#define MRI_HORIZONTAL  2
#define MRI_AXIAL       MRI_HORIZONTAL

/* vertices of an icosahedron (sp?), used by all POLV functions */
#define NVERTICES    22
extern float ic_x_vertices[]  ;
extern float ic_y_vertices[]  ;
extern float ic_z_vertices[]  ;

#define MRI_GZIPPED                  -2
#define MRI_VOLUME_TYPE_UNKNOWN      -1
#define MRI_CORONAL_SLICE_DIRECTORY   0
#define MRI_MINC_FILE                 1
#define MRI_ANALYZE_FILE              2
#define MRI_MGH_FILE                  3
#define GENESIS_FILE                  4
#define GE_LX_FILE                    5
#define SIEMENS_FILE                  6
#define BRIK_FILE                     7
#define BSHORT_FILE                   8
#define BFLOAT_FILE                   9
#define SDT_FILE                      10
#define OTL_FILE                      11
#define GDF_FILE                      12
#define RAW_FILE                      13
#define SIGNA_FILE                    14
#define DICOM_FILE                    15
#define MRI_ANALYZE4D_FILE            16
#define SIEMENS_DICOM_FILE            17
#define BRUKER_FILE                   18
#define XIMG_FILE                     19
#define NIFTI1_FILE                   20 // NIfTI-1 .img + .hdr
#define IMAGE_FILE                    21
#define MRI_GCA_FILE                  22
#define BHDR                          23 // for bshort or bfloat
#define NII_FILE                      24 // NIfTI-1 .nii (single file)
#define MRI_CURV_FILE                 25 // surface curv format
#define NRRD_FILE                     26 // NRRD .nrrd single file
#define GIFTI_FILE                    27 // GIFTI func data frames
#define VTK_FILE                      28 // VTK
#define MGH_MORPH                     29 // .m3z, .m3d
#define MGH_AUTOENCODER               30
#define ITK_MORPH                     31 // ITK (e.g., ANTs synWarp)
#define MGH_LABEL_FILE                32
#define MGH_ANNOT                     33 // surface annotation file
#define NPY_FILE                      34 // numpy npy file
#define ASCII_FILE                    35 // .asc

int        MRImatchDimensions(MRI *mri1, MRI *mri2) ;
int        MRImatch(MRI *mri1, MRI *mri2) ;
int        MRInonzeroValRange(MRI *mri, float *pmin, float *pmax) ;
int        MRIvalRange(MRI *mri, float *pmin, float *pmax) ;
double     MRImaxNorm(MRI *mri) ;

int        MRIlabelValRange(MRI *mri,
                            MRI *mri_labeled,
                            int label,
                            float *pmin,
                            float *pmax) ;
int        MRIvalRangeFrame(MRI *mri, float *pmin, float *pmax, int frame) ;
MRI        *MRIvalScale(MRI *mri_src, MRI *mri_dst, float fmin, float fmax) ;

#include "histo.h" // HISTOGRAM
double MRIfindPercentile(MRI *mri, double percentile, int frame) ;
HISTOGRAM *MRIhistogramVoxel(MRI *mri,
                             int nbins,
                             HISTOGRAM *histo,
                             int x0, int y0, int z0,
                             int wsize,
                             MRI *mri_thresh,
                             float thresh) ;
int        *MRIhistogramLabels(MRI *mri, int *counts, int max_label);
HISTOGRAM  *MRIhistogram(MRI *mri, int nbins) ;
HISTOGRAM  *MRIhistogramLabel(MRI *mri, MRI *mri_labeled,
                              int label, int nbins);
HISTOGRAM *MRIhistogramRegionWithThreshold(MRI *mri,
					   int nbins,
					   HISTOGRAM *histo,
					   MRI_REGION *region,
					   MRI *mri_thresh, 
					   float thresh,
					   int frame);
HISTOGRAM  *MRIhistogramLabelRegion(MRI *mri,
                                    MRI *mri_labeled,
                                    MRI_REGION *region,
                                    int label, int nbins);
MRI        *MRIhistogramNormalize(MRI *mri_src, MRI *mri_template, MRI *mri_dst) ;
MRI        *MRIhistoEqualize(MRI *mri_src, MRI *mri_template, MRI *mri_dst,
                             int low, int high) ;
MRI        *MRIapplyHistogram(MRI *mri_src, MRI *mri_dst, HISTOGRAM *histo) ;
MRI        *MRIcrunch(MRI *mri_src, MRI *mri_dst) ;
HISTOGRAM  *MRIgetEqualizeHisto(MRI *mri, HISTOGRAM *histo_eq, int low,
                                int high, int norm) ;

/* these are adaptive (i.e. only operate on a subregion of the whole image */
MRI_REGION *MRIclipRegion(MRI *mri, MRI_REGION *reg_src, MRI_REGION *reg_clip);
int        MRIvalRangeRegion(MRI *mri, float *pmin, float *pmax,
                             MRI_REGION *region) ;
HISTOGRAM  *MRIhistogramRegion(MRI *mri, int nbins, HISTOGRAM *histo,
                               MRI_REGION *region) ;
MRI        *MRIhistoEqualizeRegion(MRI *mri_src, MRI *mri_dst, int low,
                                   MRI_REGION *region) ;
MRI        *MRIapplyHistogramToRegion(MRI *mri_src, MRI *mri_dst,
                                      HISTOGRAM *histo, MRI_REGION *region) ;
HISTOGRAM  *MRIgetEqualizeHistoRegion(MRI *mri, HISTOGRAM *histo_eq, int low,
                                      MRI_REGION *region, int norm) ;
int        MRIfileType(char *fname) ;
int        MRIunpackFileName(const char *inFname, int *pframe, int *ptype,
                             char *outFname) ;
int        MRIisValid(MRI *mri) ;
MRI        *MRIflipByteOrder(MRI *mri_src, MRI *mri_dst) ;

MRI        *MRIregionGrow(MRI *mri_src, MRI *mri_distance,
                          float x0, float y0, float z0, int niter) ;
MRI        *MRIextractInterior(MRI *mri_src, MRI *mri_distance,  MRI *mri_dst);
MRI        *MRIbuildDistanceMap(MRI *mri_src, MRI *mri_distance,
                                float x0, float y0, float z0, float r) ;
MRI        *MRIupdateDistanceMap(MRI *mri_distance) ;
int         MRIvalueFill(MRI *mri, float value);
MRI        *MRIfill(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
                    int seed_z, int threshold, int fill_val, int max_count) ;
MRI        *MRIfillFG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
                      int seed_z, int threshold, int fill_val, int *npix) ;
MRI        *MRIfillBG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
                      int seed_z, int threshold, int fill_val, int *npix) ;

int   MRIneighborsInRange(MRI *mri,
                          int x0, int y0, int z0,
                          int frame,
                          float  low_val, float hi_val) ;
int   MRIneighbors3x3(MRI *mri, int x, int y, int z, int val) ;
int   MRIneighbors(MRI *mri, int x, int y, int z, int val) ;
int   MRIneighborsOn(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOff(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOn3x3(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOff3x3(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label) ;
int   MRIlabelsInNbhd6(MRI *mri, int x, int y, int z, int label) ;
int   MRIlabelsInPlanarNbhd(MRI *mri,
                            int x, int y, int z,
                            int whalf, int label, int which) ;

int MRIvol2Vol(MRI *src, MRI *targ, MATRIX *Vt2s,
               int InterpCode, float param);
int MRIvol2VolR(MRI *src, MRI *targ, MATRIX *Vt2s,
               int InterpCode, float param, MATRIX* RRot);

MRI *MRIresampleFill(MRI *src, MRI *template_vol,
                     int resample_type, float fill_val) ;
MRI *MRIreplaceList(MRI *seg, int *srclist, int *targlist, int nlist, MRI *mask, MRI *out);
MRI *MRIreplaceValuesOnly(MRI *mri_src, MRI *mri_dst,float in_val, float out_val) ;
MRI   *MRIreplaceValues(MRI *mri_src, MRI *mri_dst, float in_val, float out_val) ;
MRI   *MRIreplaceValueRange(MRI *mri_src, MRI *mri_dst,
                            float low_in_val, float hi_in_val, float out_val) ;
MRI   *MRIreplaceValuesUchar(MRI *mri_src, MRI *mri_dst,
                             BUFTYPE in_val, BUFTYPE out_val) ;
int MRImaskLabel(MRI *mri_src, MRI *mri_dst, MRI *mri_labeled, int label_to_mask, float out_val)  ;
MRI   *MRImask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int mask,
               float out_val) ;
MRI   *MRImaskZero(MRI *mri_src, MRI *mri_mask, MRI *mri_dst) ;
MRI   *MRImaskDifferentGeometry(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                                int mask, float out_val) ;
MRI *MRImaskInvert(MRI *mask, MRI *outmask);

int MRInMask(MRI *mask);
MRI *MRIframeBinarize(MRI *mri, double thresh, MRI *mask);


MRI   *MRImeanMask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                   int mask, int wsize) ;
MRI   *MRIthresholdMask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                        float mask_threshold, float out_val) ;

/* constants used in mri_dir of MRIoffsetDirection and for MRIminmax filter */
#define OFFSET_NEGATIVE_GRADIENT_DIRECTION    0
#define OFFSET_GRADIENT_DIRECTION             1
#define OFFSET_ZERO                           2

/* anything below this is not white matter */
#define WM_MIN_VAL                       5
#define MIN_WM_VAL                       WM_MIN_VAL
#define WM_EDITED_ON_VAL                 255
#define WM_EDITED_OFF_VAL                1

// For labeling bright voxels and neighbors of bright voxels
#define BRIGHT_LABEL 130
#define BRIGHT_BORDER_LABEL 100


MRI *MRIreduceMeanAndStd(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIreduceMeanAndStdByte(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIstdsToVariances(MRI *mri_std, MRI *mri_var, int source_frame) ;
MRI *MRIvariancesToStds(MRI *mri_var, MRI *mri_std, int dst_frame) ;
MRI *MRIconcatenateFrames(MRI *mri_frame1, MRI *mri_frame2, MRI *mri_dst);
MRI *MRIcopyFrame(MRI *mri_src, MRI *mri_dst, int src_frame, int dst_frame) ;
double MRImeanFrameNonzeroMask(MRI *mri, int frame, MRI *mri_mask) ;
MRI *MRIcopyFrames(MRI *mri_src, MRI *mri_dst, int src_start_frame, int src_end_frame, int dst_start_frame) ;
double MRImeanFrame(MRI *mri, int frame) ;

int   MRIcountPlanarAboveThreshold(MRI *mri_src, int vertex, int x, int y,
                                   int z, int wsize, int lo_lim, int hi_lim);
int   MRIcountCpolvAtVoxel(MRI *mri_src, int x, int y, int z, int wsize,
                           int *pnum, int label_to_check) ;

MRI *MRIhistoSegment(MRI *mri_src, MRI *mri_labeled, int wm_low, int wm_hi,
                     int gray_hi, int wsize, float sigma) ;
MRI *MRIhistoSegmentVoxel(MRI *mri_src, MRI *mri_labeled, int wm_low,
                          int wm_hi, int gray_hi, int wsize, int x, int y,
                          int z, HISTOGRAM *histo, HISTOGRAM *hsmooth,
                          float sigma) ;
MRI *MRIcpolvSmooth(MRI *mri_orig,MRI *mri_src, MRI *mri_dst, int wsize,
                    int low_val, int hi_val, int niter) ;
MRI *MRIextractVertexCoords(MRI *mri_src, int *px, int *py, int *pz,
                            int vertex, int x, int y,int z, int wsize) ;

MRI *MRIextractVertexPlane(MRI *mri_src, MRI *mri_dst, int vertex, int x,
                           int y, int z, int wsize) ;
int  MRIwhiteInPlane(MRI *mri_src, int x, int y, int z, int vertex, int wsize);
int  MRIneighborhoodPlanarDirection(MRI *mri_src, int xv, int yv, int zv,
                                    int nsize, int wsize) ;
int  MRIneighborhoodCpolv(MRI *mri_src, int xv, int yv, int zv,int nsize,
                          int wsize, int *pnum_white) ;
int MRIneighborsInWindow(MRI *mri, int x, int y, int z, int wsize, int val) ;
int  MRIneighborhoodBlackCpolv(MRI *mri_src, int xv, int yv, int zv,
                               int nsize, int wsize, int *pnum_black) ;

MRI *MRIorderSegment(MRI *mri_src, MRI *mri_labeled, float thresh, int wsize);
MRI *MRIthresholdLabel(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst,
                       int wm_low) ;
MRI *MRIintensitySegmentation(MRI *mri_src, MRI *MinWMMask, MRI *mri_labeled, float wm_low, float wm_hi, float gray_hi);
MRI *MRImeanLabel(MRI *mri_src, MRI *mri_label, MRI*mri_dst, int wsize) ;
int MRIvoxelsInLabel(MRI *mri, int label) ;
int MRItotalVoxelsOn(MRI *mri, int thresh) ;
int MRIcopyLabel(MRI *mri_src, MRI *mri_dst, int val) ;
int MRIsetVoxelsWithValue(MRI *mri_src,
                          MRI *mri_dst,
                          int src_val,
                          int dst_val) ;
int MRIsetDifferentVoxelsWithValue(MRI *mri1,
                                   MRI *mri2,
                                   MRI *mri_dst,
                                   int dst_val) ;
double MRIpercentThresh(MRI *mri, MRI *mask, int frame, double pct);
int MRIcopyLabeledVoxels(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst,
                         int label) ;
MRI *MRIcpolvVote(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int wsize,
                  int niter, int use_all) ;
MRI *MRIcpolvThreshold(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst,
                       int wm_low, int gray_hi,int wsize) ;
MRI *MRImaskLabels(MRI *mri_src, MRI *mri_mask, MRI *mri_dst) ;
MRI *MRIsphereMask(int ncols, int nrows, int nslices, int nframes,
                   int c0, int r0, int s0, double voxradius, double val,
                   MRI *mri);
MRI *MRIwmfilterMarked(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int wsize,
                       float pct, int onoff) ;
int  MRIcountCpolvAtVoxel(MRI *mri_src, int x, int y, int z, int wsize,
                          int *pnum, int label_to_check) ;
int  MRIcountCpolvOnAtVoxel(MRI *mri_src, int x, int y, int z, int wsize,
                            int *pnum) ;
double MRImeanAndStdInLabel(MRI *mri_src, MRI *mri_labeled, int label, double *pstd) ;
MRI *MRIcentralPlaneOfLeastVarianceNormalMarked(MRI *mri_src, MRI *mri_mask,
    MRI *mri_dst, int wsize) ;
int  MRIcountCpolvOffAtVoxel(MRI *mri_src,int x, int y, int z, int wsize,
                             int *pnum) ;
int  MRIcentralPlaneOfLeastVarianceNormalVoxel(MRI *mri_src, int wsize,
    int x, int y, int z) ;
int  MRIvertexToVector(int vertex, float *pdx, float *pdy, float *pdz) ;
int MRIcpolvMedianCurveVoxel(MRI *mri, MRI *mri_labeled, int x0, int y0,
                             int z0, int wsize, float len, float gray_hi,
                             float wm_low) ;
float MRIcpolvMedianAtVoxel(MRI *mri_src, int vertex,
                            float x, float y, float z, int wsize);
MRI   *MRIcpolvMedianCurveSegment(MRI *mri,MRI *mri_labeled, MRI *mri_dst,
                                  int wsize,float len, float gray_hi,
                                  float wm_low);
MRI   *MRImarkBorderVoxels(MRI *mri_src, MRI *mri_dst) ;
int   MRIborderClassifyVoxel(MRI *mri_src, MRI *mri_labeled, int wsize, int x,
                             int y, int z, float *ppw, float *ppg) ;
int   MRIreclassifyBorder(MRI *mri_src, MRI *mri_labeled, MRI *mri_border,
                          MRI *mri_dst, int wsize) ;

int   MRIclassifyAmbiguous(MRI *mri_src, MRI *mri_labeled, MRI *mri_border,
                           MRI *mri_dst, int wsize) ;
MRI   *MRIremoveBrightStuff(MRI *mri_src, MRI *mri_dst, int threshold) ;
int   MRIreclassify(MRI *mri_src, MRI *mri_labeled,
                    MRI *mri_dst, float wm_low, float gray_hi, int wsize) ;
MRI *MRIreclassifyWMCtxNonPar(MRI *norm, MRI *seg, int nitersmax, MRI *newseg);
MRI *MRIclassifyAmbiguousNonPar(MRI *norm, MRI *seg, MRI *statseg, double Qwm, double Qctx, int NdilWM, int NdilCtx);
MRI *MRImaskThreshold(MRI *mri_src, MRI *mri_mask, MRI *mri_dst,
                      float threshold, int out_label) ;
int MRIgrowLabel(MRI *mri, MRI *mri_bg, int in_label, int out_label) ;
int MRIturnOnFG(MRI *mri, MRI *mri_fg, MRI *mri_bg) ;
int MRIturnOffBG(MRI *mri, MRI *mri_bg) ;
/* mriprob.c */
MRI *MRIcomputeConditionalProbabilities(MRI *mri_T1, MRI *mri_mean,
                                        MRI *mri_std, MRI *mri_dst) ;
MRI *MRIapplyBayesLaw(MRI *mri_priors, MRI *mri_p1, MRI *mri_p2,MRI *mri_dst);
MRI *MRIprobabilityThresholdNeighborhoodOff(MRI *mri_src, MRI *mri_prob,
    MRI *mri_dst, float threshold,
    int nsize) ;
MRI *MRIprobabilityThresholdNeighborhoodOn(MRI *mri_src, MRI *mri_prob,
    MRI *mri_dst, float threshold,
    int nsize, int out_label) ;
MRI *MRIprobabilityThreshold(MRI *mri_src, MRI *mri_prob, MRI *mri_dst,
                             float threshold, int out_label) ;
MRI *MRIdilateLabel(MRI *mri_src, MRI *mri_dst, int label, int niter) ;
MRI *MRIdilateLabelUchar(MRI *mri_src, MRI *mri_dst, int label, int niter) ;
MRI *MRIdilateThreshLabel(MRI *mri_src, MRI *mri_val, MRI *mri_dst, int label,
                          int niter, int thresh) ;
MRI *MRIdilateInvThreshLabel(MRI *mri_src, MRI *mri_val, MRI *mri_dst,
                             int label,
                             int niter, int thresh) ;
MRI *MRIsoapBubbleLabel(MRI *mri_src, MRI *mri_label, MRI *mri_dst,
                        int label,
                        int niter);
MRI    *MRIsetLabelValues(MRI *mri_src, MRI *mri_label, MRI *mri_dst,
                          int label, float val);
int    MRIwriteImageViews(MRI *mri, const char *base_name, int target_size) ;
int MRIsetValues(MRI *mri, float val) ;
int MRIsetFrameValues(MRI *mri, int frame, float val) ;
MRI    *MRIwindow(MRI *mri_src, MRI *mri_dst, int which, float x0, float y0,
                  float z0, float parm) ;
int  MRIcomputeClassStatistics(MRI *mri_T1, MRI *mri_labeled,
                               float gray_low, float gray_hi,
                               float *pmean_wm, float *psigma_wm,
                               float *pmean_gm, float *psigma_gm) ;
int MRIcomputeClassStatisticsSeg(MRI *norm, MRI *seg, float *wmmean, float *wmstd, float *ctxmean, float *ctxstd);
int MRIcomputeNbhdMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, int label, int x, int y, int z, int nbhd, MATRIX **p_mcov, VECTOR **p_vmeans) ;
int MRIcomputeLabelMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, MATRIX **p_mcov, VECTOR **p_vmeans, int *labels, int nlabels, int nbhd);
int MRIcomputeWMMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, MATRIX **p_mcov, VECTOR **p_vmeans) ;
int MRIcomputeWMSAMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, MATRIX **p_mcov, VECTOR **p_vmeans) ;
int MRIcomputeCaudateMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, MATRIX **p_mcov, VECTOR **p_vmeans) ;
int MRIcomputeVentMeansandCovariances(MRI *mri_inputs, MRI *mri_labeled, MATRIX **p_mcov, VECTOR **p_vmeans) ;

#define WINDOW_GAUSSIAN  0
#define WINDOW_HAMMING   1
#define WINDOW_HANNING   2

#define MRI_NOT_WHITE   1
#define MRI_AMBIGUOUS   128
#define MRI_WHITE       255

#define MRI_LEFT_HEMISPHERE     255
#define MRI_RIGHT_HEMISPHERE    127
#define MRI_RIGHT_HEMISPHERE2   80

#define MRI_NONBRAIN            0
#define MRI_PIAL_INTERIOR       1
#define MRI_WHITE_INTERIOR      2

/* STATS volumes have 2 images each */
#define TRI_HI_PRIORS            0
#define TRI_HI_STATS             1
#define TRI_LOW_PRIORS           3
#define TRI_LOW_STATS            4
#define TRI_OFF_STATS            6

#define REMOVE_1D          2
#define REMOVE_WRONG_DIR   3

#define BASAL_GANGLIA_FILL   50
#define MAX_WM_VAL           (THICKEN_FILL-1)
#define WM_MAX_VAL           MAX_WM_VAL
#define THICKEN_FILL         200
#define NBHD_FILL            210
#define VENTRICLE_FILL       220
#define DIAGONAL_FILL        230
#define DEGENERATE_FILL      240
#define OFFSET_FILTER_FILL   245
#define AUTO_FILL            250
#define PRETESS_FILL         215

MRI *MRIchangeType(MRI *src, int dest_type, float f_low,
                   float f_high, int no_scale_option_flag);
MRI *MRISeqchangeType(MRI *vol, int dest_type, float f_low,
                      float f_high, int no_scale_option_flag);

MRI *MRIresample(MRI *src, MRI *template_vol, int resample_type);
MATRIX *MRIgetResampleMatrix(MRI *src, MRI *template_vol);
void MRIlimitsMultipleTimes(MRI *mri_src, float *psrc_min, float *psrc_max, int ntimes);
int MRIlimits(MRI *mri, float *min, float *max);
int MRIprintStats(MRI *mri, FILE *stream);
int MRIstats(MRI *mri, float *min, float *max, int *n_voxels,
             float *mean, float *std);

float MRIvolumeDeterminant(MRI *mri);

MRI *MRISreadCurvAsMRI(const char *curvfile, int read_volume);

int mriio_command_line(int argc, char *argv[]);
void mriio_set_gdf_crop_flag(int new_gdf_crop_flag);
int MRIgetVolumeName(const char *string, char *name_only);
MRI *MRIread(const char *fname, std::vector<MRI*> *mriVector=NULL);
MRI *MRIreadEx(const char *fname, int nthframe);
MRI *MRIreadType(const char *fname, int type);
MRI *MRIreadInfo(const char *fname);
MRI *MRIreadHeader(const char *fname, int type);
int GetSPMStartFrame(void);
int MRIwrite(MRI *mri,const  char *fname, std::vector<MRI*> *mriVector=NULL);
int MRIwriteFrame(MRI *mri,const  char *fname, int frame) ;
int MRIwriteType(MRI *mri,const  char *fname, int type);
MRI *MRIreadRaw(FILE *fp, int width, int height, int depth, int type);
int MRIreorderVox2RAS(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim);
MRI *MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim);
MRI *MRIreorder4(MRI *mri, int order[4]);
MRI *MRIsmoothParcellation(MRI *mri, int smooth_parcellation_count);
MRI *MRIsmoothLabel(MRI *mri_intensity,
                    MRI *mri_label,
                    MRI *mri_smooth,
                    int niter,
                    int label,
		    float min_change) ;
MRI *MRIsmoothLabel6Connected(MRI *mri_intensity,
			      MRI *mri_label,
			      MRI *mri_smooth,
			      int niter,
			      int label,
			      int fixed_label,
			      float min_change) ;
MRI *MRIreadGeRoi(const char *fname, int n_slices);

int decompose_b_fname(const char *fname_passed, char *directory, char *stem);

#define READ_OTL_READ_VOLUME_FLAG       0x01
#define READ_OTL_FILL_FLAG              0x02
#define READ_OTL_TRANSLATE_LABELS_FLAG  0x04
#define READ_OTL_ZERO_OUTLINES_FLAG     0x08
MRI *MRIreadOtl(const char *fname, int width, int height, int slices,
                const char *color_file_name, int flags);

int list_labels_in_otl_file(FILE *fp);

MATRIX *extract_i_to_r(const VOL_GEOM *mri);
int apply_i_to_r(VOL_GEOM *mri, MATRIX *m);

int stuff_four_by_four(MATRIX *m,
                       float m11, float m12, float m13, float m14,
                       float m21, float m22, float m23, float m24,
                       float m31, float m32, float m33, float m34,
                       float m41, float m42, float m43, float m44);

MATRIX *extract_r_to_i(const VOL_GEOM *mri) ;
#define MRIgetVoxelToRasXform   extract_i_to_r
#define MRIgetRasToVoxelXform   extract_r_to_i
int    MRIsetVoxelToRasXform(MRI *mri, MATRIX *m_vox2ras) ;
MATRIX *MRIvoxelXformToRasXform(MRI *mri_src, MRI *mri_dst,
                                MATRIX *m_voxel_xform, MATRIX *m_ras_xform);
MATRIX *MRIrasXformToVoxelXform(MRI *mri_src, MRI *mri_dst,
                                MATRIX *m_ras_xform, MATRIX *m_voxel_xform);


int MRIsincSampleVolume( const MRI *mri,
                         double x, double y, double z,
                         int hw, double *pval );
int MRIcubicSampleVolume( const MRI *mri,
                          double x, double y, double z,
                          double *pval ); /*E*/
DMATRIX *MRIgradCubicInterp(const MRI *mri, double x, double y, double z, DMATRIX *grad);
double MRIcubicCoeff(double x, int iter);
double MRIgradCubicCoeff(double x, int iter);

MRI *MRIsincTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA, int hw);
int MRIlabelOverlap(MRI *mri1, MRI *mri2, int label) ;
int MRIlabelUnion(MRI *mri1, MRI *mri2, int label) ;
int MRIeraseBorderPlanes(MRI *mri, int mask_size) ;

MRI *MRIzeroMeanTimecourse(MRI *mri_src, MRI *mri_dst);
MRI *MRImeanTimecourse(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIzeroMean(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIlog10(MRI *inmri, MRI *mask, MRI *outmri, int negflag);
MRI *MRIlog(MRI *in, MRI *mask, double a, double b, MRI *out);
MRI *MRIrandn(int ncols, int nrows, int nslices, int nframes,
              float avg, float stddev, MRI *mri);
MRI *MRIrandexp(MRI *mrimean, MRI *binmask, unsigned long int seed, int nreps, MRI *mrirandexp);
MRI *MRIrande(int ncols, int nrows, int nslices, int nframes,
              float avg, int order, MRI *mri);
MRI *MRIdrand48(int ncols, int nrows, int nslices, int nframes,
                float min, float max, MRI *mri);
MRI *MRIsampleCDF(int ncols, int nrows, int nslices, int nframes,
                  double *xCDF, double *CDF, int nCDF, MRI *mri);
MRI *MRIconst(int ncols, int nrows, int nslices, int nframes,
              float val, MRI *mri);
int MRInormalizeSequence(MRI *mri, float target) ;

int setDirectionCosine(MRI *mri, int orientation);
int getSliceDirection(MRI *mri);
int mriOKforSurface(MRI *mri); // check whether the volume is conformed or not
int mriConformed(MRI *mri) ;
int MRIfindRightSize(VOL_GEOM *mri, float conform_size);
float MRIfindMinSize(VOL_GEOM *mri, int *conform_width);

void setMRIforSurface(MRI *mri); // set c_(r,a,s) = 0 for a conformed volume
MRI *MRIremoveNaNs(MRI *mri_src, MRI *mri_dst) ;
MRI *MRImakePositive(MRI *mri_src, MRI *mri_dst);
MRI *MRIeraseNegative(MRI *mri_src, MRI *mri_dst) ;

MRI *MRImarkLabelBorderVoxels( const MRI *mri_src, MRI *mri_dst,
			       int label, int mark, int six_connected) ;
int MRIcomputeLabelNbhd( const MRI *mri_labels, const MRI *mri_vals,
			 const int x, const int y, const int z,
			 int *label_counts, float *label_means,
			 const int whalf, const int max_labels );
float MRIvoxelsInLabelWithPartialVolumeEffects( const MRI *mri,
						const MRI *mri_vals,
						int label,
						MRI *mri_mixing_coef,
						MRI *mri_nbr_labels );
MRI   *MRImakeDensityMap(MRI *mri, MRI *mri_vals, int label, MRI *mri_dst,
                         float orig_res) ;
int MRIfillBox(MRI *mri, MRI_REGION *box, float fillval) ;
int MRIcropBoundingBox(MRI *mri, MRI_REGION *box) ;
MRI *MRIapplyBiasCorrection(MRI *mri_in, MRI *mri_bias, MRI *mri_out) ;
MRI *MRIapplyBiasCorrectionSameGeometry(MRI *mri_in, MRI *mri_bias,
                                        MRI *mri_out, float target_val) ;
MATRIX *MRIgetVoxelToVoxelXform(MRI *mri_src, MRI *mri_dst) ;

/* extract the RASToVoxeMatrix from an MRI */
MATRIX *GetSurfaceRASToVoxelMatrix(MRI *mri);

// functions read/write MRI_MGH_FILE
MRI *mghRead(const char *fname, int read_volume, int frame);
int mghWrite(MRI *mri, const char *fname, int frame);

/* Zero-padding for 3d analyze (ie, spm) format */
#ifdef _MRIIO_SRC
int N_Zero_Pad_Input  = -1;
int N_Zero_Pad_Output = -1;
int MRIIO_Strip_Pound = 1;
#else
extern int N_Zero_Pad_Input;
extern int N_Zero_Pad_Output;
extern int MRIIO_Strip_Pound;
#endif

float MRIfovCol(MRI *mri);
int MRIdircosToOrientationString(VOL_GEOM *mri, char *ostr);
int MRIorientationStringToDircos(MRI *mri, const char *ostr);
char *MRIcheckOrientationString(const char *ostr);
char *MRIsliceDirectionName(MRI *mri);
MRI *MRIreverseSliceOrder(MRI *invol, MRI *outvol);
MRI *MRIconformSliceOrder(MRI *mri);

/* different modes for distance transform - signed
   (<0 in interior) unsigned from border, or
   just outside (interior is 0) */
#define DTRANS_MODE_SIGNED   1
#define DTRANS_MODE_UNSIGNED 2
#define DTRANS_MODE_OUTSIDE  3
#define DTRANS_MODE_INSIDE   4

/** This is deprecated.
    Please use MRIextractDistanceMap in fastmarching.h instead */
MRI *MRIdistanceTransform(MRI *mri_src, MRI *mri_dist,
                          int label, float max_dist, int mode, MRI *mri_mask);
int MRIaddCommandLine(MRI *mri, const std::string& cmdline);
MRI *MRInonMaxSuppress(MRI *mri_src, MRI *mri_sup,
                       float thresh, int thresh_dir) ;
MRI *MRIextractRegionAndPad(MRI *mri_src, MRI *mri_dst,
                            MRI_REGION *region, int pad) ;
MRI *MRIsetValuesOutsideRegion(MRI *mri_src,
                               MRI_REGION *region,
                               MRI *mri_dst,
                               float val) ;
double MRImeanAndVarianceInNbhd(MRI *mri, int wsize, int x, int y, int z, int frame, double *pvar) ;
int MRIcountNonzeroInNbhd(MRI *mri, int wsize, int x, int y, int z) ;
int MRImeanNonzeroInNbhd(MRI *mri, int wsize, int x, int y, int z, int frame) ;
int MRImaxInNbhd(MRI *mri, int wsize, int x, int y, int z, int frame) ;
int MRImaxInNbhd6Connected(MRI *mri, int x, int y, int z, int frame) ;
int MRIcountValInNbhd(MRI *mri, int wsize, int x, int y, int z, int val) ;
int MRIcountNonzero(MRI *mri) ;
int MRIcountThreshInNbhd(MRI *mri, int wsize, int x,int y,int z, float thresh);
MRI *MRImatchMeanIntensity(MRI *mri_source,
                           MRI *mri_target,
                           MRI *mri_source_scaled) ;
MRI *MRIsqr(MRI *mri_src, MRI *mri_dst)  ;
MRI *MRIsqrt(MRI *mri_src, MRI *mri_dst)  ;
float MRImaxInRegion(MRI *mri, int x, int y, int z, int whalf) ;
double MRImaxInLabelInRegion(MRI *mri_src,
                             MRI *mri_labeled,
                             int label,
                             int x0, int y0, int z0,
                             int whalf) ;

double MRIestimateTIV(char* theLtaFile,
                      double theScaleFactor,
                      double* theAtlasDet);

int MRInormalizeFrames(MRI *mri);
int MRInormalizeFramesMean(MRI *mri);
int MRInormalizeFramesFirst(MRI *mri);
MRI *MRIaverageFrames(MRI *mri_src, MRI *mri_dst, int start_frame, int end_frame) ;
MRI *MRIsort(MRI *in, MRI *mask, MRI *sorted);
int CompareDoubles(const void *a, const void *b);
int MRIlabeledVoxels(MRI *mri_src, int label) ;
int MRIlabelInVolume(MRI *mri_src, int label) ;
#define MRI_MEAN_MIN_DISTANCE 0
double MRIcomputeLabelAccuracy(MRI *mri_src, MRI *mri_ref,
                               int which, FILE *fp) ;
double MRIcomputeMeanMinLabelDistance(MRI *mri_src, MRI *mri_ref, int label) ;
int MRIcomputeLabelCentroid(MRI *mri_aseg, int label,
														double *pxc, double *pyc, double *pzc) ;
  int MRIkarcherMean(MRI *mri, float low_val, float hi_val, int *px, int *py, int *pz) ;
int MRIcomputeCentroid(MRI *mri, double *pxc, double *pyc, double *pzc) ;
MRI *MRIdivideAseg(MRI *mri_src, MRI *mri_dst, int label, int nunits);
int MRIgeometryMatched(MRI *mri1, MRI *mri2) ;

MRI *MRIsegmentationSurfaceNormals(MRI *mri_seg,
                                   MRI *mri_normals,
                                   int target_label,
                                   MRI **pmri_ctrl) ;
MRI *MRIbinMaskToCol(MRI *binmask, MRI *bincol);
MRI *MRIfillHoles(MRI *mri_src, MRI *mri_fill, int thresh)  ;
int  MRIfillRegion(MRI *mri, int x,int y,int z,float fill_val,int whalf) ;
MRI *MRIfloodFillRegion(MRI *mri_src, MRI *mri_dst,
                        int threshold, int fill_val, int max_count) ;
int  MRIcomputeBorderNormalAtVoxel(MRI *mri_seg, int x0, int y0, int z0,
                                   float *pnx, float *pny, float *pnz,
                                   int label) ;
MRI  *MRImatchIntensityRatio(MRI *mri_source, MRI *_target, MRI *mri_matched,
                             double min_scale, double max_scale,
                             double low_thresh, double high_thresh);

// types of MRI sequences
#define MRI_UNKNOWN          0
#define MRI_MGH_MPRAGE       1
#define MRI_ADNI_MPRAGE      2
#define MRI_WASHU_MPRAGE     3
#define MRI_MIND_MPRAGE      4


MRI *MRIcreateDistanceTransforms(MRI *mri, MRI *mri_all_dtrans,
                                 float max_dist,
                                 int *labels, int nlabels);
MRI *MRIapplyMorph(MRI *mri_source,
                   MRI *mri_morph,
                   MRI *mri_dst,
                   int sample_type);
int MRIorderIndices(MRI *mri,
                    short *x_indices,
                    short *y_indices,
                    short *z_indices) ;
int MRIcomputeVoxelPermutation(MRI *mri,
                               short *x_indices,
                               short *y_indices,
                               short *z_indices);
MRI *MRInormalizeInteriorDistanceTransform(MRI *mri_src_dist,
                                           MRI *mri_ref_dist,
                                           MRI *mri_dst_dist);

const char* MRItype2str(int type);
MRI *MRIaddNoise(MRI *mri_in, MRI *mri_out, float amp);
int MRIfindSliceWithMostStructure(MRI *mri_aseg, int slice_direction, int label) ;
int MRIcomputeVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox,
			      MRI *mri_seg, MRI *mri_fractions) ;

MRI *MRInbrThresholdLabel(MRI *mri_src, MRI *mri_dst,  int label, int out_label, int whalf,  float thresh) ;
MRI *MRIsolveLaplaceEquation(MRI *mri_interior, MRI *mri_seg, int source_label, int target_label,
			       float source_val,float target_val, float outside_val);

int MRIsampleVolumeFrameMasked(const MRI *mri, const MRI *mri_mask, double x, double y, double z, const int frame, double *pval);

int MRIclipBrightWM(MRI *mri_T1, const MRI *mri_wm);
std::vector<std::vector<double>> *MRIdiff2Vect(const MRI *mri1, const MRI *mri2, const double thresh, const MRI *mask);
MRI *MRIapplyDiffVect(MRI *mri1, const MRI *mask, std::vector<std::vector<double>> diffvec, int *napply, MRI *mri2);

#endif



