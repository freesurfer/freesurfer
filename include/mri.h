#ifndef MRI_H
#define MRI_H

#include "const.h"
#include "matrix.h"
#include "volume_io.h"
#include "box.h"

#define BUFTYPE  unsigned char

#define MRI_UCHAR   0
#define MRI_INT     1
#define MRI_LONG    2
#define MRI_FLOAT   3
#define MRI_SHORT   4

typedef struct
{
  int  x ;
  int  y ;
  int  z ;
  int  dx ;
  int  dy ;
  int  dz ;
} MRI_REGION ;

typedef struct
{
  int           width ;
  int           height ;
  int           depth ;     /* # of slices */
  int           type ;      /* data type for slices below */
  int           slice_direction ;   /* initially coronal, but may change */
  int           imnr0 ;     /* starting image # */
  int           imnr1 ;     /* ending image # */
  int           ptype ;     /* not used */
  float         fov ;
  float         thick ;
  float         ps ;   
  float         location ;  /* not used */
  float         xsize ;     /* size of a voxel in the x direction */ 
  float         ysize ;     /* size of a voxel in the y direction */ 
  float         zsize ;     /* size of a voxel in the z direction */ 
  float         xstart ;    /* start x (in xsize units) */
  float         xend ;      /* end x  (in xsize units) */
  float         ystart ;    /* start y   (in ysize units) */  
  float         yend ;      /* end y (in ysize units) */ 
  float         zstart ;    /* start z */  
  float         zend ;      /* end z */
  float         tr ;        /* time to recovery */
  float         te ;        /* time to echo */
  float         ti ;        /* time to inversion */
  char          fname[STR_LEN] ;

/* 
  each slice is an array of rows (mri->height of them) each of which is mri->width 
  long.
*/
  BUFTYPE       ***slices ;
  int           scale ;
  char          transform_fname[STR_LEN] ;
  General_transform transform ;   /* the next two are from this struct */
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
  int           free_transform ;   /* are we responsible for freeing it? */
  int           nframes ;          /* # of concatenated images */

  /* these are used to handle boundary conditions (arrays of indices) */
  int           *xi ;
  int           *yi ;
  int           *zi ;
  int           yinvert ;  /* for converting between MNC and coronal slices */
  MRI_REGION    roi ;
} MRI_IMAGE, MRI ;

/* single pixel filtering */
float MRIvoxelMean(MRI *mri, int x, int y, int z, int wsize) ;
float MRIvoxelStd(MRI *mri, int x, int y, int z, float mean, int wsize) ;
float MRIvoxelZscore(MRI *mri, int x, int y, int z, int wsize) ;
float MRIvoxelDx(MRI *mri, int x, int y, int z) ;
float MRIvoxelDy(MRI *mri, int x, int y, int z) ;
float MRIvoxelDz(MRI *mri, int x, int y, int z) ;
float MRIvoxelGradient(MRI *mri, int x, int y, int z, float *pdx, float *pdy, 
                       float *pdz) ;
float MRIvoxelDirection(MRI *mri, int x, int y, int z, int wsize) ;
/* use these constants for MRIreorder */
#define XDIM  1
#define YDIM  2
#define ZDIM  3
MRI  *MRIreorder(MRI *mri_src, MRI *mri_dst, int xdim, int ydim, int zdim);

/* I/O functions */
int    MRIwrite(MRI *mri, char *fpref) ;
int    MRIwriteInfo(MRI *mri, char *fpref) ;
MRI   *MRIread(char *fpref) ;
MRI   *MRIreadInfo(char *fpref) ;

/* memory allocation routines */
int   MRIfree(MRI **pmri) ;
MRI   *MRIalloc(int width, int height, int depth, int type) ;
MRI   *MRIallocSequence(int width, int height,int depth,int type,int nframes);
MRI   *MRIallocHeader(int width, int height, int depth, int type) ;
int   MRIallocIndices(MRI *mri) ;
int   MRIsetResolution(MRI *mri, float xres, float yres, float zres) ;
int   MRIsetTransform(MRI *mri,   General_transform *transform) ;


/* correlation routines */
MRI   *MRIxcorr(MRI *mri_ref, MRI *mri_in, MRI *mri_dst) ;
MRI   *MRIxcorrWindow(MRI *mri_ref, MRI *mri_in,MRI *mri_dst,int window_size) ;
MRI   *MRInxcorr(MRI *mri_ref, MRI *mri_in, MRI *mri_dst) ;
MRI   *MRInxcorrWindow(MRI *mri_ref,MRI *mri_in,MRI *mri_dst,int window_size) ;
long  MRIcorrelate(MRI *mri_ref, MRI *mri_in, int xoff, int yoff, int zoff) ;


int   MRIpeak(MRI *mri, int *px, int *py, int *pz) ;
int   MRIcopyHeader(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcopy(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIreslice(MRI *mri_src, MRI *mri_dst, int slice_direction) ;
int   MRIboundingBox(MRI *mri, int thresh, MRI_REGION *region) ;

/* coordinate transforms */
MRI   *MRItranslate(MRI *mri_src, MRI *mri_dst, int dx, int dy, int dz) ;
MRI   *MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;
MRI   *MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz) ;
MRI   *MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB) ;
MRI   *MRIinterpolate(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIconfThresh(MRI *mri_src, MRI *mri_probs, MRI *mri_classes, 
                     MRI *mri_dst,float thresh, int min_target,int max_target);

/* debugging */
int   MRIdump(MRI *mri, FILE *fp) ;
int   MRIdumpBuffer(MRI *mri, FILE *fp) ;

/* arithmetic operations */
MRI   *MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabsdiff(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIaverage(MRI *mri_src, int dof, MRI *mri_dst) ;
MRI   *MRIdivide(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRImultiply(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabs(MRI *mri, MRI *mri_dst) ;
MRI   *MRIscalarMul(MRI *mri_src, MRI *mri_dst, float scalar) ;

/* filtering */
MRI   *MRIfindThinWMStrands(MRI *mri_src, MRI *mri_dst, int wsize);
MRI   *MRIcentralPlaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, 
                                            int wsize);
MRI   *MRIplaneOfLeastVarianceNormal(MRI *mri_src, MRI *mri_dst, int wsize) ;
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

MRI   *MRIsobel(MRI *mri_src, MRI *mri_grad, MRI *mri_mag);
MRI   *MRIxSobel(MRI *mri_src, MRI *mri_x, int frame) ;
MRI   *MRIySobel(MRI *mri_src, MRI *mri_y, int frame) ;
MRI   *MRIzSobel(MRI *mri_src, MRI *mri_z, int frame) ;
MRI   *MRIsobelRegion(MRI *mri_src, MRI *mri_grad, int domag, 
                      MRI_REGION *region);
MRI   *MRIxSobelRegion(MRI *mri_src, MRI *mri_x, int frame,MRI_REGION *region);
MRI   *MRIySobelRegion(MRI *mri_src, MRI *mri_y, int frame,MRI_REGION *region);
MRI   *MRIzSobelRegion(MRI *mri_src, MRI *mri_z, int frame,MRI_REGION *region);

MRI   *MRIreduce(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *kernel, 
                     int len, int axis) ;
MRI   *MRIreduce1d(MRI *mri_src, MRI *mri_dst,float *kernel,int len,int axis);
MRI   *MRIdiffuse(MRI *mri_src, MRI *mri_dst, double k, 
                    int niter, int which, double slope) ;
MRI   *MRIdiffuseCurvature(MRI *mri_src, MRI *mri_dst, 
                            double A,int niter, double slope) ;
MRI   *MRIdiffusePerona(MRI *mri_src, MRI *mri_dst, 
                             double k, int niter,double slope);
MRI   *MRIdirectionMap(MRI *mri_grad, MRI *mri_direction, int wsize);

/* offset stuff */
MRI   *MRIoffsetDirection(MRI *mri_grad, int wsize, MRI *mri_direction,
                          MRI *mri_dir);
MRI   *MRIoffsetMagnitude(MRI *mri_src, MRI *mri_dst, int maxsteps) ;
MRI   *MRIapplyOffset(MRI *mri_src, MRI *mri_dst, MRI *mri_offset) ;


MRI   *MRIclone(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIcloneRoi(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIthreshold(MRI *mri_src, MRI *mri_dst, BUFTYPE threshold) ;
MRI   *MRIthresholdRangeInto(MRI *mri_src, MRI *mri_dst, 
                             BUFTYPE low_val, BUFTYPE hi_val) ;
int   MRIprincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues,
                              int *means, BUFTYPE theshold) ;
int   MRIclear(MRI *mri_src) ;

/* these routines use trilinear interpolation */
MRI   *MRIrotateX_I(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY_I(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ_I(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate_I(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;

/* extraction routines */
MRI   *MRIextract(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                  int dx, int dy, int dz) ;
MRI   *MRIextractInto(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                  int dx, int dy, int dz, int x1, int y1, int z1) ;
MRI   *MRIextractIntoRegion(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                            MRI_REGION *region) ;

MRI   *MRIextractRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *region) ;
MRI   *MRIextractPlane(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, int wsize,
                        int x, int y, int z);
MRI   *MRIextractCpolv(MRI *mri_src, MRI *mri_dst, MRI *mri_polv, 
                           int wsize, int x, int y, int z);
MRI   *MRIextractCpolvCoords(MRI *mri_src, int *px, int *py, int *pz, 
                                 MRI *mri_polv, int x, int y,int z,int wsize);
MRI   *MRIextractValues(MRI *mri_src, MRI *mri_dst, float min_val, 
                        float max_val) ;
MRI   *MRIwmfilter(MRI *mri_src, MRI *mri_cpolv, MRI *mri_dst) ;
MRI   *MRIorder(MRI *mri_src, MRI *mri_dst, int wsize, float pct) ;
MRI   *MRIremoveHoles(MRI *mri_src, MRI*mri_dst, int wsize, float pct) ;

/* morphology */
MRI   *MRImorph(MRI *mri_src, MRI *mri_dst, int which) ;
MRI   *MRIerode(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIerodeRegion(MRI *mri_src, MRI *mri_dst,int wsize,MRI_REGION *region);
MRI   *MRIdilate(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIopen(MRI *mri_src, MRI *mri_dst) ;
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
MRI   *MRIcomputeResidual(MRI *mri1, MRI *mri2, MRI *mri_dst, int t1, int t2) ;

/* filtering operations */
MRI   *MRIminmax(MRI *mri_src, MRI *mri_dst, MRI *mri_dir, int wsize) ;
MRI   *MRIgaussian1d(float sigma, int max_len) ;
MRI   *MRIconvolveGaussian(MRI *mri_src, MRI *mri_dst, MRI *mri_gaussian) ;
MRI   *MRImedian(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRImean(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRImeanByte(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRIstd(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize) ;
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
MRI   *MRIreadRaw(FILE *fp, int width, int height, int depth, int type) ;
int   MRIinitHeader(MRI *mri) ;
int   MRIvoxelToWorld(MRI *mri, Real xv, Real yv, Real zv, 
                      Real *xw, Real *yw, Real *zw) ;
int   MRIworldToVoxel(MRI *mri, Real xw, Real yw, Real zw,
                Real *pxv, Real *pyv, Real *pzv) ;
int   MRIworldToTalairachVoxel(MRI *mri, Real xw, Real yw, Real zw,
                Real *pxv, Real *pyv, Real *pzv) ;
int   MRIvoxelToTalairachVoxel(MRI *mri, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt) ;
int   MRIvoxelToVoxel(MRI *mri_src, MRI *mri_dst, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt) ;
int   MRItalairachVoxelToVoxel(MRI *mri, Real xt, Real yt, Real zt,
                               Real *pxv, Real *pyv, Real *pzv) ;
int   MRItalairachVoxelToWorld(MRI *mri, Real xt, Real yt, Real zt,
                               Real *pxw, Real *pyw, Real *pzw) ;
int   MRIvoxelToTalairach(MRI *mri, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt) ;
int   MRItalairachToVoxel(MRI *mri, Real xt, Real yt, Real zt,
                               Real *pxv, Real *pyv, Real *pzv) ;

int   MRItransformRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *src_region,
                                 MRI_REGION *dst_region) ;
MRI   *MRIextractTalairachPlane(MRI *mri_src, MRI *mri_dst, int orientation, 
                                int x, int y, int z, int size) ;
int   MRIeraseTalairachPlane(MRI *mri, MRI *mri_mask,
                              int orientation, int x, int y, int z, int size) ;

/* resampling routines */
MRI   *MRIupsample2(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIdownsample2(MRI *mri_src, MRI *mri_dst) ;


#include "image.h"

IMAGE *MRItoImage(MRI *mri, IMAGE *I, int slice) ;
IMAGE *MRItoImageView(MRI *mri, IMAGE *I, int slice, int view, int frame) ;



#define MRISvox(mri,x,y,z)  (((short *)mri->slices[z][y])[x])
#define MRIFvox(mri,x,y,z)  (((float *)(mri->slices[z][y]))[x])
#define MRIvox(mri,x,y,z)   (((BUFTYPE *)mri->slices[z][y])[x])
#define MRISCvox(mri,x,y,z) (((signed char *)mri->slices[z][y])[x])
#define MRIIvox(mri,x,y,z)  (((int *)mri->slices[z][y])[x])
#define MRILvox(mri,x,y,z)  (((long *)mri->slices[z][y])[x])

#define MRISseq_vox(mri,x,y,z,n)  (((short*)mri->slices[z+(n)*mri->depth][y])\
                                   [x])
#define MRISCseq_vox(mri,x,y,z,n)  (((signed char*)mri->slices[z+(n)*mri->depth][y])\
                                   [x])
#define MRIFseq_vox(mri,x,y,z,n)  (((float*)(mri->slices[z+((n)*mri->depth)]\
                                                        [y]))[x])
#define MRIseq_vox(mri,x,y,z,n)   (((BUFTYPE *)\
                                    mri->slices[z+(n)*mri->depth][y])[x])
#define MRIIseq_vox(mri,x,y,z,n)  (((int *)mri->slices[z+(n)*mri->depth][y])\
                                   [x])
#define MRILseq_vox(mri,x,y,z,n)  (((long *)mri->slices[z+(n)*mri->depth][y])\
                                   [x])

#define MRI_HEIGHT      0
#define MRI_WIDTH       1
#define MRI_DEPTH       2

#define MRI_CORONAL     0
#define MRI_SAGITTAL    1
#define MRI_HORIZONTAL  2
#define MRI_AXIAL       MRI_HORIZONTAL

/* vertices of an icosahedron (sp?), used by all POLV functions */
#define NVERTICES    22
extern float ic_x_vertices[]  ;
extern float ic_y_vertices[]  ;
extern float ic_z_vertices[]  ;


#include "histo.h"

#define MRI_CORONAL_SLICE_DIRECTORY   0
#define MRI_MINC_FILE                 1
#define MRI_ANALYZE_FILE              2

int        MRImatch(MRI *mri1, MRI *mri2) ;
int        MRIvalRange(MRI *mri, float *pmin, float *pmax) ;
MRI        *MRIvalScale(MRI *mri_src, MRI *mri_dst, float fmin, float fmax) ;
HISTOGRAM  *MRIhistogram(MRI *mri, int nbins) ;
MRI        *MRIhistoEqualize(MRI *mri_src, MRI *mri_dst, int low) ;
MRI        *MRIapplyHistogram(MRI *mri_src, MRI *mri_dst, HISTOGRAM *histo) ;
MRI        *MRIcrunch(MRI *mri_src, MRI *mri_dst) ;
HISTOGRAM  *MRIgetEqualizeHisto(MRI *mri, HISTOGRAM *histo_eq, int low, 
                                int norm) ;

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
int        MRIunpackFileName(char *inFname, int *pframe, int *ptype, 
                             char *outFname) ;
Volume     MRItoVolume(MRI *mri) ;
MRI        *MRIfromVolume(Volume volume, int start_frame, int end_frame) ;
int        MRIisValid(MRI *mri) ;
MRI        *MRIflipByteOrder(MRI *mri_src, MRI *mri_dst) ;

MRI        *MRIregionGrow(MRI *mri_src, MRI *mri_distance, 
                    float x0, float y0, float z0, int niter) ;
MRI        *MRIextractInterior(MRI *mri_src, MRI *mri_distance,  MRI *mri_dst);
MRI        *MRIbuildDistanceMap(MRI *mri_src, MRI *mri_distance,
                                float x0, float y0, float z0, float r) ;
MRI        *MRIupdateDistanceMap(MRI *mri_distance) ;
MRI        *MRIfill(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
                    int seed_z, int threshold, int fill_val) ;
MRI        *MRIfillFG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
                    int seed_z, int threshold, int fill_val, int *npix) ;
                    

/* constants used in mri_dir of MRIoffsetDirection and for MRIminmax filter */
#define OFFSET_NEGATIVE_GRADIENT_DIRECTION    0
#define OFFSET_GRADIENT_DIRECTION             1
#define OFFSET_ZERO                           2

#endif
