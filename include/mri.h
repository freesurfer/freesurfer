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
#define MRI_BITMAP  5

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
  int           xdir ;      /* these are the actual directions of the structure axes, */
  int           ydir ;      /* set when the volume is read, to be used to orient the  */
  int           zdir ;      /* structure (if needed) with                             */
                            /* MRIreorder(src, dest, src->xdir, src->ydir, src->zdir) */

/* 
  each slice is an array of rows (mri->height of them) each of which is 
  mri->width long.
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
  int           dof ;
  double        mean ;   
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
int    MRIappend(MRI *mri, char *fpref) ;
int    MRIwriteInfo(MRI *mri, char *fpref) ;
MRI   *MRIread(char *fpref) ;
MRI   *MRIreadInfo(char *fpref) ;

/* memory allocation routines */
int   MRIfree(MRI **pmri) ;
int   MRIfreeFrames(MRI *mri, int start_frame) ;
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
MRI   *MRItranslate(MRI *mri_src, MRI *mri_dst, double dx, double dy, double dz) ;
MRI   *MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;
MRI   *MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz) ;
MRI   *MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB) ;
MRI   *MRIinverseLinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA) ;
MRI   *MRIlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA) ;
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
MRI   *MRIremoveIslands(MRI *mri_src, MRI*mri_dst, int wsize, int thresh) ;
MRI   *MRIresegmentThinWMStrands(MRI *mri_src, MRI *mri_dst, int thickness);
MRI   *MRIthickenThinWMStrands(MRI *mri_src, MRI *mri_dst, int thickness,
                               int nvoxels) ;
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
MRI   *MRIreduceByte(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIconvolve1dByte(MRI *mri_src, MRI *mri_dst, float *k, int len, 
                         int axis, int src_frame, int dst_frame) ;
MRI   *MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *kernel, 
                     int len, int axis, int src_frame, int dst_frame) ;
MRI   *MRIreduce1d(MRI *mri_src, MRI *mri_dst,float *kernel,int len,int axis);
MRI   *MRIreduce1dByte(MRI *mri_src, MRI *mri_dst,float *kernel,int len,
                       int axis);
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
MRI   *MRIbinarize(MRI *mri_src, MRI *mri_dst, BUFTYPE threshold,
                   BUFTYPE low_val, BUFTYPE hi_val) ;
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
MRI   *MRIconvolveGaussianMeanAndStdByte(MRI *mri_src, MRI *mri_dst,
                                          MRI *mri_gaussian) ;
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
int   MRIeraseTalairachPlane(MRI *mri, MRI *mri_mask, int orientation, 
                             int x, int y, int z,int size,int fill_val);
int   MRIeraseTalairachPlaneNew(MRI *mri, MRI *mri_mask, int orientation, 
                             int x, int y, int z,int size,int fill_val);

MRI   *MRIextractPlane(MRI *mri_src, MRI *mri_dst, int orientation,int where);
int   MRIerasePlane(MRI *mri, float x0, float y0, float z0,
                    float dx, float dy, float dz, int fill_val);

int   MRIeraseBorders(MRI *mri, int width) ;
int   MRIsampleVolume(MRI *mri, Real x, Real y, Real z, Real *pval) ;
int   MRIsampleVolumeFrame(MRI *mri,Real x,Real y,Real z,int frame,Real *pval);
int   MRIsampleVolumeGradient(MRI *mri, Real x, Real y, Real z, 
                              Real *pdx, Real *pdy, Real *pdz) ;
int   MRIsampleVolumeDerivative(MRI *mri, Real x, Real y, Real z,
                                Real dx, Real dy, Real dz, Real *pmag) ;
float MRIsampleCardinalDerivative(MRI *mri, int x, int y, int z,
                                  int xk, int yk, int zk) ;
float MRIsampleXDerivative(MRI *mri, int x, int y, int z, int dir) ;
float MRIsampleYDerivative(MRI *mri, int x, int y, int z, int dir) ;
float MRIsampleZDerivative(MRI *mri, int x, int y, int z, int dir) ;

/* resampling routines */
MRI   *MRIupsample2(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIdownsample2(MRI *mri_src, MRI *mri_dst) ;



#include "image.h"

IMAGE *MRItoImage(MRI *mri, IMAGE *I, int slice) ;
IMAGE *MRItoImageView(MRI *mri, IMAGE *I, int slice, int view, int frame) ;

/* bitmap image access macros */
#define MRIset_bit(mri,x,y,z)    MRIvox(mri,(x)/8,y,z) |= (0x001 << ((x)%8))
#define MRItest_bit(mri,x,y,z)   (MRIvox(mri,(x)/8,(y),(z))&(0x001 << ((x)%8)))
#define MRIclear_bit(mri,x,y,z)  MRIvox(mri,(x)/8,y,z) &= ~(0x001 << ((x)%8))

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


#include "histo.h"

#define MRI_CORONAL_SLICE_DIRECTORY   0
#define MRI_MINC_FILE                 1
#define MRI_ANALYZE_FILE              2
#define MRI_MGH_FILE                  3
#define GENESIS_FILE                  4
#define GE_LX_FILE                    5
#define SIEMENS_FILE                  6
#define BRIK_FILE                     7
#define BSHORT_FILE                   8
#define SDT_FILE                      9

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
                    int seed_z, int threshold, int fill_val, int max_count) ;
MRI        *MRIfillFG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y, 
                    int seed_z, int threshold, int fill_val, int *npix) ;
                    

int   MRIneighborsOn(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOff(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOn3x3(MRI *mri, int x0, int y0, int z0, int min_val) ;
int   MRIneighborsOff3x3(MRI *mri, int x0, int y0, int z0, int min_val) ;

MRI   *MRIreplaceValues(MRI *mri_src, MRI *mri_dst, 
                       BUFTYPE in_val, BUFTYPE out_val) ;
MRI   *MRImask(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, BUFTYPE mask,
               BUFTYPE out_val) ;

/* constants used in mri_dir of MRIoffsetDirection and for MRIminmax filter */
#define OFFSET_NEGATIVE_GRADIENT_DIRECTION    0
#define OFFSET_GRADIENT_DIRECTION             1
#define OFFSET_ZERO                           2

/* anything below this is not white matter */
#define WM_MIN_VAL                       2 
#define WM_EDITED_ON_VAL                 255

MRI *MRIreduceMeanAndStdByte(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIstdsToVariances(MRI *mri_std, MRI *mri_var, int source_frame) ;
MRI *MRIvariancesToStds(MRI *mri_var, MRI *mri_std, int dst_frame) ;
MRI *MRIconcatenateFrames(MRI *mri_frame1, MRI *mri_frame2, MRI *mri_dst);
MRI *MRIcopyFrame(MRI *mri_src, MRI *mri_dst, int src_frame, int dst_frame) ;
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
int  MRIneighborhoodBlackCpolv(MRI *mri_src, int xv, int yv, int zv, 
                               int nsize, int wsize, int *pnum_black) ;

MRI *MRIorderSegment(MRI *mri_src, MRI *mri_labeled, float thresh, int wsize);
MRI *MRIthresholdLabel(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, 
                       int wm_low) ;
MRI *MRIintensitySegmentation(MRI *mri_src, MRI *mri_labeled,
                              float wm_low, float wm_hi, float gray_hi) ;
MRI *MRImeanLabel(MRI *mri_src, MRI *mri_label, MRI*mri_dst, int wsize) ;
MRI *MRIcpolvVote(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst, int wsize, 
                  int niter, int use_all) ;
MRI *MRIcpolvThreshold(MRI *mri_src, MRI *mri_labeled, MRI *mri_dst,
                       int wm_low, int gray_hi,int wsize) ;
MRI *MRImaskLabels(MRI *mri_src, MRI *mri_mask, MRI *mri_dst) ;


MRI *MRIwmfilterMarked(MRI *mri_src, MRI *mri_mask, MRI *mri_dst, int wsize, 
                       float pct, int onoff) ;
int  MRIcountCpolvAtVoxel(MRI *mri_src, int x, int y, int z, int wsize, 
                          int *pnum, int label_to_check) ;
int  MRIcountCpolvOnAtVoxel(MRI *mri_src, int x, int y, int z, int wsize, 
                            int *pnum) ;
MRI *MRIcentralPlaneOfLeastVarianceNormalMarked(MRI *mri_src, MRI *mri_mask,
                                                MRI *mri_dst, int wsize) ;
int  MRIcountCpolvOffAtVoxel(MRI *mri_src,int x, int y, int z, int wsize, 
                             int *pnum) ;
int  MRIcentralPlaneOfLeastVarianceNormalVoxel(MRI *mri_src, int wsize,
                                               int x, int y, int z) ;
int  MRIvertexToVector(int vertex, float *pdx, float *pdy, float *pdz) ;
int MRIcpolvMedianCurveVoxel(MRI *mri, MRI *mri_labeled, int x0, int y0, 
                             int z0, int wsize, float len) ;
float MRIcpolvMedianAtVoxel(MRI *mri_src, int vertex, 
                             float x, float y, float z, int wsize);
MRI   *MRIcpolvMedianCurveSegment(MRI *mri,MRI *mri_labeled, MRI *mri_dst,
                                int wsize,float len);

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
int MRIsetValues(MRI *mri, int val) ;

#define MRI_NOT_WHITE   1
#define MRI_AMBIGUOUS   128
#define MRI_WHITE       255

#define MRI_LEFT_HEMISPHERE     255
#define MRI_RIGHT_HEMISPHERE    127

/* STATS volumes have 2 images each */
#define TRI_HI_PRIORS            0
#define TRI_HI_STATS             1
#define TRI_LOW_PRIORS           3
#define TRI_LOW_STATS            4
#define TRI_OFF_STATS            6

#endif
