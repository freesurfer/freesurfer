#ifndef MRI_H
#define MRI_H

#include "const.h"
#include "matrix.h"
#include "volume_io.h"

#define BUFTYPE  unsigned char

#define MRI_UCHAR   0
#define MRI_INT     1
#define MRI_LONG    2
#define MRI_FLOAT   3

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
  int           xnum ;
  int           ynum ;
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
} MRI_IMAGE, MRI ;

/* I/O functions */
int    MRIwrite(MRI *mri, char *fpref) ;
int    MRIwriteInfo(MRI *mri, char *fpref) ;
MRI   *MRIread(char *fpref) ;
MRI   *MRIreadInfo(char *fpref) ;

/* memory allocation routines */
int   MRIfree(MRI **pmri) ;
MRI   *MRIalloc(int width, int height, int depth, int type) ;


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

/* coordinate transforms */
MRI   *MRItranslate(MRI *mri_src, MRI *mri_dst, int dx, int dy, int dz) ;
MRI   *MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle) ;
MRI   *MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle) ;
MRI   *MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle) ;
MRI   *MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO) ;
MRI   *MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz) ;
MRI   *MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB) ;

/* debugging */
int   MRIdump(MRI *mri, FILE *fp) ;
int   MRIdumpBuffer(MRI *mri, FILE *fp) ;

/* arithmetic operations */
MRI   *MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabsdiff(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIadd(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIdivide(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRImultiply(MRI *mri1, MRI *mri2, MRI *mri_dst) ;
MRI   *MRIabs(MRI *mri, MRI *mri_dst) ;

/* filtering */
MRI   *MRIreduce(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIconvolve1d(MRI *mri_src, MRI *mri_dst, float *kernel, 
                     int len, int axis) ;
MRI   *MRIreduce1d(MRI *mri_src, MRI *mri_dst, float *kernel, 

                     int len, int axis) ;
MRI   *MRIclone(MRI *mri_src, MRI *mri_dst) ;
MRI   *MRIthreshold(MRI *mri_src, MRI *mri_dst, BUFTYPE threshold) ;
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


MRI   *MRImean(MRI *mri_src, MRI *mri_dst, int wsize) ;
MRI   *MRIstd(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, int wsize) ;
MRI   *MRInorm(MRI *mri_src, MRI*mri_dst, MRI *mri_mean, MRI *mri_std) ;
MRI   *MRIreadRaw(FILE *fp, int width, int height, int depth, int type) ;
int   MRIinitHeader(MRI *mri) ;
int   MRIvoxelToWorld(MRI *mri, Real xv, Real yv, Real zv, 
                      Real *xw, Real *yw, Real *zw) ;
int   MRIworldToVoxel(MRI *mri, Real xw, Real yw, Real zw,
                Real *pxv, Real *pyv, Real *pzv) ;
int   MRIvoxelToTalairachVoxel(MRI *mri, Real xv, Real yv, Real zv,
                               Real *pxt, Real *pyt, Real *pzt) ;
int   MRItalairachVoxelToVoxel(MRI *mri, Real xt, Real yt, Real zt,
                               Real *pxv, Real *pyv, Real *pzv) ;


#include "image.h"

IMAGE *MRItoImage(MRI *mri, IMAGE *I, int slice) ;
IMAGE *MRItoImageView(MRI *mri, IMAGE *I, int slice, int view) ;



#define MRIFvox(mri,x,y,z)  (((float *)mri->slices[z][y])[x])
#define MRIvox(mri,x,y,z)   (((BUFTYPE *)mri->slices[z][y])[x])
#define MRIIvox(mri,x,y,z)  (((int *)mri->slices[z][y])[x])
#define MRILvox(mri,x,y,z)  (((long *)mri->slices[z][y])[x])

#define MRI_HEIGHT      0
#define MRI_WIDTH       1
#define MRI_DEPTH       2

#define MRI_CORONAL     0
#define MRI_SAGITAL     1
#define MRI_HORIZONTAL  2
#define MRI_AXIAL       MRI_HORIZONTAL

#endif
