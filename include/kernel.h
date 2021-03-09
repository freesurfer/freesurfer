/*
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


#ifndef KERNEL_H
#define KERNEL_H

/*
  adding row0, col0 to a point in kernel space maps it into a point
  int image space, conversely, subtracting row0,col0 from a point in
  image space yields a point in kernel space.
*/

typedef struct
{
  int    row0, col0 ;      /* top-left hand point of this kernel */
  int    row,col ;         /* location of this pixel */
  int    rows, cols ;      /* size of this kernel */
  float  **weights ;       /* array of rows */
}
KERNEL ;

typedef struct
{
  int    rows, cols ;      /* size of the image */
  int    krows, kcols ;    /* size of the kernels */
  KERNEL kernels[1] ;      /* actually allocated to be width x height */
  char   *fname ;          /* name of image that kernels were derived from */
}
KERNEL_IMAGE, KIMAGE ;

/* include this after definition of KIMAGE as image.h uses it */
#include "image.h"

extern KIMAGE *KernelImageAlloc(int rows, int cols,int krows,int kcols);
extern KIMAGE *KernelImageClone(KIMAGE *kimage) ;
extern void   KernelImageCopy(KIMAGE *ksrc, KIMAGE *kdst) ;
extern void   KernelImageFree(KIMAGE *kimage) ;
extern void   KernelInit(KIMAGE *kimage, int row, int col) ;
extern void   KernelFree(KERNEL *kernel) ;
extern void   KernelImageDump(KIMAGE *kimage, FILE *fp) ;
extern void   KernelDiscount(KIMAGE *kimage, int row, int col, float weight) ;
extern void   KernelUpdate(KIMAGE *ksrc, KIMAGE *kdst, int dst_row,
                             int dst_col, int src_row,int src_col,float weight);
extern void   KernelImageNormalize(KIMAGE *kimage) ;
extern void   KernelNormalize(KIMAGE *kimage, int row, int col) ;
extern void   KernelCopy(KIMAGE *ksrc, KIMAGE *kdst, int src_row, int src_col,
                           int dst_row, int dst_col) ;
extern void   KernelImageConvolve(KIMAGE *kimage, IMAGE *src_image,
                                    IMAGE *dst_image) ;
extern void   KernelImageWrite(KIMAGE *kimage, char *fname, int arc,
                                 char *argv[]) ;
extern KIMAGE    *KernelImageRead(char *fname) ;
extern IMAGE *KernelImageToSeq(KIMAGE *kimage) ;
extern KIMAGE    *KernelImageFromSeq(IMAGE *image) ;

#define KIMAGEpix(k,row,col)  (k->kernels + row * k->cols + col)

#endif
