/**
 * @brief structures and prototypes for discrete cosine transform
 *
 * 
 *  structures and prototypes for discrete cosine transform
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


//
// dct.h
//
// written by Bruce Fischl
// Nov. 9th ,2000
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "transform.h"
#include "matrix.h"
#include "voxlist.h"

#ifndef DCT_H
#define DCT_H

typedef struct
{
  int      ncoef ;
  VOL_GEOM *VG ;
  LTA      *lta ;
  MRI      *mri_source ;
  MATRIX   *m_x_basis ;
  MATRIX   *m_y_basis ;
  MATRIX   *m_z_basis ;
  VECTOR   *v_xk ;   // vector of x coefs
  VECTOR   *v_yk ;   // vector of y coefs
  VECTOR   *v_zk ;   // vector of z coefs
  double   *x ;      // forward mapping of x
  double   *y ;      // forward mapping of y
  double   *z ;      // forward mapping of z
  double   *x_inv ;  // inverse of x indices
  double   *y_inv ;  // inverse of y indices
  double   *z_inv ;  // inverse of z indices
  double   res ;     // resolution of inverse transform tables
  double   b ;       // intensity scale factor
} DISCRETE_COSING_TRANSFORM, DCT ;

DCT   *DCTalloc(int ncoef, MRI *mri_source) ;
int   DCTfree(DCT **pdct) ;
int   DCTcreateMatrix(DCT *dct, MRI *mri, int skip) ;
DCT   *DCTcopy(DCT *dct_src, DCT *dct_dst) ;
MRI   *DCTapply(DCT *dct, MRI *mri_src, MRI *mri_template, MRI *mri_dst, int sample_type) ;
MRI   *DCTapplyInverse(DCT *dct, MRI *mri_src, MRI *mri_dst, int sample_type) ;
int   DCTtransformVoxlist(DCT *dct, VOXEL_LIST *vl) ;
int   DCTinverseTransformVoxlist(DCT *dct, VOXEL_LIST *vl) ;
int   DCTtransformPoint(DCT *dct, int x, int y, int z, double *px, double *py,  double *pz) ;
int   DCTdump(DCT *dct, FILE *fp) ;
int   DCTupdate(DCT *dct) ;
int   DCTinverseTransformPoint(DCT *dct, double x, double y, double z, double *px, double *py,  double *pz) ;

#endif
