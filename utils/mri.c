/**
 * @File  mri.c
 * @brief utilities for MRI data structure
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/01/08 19:49:48 $
 *    $Revision: 1.520 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

extern const char* Progname;
const char *MRI_C_VERSION = "$Revision: 1.520 $";


/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#define USE_ELECTRIC_FENCE  1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <errno.h>
#include <ctype.h>
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "minc_volume_io.h"
#include "filter.h"
#include "box.h"
#include "region.h"
#include "mri_transform.h"
#include "utils.h"
#include "matrix.h"
#include "pdf.h"
#include "cma.h"
#include "talairachex.h"
#include "voxlist.h"
#include "fastmarching.h"
#include "mriBSpline.h"

extern int errno;

/*-----------------------------------------------------
  MACROS AND CONSTANTS
  -------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x==8&&y==9) || (x==9&&y==8)) &&((z)==15))


#ifndef UCHAR_MIN
#define UCHAR_MIN  0.0
#endif
#ifndef UCHAR_MAX
#define UCHAR_MAX  255.0
#endif
#ifndef SHORT_MIN
#define SHORT_MIN  -32768.0
#endif
#ifndef SHORT_MAX
#define SHORT_MAX  32767.0
#endif
#ifndef INT_MIN
#define INT_MIN    -2147483648.0
#endif
#ifndef INT_MAX
#define INT_MAX    2147483647.0
#endif
#ifndef LONG_MIN
#define LONG_MIN   -2147483648.0
#endif
#ifndef LONG_MAX
#define LONG_MAX   2147483647.0
#endif

#define N_HIST_BINS 1000



/*-----------------------------------------------------
  STATIC DATA
  -------------------------------------------------------*/

static long mris_alloced = 0 ;

/*-----------------------------------------------------
  STATIC PROTOTYPES
  -------------------------------------------------------*/

/*-----------------------------------------------------
  GLOBAL FUNCTIONS
  -------------------------------------------------------*/

/*----------------------------------------------------------
  MRIxfmCRS2XYZ() - computes the matrix needed to compute the
  XYZ of the center of a voxel at a given Col, Row, and Slice
  from the native geometry of the volume (ie, native or scanner
  Vox2RAS matrix).

  x         col
  y  = T *  row
  z        slice
  1          1

  T = [Mdc*D Pxyz0]
  [0 0 0   1  ]

  Mdc = [Vcol Vrow Vslice]
  V<dim> = the direction cosine pointing from the center of one voxel
  to the center of an adjacent voxel in the next dim, where
  dim is either colum, row, or slice. Vcol = [x_r x_a x_s],
  Vrow = [y_r y_a y_s], Vslice = [z_r z_a z_s]. Vcol can also
  be described as the vector normal to the plane formed by
  the rows and slices of a given column (ie, the column normal).

  D = diag([colres rowres sliceres])
  dimres = the distance between adjacent dim, where colres = mri->xsize,
  rowres = mri->ysize, and sliceres = mri->zsize.

  Pxyz0 = the XYZ location at CRS=0. This number is not part of the
  mri structure, so it is computed here according to the formula:
  Pxyz0 = PxyzCenter - Mdc*D*PcrsCenter

  PcrsCenter = the col, row, and slice at the center of the volume,
  = [ ncols/2 nrows/2 nslices/2 ]

  PxyzCenter = the X, Y, and Z at the center of the volume and does
  exist in the header as mri->c_r, mri->c_a, and mri->c_s,
  respectively.

  Note: coordinates are at the center of the voxel.

  Note: to compute the matrix with respect to the first voxel being
  at CRS 1,1,1 instead of 0,0,0, then set base = 1. This is
  necessary with SPM matrices.

  See also: MRIxfmCRS2XYZtkreg, MRItkReg2Native, extract_i_to_r().
  surfaceRASFromVoxel_(MRI *mri), voxelFromSurfaceRAS_().

  Note: MRIgetVoxelToRasXform is #defined to be extract_i_to_r().
  ----------------------------------------------------------------*/
MATRIX *MRIxfmCRS2XYZ( const MRI *mri, int base)
{
  MATRIX *m;
  MATRIX *Pcrs, *PxyzOffset;

  m = MatrixAlloc(4, 4, MATRIX_REAL);

  /* direction cosine between columns scaled by
     distance between colums */
  *MATRIX_RELT(m, 1, 1) = mri->x_r * mri->xsize;
  *MATRIX_RELT(m, 2, 1) = mri->x_a * mri->xsize;
  *MATRIX_RELT(m, 3, 1) = mri->x_s * mri->xsize;

  /* direction cosine between rows scaled by
     distance between rows */
  *MATRIX_RELT(m, 1, 2) = mri->y_r * mri->ysize;
  *MATRIX_RELT(m, 2, 2) = mri->y_a * mri->ysize;
  *MATRIX_RELT(m, 3, 2) = mri->y_s * mri->ysize;

  /* direction cosine between slices scaled by
     distance between slices */
  *MATRIX_RELT(m, 1, 3) = mri->z_r * mri->zsize;
  *MATRIX_RELT(m, 2, 3) = mri->z_a * mri->zsize;
  *MATRIX_RELT(m, 3, 3) = mri->z_s * mri->zsize;

  /* Preset the offsets to 0 */
  *MATRIX_RELT(m, 1, 4) = 0.0;
  *MATRIX_RELT(m, 2, 4) = 0.0;
  *MATRIX_RELT(m, 3, 4) = 0.0;

  /* Last row of matrix */
  *MATRIX_RELT(m, 4, 1) = 0.0;
  *MATRIX_RELT(m, 4, 2) = 0.0;
  *MATRIX_RELT(m, 4, 3) = 0.0;
  *MATRIX_RELT(m, 4, 4) = 1.0;

  /* At this point, m = Mdc * D */

  /* Col, Row, Slice at the Center of the Volume */
  Pcrs = MatrixAlloc(4, 1, MATRIX_REAL);
  *MATRIX_RELT(Pcrs, 1, 1) = mri->width/2.0  + base;
  *MATRIX_RELT(Pcrs, 2, 1) = mri->height/2.0 + base;
  *MATRIX_RELT(Pcrs, 3, 1) = mri->depth/2.0  + base;
  *MATRIX_RELT(Pcrs, 4, 1) = 1.0;

  /* XYZ offset the first Col, Row, and Slice from Center */
  /* PxyzOffset = Mdc*D*PcrsCenter */
  PxyzOffset = MatrixMultiply(m,Pcrs,NULL);

  /* XYZ at the Center of the Volume is mri->c_r, c_a, c_s  */

  /* The location of the center of the voxel at CRS = (0,0,0)*/
  *MATRIX_RELT(m, 1, 4) = mri->c_r - PxyzOffset->rptr[1][1];
  *MATRIX_RELT(m, 2, 4) = mri->c_a - PxyzOffset->rptr[2][1];
  *MATRIX_RELT(m, 3, 4) = mri->c_s - PxyzOffset->rptr[3][1];

  MatrixFree(&Pcrs);
  MatrixFree(&PxyzOffset);

  return(m);
}
/*!
  \fn MATRIX *MRImatrixOfDirectionCosines(MRI *mri, MATRIX *Mdc)
  \brief Fills Mdc with direction cosines
*/
MATRIX *MRImatrixOfDirectionCosines(MRI *mri, MATRIX *Mdc)
{
  if(Mdc == NULL) Mdc = MatrixZero(4,4,NULL);
  Mdc->rptr[1][1] = mri->x_r;
  Mdc->rptr[2][1] = mri->x_a;
  Mdc->rptr[3][1] = mri->x_s;
  Mdc->rptr[1][2] = mri->y_r;
  Mdc->rptr[2][2] = mri->y_a;
  Mdc->rptr[3][2] = mri->y_s;
  Mdc->rptr[1][3] = mri->z_r;
  Mdc->rptr[2][3] = mri->z_a;
  Mdc->rptr[3][3] = mri->z_s;
  Mdc->rptr[4][4] = 1;
  return(Mdc);
}
/*!
  \fn MATRIX *MRImatrixOfVoxelSizes(MRI *mri, MATRIX *D)
  \brief Creaetes diagonal matrix D with voxel sizes on diagonal
*/
MATRIX *MRImatrixOfVoxelSizes(MRI *mri, MATRIX *D)
{
  if(D == NULL) D = MatrixZero(4,4,NULL);
  D->rptr[1][1] = mri->xsize;
  D->rptr[2][2] = mri->ysize;
  D->rptr[3][3] = mri->zsize;
  D->rptr[4][4] = 1;
  return(D);
}
/*!
  \fn MATRIX *MRImatrixOfTranslations(MRI *mri, MATRIX *P0)
  \brief Creates 4x4 matrix which implements the translation
*/
MATRIX *MRImatrixOfTranslations(MRI *mri, MATRIX *P0)
{
  MATRIX *Vox2ScannerRAS;
  int k;
  if(P0 == NULL) P0 = MatrixZero(4,4,NULL);
  Vox2ScannerRAS = MRIxfmCRS2XYZ(mri,0);
  for(k=1; k<=4; k++) {
    P0->rptr[k][4] = Vox2ScannerRAS->rptr[k][4];
    P0->rptr[k][k] = 1;
  }
  MatrixFree(&Vox2ScannerRAS);
  return(P0);
}
/*--------------------------------------------------------------------------
  extract_i_to_r() - computes scanner vox2ras. On 2/27/06, this was replaced
  with a simple call to MRIxfmCRS2XYZ(). The original code is below (but
  removed with #defines). MRIxfmCRS2XYZ() is more general in that it handles
  non-zero voxels base correctly (eg, SPM expects vox2ras to be 1-based).
  Note: MRIgetVoxelToRasXform is #defined to be extract_i_to_r().
  ---------------------------------------------------------------------------*/
MATRIX *extract_i_to_r(const MRI *mri)
{
  MATRIX *m;
  m = MRIxfmCRS2XYZ(mri, 0);
  return(m);
}
/*---------------------------------------------------------------------
  extract_r_to_i() - computes scanner ras2vox. See also extract_i_to_r()
  and MRIxfmCRS2XYZ()
  ---------------------------------------------------------------------*/
MATRIX *extract_r_to_i(const MRI *mri)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras ;
  m_voxel_to_ras = extract_i_to_r(mri) ;
  m_ras_to_voxel = MatrixInverse(m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ;
  return(m_ras_to_voxel) ;
}

/*!
  \fn int MRIsetVox2RASFromMatrix(MRI *mri, MATRIX *m_vox2ras)
  \brief Takes a vox2ras matrix and assigns the MRI structure
  geometry fields such that it will realize this matrix. WARNING:
  this matrix can only be 9 DOF. It cannot have shear because
  shear is not represented in the MRI struct. See also 
  MRIsetVox2RASFromMatrixUnitTest(). See also niftiSformToMri()
  int mriio.c.
*/
int MRIsetVox2RASFromMatrix(MRI *mri, MATRIX *m_vox2ras)
{
  double rx,ax,sx,ry,ay,sy,rz,az,sz;
  double P0r, P0a, P0s;
  double xsize, ysize, zsize;

  rx = m_vox2ras->rptr[1][1];  ry = m_vox2ras->rptr[1][2];  rz = m_vox2ras->rptr[1][3];
  ax = m_vox2ras->rptr[2][1];  ay = m_vox2ras->rptr[2][2];  az = m_vox2ras->rptr[2][3];
  sx = m_vox2ras->rptr[3][1];  sy = m_vox2ras->rptr[3][2];  sz = m_vox2ras->rptr[3][3];
  P0r = m_vox2ras->rptr[1][4];  
  P0a = m_vox2ras->rptr[2][4];  
  P0s = m_vox2ras->rptr[3][4];

  xsize = sqrt(rx*rx + ax*ax + sx*sx);
  ysize = sqrt(ry*ry + ay*ay + sy*sy);
  zsize = sqrt(rz*rz + az*az + sz*sz);
  if(fabs(xsize-mri->xsize) > .001 ||
     fabs(ysize-mri->ysize) > .001 || 
     fabs(zsize-mri->zsize) > .001){
    printf("WARNING: MRIsetRas2VoxFromMatrix(): voxels sizes are inconsistent\n");
    printf("   (%g,%g) (%g,%g) (%g,%g) \n",mri->xsize,xsize, mri->ysize,ysize,
	   mri->zsize,zsize);
    printf("This is probably due to shear in the vox2ras matrix\n");
    printf("Input Vox2RAS ------\n");
    MatrixPrint(stdout,m_vox2ras);
  }
  mri->x_r = rx/xsize;
  mri->x_a = ax/xsize;
  mri->x_s = sx/xsize;

  mri->y_r = ry/ysize;
  mri->y_a = ay/ysize;
  mri->y_s = sy/ysize;

  mri->z_r = rz/zsize;
  mri->z_a = az/zsize;
  mri->z_s = sz/zsize;

  MRIp0ToCRAS(mri, P0r, P0a, P0s);
  return(NO_ERROR);
} 
/*!
  \fn int MRIsetVox2RASFromMatrixUnitTest(MRI *mri)
  \brief Unit test for MRIsetVox2RASFromMatrix(). Note
  the geometry params of the mri struct might be changed.
*/
int MRIsetVox2RASFromMatrixUnitTest(MRI *mri)
{
  MATRIX *Vox2RAS, *Vox2RASb;
  int c,r;
  double err;

  Vox2RAS = MRIxfmCRS2XYZ(mri,0);
  MRIsetVox2RASFromMatrix(mri,Vox2RAS);
  Vox2RASb = MRIxfmCRS2XYZ(mri,0);

  for(r=0; r<4; r++){
    for(c=0; c<4; c++){
      err = fabs(Vox2RAS->rptr[r+1][c+1]-Vox2RASb->rptr[r+1][c+1]);
      if(err > .00001){
	printf("MRIsetRas2VoxFromMatrixUnitTest() failed on element r=%d c=%d\n",r+1,c+1);
	printf("%g %g\n",Vox2RAS->rptr[r+1][c+1],Vox2RASb->rptr[r+1][c+1]);
	printf("Original Vox2RAS ------\n");
	MatrixPrint(stdout,Vox2RAS);
	printf("New Vox2RAS ------\n");
	MatrixPrint(stdout,Vox2RASb);
	return(1);
      }
    }
  }
  return(0);
}

/*-------------------------------------------------------------
  MRIxfmCRS2XYZtkreg() - computes the TkReg vox2ras, ie, linear
  transform between the column, row, and slice of a voxel and the x,
  y, z of that voxel as expected by tkregister (or for when using a
  tkregister compatible matrix). For tkregister, the column DC points
  in the "x" direction, the row DC points in the "z" direction, and
  the slice DC points in the "y" direction. The center of the
  coordinates is set to the center of the FOV. These definitions are
  arbitrary (and more than a little confusing). Since they are
  arbitrary, they must be applied consistently. See also:
  surfaceRASFromVoxel_ and voxelFromSurfaceRAS_.
  -------------------------------------------------------------*/
MATRIX *MRIxfmCRS2XYZtkreg( const MRI *mri)
{
  MRI *tmp;
  MATRIX *K;

  tmp = MRIallocHeader(mri->width, mri->height, mri->depth, mri->type, 1);

  /* Set tkregister defaults */
  /* column         row           slice          center      */
  tmp->x_r = -1;
  tmp->y_r =  0;
  tmp->z_r =  0;
  tmp->c_r = 0.0;
  tmp->x_a =  0;
  tmp->y_a =  0;
  tmp->z_a =  1;
  tmp->c_a = 0.0;
  tmp->x_s =  0;
  tmp->y_s = -1;
  tmp->z_s =  0;
  tmp->c_s = 0.0;

  /* Copy the voxel resolutions */
  tmp->xsize = mri->xsize;
  tmp->ysize = mri->ysize;
  tmp->zsize = mri->zsize;

  K = MRIxfmCRS2XYZ(tmp,0);

  MRIfree(&tmp);

  return(K);
}
/*-------------------------------------------------------------
  MRIxfmCRS2XYZfsl() - computes the FSL vox2ras, ie, linear
  transform between the column, row, and slice of a voxel and the x,
  y, z of that voxel as expected by FSL/FLIRT.
  -------------------------------------------------------------*/
MATRIX *MRIxfmCRS2XYZfsl(MRI *mri)
{
  MATRIX *vox2ras;
  vox2ras = MatrixAlloc(4,4,MATRIX_REAL);
  vox2ras->rptr[1][1] = mri->xsize;
  vox2ras->rptr[2][2] = mri->ysize;
  vox2ras->rptr[3][3] = mri->zsize;
  vox2ras->rptr[4][4] = 1.0;
  return(vox2ras);
}
/*-------------------------------------------------------------------
  MRItkReg2Native() - converts a tkregister-compatible registration
  matrix R to one that works with the native/scanner/xfm geometry. Use
  this function to convert from a MINC XFM to a TkReg. Note that there
  is a reversal in direction here as the tkreg R maps from Ref RAS to
  the Mov RAS, whereas the native D goes from Mov to Ref.  If R is
  null, it is assumed to be the identity.  In a typical application,
  ref is the anatomical volume and mov is the functional volume.

  See also: MRItkRegMtx, MRIxfmCRS2XYZtkreg, MRIxfmCRS2XYZ
  -------------------------------------------------------------------*/
MATRIX *MRItkReg2Native(MRI *ref, MRI *mov, MATRIX *R)
{
  MATRIX *Kref, *Kmov;
  MATRIX *Tref, *Tmov, *D;
  MATRIX *invKmov, *invTref;

  Tref = MRIxfmCRS2XYZ(ref,0); // Ref Native Vox2RAS
  Tmov = MRIxfmCRS2XYZ(mov,0); // Mov Native Vox2RAS

  Kref = MRIxfmCRS2XYZtkreg(ref); // Ref TkReg Vox2RAS
  Kmov = MRIxfmCRS2XYZtkreg(mov); // Mov TkReg Vox2RAS

  // D = Tref * inv(Kref) * inv(R) * Kmov * inv(Tmov)
  //   = inv( Tmov * inv(Kmov) * R * Kref * inv(Tref) )
  invKmov = MatrixInverse(Kmov,NULL);
  invTref = MatrixInverse(Tref,NULL);

  D = MatrixMultiply(Tmov,invKmov,NULL);
  if (R!=NULL) MatrixMultiply(D,R,D);
  MatrixMultiply(D,Kref,D);
  MatrixMultiply(D,invTref,D);
  MatrixInverse(D,D);

  if (0)
  {
    printf("MRITkReg2Native -----------------------------\n");
    printf("Tref ----------------\n");
    MatrixPrint(stdout,Tref);
    printf("Tmov ----------------\n");
    MatrixPrint(stdout,Tmov);
    printf("Kref ----------------\n");
    MatrixPrint(stdout,Kref);
    printf("Kmov ----------------\n");
    MatrixPrint(stdout,Kmov);
    printf("------------------------------------------\n");
  }

  MatrixFree(&Kref);
  MatrixFree(&Tref);
  MatrixFree(&Kmov);
  MatrixFree(&Tmov);
  MatrixFree(&invKmov);
  MatrixFree(&invTref);

  return(D);
}
/*----------------------------------------------------------------
  MRItkRegMtx() - creates a tkregsiter-compatible matrix from the
  matrix D that aligns the two volumes assuming the native/scanner/xfm
  geometry.  Use this function to convert from a TkReg to a MINC XFM.
  Note that there is a reversal in direction here as the tkreg R maps
  from Ref RAS to the Mov RAS, whereas the native D goes from Mov to
  Ref.  If D was created from running minc_tracc, then Ref should be
  the target volume and Mov should be the source volume. This is the
  counterpart to MRITkReg2Native(). If D is null, it is assumed to be
  the identity. This function could have been called
  MRInative2TkReg().
  ---------------------------------------------------------------*/
MATRIX *MRItkRegMtx(MRI *ref, MRI *mov, MATRIX *D)
{
  MATRIX *Kref, *Kmov;
  MATRIX *Tref, *Tmov, *R;
  MATRIX *invTmov, *invKref, *invD;

  /* Native Goemetry */
  Tref = MRIxfmCRS2XYZ(ref,0);  // Ref Native/Scanner/XFM vox2ras
  Tmov = MRIxfmCRS2XYZ(mov,0);  // Mov Native/Scanner/XFM vox2ras

  /* TkReg Goemetry */
  Kref = MRIxfmCRS2XYZtkreg(ref); // Ref TkReg vox2ras
  Kmov = MRIxfmCRS2XYZtkreg(mov); // Mov TkReg vox2ras

  invTmov = MatrixInverse(Tmov,NULL);
  invKref = MatrixInverse(Kref,NULL);

  // R = Kmov * inv(Tmov) * inv(D) * Tref * inv(Kref)
  R = MatrixMultiply(Kmov,invTmov,NULL);
  if (D != NULL)
  {
    invD = MatrixInverse(D,NULL);
    MatrixMultiply(R,invD,R);
    MatrixFree(&invD);
  }
  MatrixMultiply(R,Tref,R);
  MatrixMultiply(R,invKref,R);

  MatrixFree(&Kref);
  MatrixFree(&Tref);
  MatrixFree(&Kmov);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&invKref);

  return(R);
}
/*!
  \fn MATRIX *MRItkRegMtxFromVox2Vox(MRI *ref, MRI *mov, MATRIX *vox2vox)
  \brief Creates a tkregsiter-compatible ras2ras matrix from the
  vox2vox matrix. It is assumed that vox2vox maps the indices from the
  ref/target volume to that of the mov volume. If vox2vox is NULL, it
  is assumed to be the identity (which will return the same thing as
  MRItkRegMtx()). See also MRIvoxToVoxFromTkRegMtx().
*/
MATRIX *MRItkRegMtxFromVox2Vox(MRI *ref, MRI *mov, MATRIX *vox2vox)
{
  MATRIX *Kref, *Kmov;
  MATRIX *R, *invKref;

  if(vox2vox == NULL)  return(MRItkRegMtx(ref, mov, NULL));

  /* TkReg Goemetry */
  Kref = MRIxfmCRS2XYZtkreg(ref); // Ref TkReg vox2ras
  Kmov = MRIxfmCRS2XYZtkreg(mov); // Mov TkReg vox2ras

  invKref = MatrixInverse(Kref,NULL);

  // R = Kmov * vox2vox * inv(Kref)
  R = MatrixMultiply(Kmov,vox2vox,NULL);
  MatrixMultiply(R,invKref,R);

  MatrixFree(&Kref);
  MatrixFree(&Kmov);
  MatrixFree(&invKref);

  return(R);
}
/*!
  \fn MATRIX *MRIvoxToVoxFromTkRegMtx(MRI *mov, MRI *targ, MATRIX *tkR)
  \brief Creates a vox2vox from a tkregsiter-compatible ras2ras matrix
  The vox2vox maps the indices from the target volume to that of the
  mov volume. If tkR is NULL, it is assumed to be the identity.
  See also MRItkRegMtxFromVox2Vox().
*/
MATRIX *MRIvoxToVoxFromTkRegMtx(MRI *mov, MRI *targ, MATRIX *tkR)
{
  MATRIX *ras2vox_mov, *vox2ras_mov, *vox2ras_targ, *vox2vox ;
  vox2ras_mov  = MRIxfmCRS2XYZtkreg(mov);
  ras2vox_mov  = MatrixInverse(vox2ras_mov,NULL);
  vox2ras_targ = MRIxfmCRS2XYZtkreg(targ);
  if(tkR){
    vox2vox = MatrixMultiply(ras2vox_mov, tkR, NULL) ;
    MatrixMultiply(vox2vox,vox2ras_targ,vox2vox);
  }
  else 
    vox2vox = MatrixMultiply(ras2vox_mov, vox2ras_targ, NULL) ;
  MatrixFree(&vox2ras_mov);
  MatrixFree(&ras2vox_mov);
  MatrixFree(&vox2ras_targ);
  return(vox2vox) ;
}

/*-------------------------------------------------------------
  MRIfixTkReg() - this routine will adjust a matrix created by the
  "old" tkregister program. The old program had a problem in the way
  it chose the CRS of a voxel in the functional volume based on a
  point in the anatomical volume. The functional CRS of a point in
  anatomical space rarely (if ever) falls directly on a functional
  voxel, so it's necessary to choose a functional voxel given that the
  point falls between functional voxels (or it can be interpolated).
  The old tkregister program did not interpolate, rather it would
  choose the CRS in the following way: iC = floor(fC), iR = ceil(fR),
  and iS = floor(fS), where iC is the integer column number and fC is
  the floating point column, etc. Unfortunately, this is not nearest
  neighbor and it's not invertible. The right way to do it is to do
  nearest neighbor (ie, round to the closest integer). Unfortunately,
  there are a lot of data sets out there that have been regsitered
  with the old program, and we don't want to force poeple to
  reregister with the "new" program. This routine attempts to adjust
  the matrix created with the old program so that it will work with
  code that assumes that pure nearest neighbor was used.

  It does this by randomly sampling the anatomical volume in xyz
  and computing the tkreg'ed CRS for each point.

  Pcrs = inv(Tmov)*R*Pxyz
  PcrsTkReg = fcf(Pcrs) -- fcf is floor ceiling floor

  We seek a new R (Rfix) define with

  PcrsFix = inv(Tmov)*Rfix*Pxyz

  such that that the difference between PcrsFix and PcrsTkReg are
  minimized. To do this, we set

  PcrsFix = PcrsTkReg = inv(Tmov)*Rfix*Pxyz

  and solve for Rfix (this is an LMS solution):

  Rfix = Tmov*(PcrsTkReg*Pxyz')*inv(Pxyz*Pxyz');

  Applications that read in the registration matrix should detect the
  truncation method used (see below), fix the matrix if necessary, and
  proceed as if nearest neighbor/rounding was used.  The type of
  truncation can be determined from the last line of the registration
  file (after the matrix itself). If there is nothing there or the
  string "tkregister" is there, then the matrix should be
  converted. Otherwise, the string "round" should be there. The
  function regio_read_register (from registerio.c) will return the
  type of matrix in the float2int variable. It will be either
  FLT2INT_TKREG or FLT2INT_ROUND (constants defined in resample.h).
  ---------------------------------------------------------------*/
MATRIX *MRIfixTkReg(MRI *mov, MATRIX *R)
{
  int n, ntest = 1000;
  MATRIX *Pxyz, *Pcrs, *PcrsTkReg;
  MATRIX *PxyzT, *PxyzPxyzT, *invPxyzPxyzT;
  MATRIX *Tmov,*invTmov,*Rfix;
  MATRIX *tmp;
  float xrange, yrange, zrange;

  /* Assume a COR reference image */
  xrange = 256.0;
  yrange = 256.0;
  zrange = 256.0;

  Tmov = MRIxfmCRS2XYZtkreg(mov);
  invTmov = MatrixInverse(Tmov,NULL);

  Pxyz      = MatrixAlloc(4,ntest,MATRIX_REAL);
  PcrsTkReg = MatrixAlloc(4,ntest,MATRIX_REAL);

  /* Fill xyz with rand within the reference volume range */
  for (n=0;n<ntest;n++)
  {
    Pxyz->rptr[1][n+1] = xrange * (drand48()-0.5);
    Pxyz->rptr[2][n+1] = yrange * (drand48()-0.5);
    Pxyz->rptr[3][n+1] = zrange * (drand48()-0.5);
    Pxyz->rptr[4][n+1] = 1;
  }

  /* Compute floating mov CRS from targ XYZ */
  /* Pcrs = inv(Tmov)*R*Pxyz */
  tmp = MatrixMultiply(R,Pxyz,NULL);
  Pcrs = MatrixMultiply(invTmov,tmp,NULL);
  MatrixFree(&tmp);

  /* Truncate floating mov CRS using tkregister method*/
  for (n=0;n<ntest;n++)
  {
    PcrsTkReg->rptr[1][n+1] = floor(Pcrs->rptr[1][n+1]);
    PcrsTkReg->rptr[2][n+1] =  ceil(Pcrs->rptr[2][n+1]);
    PcrsTkReg->rptr[3][n+1] = floor(Pcrs->rptr[3][n+1]);
    PcrsTkReg->rptr[4][n+1] = 1;
  }
  MatrixFree(&Pcrs);

  //Rfix = Tmov*(PcrsTkreg*Pxyz')*inv(Pxyz*Pxyz');
  PxyzT = MatrixTranspose(Pxyz,NULL);
  PxyzPxyzT = MatrixMultiply(Pxyz,PxyzT,NULL);
  invPxyzPxyzT = MatrixInverse(PxyzPxyzT,NULL);
  tmp = MatrixMultiply(PcrsTkReg,PxyzT,NULL);
  MatrixMultiply(Tmov,tmp,tmp);
  Rfix = MatrixMultiply(tmp,invPxyzPxyzT,NULL);


  MatrixFree(&Pxyz);
  MatrixFree(&PcrsTkReg);
  MatrixFree(&PxyzT);
  MatrixFree(&PxyzPxyzT);
  MatrixFree(&invPxyzPxyzT);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&tmp);

  return(Rfix);
}

/*-------------------------------------------------------------------
  MRIfsl2TkReg() - converts an FSL registration matrix to one
  compatible with tkregister.  Note: the FSL matrix is assumed to map
  from the mov to the ref whereas the tkreg matrix maps from the ref
  to the mov.
  -------------------------------------------------------------------*/
MATRIX *MRIfsl2TkReg(MRI *ref, MRI *mov, MATRIX *FSLRegMat)
{
  MATRIX *RegMat=NULL, *invDmov, *Tmov, *Dref;
  MATRIX *invFSLRegMat, *invTref, *Tref;
  MATRIX *Qmov,*Qref;
  char *FSLOUTPUTTYPE=NULL;

  /* R = Tmov * inv(Dmov) * inv(Mfsl) * Dref * inv(Tref) */
  invDmov = MatrixAlloc(4,4,MATRIX_REAL);
  invDmov->rptr[1][1] = 1.0/mov->xsize;
  invDmov->rptr[2][2] = 1.0/mov->ysize;
  invDmov->rptr[3][3] = 1.0/mov->zsize;
  invDmov->rptr[4][4] = 1.0;
  Dref = MatrixAlloc(4,4,MATRIX_REAL);
  Dref->rptr[1][1] = ref->xsize;
  Dref->rptr[2][2] = ref->ysize;
  Dref->rptr[3][3] = ref->zsize;
  Dref->rptr[4][4] = 1.0;

  /*-------------------------------------------------------------------
    This next section of code alters the way the reg mat is computed
    based on the current FSL format and the determinant of the vox2ras
    matrices. When using ANALYZE, FSL did not know what the true
    vox2ras was. However, with NIFTI, FSL does know, and it treats the
    data differently depending upon whether the vox2ras has a negative
    or positive determinant. If positive, FSL (flirt) will flip the
    data left-right internally to make it negative, and this has to be
    taken into account when computing the tkreg matrix. It is not
    possible to stop FSL/flirt from flipping the data, even when
    supplying an initial matrix that has a flip in it.  With ANALYZE,
    FSL always assumes negative determinant, so no flipping  occurs.
    This change was also implemented in MRItkreg2FSL().
    -------------------------------------------------------------------*/
  FSLOUTPUTTYPE = getenv("FSLOUTPUTTYPE");
  if (FSLOUTPUTTYPE == NULL)
  {
    printf("ERROR: trying to convert FSL registration "
           "matrix to FreeSurfer,\n");
    printf("       but FSLOUTPUTTYPE variable is not set.\n");
    exit(1);
  }
  printf("FSLOUTPUTTYPE %s \n",FSLOUTPUTTYPE);
  if (!strcmp(FSLOUTPUTTYPE,"NIFTI") ||
      !strcmp(FSLOUTPUTTYPE,"NIFTI_PAIR") ||
      !strcmp(FSLOUTPUTTYPE,"NIFTI_GZ") )
  {
    Qref = MRIxfmCRS2XYZ(ref,0);
    Qmov = MRIxfmCRS2XYZ(mov,0);
    printf("fsl2TkReg: mov det = %g, ref det = %g\n",
           MatrixDeterminant(Qmov),MatrixDeterminant(Qref));
    if (MatrixDeterminant(Qmov) > 0)
    {
      printf("INFO: FSL2FreeSurfer: Mov volume is NIFTI with positive det,\n");
      printf("      applying LR flip to registration matrix.\n");
      invDmov->rptr[1][1] *= -1;
      invDmov->rptr[1][4] = mov->width-1;
      // Note: this is diff than for Dref because this is an inverse
    }
    if (MatrixDeterminant(Qref) > 0)
    {
      printf("INFO: FSL2FreeSurfer: Ref volume is NIFTI with positive det,\n");
      printf("      applying LR flip to registration matrix.\n");
      Dref->rptr[1][1] *= -1;
      Dref->rptr[1][4] = ref->xsize*(ref->width-1);
    }
    MatrixFree(&Qref);
    MatrixFree(&Qmov);
  }
  /*-------------------------------------------------------------------*/

  Tmov = MRIxfmCRS2XYZtkreg(mov);
  Tref = MRIxfmCRS2XYZtkreg(ref);
  invTref = MatrixInverse(Tref,NULL);

  invFSLRegMat = MatrixInverse(FSLRegMat,NULL);

  RegMat = MatrixMultiply(Tmov,invDmov,RegMat);
  RegMat = MatrixMultiply(RegMat,invFSLRegMat,RegMat);
  RegMat = MatrixMultiply(RegMat,Dref,RegMat);
  RegMat = MatrixMultiply(RegMat,invTref,RegMat);

  MatrixFree(&invDmov);
  MatrixFree(&invFSLRegMat);
  MatrixFree(&Tmov);
  MatrixFree(&Tref);
  MatrixFree(&invTref);
  MatrixFree(&Dref);

  return(RegMat);
}
/*-------------------------------------------------------------------
  MRItkreg2FSL() - converts tkregister registration matrix to one
  compatible with FSL. Note: the FSL matrix is assumed to map from the
  mov to the ref whereas the tkreg matrix maps from the ref to the
  mov.
  -------------------------------------------------------------------*/
MATRIX *MRItkreg2FSL(MRI *ref, MRI *mov, MATRIX *tkRegMat)
{
  MATRIX *FSLRegMat=NULL, *Dmov, *Tmov, *invTmov, *Tref, *Dref, *invDref;
  MATRIX *Qmov,*Qref;
  char *FSLOUTPUTTYPE=NULL;

  // R = Tmov * inv(Dmov) * inv(Mfsl) * Dref * inv(Tref)
  // Mfsl = Dref * inv(Tref) * inv(R) * Tmov * inv(Dmov)
  //      = inv( Dmov * inv(Tmov) * R * Tref * inv(Dref) )
  Dmov = MatrixAlloc(4,4,MATRIX_REAL);
  Dmov->rptr[1][1] = mov->xsize;
  Dmov->rptr[2][2] = mov->ysize;
  Dmov->rptr[3][3] = mov->zsize;
  Dmov->rptr[4][4] = 1.0;

  Dref = MatrixAlloc(4,4,MATRIX_REAL);
  Dref->rptr[1][1] = ref->xsize;
  Dref->rptr[2][2] = ref->ysize;
  Dref->rptr[3][3] = ref->zsize;
  Dref->rptr[4][4] = 1.0;

  /*-------------------------------------------------------------------
    This next section of code alters the way the reg mat is computed
    based on the current FSL format and the determinant of the vox2ras
    matrices. See MRIfsl2TkReg() above for a fuller explanation.
    -------------------------------------------------------------------*/
  FSLOUTPUTTYPE = getenv("FSLOUTPUTTYPE");
  if (FSLOUTPUTTYPE == NULL)
  {
    printf("ERROR: trying to convert FSL registration "
           "matrix to FreeSurfer,\n");
    printf("       but FSLOUTPUTTYPE variable is not set.\n");
    exit(1);
  }
  printf("FSLOUTPUTTYPE %s \n",FSLOUTPUTTYPE);
  if (!strcmp(FSLOUTPUTTYPE,"NIFTI") ||
      !strcmp(FSLOUTPUTTYPE,"NIFTI_PAIR") ||
      !strcmp(FSLOUTPUTTYPE,"NIFTI_GZ") )
  {
    Qref = MRIxfmCRS2XYZ(ref,0);
    Qmov = MRIxfmCRS2XYZ(mov,0);
    printf("tkreg2FSL: mov det = %g, ref det = %g\n",
           MatrixDeterminant(Qmov),MatrixDeterminant(Qref));
    if (MatrixDeterminant(Qmov) > 0)
    {
      printf("INFO: FSL2FreeSurfer: Mov volume is NIFTI with positive det,\n");
      printf("      applying LR flip to registration matrix.\n");
      Dmov->rptr[1][1] *= -1;
      Dmov->rptr[1][4] = mov->xsize*(mov->width-1);
    }
    if (MatrixDeterminant(Qref) > 0)
    {
      printf("INFO: FSL2FreeSurfer: Ref volume is NIFTI with positive det,\n");
      printf("      applying LR flip to registration matrix.\n");
      Dref->rptr[1][1] *= -1;
      Dref->rptr[1][4] = ref->xsize*(ref->width-1);
    }
    MatrixFree(&Qref);
    MatrixFree(&Qmov);
  }
  /*-------------------------------------------------------------------*/

  invDref = MatrixInverse(Dref,NULL);
  Tmov = MRIxfmCRS2XYZtkreg(mov);
  invTmov = MatrixInverse(Tmov,NULL);
  Tref = MRIxfmCRS2XYZtkreg(ref);

  FSLRegMat = MatrixMultiply(Dmov,invTmov,FSLRegMat);
  FSLRegMat = MatrixMultiply(FSLRegMat,tkRegMat,FSLRegMat);
  FSLRegMat = MatrixMultiply(FSLRegMat,Tref,FSLRegMat);
  FSLRegMat = MatrixMultiply(FSLRegMat,invDref,FSLRegMat);
  FSLRegMat = MatrixInverse(FSLRegMat,FSLRegMat);

  if (0)
  {
    printf("--- Dmov ---------------------\n");
    MatrixPrint(stdout,Dmov);
    printf("--- Tmov ---------------------\n");
    MatrixPrint(stdout,Tmov);
    printf("--- R ---------------------\n");
    MatrixPrint(stdout,tkRegMat);
    printf("--- Tref ---------------------\n");
    MatrixPrint(stdout,Tref);
    printf("--- Dref ---------------------\n");
    MatrixPrint(stdout,Dref);
    printf("--- Rfsl ---------------------\n");
    MatrixPrint(stdout,FSLRegMat);
    printf("--- R (from Rfsl) ------------\n");
    tkRegMat = MRIfsl2TkReg(ref,mov,FSLRegMat);
    MatrixPrint(stdout,tkRegMat);
  }

  MatrixFree(&Dmov);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&Tref);
  MatrixFree(&Dref);
  MatrixFree(&invDref);

  return(FSLRegMat);
}

/*--------------------------------------------------------------------------
  MtxCRS1toCRS0() - generates a matrix that will convert 1-based CRS (as
  found in SPM matrices) to 0-based CRS, ie, CRS0 = Q*CRS1 (assuming that
  CRS1 has been packed with a 1 in the 4th component.
  --------------------------------------------------------------------------*/
MATRIX *MtxCRS1toCRS0(MATRIX *Q)
{
  int r,c;

  if (Q == NULL) Q = MatrixAlloc(4,4,MATRIX_REAL);
  else
  {
    if (Q->rows != 4 || Q->cols != 4)
    {
      printf("ERROR: MtxCRS1toCRS0(): input matrix is not 4x4\n");
      return(NULL);
    }
  }

  for (r=1;r<=4;r++)
  {
    for (c=1;c<=4;c++)
    {
      if (r==c || c==4) Q->rptr[r][c] = 1.0;
      else             Q->rptr[r][c] = 0.0;
    }
  }

  return(Q);
}
/*------------------------------------------------------------------
  MRIp0ToCRAS() - Computes and sets the c_{ras} values in an MRI
  struct given P0 and assuming that the direction cosine and voxel
  resolution to be correct. P0 is the RAS coordinate at the "first"
  voxel (ie, col,row,slice=0). Typically, P0 is known from some
  source (eg, dicom header) and the c_{ras} needs to be computed.
  -----------------------------------------------------------------*/
int MRIp0ToCRAS(MRI *mri, double r0, double a0, double s0)
{
  MATRIX *vox2ras,*CRScenter,*RAScenter;

  // Get the vox2ras matrix.
  vox2ras = MRIxfmCRS2XYZ(mri, 0);
  // The last column will be wrong because the c_{ras} has not been
  // set properly (which is why you are calling this function).
  // So replace the last col with the RAS at the first voxel. This
  // makes the vox2ras correct.
  vox2ras->rptr[1][4] = r0;
  vox2ras->rptr[2][4] = a0;
  vox2ras->rptr[3][4] = s0;
  // Now set up a vector with the indices at the "center" of the
  // volume (not quite the center though).
  CRScenter = MatrixZero(4,1,NULL);
  CRScenter->rptr[1][1] = mri->width/2.0;
  CRScenter->rptr[2][1] = mri->height/2.0;
  CRScenter->rptr[3][1] = mri->depth/2.0;
  CRScenter->rptr[4][1] = 1;
  // Compute the RAS at the center as vox2ras*CRScenter
  RAScenter = MatrixMultiply(vox2ras,CRScenter,NULL);
  // Set the values in the MRI struct.
  mri->c_r = RAScenter->rptr[1][1];
  mri->c_a = RAScenter->rptr[2][1];
  mri->c_s = RAScenter->rptr[3][1];

  // Recompute matrix
  MATRIX *tmp;

  tmp = extract_i_to_r( mri );
  AffineMatrixAlloc( &(mri->i_to_r__ ) );
  SetAffineMatrix( mri->i_to_r__, tmp );
  MatrixFree( &tmp );

  if( mri->r_to_i__ ) {
    MatrixFree(&mri->r_to_i__);
  }
  mri->r_to_i__ = extract_r_to_i(mri);


  // Clean up
  MatrixFree(&vox2ras);
  MatrixFree(&CRScenter);
  MatrixFree(&RAScenter);
  // Get out of town
  return(0);
}
/*---------------------------------------------------------------
  MRIhfs2Sphinx() - reorient to sphinx the position. This function is
  applicable when the input geometry information is correct but the
  subject was in the scanner in the "sphinx" position (ie, AP in line
  with the bore) instead of head-first-supine (HFS). This is often the
  case with monkeys.
  ---------------------------------------------------------------*/
int MRIhfs2Sphinx(MRI *mri)
{
  double tmpxa,tmpya,tmpza,tmpca;

  // Negate right to make it left
  mri->x_r *= -1.0;
  mri->y_r *= -1.0;
  mri->z_r *= -1.0;
  mri->c_r *= -1.0;


  // Swap ant and sup
  tmpxa = mri->x_a;
  tmpya = mri->y_a;
  tmpza = mri->z_a;
  tmpca = mri->c_a;

  mri->x_a = mri->x_s;
  mri->y_a = mri->y_s;
  mri->z_a = mri->z_s;
  mri->c_a = mri->c_s;

  mri->x_s = tmpxa;
  mri->y_s = tmpya;
  mri->z_s = tmpza;
  mri->c_s = tmpca;

  return(0);
}
/*------------------------------------------------------*/
/*!
 \fn size_t MRIsizeof(int mritype)
 \brief Returns the sizeof() the MRI data type
*/
size_t MRIsizeof(int mritype)
{
  switch (mritype)
  {
  case MRI_UCHAR:
    return(sizeof(char));
    break;
  case MRI_SHORT:
    return(sizeof(short));
    break;
  case MRI_INT:
    return(sizeof(int));
    break;
  case MRI_LONG:
    return(sizeof(long));
    break;
  case MRI_FLOAT:
    return(sizeof(float));
    break;
  }
  // should never get here
  return(-1);
}

/*--------------------------------------------------------*/
/*!
 \fn double MRIptr2dbl(void *pmric, int mritype)
 \brief Derefences the pointer to a double.
 \param pmric - pointer to a column (eg, mri->slices[n][r])
 \param mritype - mri->type
 This is somewhat like MRIgetVoxVal(). It uses the pointer
 to the pixel data (instead of an MRI struct), which makes
 it several times faster but less general.
*/
inline double MRIptr2dbl(void *pmric, int mritype)
{
  double v=0;
  switch (mritype)
  {
  case MRI_UCHAR:
    v = (double)(*((char *)pmric));
    break;
  case MRI_SHORT:
    v = (double)(*((short*)pmric));
    break;
  case MRI_INT:
    v = (double)(*((int  *)pmric));
    break;
  case MRI_LONG:
    v = (double)(*((long *)pmric));
    break;
  case MRI_FLOAT:
    v = (double)(*((float*)pmric));
    break;
  }
  return(v);
}
/*--------------------------------------------------------*/
/*!
 \fn void MRIdbl2ptr(double v, void *pmric, int mritype)
 \brief Sets the derefenced pointer to the double.
 \param v - double value to use
 \param pmric - pointer to a column (eg, mri->slices[n][r])
 \param mritype - mri->type
 This is somewhat like MRIsetVoxVal(). It uses the pointer
 to the pixel data (instead of an MRI struct), which makes
 it several times faster but less general.
*/
inline void MRIdbl2ptr(double v, void *pmric, int mritype)
{
  switch (mritype)
  {
  case MRI_UCHAR:
    *((char*)pmric)  = nint(v);
    break;
  case MRI_SHORT:
    *((short*)pmric) = nint(v);
    break;
  case MRI_INT:
    *((int*)pmric)   = nint(v);
    break;
  case MRI_LONG:
    *((long*)pmric)  = nint(v);
    break;
  case MRI_FLOAT:
    *((float*)pmric) = v;
    break;
  }
}

/*-------------------------------------------------------------------*/
/*!
  \fn float MRIgetVoxDx(MRI *mri, int c, int r, int s, int f)
  \brief Returns voxel x derivative as a float regardless of the underlying data type.
  \param MRI *mri - input MRI
  \param int c - column
  \param int r - row
  \param int s - slice
  \param int f - frame
  \return float intensity x derivative at the given col, row, slice, frame
*/
float
MRIgetVoxDx(MRI *mri, int c, int r, int s, int f)
{
  float Ip1, Im1 ;

  Ip1 = MRIgetVoxVal(mri, c+1, r, s, f) ;
  Im1 = MRIgetVoxVal(mri, c-1, r, s, f) ;
  return((Ip1-Im1)/(2.0*mri->xsize)) ;
}
/*-------------------------------------------------------------------*/
/*!
  \fn float MRIgetVoxDy(MRI *mri, int c, int r, int s, int f)
  \brief Returns voxel y derivative as a float regardless of the underlying data type.
  \param MRI *mri - input MRI
  \param int c - column
  \param int r - row
  \param int s - slice
  \param int f - frame
  \return float intensity y derivative at the given col, row, slice, frame
  This function is general but slow. See also MRIptr2dbl().
*/
float
MRIgetVoxDy(MRI *mri, int c, int r, int s, int f)
{
  float Ip1, Im1 ;

  Ip1 = MRIgetVoxVal(mri, c, r+1, s, f) ;
  Im1 = MRIgetVoxVal(mri, c, r-1, s, f) ;
  return((Ip1-Im1)/(2.0*mri->ysize)) ;
}
/*-------------------------------------------------------------------*/
/*!
  \fn float MRIgetVoxDz(MRI *mri, int c, int r, int s, int f)
  \brief Returns voxel z derivative as a float regardless of the underlying data type.
  \param MRI *mri - input MRI
  \param int c - column
  \param int r - row
  \param int s - slice
  \param int f - frame
  \return float intensity z derivative at the given col, row, slice, frame
*/
float
MRIgetVoxDz(MRI *mri, int c, int r, int s, int f)
{
  float Ip1, Im1 ;

  Ip1 = MRIgetVoxVal(mri, c, r, s+1, f) ;
  Im1 = MRIgetVoxVal(mri, c, r, s-1, f) ;
  return((Ip1-Im1)/(2.0*mri->zsize)) ;
}
/*-------------------------------------------------------------------*/
/*!
  \fn float MRIgetVoxVal(MRI *mri, int c, int r, int s, int f)
  \brief Returns voxel value as a float regardless of the underlying data type.
  \param MRI *mri - input MRI
  \param int c - column
  \param int r - row
  \param int s - slice
  \param int f - frame
  \return float value at the given col, row, slice, frame
  This function is general but slow. See also MRIptr2dbl().
*/
inline float MRIgetVoxVal(const MRI *mri, int c, int r, int s, int f)
{
  static void *p=NULL;

  // bounds checks:
  if (c < 0) return mri->outside_val;
  if (r < 0) return mri->outside_val;
  if (s < 0) return mri->outside_val;

  if (mri->ischunked)
  {
    p = mri->chunk + c + r*mri->bytes_per_row +
        s*mri->bytes_per_slice + f*mri->bytes_per_vol;
    switch (mri->type)
    {
    case MRI_UCHAR:
      return( (float) *(unsigned char *)p);
      break;
    case MRI_SHORT:
      return( (float) *(short *)p);
      break;
    case MRI_INT:
      return( (float) *(int *)p);
      break;
    case MRI_LONG:
      return( (float) *(long *)p);
      break;
    case MRI_FLOAT:
      return( (float) *(float *)p);
      break;
    }
  }

  switch (mri->type)
  {
  case MRI_UCHAR:
    return( (float) MRIseq_vox(mri,c,r,s,f));
    break;
  case MRI_SHORT:
    return( (float) MRISseq_vox(mri,c,r,s,f));
    break;
  case MRI_INT:
    return( (float) MRIIseq_vox(mri,c,r,s,f));
    break;
  case MRI_LONG:
    return( (float) MRILseq_vox(mri,c,r,s,f));
    break;
  case MRI_FLOAT:
    return( (float) MRIFseq_vox(mri,c,r,s,f));
    break;
  }
  return(-10000000000.9);
}
/*-------------------------------------------------------------------*/
/*!
  \fn int MRIsetVoxVal(MRI *mri, int c, int r, int s, int f, float voxval)
  \brief Sets a voxel value regardless of the underlying data type.
  \param MRI *mri - input MRI
  \param int c - column
  \param int r - row
  \param int s - slice
  \param int f - frame
  \return int - 0 if ok, 1 if mri->type is unrecognized.
  This function is general but slow. See also MRIdbl2ptr().
*/
inline int MRIsetVoxVal(MRI *mri, int c, int r, int s, int f, float voxval)
{
  static void *p=NULL;
  
  //clipping
  switch (mri->type)
  {
  case MRI_UCHAR:
    if (voxval < UCHAR_MIN) voxval = UCHAR_MIN;
    if (voxval > UCHAR_MAX) voxval = UCHAR_MAX;
    break;
  case MRI_SHORT:
    if (voxval < SHORT_MIN) voxval = SHORT_MIN;
    if (voxval > SHORT_MAX) voxval = SHORT_MAX;
    break;
  case MRI_INT:
    if (voxval < INT_MIN) voxval = INT_MIN;
    if (voxval > INT_MAX) voxval = INT_MAX;
    break;
  case MRI_LONG:
    if (voxval < LONG_MIN) voxval = LONG_MIN;
    if (voxval > LONG_MAX) voxval = LONG_MAX;
    break;
  }  
  
  if (mri->ischunked)
  {
    p = mri->chunk + c + r*mri->bytes_per_row +
        s*mri->bytes_per_slice + f*mri->bytes_per_vol;
    switch (mri->type)
    {
    case MRI_UCHAR:
      *((unsigned char *)p) = nint(voxval);
      break;
    case MRI_SHORT:
      *((short *)p) = nint(voxval);
      break;
    case MRI_INT:
      *((int *)p)   = nint(voxval);
      break;
    case MRI_LONG:
      *((long *)p)  = nint(voxval);
      break;
    case MRI_FLOAT:
      *((float *)p) = voxval;
      break;
    }
  }

  switch (mri->type)
  {
  case MRI_UCHAR:
    MRIseq_vox(mri,c,r,s,f)  = nint(voxval);
    break;
  case MRI_SHORT:
    MRISseq_vox(mri,c,r,s,f) = nint(voxval);
    break;
  case MRI_INT:
    MRIIseq_vox(mri,c,r,s,f) = nint(voxval);
    break;
  case MRI_LONG:
    MRILseq_vox(mri,c,r,s,f) = nint(voxval);
    break;
  case MRI_FLOAT:
    MRIFseq_vox(mri,c,r,s,f) = voxval;
    break;
  default:
    return(1);
    break;
  }
  return(0);
}

/*------------------------------------------------------------------
  MRIinterpCode() - returns the numeric interpolation code given the
  name of the interpolation method.
  -----------------------------------------------------------------*/
int MRIinterpCode(char *InterpString)
{
  if (!strncasecmp(InterpString,"nearest",3))
    return(SAMPLE_NEAREST);
  if (!strncasecmp(InterpString,"trilinear",3))
    return(SAMPLE_TRILINEAR);
  if (!strncasecmp(InterpString,"tli",3))
    return(SAMPLE_TRILINEAR);
  if (!strncasecmp(InterpString,"sinc",3))
    return(SAMPLE_SINC);
  if (!strncasecmp(InterpString,"cubic",3))
    return(SAMPLE_CUBIC_BSPLINE);

  return(-1);
}

/*------------------------------------------------------------------
  MRIinterpString() - returns the the name of the interpolation method
  given numeric interpolation code
  -----------------------------------------------------------------*/
char * MRIinterpString(int InterpCode)
{
  switch (InterpCode)
  {
  case SAMPLE_NEAREST:
    return("nearest");
    break ;
  case SAMPLE_TRILINEAR:
    return("trilinear");
    break ;
  case SAMPLE_SINC:
    return("sinc");
    break ;
  case SAMPLE_CUBIC_BSPLINE:
    return("cubic");
    break ;
  }
  return(NULL);
}
/*------------------------------------------------------------------
  MRIprecisionCode() - returns the numeric code given the
  name of the precision. This corresponds to the value of the type
  field in the MRI structure.
  -----------------------------------------------------------------*/
int MRIprecisionCode(char *PrecisionString)
{
  if (!strcasecmp(PrecisionString,"uchar"))
    return(MRI_UCHAR);
  if (!strcasecmp(PrecisionString,"short"))
    return(MRI_SHORT);
  if (!strcasecmp(PrecisionString,"int"))
    return(MRI_INT);
  if (!strcasecmp(PrecisionString,"long"))
    return(MRI_LONG);
  if (!strcasecmp(PrecisionString,"float"))
    return(MRI_FLOAT);

  return(-1);
}

/*------------------------------------------------------------------
  MRIprecisionString() - returns the the name of the precision given
  numeric precision code. The code corresponds to the value of the
  type field in the MRI structure.
  -----------------------------------------------------------------*/
char * MRIprecisionString(int PrecisionCode)
{
  switch (PrecisionCode)
  {
  case MRI_UCHAR:
    return("uchar");
    break ;
  case MRI_SHORT:
    return("short");
    break ;
  case MRI_INT:
    return("int");
    break ;
  case MRI_LONG:
    return("long");
    break ;
  case MRI_FLOAT:
    return("float");
    break ;
  }
  return(NULL);
}

/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRImatch(MRI *mri1, MRI *mri2)
{
  return(
          (mri1->width == mri2->width) &&
          (mri1->height == mri2->height) &&
          (mri1->depth == mri2->depth) &&
          (mri1->type == mri2->type)
        ) ;
}
MRI *
MRIlinearScale(MRI *mri_src, MRI *mri_dst, float scale, float offset, 
               int only_nonzero) 
{
  int     width, height, depth, x, y, z, frame ;
  float   val ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, frame) ;
          if (!only_nonzero || !DZERO(val))
            val = val*scale+offset ;
          if (mri_dst->type == MRI_UCHAR)
          {
            if (val > 255)
              val = 255 ;
            else if (val < 0) 
              val = 0 ;
          }
          /*printf(" %d %d %d %d %f \n",x,y,z,frame,val);*/
          MRIsetVoxVal(mri_dst, x, y, z, frame, val) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
double 
MRIrmsDifferenceNonzero(MRI *mri1, MRI *mri2) 
{
  int     width, height, depth, x, y, z, frame, nvox ;
  float   val1, val2 ;
  double  sse ;

  width = mri1->width ; height = mri1->height ; depth = mri1->depth ;

  for (sse = 0.0, nvox = frame = 0 ; frame < mri1->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < height ; x++)
        {
          val1 = MRIgetVoxVal(mri1, x, y, z, frame) ;
          val2 = MRIgetVoxVal(mri2, x, y, z, frame) ;
          if (!DZERO(val1) && !DZERO(val2))
          {
            sse += (val1-val2)*(val1-val2) ;
            nvox++ ;
          }
        }
      }
    }
  }
  if (nvox == 0)
    return(0) ;
  return(sqrt(sse/nvox)) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIscalarMul(MRI *mri_src, MRI *mri_dst, float scalar)
{
  int     width, height, depth, x, y, z, frame ;
  BUFTYPE *psrc, *pdst ;
  float   *pfsrc, *pfdst, dval ;
  short   *pssrc, *psdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        if (1)
        {
          for (x = 0 ; x < width ; x++)
          {
            dval = MRIgetVoxVal(mri_src, x, y, z, frame) ;
            MRIsetVoxVal(mri_dst, x, y, z, frame, dval*scalar) ;
          }
        }
        else switch (mri_src->type)
        {
        case MRI_UCHAR:
          psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            dval = *psrc++ * scalar ;
            if (dval < 0)
              dval = 0 ;
            if (dval > 255)
              dval = 255 ;
            *pdst++ = dval ;
          }
          break ;
        case MRI_FLOAT:
          pfsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
          pfdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *pfdst++ = *pfsrc++ * scalar ;
          break ;
        case MRI_SHORT:
          pssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
          psdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
            *psdst++ = (short)nint((float)*pssrc++ * scalar) ;
          break ;
        default:
          ErrorReturn
          (NULL,
           (ERROR_UNSUPPORTED,
            "MRIscalarMul: unsupported type %d", mri_src->type)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIscalarMulFrame(MRI *mri_src, MRI *mri_dst, float scalar, int frame)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst ;
  float   *pfsrc, *pfdst, dval ;
  short   *pssrc, *psdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      switch (mri_src->type)
      {
      case MRI_UCHAR:
        psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
        pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
        {
          dval = *psrc++ * scalar ;
          if (dval < 0)
            dval = 0 ;
          if (dval > 255)
            dval = 255 ;
          *pdst++ = dval ;
        }
        break ;
      case MRI_FLOAT:
        pfsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
        pfdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
          *pfdst++ = *pfsrc++ * scalar ;
        break ;
      case MRI_SHORT:
        pssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
        psdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
          *psdst++ = (short)nint((float)*pssrc++ * scalar) ;
        break ;
      default:
        ErrorReturn
        (NULL,
         (ERROR_UNSUPPORTED,
          "MRIscalarMulFrame: unsupported type %d", mri_src->type)) ;
      }
    }
  }
  return(mri_dst) ;
}
int
MRIlabelValRange(MRI *mri,
                 MRI *mri_labeled,
                 int label,
                 float *pmin,
                 float *pmax)
{
  int      width, height, depth, x, y, z, frame, l, first = 1 ;
  float    fmin, fmax, val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri, x, y, z, frame) ;
          l = nint(MRIgetVoxVal(mri_labeled, x, y, z, frame)) ;
          if (l != label)
            continue ;
          if (first)
          {
            first = 0 ;
            fmax = fmin = val ;
          }
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }
  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvalRange(MRI *mri, float *pmin, float *pmax)
{
  int      width, height, depth, x, y, z, frame ;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  switch (mri->type)
  {
  case MRI_FLOAT:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          pf = &MRIFseq_vox(mri, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            val = *pf++ ;
            if (val < fmin)
              fmin = val ;
            if (val > fmax)
              fmax = val ;
          }
        }
      }
    }
    break ;
  case MRI_INT:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            val = (float)MRIIseq_vox(mri, x, y, z, frame) ;
            if (val < fmin)
              fmin = val ;
            if (val > fmax)
              fmax = val ;
          }
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            val = (float)MRISseq_vox(mri, x, y, z, frame) ;
            if (val < fmin)
              fmin = val ;
            if (val > fmax)
              fmax = val ;
          }
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          pb = &MRIseq_vox(mri, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            val = (float)*pb++ ;
            if (val < fmin)
              fmin = val ;
            if (val > fmax)
              fmax = val ;
          }
        }
      }
    }
    break ;
  default:
    for (frame = 0 ; frame < mri->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            val = (float)MRIgetVoxVal(mri, x, y, z, frame) ;
            if (val < fmin)
              fmin = val ;
            if (val > fmax)
              fmax = val ;
          }
        }
      }
    }
    break ;
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
int
MRInonzeroValRange(MRI *mri, float *pmin, float *pmax)
{
  int      width, height, depth, x, y, z, frame ;
  float    fmin, fmax, val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  for (frame = 0 ; frame < mri->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri, x, y, z, 0) ;
          if (FZERO(val))
            continue ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvalRangeFrame(MRI *mri, float *pmin, float *pmax, int frame)
{
  int      width, height, depth, x, y, z ;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  switch (mri->type)
  {
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf = &MRIFseq_vox(mri, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_INT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRIIseq_vox(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = (float)MRISseq_vox(mri, x, y, z, frame) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb = &MRIseq_vox(mri, 0, y, z, frame) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, "MRIvalRange: unsupported type %d",
                 mri->type)) ;
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvalRangeRegion(MRI *mri, float *pmin, float *pmax, MRI_REGION *region)
{
  int      width, height, depth, x, y, z, x0, y0, z0 ;
  float    fmin, fmax, *pf, val ;
  BUFTYPE  *pb ;

  width = region->x + region->dx ;
  if (width > mri->width)
    width = mri->width ;
  height = region->y + region->dy ;
  if (height > mri->height)
    height = mri->height ;
  depth = region->z + region->dz ;
  if (depth > mri->depth)
    depth = mri->depth ;
  x0 = region->x ;
  if (x0 < 0)
    x0 = 0 ;
  y0 = region->y ;
  if (y0 < 0)
    y0 = 0 ;
  z0 = region->z ;
  if (z0 < 0)
    z0 = 0 ;

  fmin = 10000.0f ;
  fmax = -10000.0f ;
  switch (mri->type)
  {
  default:
    for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
      {
        pf = &MRIFvox(mri, x0, y, z) ;
        for (x = x0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri, x, y, z, 0) ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;

  case MRI_FLOAT:
    for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
      {
        pf = &MRIFvox(mri, x0, y, z) ;
        for (x = x0 ; x < width ; x++)
        {
          val = *pf++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = z0 ; z < depth ; z++)
    {
      for (y = y0 ; y < height ; y++)
      {
        pb = &MRIvox(mri, x0, y, z) ;
        for (x = x0 ; x < width ; x++)
        {
          val = (float)*pb++ ;
          if (val < fmin)
            fmin = val ;
          if (val > fmax)
            fmax = val ;
        }
      }
    }
    break ;
  }

  *pmin = fmin ;
  *pmax = fmax ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI_REGION *
MRIclipRegion(MRI *mri, MRI_REGION *reg_src, MRI_REGION *reg_clip)
{
  int  x2, y2, z2 ;

  x2 = MIN(mri->width-1, reg_src->x + reg_src->dx - 1) ;
  y2 = MIN(mri->height-1, reg_src->y + reg_src->dy - 1) ;
  z2 = MIN(mri->depth-1, reg_src->z + reg_src->dz - 1) ;
  reg_clip->x = MAX(0, reg_src->x) ;
  reg_clip->y = MAX(0, reg_src->y) ;
  reg_clip->z = MAX(0, reg_src->z) ;
  reg_clip->dx = x2 - reg_clip->x + 1 ;
  reg_clip->dy = y2 - reg_clip->y + 1 ;
  reg_clip->dz = z2 - reg_clip->z + 1 ;
  return(reg_clip) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIvalScale(MRI *mri_src, MRI *mri_dst, float flo, float fhi)
{
  int      width, height, depth, x, y, z ;
  float    fmin, fmax, *pf_src, *pf_dst, val, scale ;
  short    *ps_src, *ps_dst ;
  BUFTYPE  *pb_src, *pb_dst ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRIvalRange(mri_src, &fmin, &fmax) ;
  scale = (fhi - flo) / (fmax - fmin) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  switch (mri_src->type)
  {
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pf_src = &MRIFvox(mri_src, 0, y, z) ;
        pf_dst = &MRIFvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = *pf_src++ ;
          *pf_dst++ = (val - fmin) * scale + flo ;
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        ps_src = &MRISvox(mri_src, 0, y, z) ;
        ps_dst = &MRISvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)(*ps_src++) ;
          *ps_dst++ = (short)nint((val - fmin) * scale + flo) ;
        }
      }
    }
    break ;
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pb_src = &MRIvox(mri_src, 0, y, z) ;
        pb_dst = &MRIvox(mri_dst, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          val = (float)*pb_src++ ;
          *pb_dst++ = (BUFTYPE)nint((val - fmin) * scale + flo) ;
        }
      }
    }
    break ;
  default:
    ErrorReturn(mri_dst,
                (ERROR_UNSUPPORTED, "MRIvalScale: unsupported type %d",
                 mri_src->type)) ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIconfThresh(MRI *mri_src, MRI *mri_probs, MRI *mri_classes, MRI *mri_dst,
              float thresh, int min_target, int max_target)
{
  int      x, y, z, width, height, depth, class ;
  float    *pprobs, prob ;
  BUFTYPE  *pclasses, *pdst, *psrc, src ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_classes, NULL) ;

  width = mri_classes->width ;
  height = mri_classes->height ;
  depth = mri_classes->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pprobs = &MRIFvox(mri_probs, 0, y, z) ;
      pclasses = &MRIvox(mri_classes, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        src = *psrc++ ;
        prob = *pprobs++ ;
        class = (int)*pclasses++ ;
        if (prob >= thresh &&
            ((class >= min_target) && (class <= max_target)))
          *pdst++ = src ;
        else if ((class >= min_target) && (class <= max_target))
          *pdst++ = 25 ;
        else
          *pdst++ = 0 ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIboundingBoxNbhd(MRI *mri, int thresh, int wsize,MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1, xi, yi, zi, xk, yk, zk,
  whalf, in_brain ;
  BUFTYPE  *psrc ;
  float    *pfsrc ;
  short    *pssrc ;

  whalf = (wsize-1)/2 ;
  box->dx = width = mri->width ;
  box->dy = height = mri->height ;
  box->dz = depth = mri->depth ;

  x1 = y1 = z1 = 0 ;
  box->x = width-1 ;
  box->y = height-1 ;
  box->z = depth-1 ;
  switch (mri->type)
  {
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*psrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ;
                     in_brain && xk <= whalf ;
                     xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRIvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
              if (x < box->x)
                box->x = x ;
              if (y < box->y)
                box->y = y ;
              if (z < box->z)
                box->z = z ;
              if (x > x1)
                x1 = x ;
              if (y > y1)
                y1 = y ;
              if (z > z1)
                z1 = z ;
            }
          }
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pssrc = &MRISvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pssrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ;
                     in_brain && xk <= whalf ;
                     xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRISvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
              if (x < box->x)
                box->x = x ;
              if (y < box->y)
                box->y = y ;
              if (z < box->z)
                box->z = z ;
              if (x > x1)
                x1 = x ;
              if (y > y1)
                y1 = y ;
              if (z > z1)
                z1 = z ;
            }
          }
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pfsrc = &MRIFvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pfsrc++ > thresh)
          {
            in_brain = 1 ;
            for (zk = -whalf ; in_brain && zk <= whalf ; zk++)
            {
              zi = mri->zi[z+zk] ;
              for (yk = -whalf ; in_brain && yk <= whalf ; yk++)
              {
                yi = mri->yi[y+yk] ;
                for (xk = -whalf ;
                     in_brain && xk <= whalf ;
                     xk++)
                {
                  xi = mri->xi[x+xk] ;
                  if (MRIFvox(mri, xi, yi, zi) < thresh)
                    in_brain = 0 ;
                }
              }
            }
            if (in_brain)
            {
              if (x < box->x)
                box->x = x ;
              if (y < box->y)
                box->y = y ;
              if (z < box->z)
                box->z = z ;
              if (x > x1)
                x1 = x ;
              if (y > y1)
                y1 = y ;
              if (z > z1)
                z1 = z ;
            }
          }
        }
      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, "MRIboundingBoxNbd: unsupported type %d",
                 mri->type)) ;
    break ;
  }
  box->x -= whalf+1 ;
  box->y -= whalf+1 ;
  box->z -= whalf+1 ;
  x1 += whalf+1 ;
  y1 += whalf+1 ;
  z1 += whalf+1 ;
  box->x = MAX(0,box->x) ;
  box->y = MAX(0,box->y) ;
  box->z = MAX(0,box->z) ;
  x1 = MIN(width-1,x1) ;
  y1 = MIN(height-1, y1) ;
  z1 = MIN(depth-1, z1) ;
  box->dx = x1 - box->x + 1 ;
  box->dy = y1 - box->y + 1 ;
  box->dz = z1 - box->z + 1 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
#define MIN_DARK 10

int
MRIfindApproximateSkullBoundingBox(MRI *mri, int thresh,MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1;
  int      ndark, max_dark, start, nlight, max_light, done = 0 ;
  double   means[3] ;
  double     val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  MRIcenterOfMass(mri, means, thresh) ;


#define MAX_LIGHT 30   /* don't let there by 3 cm of
  bright stuff 'outside' of brain */

  /* search for left edge */
  do
  {
    nlight = ndark = max_dark = 0 ;
    y = nint(means[1]) ;
    z = nint(means[2]) ;
    for (start = x1 = x = nint(means[0]) ; x >= 0 ; x--)
    {
      MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
      
      if (val < thresh)
      {
	if (!ndark)
	  start = x ;
	ndark++ ;
	nlight = 0  ;
      }
      else
      {
	if (++nlight > MAX_LIGHT)
	  max_dark = 0 ;
	if (ndark > max_dark)
	{
	  max_dark = ndark ;
	  x1 = start ;
	}
	ndark = 0 ;
      }
    }
    if (ndark > max_dark)
    {
      max_dark = ndark ;
      x1 = start ;
    }
    if (max_dark < MIN_DARK)
      x1 = 0 ;
    box->x = x1 ;

    /* search for right edge */
    nlight = ndark = max_dark = 0 ;
    y = nint(means[1]) ;
    z = nint(means[2]) ;
    for (start = x1 = x = nint(means[0]) ; x < width ; x++)
    {
      MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
      if (val < thresh)
      {
	if (!ndark)
	  start = x ;
	ndark++ ;
	nlight = 0 ;
      }
      else
      {
	if (++nlight > MAX_LIGHT)
	  max_dark = 0 ;
	if (ndark >= max_dark)
	{
	  max_dark = ndark ;
	  x1 = start ;
	}
	ndark = 0 ;
      }
    }
    if (ndark > max_dark)
    {
      max_dark = ndark ;
      x1 = start ;
    }
    if (max_dark < MIN_DARK)
      x1 = mri->width-1 ;
    box->dx = x1 - box->x + 1 ;
    if (box->dx <= 10)  // too small
    {
      done = 0 ;
      means[1] -= 10 ;
      printf("left/right detection failed, moving y coord to %d from %d\n", nint(means[1])+10, nint(means[1])) ;
    }
    else
      done = 1 ;
  } while (!done) ;

  /* search for superior edge */
  nlight = ndark = max_dark = max_light = 0 ;
  x = MAX(0,nint(means[0])-20) ; // avoid inter-hemispheric fissure
  z = nint(means[2]) ;
  for (start = y1 = y = nint(means[1]) ; y >= 0 ; y--)
  {
    MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
    if (val < thresh)
    {
      if (nlight > max_light)
        max_light = nlight ;
      if (!ndark)
        start = y ;
      ndark++ ;
      nlight = 0 ;
    }
    else
    {
      if (++nlight > MAX_LIGHT)
        max_dark = 0 ;
      if (ndark >= max_dark)
      {
        max_dark = ndark ;
        y1 = start ;
        max_light = 0 ;  // max_light is max in a row light above dark
      }
      ndark = 0 ;
    }
  }

  /* if we ended on a string of dark voxels, check two things:
     1. the string was longer than the previous longest
     2. the strong was longer than 1/2 the previous longest, and there
     was an intervening string of light voxels indicated it was still in
     brain.
  */
  if ((ndark > max_dark) || (y < 0 &&
                             (ndark > max_dark/2) && max_light > MAX_LIGHT/2))
  {
    max_dark = ndark ;
    y1 = start ;
  }
  if (max_dark < MIN_DARK)
    y1 = 0 ;
  box->y = y1 ;

  /* search for inferior edge */
  nlight = ndark = max_dark = 0 ;
  x = nint(means[0]) ;
  z = nint(means[2]) ;
  for (start = y = y1 = nint(means[1]) ; y < height ; y++)
  {
    MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
    if (val < thresh)
    {
      if (!ndark)
        start = y ;
      ndark++ ;
      nlight = 0 ;
    }
    else
    {
      if (++nlight > MAX_LIGHT)
        max_dark = 0 ;
      if (ndark >= max_dark)
      {
        max_dark = ndark ;
        y1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    y1 = start ;
  }
  if (max_dark < MIN_DARK)
    y1 = mri->height-1 ;
  box->dy = y1 - box->y + 1 ;

  /* search for posterior edge */
  nlight = ndark = max_dark = 0 ;
  x = nint(means[0]) ;
  y = nint(means[1]) ;
  for (z1 = start = z = nint(means[2]) ; z >= 0 ; z--)
  {
    MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
    if (val < thresh)
    {
      if (!ndark)
        start = z ;
      ndark++ ;
      nlight = 0 ;
    }
    else
    {
      if (++nlight > MAX_LIGHT)
        max_dark = 0 ;
      if (ndark >= max_dark)
      {
        max_dark = ndark ;
        z1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    z1 = start ;
  }
  if (max_dark < MIN_DARK)
    z1 = 0 ;
  box->z = z1 ;

  /* search for anterior edge */
  nlight = ndark = max_dark = 0 ;
  x = nint(means[0]) ;
  y = nint(means[1]) ;
  for (start = z = nint(means[2]) ; z < depth ; z++)
  {
    MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
    if (val < thresh)
    {
      if (!ndark)
        start = z ;
      ndark++ ;
      nlight = 0 ;
    }
    else
    {
      if (++nlight > MAX_LIGHT)
        max_dark = 0 ;
      if (ndark >= max_dark)
      {
        max_dark = ndark ;
        z1 = start ;
      }
      ndark = 0 ;
    }
  }
  if (ndark > max_dark)
  {
    max_dark = ndark ;
    z1 = start ;
  }
  if (max_dark < MIN_DARK)
    z1 = mri->depth-1 ;
  box->dz = z1 - box->z + 1 ;
  if (box->dz < 5) // have to avoid being right on the midline
  {
    nlight = ndark = max_dark = 0 ;
    x -= 15 ;
    y = nint(means[1]) ;
    for (z1 = start = z = nint(means[2]) ; z >= 0 ; z--)
    {
      MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
      if (val < thresh)
      {
        if (!ndark)
          start = z ;
        ndark++ ;
        nlight = 0 ;
      }
      else
      {
        if (++nlight > MAX_LIGHT)
          max_dark = 0 ;
        if (ndark >= max_dark)
        {
          max_dark = ndark ;
          z1 = start ;
        }
        ndark = 0 ;
      }
    }
    if (ndark > max_dark)
    {
      max_dark = ndark ;
      z1 = start ;
    }
    if (max_dark < MIN_DARK)
      z1 = 0 ;
    box->z = z1 ;
    nlight = ndark = max_dark = 0 ;
    for (start = z = nint(means[2]) ; z < depth ; z++)
    {
      MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
      if (val < thresh)
      {
        if (!ndark)
          start = z ;
        ndark++ ;
        nlight = 0 ;
      }
      else
      {
        if (++nlight > MAX_LIGHT)
          max_dark = 0 ;
        if (ndark >= max_dark)
        {
          max_dark = ndark ;
          z1 = start ;
        }
        ndark = 0 ;
      }
    }
    if (ndark > max_dark)
    {
      max_dark = ndark ;
      z1 = start ;
    }
    if (max_dark < MIN_DARK)
      z1 = mri->depth-1 ;
    box->dz = z1 - box->z + 1 ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIboundingBox(MRI *mri, int thresh, MRI_REGION *box)
{
  int      width, height, depth, x, y, z, x1, y1, z1 ;
  BUFTYPE  *psrc ;
  float    *pfsrc ;
  short    *pssrc ;
  int      *pisrc ;

  box->dx = width = mri->width ;
  box->dy = height = mri->height ;
  box->dz = depth = mri->depth ;

  x1 = y1 = z1 = 0 ;
  box->x = width-1 ;
  box->y = height-1 ;
  box->z = depth-1 ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        psrc = &MRIvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*psrc++ > thresh)
          {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
          }
        }
      }
    }
    break ;
  case MRI_FLOAT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pfsrc = &MRIFvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pfsrc++ > thresh)
          {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
          }
        }
      }
    }
    break ;
  case MRI_SHORT:
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        pssrc = &MRISvox(mri, 0, y, z) ;
        for (x = 0 ; x < width ; x++)
        {
          if (*pssrc++ > thresh)
          {
            if (x < box->x)
              box->x = x ;
            if (y < box->y)
              box->y = y ;
            if (z < box->z)
              box->z = z ;
            if (x > x1)
              x1 = x ;
            if (y > y1)
              y1 = y ;
            if (z > z1)
              z1 = z ;
          }
        }
      }
    }
    break ;
  case MRI_INT:
    for (z = 0 ; z < depth ; z++){
      for (y = 0 ; y < height ; y++){

        pisrc = &MRIIvox(mri, 0, y, z);
        for (x = 0 ; x < width ; x++)  {
          if(*pisrc++ > thresh) {
            if (x < box->x)          box->x = x ;
            if (y < box->y)          box->y = y ;
            if (z < box->z)          box->z = z ;
            if (x > x1)              x1 = x ;
            if (y > y1)              y1 = y ;
            if (z > z1)              z1 = z ;
          }
        }

      }
    }
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED, "MRIboundingBox: unsupported type %d",
                 mri->type)) ;
    break ;
  }
  box->dx = x1 - box->x + 1 ;
  box->dy = y1 - box->y + 1 ;
  box->dz = z1 - box->z + 1 ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIcheckSize(MRI *mri_src, MRI *mri_check, int width, int height, int depth)
{
  if (!mri_check)
    return(0) ;

  if (!width)
    width = mri_src->width ;
  if (!height)
    height = mri_src->height ;
  if (!depth)
    depth = mri_src->depth ;

  if (width != mri_check->width ||
      height != mri_check->height ||
      depth != mri_check->depth)
    return(0) ;

  return(1) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRItransformRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *src_region,
                   MRI_REGION *dst_region)
{
  double  xw, yw, zw, xt, yt, zt, xv, yv, zv ;

  if (getSliceDirection(mri_src) != getSliceDirection(mri_dst))
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItransformRegion(%s): slice directions must match",
                 mri_src->fname)) ;

  if (!mri_src->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItransformRegion(%s): no transform loaded",
                 mri_src->fname)) ;
  if (!mri_dst->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItransformRegion(%s): no transform loaded",
                 mri_dst->fname)) ;
  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri_src))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIregionToTalairachRegion: unsupported slice direction %d",
      getSliceDirection(mri_src))) ;
  }

  xv = (double)src_region->x ;
  yv = (double)src_region->y ;
  zv = (double)src_region->z ;
  MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
  transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  transform_point(mri_dst->inverse_linear_transform, xt, yt, zt, &xw,&yw,&zw);
  MRIworldToVoxel(mri_dst, xw, yw, zw, &xv, &yv, &zv) ;
  dst_region->x = nint(xv) ;
  dst_region->y = nint(yv) ;
  dst_region->z = nint(zv) ;

  xv = (double)(src_region->x + src_region->dx - 1) ;
  yv = (double)(src_region->y + src_region->dy - 1) ;
  zv = (double)(src_region->z + src_region->dz - 1) ;
  MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
  transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  transform_point(mri_dst->inverse_linear_transform, xt, yt, zt, &xw,&yw,&zw);
  MRIworldToVoxel(mri_dst, xw, yw, zw, &xv, &yv, &zv) ;
  dst_region->dx = nint(xv - (double)dst_region->x) + 1 ;
  dst_region->dy = nint(yv - (double)dst_region->y) + 1 ;
  dst_region->dz = nint(zv - (double)dst_region->z) + 1 ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvoxelToVoxel(MRI *mri_src, MRI *mri_dst, double xv, double yv, double zv,
                double *pxt, double *pyt, double *pzt)
{
  double  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri_src->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToVoxel(%s): no transform loaded", mri_src->fname));
#endif

  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri_src))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToVoxel: unsupported slice direction %d",
                 getSliceDirection(mri_src))) ;
  }

  if (!mri_src->linear_transform || !mri_dst->inverse_linear_transform)
  {
    /*
      if either doesn't have a transform defined, assume they are in
      the same coordinate system.
    */
    *pxt = xv ;
    *pyt = yv ;
    *pzt = zv ;
  }
  else
  {
    MRIvoxelToWorld(mri_src, xv, yv, zv, &xw, &yw, &zw) ;
    if (mri_src->linear_transform)
      transform_point(mri_src->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
    else
    {
      xt = xw ;
      yt = yw ;
      zt = zw ;
    }
    if (mri_dst->inverse_linear_transform)
      transform_point
      (mri_dst->inverse_linear_transform, xt,yt,zt,&xw,&yw,&zw);
    else
    {
      xw = xt ;
      yw = yt ;
      zw = zt ;
    }
    MRIworldToVoxel(mri_dst, xw, yw, zw, pxt, pyt, pzt) ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvoxelToTalairachVoxel(MRI *mri, double xv, double yv, double zv,
                         double *pxt, double *pyt, double *pzt)
{
  double  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif
  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 getSliceDirection(mri))) ;
  }

  MRIvoxelToWorld(mri, xv, yv, zv, &xw, &yw, &zw) ;
  if (mri->linear_transform)
    transform_point(mri->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  else
  {
    xt = xw ;
    yt = yw ;
    zt = zw ;
  }
  MRIworldToVoxel(mri, xt, yt, zt, pxt, pyt, pzt) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIvoxelToTalairach(MRI *mri, double xv, double yv, double zv,
                    double *pxt, double *pyt, double *pzt)
{
  double  xw, yw, zw ;

#if 0
  if (!mri->linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 getSliceDirection(mri))) ;
  }

  MRIvoxelToWorld(mri, xv, yv, zv, &xw, &yw, &zw) ;
  if (mri->linear_transform)
    transform_point(mri->linear_transform, xw, yw, zw, pxt, pyt, pzt) ;
  else
  {
    *pxt = xw ;
    *pyt = yw ;
    *pzt = zw ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRItalairachToVoxel(MRI *mri, double xt, double yt, double zt,
                    double *pxv, double *pyv, double *pzv)
{
  double  xw, yw, zw ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItalairachToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 getSliceDirection(mri))) ;
  }

  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  {
    xw = xt ;
    yw = yt ;
    zw = zt ;
  }
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRItalairachVoxelToVoxel(MRI *mri, double xtv, double ytv, double ztv,
                         double *pxv, double *pyv, double *pzv)
{
  double  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItalairachVoxelToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 getSliceDirection(mri))) ;
  }

  MRIvoxelToWorld(mri, xtv, ytv, ztv, &xt, &yt, &zt) ;
  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  {
    xw = xt ;
    yw = yt ;
    zw = zt ;
  }
  MRIworldToVoxel(mri, xw, yw, zw, pxv, pyv, pzv) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRItalairachVoxelToWorld(MRI *mri, double xtv, double ytv, double ztv,
                         double *pxw, double *pyw, double *pzw)
{
  double  xw, yw, zw, xt, yt, zt ;

#if 0
  if (!mri->inverse_linear_transform)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRItalairachVoxelToVoxel(%s): no transform loaded",
                 mri->fname)) ;
#endif

  /*
    The convention  is  that  positive xspace coordinates run
    from the patient's  left  side  to  right  side,  positive
    yspace  coordinates run from patient posterior to anterior
    and positive zspace coordinates run from inferior to superior.
  */
  switch (getSliceDirection(mri))
  {
  case MRI_CORONAL:
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIvoxelToTalairachVoxel: unsupported slice direction %d",
                 getSliceDirection(mri))) ;
  }

  MRIvoxelToWorld(mri, xtv, ytv, ztv, &xt, &yt, &zt) ;
  if (mri->inverse_linear_transform)
    transform_point(mri->inverse_linear_transform, xt, yt, zt, &xw, &yw,&zw);
  else
  {
    xw = xt ;
    yw = yt ;
    zw = zt ;
  }
  *pxw = xw ;
  *pyw = yw ;
  *pzw = zw ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
#define V4_LOAD(v, x, y, z, r)  (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z, VECTOR_ELT(v,4)=r) ;

int MRIvoxelToWorld(MRI *mri, double xv, double yv, double zv,
                    double *pxw, double *pyw, double *pzw)
{

  AffineVector vw, vv;


  // if the transform is not cached yet, then
  if( !mri->i_to_r__ ) {
    AffineMatrixAlloc( &(mri->i_to_r__) );
    MATRIX *tmp = extract_i_to_r(mri);
    SetAffineMatrix( mri->i_to_r__, tmp );
    MatrixFree( &tmp );
  }

  if (!mri->r_to_i__) { 
    mri->r_to_i__ = extract_r_to_i(mri);
  }

  // Do matrix-vector multiply
  SetAffineVector( &vv, xv, yv, zv );
  AffineMV( &vw, mri->i_to_r__, &vv );

  // Extract the results
  float xwf, ywf, zwf;
  GetAffineVector( &vw, &xwf, &ywf, &zwf );

  *pxw = xwf;
  *pyw = ywf;
  *pzw = zwf;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIworldToTalairachVoxel(MRI *mri, double xw, double yw, double zw,
                         double *pxv, double *pyv, double *pzv)
{
  double  xt, yt, zt ;

  transform_point(mri->linear_transform, xw, yw, zw, &xt, &yt, &zt) ;
  MRIworldToVoxel(mri, xt, yt, zt, pxv, pyv, pzv) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int MRIworldToVoxelIndex(MRI *mri, double xw, double yw, double zw,
                         int *pxv, int *pyv, int *pzv)
{
  double xv, yv, zv;
  MRIworldToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
  // Changed on 7/25/08 so that it rounds intead of truncates
  *pxv = (int)rint(xv);
  *pyv = (int)rint(yv);
  *pzv = (int)rint(zv);
  return(NO_ERROR) ;
}

/*
  Tosa: MRIvoxelToSurfaceRAS and MRIsurfaceRASToVoxel get the
  surfaceRAS values from original voxel Note that this is different
  from MRIvoxelToWorld().  Note that currently MATRIX uses float** to
  store data.  Going around the circle of transform causes error
  accumulation quickly.  I noticed that 7 x 10^(-6) is very common.

  Doug: Tosa's code has been modified somewhat extensively
  (2/27/06). What he calls "surface" RAS is really supposed to be
  "tkregister" RAS.  They are the same with conformed volumes, but
  surface RAS is wrong otherwise.
*/

MATRIX *surfaceRASFromVoxel_(MRI *mri)
{
  MATRIX *vox2ras;

  // Compute i_to_r and r_to_i if it has not been done yet. This is
  // not necessary for this function, but it was in Tosa's original
  // code, and I don't know what else might be using it.
  if (!mri->i_to_r__) {
    MATRIX *tmp = extract_i_to_r(mri);
    AffineMatrixAlloc( &(mri->i_to_r__) );
    SetAffineMatrix( mri->i_to_r__, tmp );
    MatrixFree( &tmp );
  }

  if (!mri->r_to_i__) {
    mri->r_to_i__ = extract_r_to_i(mri);
  }

  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  if (Gdiag_no > 0 && DIAG_VERBOSE_ON)
  {
    printf("surfaceRASFromVoxel_() vox2ras --------------------\n");
    MatrixPrint(stdout,vox2ras);
  }
  return(vox2ras);

  /*-----------------------------------------------------------------
    This was tosa's old code. It is broken in that it only
    works for COR-oriented volumes.
    rasFromVoxel = mri->i_to_r__; // extract_i_to_r(mri);
    sRASFromVoxel = MatrixCopy(rasFromVoxel, NULL);
    // MatrixFree(&rasFromVoxel);
    // modify
    m14 = *MATRIX_RELT(sRASFromVoxel, 1,4);
    *MATRIX_RELT(sRASFromVoxel, 1,4) = m14 - mri->c_r;
    m24 = *MATRIX_RELT(sRASFromVoxel, 2,4);
    *MATRIX_RELT(sRASFromVoxel, 2,4) = m24 - mri->c_a;
    m34 = *MATRIX_RELT(sRASFromVoxel, 3,4);
    *MATRIX_RELT(sRASFromVoxel, 3,4) = m34 - mri->c_s;
    ---------------------------------------------------*/
}
/*------------------------------------------------------------------
  voxelFromSurfaceRAS_() - this is a Tosa function that is supposed
  to compute the voxel index from "surface" RAS. The problem is that
  it only worked for COR-oriented volumes. What it is supposed to do
  is compute the tkregister-style Vox2RAS/RAS2Vox, so now it simply
  invertes the matrix from MRIxfmCRS2XYZtkreg(mri).
  *-------------------------------------------------------------------*/
MATRIX *voxelFromSurfaceRAS_(MRI *mri)
{
  MATRIX *vox2ras, *ras2vox;
  // Compute i_to_r and r_to_i if it has not been done yet. This is
  // not necessary for this function, but it was in Tosa's original
  // code, and I don't know what else might be using it.
  if (!mri->i_to_r__) {
    MATRIX *tmp = extract_i_to_r(mri);
    AffineMatrixAlloc( &(mri->i_to_r__) );
    SetAffineMatrix( mri->i_to_r__, tmp );
    MatrixFree( &tmp );
  }

  if (!mri->r_to_i__) {
    mri->r_to_i__ = extract_r_to_i(mri);
  }
  vox2ras = MRIxfmCRS2XYZtkreg(mri);
  ras2vox = MatrixInverse(vox2ras,NULL);
  if (Gdiag_no > 0 && DIAG_VERBOSE_ON)
  {
    printf("voxelFromSurfaceRAS_() ras2vox --------------------\n");
    MatrixPrint(stdout,ras2vox);
  }
  MatrixFree(&vox2ras);
  return(ras2vox);
  /*----------------------------------------------------------
    This was tosa's old code. It is broken in that it only
    works for COR-oriented volumes.
    voxelFromSRAS = MatrixCopy(mri->r_to_i__, NULL);
    // modify translation part
    *MATRIX_RELT(voxelFromSRAS, 1,4) = (double)mri->width/2.0;
    *MATRIX_RELT(voxelFromSRAS, 2,4) = (double)mri->height/2.0;
    *MATRIX_RELT(voxelFromSRAS, 3,4) = (double)mri->depth/2.0;
    ---------------------------------------------------*/
}

/*------------------------------------------------------------------
  surfaceRASFromRAS_(MRI *mri): creates a matrix that converts from
  Scanner RAS to TkRAS (what Tosa called "surface"  RAS). Note:
  intermediate matrices are alloced, inverted, and dealloced, so
  it might not be a good thing to have inside a loop.
  --------------------------------------------------------------*/
MATRIX *surfaceRASFromRAS_(MRI *mri)
{
  MATRIX *sRASFromRAS;
  MATRIX *Vox2TkRAS, *Vox2RAS;

  Vox2RAS = MRIxfmCRS2XYZ(mri,0); // scanner vox2ras
  Vox2TkRAS = MRIxfmCRS2XYZtkreg(mri); // tkreg vox2ras
  // sRASFromRAS = Vox2TkRAS * inv(Vox2RAS)
  sRASFromRAS = MatrixInverse(Vox2RAS,NULL);
  sRASFromRAS = MatrixMultiply(Vox2TkRAS,sRASFromRAS,sRASFromRAS);
  MatrixFree(&Vox2RAS);
  MatrixFree(&Vox2TkRAS);
  return(sRASFromRAS);

  /*-----------------------------------------
  // Tosa's code: only works for conformed vols
  sRASFromRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, sRASFromRAS);
  *MATRIX_RELT(sRASFromRAS, 1,4) = - mri->c_r;
  *MATRIX_RELT(sRASFromRAS, 2,4) = - mri->c_a;
  *MATRIX_RELT(sRASFromRAS, 3,4) = - mri->c_s;
  return sRASFromRAS;
  *---------------------------------------------*/
}
/*------------------------------------------------------------------
  RASFromSurfaceRAS_(MRI *mri): creates a matrix that converts from
  TkRAS to Scanner RAS (what Tosa called "surface"  RAS). Note:
  intermediate matrices are alloced, inverted, and dealloced, so
  it might not be a good thing to have inside a loop.
  --------------------------------------------------------------*/
MATRIX *RASFromSurfaceRAS_(MRI *mri)
{
  MATRIX *RASFromSRAS;

  MATRIX *Vox2TkRAS, *Vox2RAS;

  Vox2RAS = MRIxfmCRS2XYZ(mri,0); // scanner vox2ras
  Vox2TkRAS = MRIxfmCRS2XYZtkreg(mri); // tkreg vox2ras
  // RASFromSRAS = Vox2RAS * inv(Vox2TkRAS)
  RASFromSRAS = MatrixInverse(Vox2TkRAS,NULL);
  RASFromSRAS = MatrixMultiply(Vox2RAS,RASFromSRAS,RASFromSRAS);
  MatrixFree(&Vox2RAS);
  MatrixFree(&Vox2TkRAS);
  return(RASFromSRAS);
  /*-----------------------------------------
  // Tosa's code: only works for conformed vols
  RASFromSRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, RASFromSRAS);
  *MATRIX_RELT(RASFromSRAS, 1,4) = mri->c_r;
  *MATRIX_RELT(RASFromSRAS, 2,4) = mri->c_a;
  *MATRIX_RELT(RASFromSRAS, 3,4) = mri->c_s;
  return RASFromSRAS;
  *---------------------------------------------*/
}
/*--------------------------------------------------------------
  MRIRASToSurfaceRAS() - convert from a scanner RAS to a TkReg (or
  "surface") RAS. Note: intermediate matrices are alloced, inverted,
  and dealloced, so it might not be a good thing to have inside a
  loop.
  -------------------------------------------------------------*/
int MRIRASToSurfaceRAS(MRI *mri, double xr, double yr, double zr,
                       double *xsr, double *ysr, double *zsr)
{
  MATRIX *surfaceRASFromRAS=0;
  VECTOR *v, *sr;
  v = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(v, xr, yr, zr, 1.);
  surfaceRASFromRAS = surfaceRASFromRAS_(mri);
  sr = MatrixMultiply(surfaceRASFromRAS, v, NULL);
  *xsr = V3_X(sr);
  *ysr = V3_Y(sr);
  *zsr = V3_Z(sr);
  MatrixFree(&surfaceRASFromRAS);
  VectorFree(&v);
  VectorFree(&sr);
  return (NO_ERROR);
}
/*--------------------------------------------------------------
  MRIRASToSurfaceRAS() - convert from a TkReg (or "surface") RAS to a
  scanner RAS. Note: intermediate matrices are alloced, inverted, and
  dealloced, so it might not be a good thing to have inside a loop.
  -------------------------------------------------------------*/
int MRIsurfaceRASToRAS(MRI *mri, double xsr, double ysr, double zsr,
                       double *xr, double *yr, double *zr)
{
  MATRIX *RASFromSurfaceRAS=0;
  VECTOR *v, *r;
  v = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(v, xsr, ysr, zsr, 1.);
  RASFromSurfaceRAS = RASFromSurfaceRAS_(mri);
  r = MatrixMultiply(RASFromSurfaceRAS, v, NULL);
  *xr = V3_X(r);
  *yr = V3_Y(r);
  *zr = V3_Z(r);
  MatrixFree(&RASFromSurfaceRAS);
  VectorFree(&v);
  VectorFree(&r);
  return (NO_ERROR);
}

//--------------------------------------------------------------
int MRIvoxelToSurfaceRAS(MRI *mri, double xv, double yv, double zv,
                         double *xs, double *ys, double *zs)
{
  MATRIX *sRASFromVoxel;
  VECTOR *vv, *sr;

  sRASFromVoxel = surfaceRASFromVoxel_(mri);
  // calculate the surface ras value
  vv = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(vv, xv, yv, zv, 1.);
  sr = MatrixMultiply(sRASFromVoxel, vv, NULL);
  *xs = V3_X(sr);
  *ys = V3_Y(sr);
  *zs = V3_Z(sr);

  MatrixFree(&sRASFromVoxel);
  VectorFree(&vv);
  VectorFree(&sr);

  return (NO_ERROR);
}

/* extract the RASToVoxel Matrix */
MATRIX *GetSurfaceRASToVoxelMatrix(MRI *mri)
{
  return voxelFromSurfaceRAS_(mri);
}

int MRIsurfaceRASToVoxel(MRI *mri, double xr, double yr, double zr,
                         double *xv, double *yv, double *zv)
{
  MATRIX *voxelFromSRAS;
  static VECTOR *sr = NULL, *vv = NULL;

  voxelFromSRAS=voxelFromSurfaceRAS_(mri);
  if (sr == NULL)
    sr = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(sr, xr, yr, zr, 1.);
  vv = MatrixMultiply(voxelFromSRAS, sr, vv);
  *xv = V3_X(vv);
  *yv = V3_Y(vv);
  *zv = V3_Z(vv);

  MatrixFree(&voxelFromSRAS);
  //  VectorFree(&sr);
  //  VectorFree(&vv);

  return (NO_ERROR);
}

// same as above, but don't free matrix. Won't work if mri is changing
int
MRIsurfaceRASToVoxelCached(MRI *mri, double xr, double yr, double zr,
                           double *xv, double *yv, double *zv)
{
  static MATRIX *voxelFromSRAS = NULL;
  static VECTOR *sr = NULL, *vv = NULL;

  if (voxelFromSRAS == NULL)
  {
    voxelFromSRAS=voxelFromSurfaceRAS_(mri);
    sr = VectorAlloc(4, MATRIX_REAL);
  }
  V4_LOAD(sr, xr, yr, zr, 1.);
  vv = MatrixMultiply(voxelFromSRAS, sr, vv);
  *xv = V3_X(vv);
  *yv = V3_Y(vv);
  *zv = V3_Z(vv);

  return (NO_ERROR);
}

/*------------------------------------------------------*/
int
MRIworldToVoxel(MRI *mri, double xw, double yw, double zw,
                double *pxv, double *pyv, double *pzv)
{
  /*
    These internal workspaces are now static.
    They will 'leak' in that they will not be freed at exit
    They also contribute to the further destruction and
    annihilation of thread safety
  */
  static VECTOR *vv = NULL;
  static VECTOR *vw = NULL;
  MATRIX *IfromR;

  if( vw == NULL ) {
    vw = VectorAlloc(4, MATRIX_REAL);
  }
  if( vv == NULL ) {
    vv = VectorAlloc(4, MATRIX_REAL);
  }

  // if transform is not cached yet, then

  if (!mri->r_to_i__) {
    mri->r_to_i__ = extract_r_to_i(mri);
  }

  if( !mri->i_to_r__ ) {
    MATRIX *tmp = extract_i_to_r(mri);
    AffineMatrixAlloc( &(mri->i_to_r__) );
    SetAffineMatrix( mri->i_to_r__, tmp );
    MatrixFree( &tmp );
  }

  IfromR = mri->r_to_i__;
 
  V4_LOAD(vw, xw, yw, zw, 1.) ;
  MatrixMultiply(IfromR, vw, vv) ;
  *pxv = V3_X(vv);
  *pyv = V3_Y(vv);
  *pzv = V3_Z(vv);

#if 0
  // Memory leaked here
  VectorFree(&vv);
  VectorFree(&vw);
#endif

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIinitHeader(MRI *mri)
{
  mri->ptype = 2 ;

  /* most of these are in mm */
  mri->imnr0 = 1 ;
  mri->imnr1 = mri->depth ;
  mri->fov = mri->width ;
  mri->thick = 1.0 ;
  mri->xsize = 1.0 ;
  mri->ysize = 1.0 ;
  mri->zsize = 1.0 ;
  mri->ps = 1 ;
  mri->xstart = -mri->width/2.0 ;
  mri->xend = mri->width/2.0 ;
  mri->ystart = -mri->height/2.0 ;
  mri->yend = mri->height/2.0 ;
  mri->zstart = -mri->depth/2.0 ;
  mri->zend = mri->depth/2 ;
  // coronal
  mri->x_r = -1;
  mri->x_a = 0.;
  mri->x_s = 0.;
  //
  mri->y_r = 0.;
  mri->y_a = 0.;
  mri->y_s = -1;
  //
  mri->z_r = 0.;
  mri->z_a = 1;
  mri->z_s = 0.;
  //
  mri->c_r = mri->c_a = mri->c_s = 0.0;

  mri->ras_good_flag = 1;
  mri->gdf_image_stem[0] = '\0';
  mri->tag_data = NULL;
  mri->tag_data_size = 0;

  if (!mri->i_to_r__) {
    MATRIX *tmp = extract_i_to_r(mri);
    AffineMatrixAlloc( &(mri->i_to_r__) );
    SetAffineMatrix( mri->i_to_r__, tmp );
    MatrixFree( &tmp );
  }
    
  if (!mri->r_to_i__)
    mri->r_to_i__ = extract_r_to_i(mri);

  return(NO_ERROR) ;
}

/**
 * MRIreInitCache
 *
 * @param mri MRI* whose header information was modified
 *
 * @return NO_ERROR
 */
int MRIreInitCache(MRI *mri)
{
  MATRIX *tmp;

  AffineMatrixFree( &(mri->i_to_r__) );
  AffineMatrixAlloc( &(mri->i_to_r__) );
  tmp = extract_i_to_r(mri);
  SetAffineMatrix( mri->i_to_r__, tmp );
  MatrixFree( &tmp );

  if (mri->r_to_i__)
  {
    MatrixFree(&mri->r_to_i__);
    mri->r_to_i__ = 0;
  }
  mri->r_to_i__ = extract_r_to_i(mri);

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  change the direction of slices
  ------------------------------------------------------*/
MRI *
MRIextract(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
           int dx, int dy, int dz)
{
  return(MRIextractInto(mri_src, mri_dst, x0, y0, z0, dx, dy, dz, 0, 0, 0)) ;
}
MRI *
MRIcopyFrames(MRI *mri_src, MRI *mri_dst, int src_start_frame, int src_end_frame, int dst_start_frame)
{
  int  fno, offset ;

  if (mri_dst == NULL)
  {
    mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, mri_src->type,
          src_end_frame-src_start_frame+1) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  offset = dst_start_frame-src_start_frame ;

  for (fno = src_start_frame ; fno <= src_end_frame ; fno++)
  {
    mri_dst = MRIcopyFrame(mri_src, mri_dst, fno, fno+offset) ;
    if(mri_dst == NULL) return(NULL);
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Extract a cubic region of an MR image and return it to the caller
  ------------------------------------------------------*/
MRI *
MRIextractRegion(MRI *mri_src, MRI *mri_dst, MRI_REGION *region)
{
  MRI_REGION box ;

  if (region == NULL)
  {
    region = & box ;
    box.x = box.y = box.z = 0 ;
    box.dx = mri_src->width ;
    box.dy = mri_src->height ;
    box.dz = mri_src->depth ;
  }
  return(MRIextractInto(mri_src, mri_dst, region->x, region->y, region->z,
                        region->dx, region->dy, region->dz, 0, 0, 0)) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Extract a cubic region of an MR image and return it to the caller
  ------------------------------------------------------*/
MRI *
MRIextractIntoRegion(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
                     MRI_REGION *region)
{
  return(MRIextractInto(mri_src, mri_dst, x0, y0, z0, region->dx, region->dy,
                        region->dz, region->x, region->y, region->z)) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Extract a cubic region of an MR image and return it to the caller
  ------------------------------------------------------*/
MRI *
MRIextractInto(MRI *mri_src, MRI *mri_dst, int x0, int y0, int z0,
               int dx, int dy, int dz, int x1, int y1, int z1)
{
  int  width, height, depth, ys, zs, yd, zd, bytes, frame, xsize,ysize,zsize,
  dst_alloced = 0 ;
  double c_r, c_a, c_s;

  width = mri_src->width ;
  depth = mri_src->depth ;
  height = mri_src->height ;

  if (z0 >= depth || y0 >= height || x0 >= width)
    ErrorPrintf(ERROR_BADPARM,
                "MRIextractInto: bad src location (%d, %d, %d)", x0,y0,z0);
  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(dx, dy, dz, mri_src->type, mri_src->nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->imnr0 = z0 + mri_src->imnr0 - z1 ;
    mri_dst->imnr1 = mri_dst->imnr0 + dz - 1 ;
    dst_alloced = 1 ;
  }

  if (mri_src->type != mri_dst->type)
  {
    MRIfree(&mri_dst) ;
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MRIextractInto: src and dst types must match"));
  }

  // validation
  if (x0 < 0)
    x0 = 0 ;
  if (y0 < 0)
    y0 = 0 ;
  if (z0 < 0)
    z0 = 0 ;
  if (x0+dx > width)
    dx = (width - x0) ;
  if (y0+dy > height)
    dy = (height - y0) ;
  if (z0+dz > depth)
    dz = (depth - z0) ;
  if (x1 < 0)
    x1 = 0 ;
  if (y1 < 0)
    y1 = 0 ;
  if (z1 < 0)
    z1 = 0 ;

  if (x1+dx > mri_dst->width)
    dx = (mri_dst->width - x1) ;
  if (y1+dy > mri_dst->height)
    dy = (mri_dst->height - y1) ;
  if (z1+dz > mri_dst->depth)
    dz = (mri_dst->depth - z1) ;

  xsize = mri_src->xsize ;
  ysize = mri_src->ysize ;
  zsize = mri_src->zsize ;

  if (dst_alloced)
  {
    mri_dst->xstart += x0*xsize ;
    mri_dst->xend = mri_dst->xstart + dx*xsize ;
    mri_dst->ystart += y0*ysize ;
    mri_dst->yend = mri_dst->ystart + dy*ysize ;
    mri_dst->zstart += z0*zsize ;
    mri_dst->zend = mri_dst->zstart + dz*zsize  ;
  }

  bytes = dx ;
  switch (mri_src->type)
  {
  default:
    ErrorExit(ERROR_UNSUPPORTED, "MRIextractInto: unsupported source type %d", mri_src->type) ;
    break ;
  case MRI_FLOAT:
    bytes *= sizeof(float) ;
    break ;
  case MRI_LONG:
    bytes *= sizeof(long) ;
    break ;
  case MRI_INT:
    bytes *= sizeof(int) ;
    break ;
  case MRI_SHORT:
    bytes *= sizeof(short) ;
    break ;
  case MRI_UCHAR: break ;
  }

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (zd = z1, zs = z0 ; zs < z0+dz ; zs++, zd++)
    {
      for (yd = y1, ys = y0 ; ys < y0+dy ; ys++, yd++)
      {
        switch (mri_src->type)
        {
        default:
          ErrorExit(ERROR_UNSUPPORTED, "MRIextractInto: unsupported source type %d", mri_src->type) ;
          break ;
        case MRI_UCHAR:
          memmove(&MRIseq_vox(mri_dst, x1, yd, zd,frame),
                 &MRIseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        case MRI_FLOAT:
          memmove(&MRIFseq_vox(mri_dst, x1, yd, zd,frame),
                 &MRIFseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        case MRI_SHORT:
          memmove(&MRISseq_vox(mri_dst, x1, yd, zd,frame),
                 &MRISseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        case MRI_LONG:
          memmove(&MRILseq_vox(mri_dst, x1, yd, zd,frame),
                 &MRILseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        case MRI_INT:
          memmove(&MRIIseq_vox(mri_dst, x1, yd, zd,frame),
                 &MRIIseq_vox(mri_src,x0,ys,zs,frame), bytes);
          break ;
        }
      }
    }
  }
  mri_dst->xsize = mri_src->xsize ;
  mri_dst->ysize = mri_src->ysize ;
  mri_dst->zsize = mri_src->zsize ;
  mri_dst->thick = mri_src->thick ;
  mri_dst->ps = mri_src->ps ;
  MRIcopyPulseParameters(mri_src, mri_dst) ;
  // calculate c_ras
  MRIcalcCRASforExtractedVolume
  (mri_src, mri_dst, x0, y0, z0, x1, y1, z1, &c_r, &c_a, &c_s);
  mri_dst->c_r = c_r;
  mri_dst->c_a = c_a;
  mri_dst->c_s = c_s;
  // initialize cached transform
  MRIreInitCache(mri_dst);

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  change the direction of slices
  ------------------------------------------------------*/
MRI *
MRIreslice(MRI *mri_src, MRI *mri_dst, int slice_direction)
{
  int     width, height, depth, x1, x2, x3 ;
  BUFTYPE *psrc, val, *pdst ;

  int src_slice_direction = getSliceDirection(mri_src);
  if (slice_direction == src_slice_direction)
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
    return(mri_dst) ;
  }

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;


  if ((src_slice_direction == MRI_SAGITTAL &&
       slice_direction == MRI_CORONAL) ||
      (src_slice_direction == MRI_CORONAL &&
       slice_direction == MRI_SAGITTAL))
  {
    /*
     coronal images are back to front of the head, thus the depth axis
     points from the nose to the back of the head, with x from neck to
     crown, and y from ear to ear.
    */
    /* x1 --> x3
      x2 --> x2
      x3 --> x1
    */
    if (!mri_dst)
    {
      mri_dst = MRIalloc(depth, height, width, mri_src->type) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else if (depth != mri_dst->width || height != mri_dst->height ||
             width != mri_dst->depth)
    {
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "MRIreslice: invalid destination size (%d, %d, %d)",
                   mri_dst->width, mri_dst->height, mri_dst->depth)) ;
    }

    for (x3 = 0 ; x3 < depth ; x3++)
    {
      for (x2 = 0 ; x2 < height ; x2++)
      {
        psrc = &MRIvox(mri_src, 0, x2, x3) ;
        for (x1 = 0 ; x1 < width ; x1++)
        {
          /* swap so in place transformations are possible */
          mri_dst->slices[x1][x2][x3] = *psrc++ ;
        }
      }
    }
  }
  else
    if ((src_slice_direction == MRI_HORIZONTAL &&
         slice_direction == MRI_CORONAL) ||
        (src_slice_direction == MRI_CORONAL &&
         slice_direction == MRI_HORIZONTAL))
    {
      /*
       horizontal images are top to bottom of the head, thus the depth axis
       points from the top of the head to the neck, with x from ear to ear
       and y from nose to back of head.
      */
      /* x3 --> x2
        x2 --> x3
        x1 --> x1
      */
      if (!mri_dst)
      {
        mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
        MRIcopyHeader(mri_src, mri_dst) ;
      }
      else if (depth != mri_dst->height || height != mri_dst->depth ||
               width != mri_dst->width)
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIreslice: invalid destination size (%d, %d, %d)",
                     mri_dst->width, mri_dst->height, mri_dst->depth)) ;

      for (x3 = 0 ; x3 < depth ; x3++)
      {
        for (x2 = 0 ; x2 < height ; x2++)
        {
          psrc = &MRIvox(mri_src, 0, x2, x3) ;
          pdst = &MRIvox(mri_dst, 0, x3, x2) ;
          for (x1 = 0 ; x1 < width ; x1++)
          {
            /* swap so in place transformations are possible */
            *pdst++ = *psrc++ ;
          }
        }
      }
    }
    else
      if ((src_slice_direction == MRI_SAGITTAL &&
           slice_direction == MRI_HORIZONTAL))
      {
        /*
         horizontal images are top to bottom of the head,
         thus the depth axis
         points from the top of the head to the neck, with x from
         ear to ear
         and y from nose to back of head.
        */
        /* x3 --> x2
          x1 --> x3
          x2 --> x1
        */
        if (!mri_dst)
        {
          mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
          MRIcopyHeader(mri_src, mri_dst) ;
        }
        else if (depth != mri_dst->height || height != mri_dst->depth ||
                 width != mri_dst->width)
          ErrorReturn(NULL,
                      (ERROR_BADPARM,
                       "MRIreslice: invalid destination size (%d, %d, %d)",
                       mri_dst->width, mri_dst->height, mri_dst->depth)) ;

        for (x3 = 0 ; x3 < depth ; x3++)
        {
          for (x2 = 0 ; x2 < height ; x2++)
          {
            psrc = &MRIvox(mri_src, 0, x2, x3) ;
            for (x1 = 0 ; x1 < width ; x1++)
            {
              /* swap so in place transformations are possible */
              mri_dst->slices[x2][x1][x3] = *psrc++ ;
            }
          }
        }
      }
      else
        if (src_slice_direction == MRI_HORIZONTAL &&
            slice_direction == MRI_SAGITTAL)
        {
          /*
           horizontal images are top to bottom of the head,
           thus the depth axis
           points from the top of the head to the neck,
           with x from ear to ear
           and y from nose to back of head.
          */
          /* x2 --> x3
            x3 --> x1
            x1 --> x2
          */
          if (!mri_dst)
          {
            mri_dst = MRIalloc(width, depth, height, mri_src->type) ;
            MRIcopyHeader(mri_src, mri_dst) ;
          }
          else if (depth != mri_dst->height || height != mri_dst->depth ||
                   width != mri_dst->width)
            ErrorReturn(NULL,
                        (ERROR_BADPARM,
                         "MRIreslice: invalid destination size (%d, %d, %d)",
                         mri_dst->width, mri_dst->height, mri_dst->depth)) ;

          for (x3 = 0 ; x3 < depth ; x3++)
          {
            for (x2 = 0 ; x2 < height ; x2++)
            {
              psrc = &MRIvox(mri_src, 0, x2, x3) ;
              for (x1 = 0 ; x1 < width ; x1++)
              {
                /* swap so in place transformations are possible */
                mri_dst->slices[x1][x3][x2] = *psrc++ ;
              }
            }
          }
        }
        else
          switch (src_slice_direction)
          {
          default:
            MRIfree(&mri_dst) ;
            ErrorReturn(NULL,
                        (ERROR_BADPARM,
                         "MRIreslice: mri_src unknown slice direction %d",
                         src_slice_direction)) ;
            break ;
          case MRI_CORONAL:
            /*
             coronal images are back to front of the head,
             thus the depth axis
             points from the nose to the back of the head,
             with x from neck to
             crown, and y from ear to ear.
            */
            switch (slice_direction)
            {
            case MRI_SAGITTAL:
              /* x1 --> x3
                x2 --> x2
                x3 --> x1
              */
              if (!mri_dst)
              {
                mri_dst = MRIalloc(depth, height, width, mri_src->type) ;
                MRIcopyHeader(mri_src, mri_dst) ;
              }
              else if (depth != mri_dst->width ||
                       height != mri_dst->height ||
                       width != mri_dst->depth)
                ErrorReturn
                (NULL,
                 (ERROR_BADPARM,
                  "MRIreslice: invalid destination size (%d, %d, %d)",
                  mri_dst->width, mri_dst->height, mri_dst->depth)) ;

              for (x3 = 0 ; x3 < depth ; x3++)
              {
                for (x2 = 0 ; x2 < height ; x2++)
                {
                  psrc = &MRIvox(mri_src, 0, x2, x3) ;
                  for (x1 = 0 ; x1 < width ; x1++)
                  {
                    /* swap so in place transformations
                      are possible */
                    val = *psrc++ ;
#if 0
                    mri_dst->slices[x3][x2][x1] =
                      mri_src->slices[x1][x2][x3] ;
#endif
                    mri_dst->slices[x1][x2][x3] = val ;
                  }
                }
              }
              break ;
            case MRI_HORIZONTAL:
              break ;
            }
            break ;
          case MRI_SAGITTAL:
            /*
             sagittal images are slices in
             the plane of the nose, with depth going
             from ear to ear.
            */
            break ;
          }
  setDirectionCosine(mri_dst, slice_direction);
  mri_dst->ras_good_flag = 0;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Set an MRI intensity values to 0
  ------------------------------------------------------*/
int
MRIclear(MRI *mri)
{
  int   width, depth, height, bytes, y, z, frame, nframes ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  nframes = mri->nframes ;
  bytes = width ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    bytes *= sizeof(unsigned char) ;
    break ;
  case MRI_BITMAP:
    bytes /= 8 ;
    break ;
  case MRI_FLOAT:
    bytes *= sizeof(float) ;
    break ;
  case MRI_LONG:
    bytes *= sizeof(long) ;
    break ;
  case MRI_INT:
    bytes *= sizeof(int) ;
    break ;
  case MRI_SHORT:
    bytes *= sizeof(short) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIclear: unsupported input type %d", mri->type)) ;
    break ;
  }

  for (frame = 0 ; frame < nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
        memset(mri->slices[z+frame*depth][y], 0, bytes) ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  find the principle components of a (binary) MRI. The
  eigenvectors are the columns of the matrix mEvectors, the
  eigenvalues are returned in the array evalues and the means
  in means (these last two must be three elements long.)
  ------------------------------------------------------*/
int
MRIcenterOfMass(MRI *mri,double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  long    npoints ;
  double  mx, my, mz, weight ;
  double    val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ;
  npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        MRIsampleVolumeType(mri, x,  y, z, &val, SAMPLE_NEAREST) ;
        if (val > threshold)
        {
          weight += val ;
          mx += (float)x*val ;
          my += (float)y*val ;
          mz += (float)z*val ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;


  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  find the principle components of a (binary) MRI. The
  eigenvectors are the columns of the matrix mEvectors, the
  eigenvalues are returned in the array evalues and the means
  in means (these last two must be three elements long.)
  ------------------------------------------------------*/
int
MRIprincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues,
                       double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, val ;
  long    npoints ;
  MATRIX  *mCov, *mX, *mXT, *mTmp ;
  double  mx, my, mz, weight ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIprincipleComponents: unsupported input type %d",
                 mri->type)) ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ;
  npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          weight += val ;
          mx += (float)x*val ;
          my += (float)y*val ;
          mz += (float)z*val ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;     /* zero-mean coordinate vector */
  mXT = NULL ;                              /* transpose of above */
  mTmp = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* tmp matrix for covariance */
  mCov = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* covariance matrix */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          mX->rptr[1][1] = ((float)x - mx)*val ;
          mX->rptr[2][1] = ((float)y - my)*val ;
          mX->rptr[3][1] = ((float)z - mz)*val ;
          mXT = MatrixTranspose(mX, mXT) ;
          mTmp = MatrixMultiply(mX, mXT, mTmp) ;
          mCov = MatrixAdd(mTmp, mCov, mCov) ;
        }
      }
    }
  }

  if (weight > 0)
    MatrixScalarMul(mCov, 1.0f/weight, mCov) ;

  MatrixEigenSystem(mCov, evalues, mEvectors) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  find the principle components of a (binary) MRI. The
  eigenvectors are the columns of the matrix mEvectors, the
  eigenvalues are returned in the array evalues and the means
  in means (these last two must be three elements long) of values
  low_thresh <= val <= hi_thresh
  ------------------------------------------------------*/
int
MRIprincipleComponentsRange(MRI *mri, MATRIX *mEvectors, float *evalues,
                            double *means, float low_thresh, float hi_thresh)
{
  int     width, height, depth, x, y, z ;
  long    npoints ;
  MATRIX  *mCov, *mX, *mXT, *mTmp ;
  double  mx, my, mz, weight ;
  float   val ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;

  mx = my = mz = weight = 0.0f ;
  npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri, x, y, z, 0) ;
        if (val >= low_thresh && val <= hi_thresh)
        {
          weight += val ;
          mx += (float)x*val ;
          my += (float)y*val ;
          mz += (float)z*val ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;     /* zero-mean coordinate vector */
  mXT = NULL ;                              /* transpose of above */
  mTmp = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* tmp matrix for covariance */
  mCov = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* covariance matrix */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIgetVoxVal(mri, x, y,z, 0) ;
        if (val >= low_thresh && val <= hi_thresh)
        {
          mX->rptr[1][1] = ((float)x - mx)*val ;
          mX->rptr[2][1] = ((float)y - my)*val ;
          mX->rptr[3][1] = ((float)z - mz)*val ;
          mXT = MatrixTranspose(mX, mXT) ;
          mTmp = MatrixMultiply(mX, mXT, mTmp) ;
          mCov = MatrixAdd(mTmp, mCov, mCov) ;
        }
      }
    }
  }

  if (weight > 0)
    MatrixScalarMul(mCov, 1.0f/weight, mCov) ;

  MatrixEigenSystem(mCov, evalues, mEvectors) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  find the principle components of a (binary) MRI. The
  eigenvectors are the columns of the matrix mEvectors, the
  eigenvalues are returned in the array evalues and the means
  in means (these last two must be three elements long.)
  ------------------------------------------------------*/
int
MRIbinaryPrincipleComponents(MRI *mri, MATRIX *mEvectors, float *evalues,
                             double *means, BUFTYPE threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, val ;
  long    npoints ;
  MATRIX  *mCov, *mX, *mXT, *mTmp ;
  double  mx, my, mz, weight ;

  if (mri->type != MRI_UCHAR)
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIprincipleComponents: unsupported input type %d",
                 mri->type)) ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  mx = my = mz = weight = 0.0f ;
  npoints = 0L ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          weight++ ;
          mx += (float)x ;
          my += (float)y ;
          mz += (float)z ;
          npoints++ ;
        }
      }
    }
  }

  if (weight > 0.0)
  {
    mx /= weight ;
    my /= weight  ;
    mz /= weight ;
    means[0] = mx ;
    means[1] = my ;
    means[2] = mz ;
  }
  else
    means[0] = means[1] = means[2] = 0.0f ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;     /* zero-mean coordinate vector */
  mXT = NULL ;                              /* transpose of above */
  mTmp = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* tmp matrix for covariance */
  mCov = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* covariance matrix */

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
        {
          mX->rptr[1][1] = ((float)x - mx) ;
          mX->rptr[2][1] = ((float)y - my) ;
          mX->rptr[3][1] = ((float)z - mz) ;
          mXT = MatrixTranspose(mX, mXT) ;
          mTmp = MatrixMultiply(mX, mXT, mTmp) ;
          mCov = MatrixAdd(mTmp, mCov, mCov) ;
        }
      }
    }
  }

  if (weight > 0)
    MatrixScalarMul(mCov, 1.0f/weight, mCov) ;

  MatrixEigenSystem(mCov, evalues, mEvectors) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  threshold an MRI.
  ------------------------------------------------------*/
MRI *
MRIthresholdRangeInto(MRI *mri_src,MRI *mri_dst,BUFTYPE low_val,BUFTYPE hi_val)
{
  int     width, height, depth, x, y, z, f ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (f = 0 ; f < mri_src->nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          if (val <  low_val || val > hi_val)
            val = 0 ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  threshold an MRI -- only 1st frame!
  ------------------------------------------------------*/
MRI *
MRIthreshold(MRI *mri_src, MRI *mri_dst, float threshold)
{
  int     width, height, depth, x, y, z, f ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (f = 0 ; f < mri_src->nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          if (val < threshold)
            val = 0 ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  threshold an MRI considering all the frames
  ------------------------------------------------------*/
MRI *
MRIthresholdAllFrames(MRI *mri_src, MRI *mri_dst, float threshold)
{
  int     frame, width, height, depth, x, y, z, f, n ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  frame = mri_src->nframes ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
	    {
	      n = 0;
	      for (f = 0 ; f < frame ; f++) 
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          if (val < threshold) n++;
        }
	      if(n > 0) 
        { 
          val = 0;
          for (f = 0 ; f < frame ; f++) 
            MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
        }
	      else
        { 
          for (f = 0 ; f < frame ; f++) 
          {
            val = MRIgetVoxVal(mri_src, x, y, z, f);
            MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
          }
        }
	    }
    }
  }
  
  return(mri_dst) ;
}
MRI *
MRIupperthresholdAllFrames(MRI *mri_src, MRI *mri_dst, float threshold)
{
  int     frame, width, height, depth, x, y, z, f, n ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  frame = mri_src->nframes ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
	    {
	      n = 0;
	      for (f = 0 ; f < frame ; f++) 
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          if (val > threshold) n++;
        }
	      if(n > 0) 
        { 
          val = threshold;
          for (f = 0 ; f < frame ; f++) 
            MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
        }
	      else
        { 
          for (f = 0 ; f < frame ; f++) 
          {
            val = MRIgetVoxVal(mri_src, x, y, z, f);
            MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
          }
        }
	    }
    }
  }
  
  return(mri_dst) ;
}
// Only threshold the specified frame
MRI *
MRIthresholdFrame(MRI *mri_src, MRI *mri_dst, float threshold, int frame)
{
  int     f, width, height, depth, x, y, z, w ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  f = mri_src->nframes ;

  if (frame > f+1) 
    ErrorReturn(NULL,
                (ERROR_BADPARM,"MRIthreshold: invalid frame number"));
  
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for ( w = 0; w < f ; w++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, w) ;
          if (w == frame-1)
          {	    
            if (val < threshold)
              val = 0 ;
            MRIsetVoxVal(mri_dst, x, y, z, frame-1, val) ;
          }
          else 
            MRIsetVoxVal(mri_dst, x, y, z, w, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}
MRI *
MRIupperthresholdFrame(MRI *mri_src, MRI *mri_dst, float threshold, int frame)
{
  int     f, width, height, depth, x, y, z, w ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  f = mri_src->nframes ;

  if (frame > f+1) 
    ErrorReturn(NULL,
                (ERROR_BADPARM,"MRIthreshold: invalid frame number"));
  
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for ( w = 0; w < f ; w++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, w) ;
          if (w == frame-1)
          {	    
            if (val > threshold)
              val = threshold ;
            MRIsetVoxVal(mri_dst, x, y, z, frame-1, val) ;
          }
          else 
            MRIsetVoxVal(mri_dst, x, y, z, w, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}

// Threshold all frames using one frame
MRI *
MRIthresholdByFrame(MRI *mri_src, MRI *mri_dst, float threshold, int frame)
{
  int     f, width, height, depth, x, y, z, w ;
  float   val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  f = mri_src->nframes ;

  if (frame > f+1) 
    ErrorReturn(NULL,
                (ERROR_BADPARM,"MRIthreshold: invalid frame number"));
  
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for ( w = 0; w < f ; w++)
        {
          //	    val = MRIgetVoxVal(mri_src, x, y, z, w) ;
          val = MRIgetVoxVal(mri_src, x, y, z, frame-1) ;
          if (val < threshold)
            val = 0 ;
          if (w == frame-1)
          {	    
		
            MRIsetVoxVal(mri_dst, x, y, z, frame-1, val) ;
          }
          else 
            MRIsetVoxVal(mri_dst, x, y, z, w, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  threshold an MRI.
  ------------------------------------------------------*/
MRI  *
MRIinvertContrast(MRI *mri_src, MRI *mri_dst, float threshold)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *psrc, *pdst, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        val = *psrc++ ;
        if (val > threshold)
          val = 255 - val ;
        *pdst++ = val ;
      }
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  threshold an MRI.
  ------------------------------------------------------*/
MRI *
MRIbinarize(MRI *mri_src, MRI *mri_dst, float threshold, float low_val,
            float hi_val)
{
  int     width, height, depth, x, y, z, f ;
  double    val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  for (f = 0 ; f < mri_src->nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
#if 0
          MRIsampleVolumeFrameType
          (mri_src, x, y, z, f, SAMPLE_NEAREST, &val) ;
#endif
          if (val < threshold)
            val = low_val ;
          else
            val = hi_val ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}

/*-----------------------------------------------------*/
MRI *MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     nframes, width, height, depth, x, y, z, f, s ;
  float v1, v2, v=0.0;
  BUFTYPE *p1=NULL, *p2=NULL, *pdst=NULL ;
  short *ps1=NULL, *ps2=NULL, *psdst=NULL ;
  int   *pi1=NULL, *pi2=NULL, *pidst=NULL ;
  long  *pl1=NULL, *pl2=NULL, *pldst=NULL ;
  float *pf1=NULL, *pf2=NULL, *pfdst=NULL ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;
  nframes = mri1->nframes ;
  if (nframes == 0) nframes = 1;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width, height, depth, mri1->type,nframes) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  if (mri1->type != mri2->type)
  {
    /* Generic but slow */
    for (f = 0 ; f < nframes ; f++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            v1 = MRIgetVoxVal(mri1,x,y,z,f);
            v2 = MRIgetVoxVal(mri2,x,y,z,f);
            v = v1-v2;
            MRIsetVoxVal(mri_dst,x,y,z,f,v);
          }
        }
      }
    }
    return(mri_dst) ;
  }


  s = 0;
  for (f = 0 ; f < nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {

        switch (mri_dst->type)
        {
        case MRI_UCHAR:
          pdst  =           mri_dst->slices[s][y] ;
          break;
        case MRI_SHORT:
          psdst = (short *) mri_dst->slices[s][y] ;
          break;
        case MRI_INT:
          pidst = (int *)   mri_dst->slices[s][y] ;
          break;
        case MRI_LONG:
          pldst = (long *)  mri_dst->slices[s][y] ;
          break;
        case MRI_FLOAT:
          pfdst = (float *) mri_dst->slices[s][y] ;
          break;
        }
        switch (mri1->type)
        {
        case MRI_UCHAR:
          p1 = mri1->slices[s][y] ;
          p2 = mri2->slices[s][y] ;
          break;
        case MRI_SHORT:
          ps1 = (short *) mri1->slices[s][y] ;
          ps2 = (short *) mri2->slices[s][y] ;
          break;
        case MRI_INT:
          pi1 = (int *) mri1->slices[s][y] ;
          pi2 = (int *) mri2->slices[s][y] ;
          break;
        case MRI_LONG:
          pl1 = (long *) mri1->slices[s][y] ;
          pl2 = (long *) mri2->slices[s][y] ;
          break;
        case MRI_FLOAT:
          pf1 = (float *) mri1->slices[s][y] ;
          pf2 = (float *) mri2->slices[s][y] ;
          break;
        }

        for (x = 0 ; x < width ; x++)
        {

          switch (mri1->type)
          {
          case MRI_UCHAR:
            v = (float)(*p1++)  - (float)(*p2++);
            break;
          case MRI_SHORT:
            v = (float)(*ps1++) - (float)(*ps2++);
            break;
          case MRI_INT:
            v = (float)(*pi1++) - (float)(*pi2++);
            break;
          case MRI_LONG:
            v = (float)(*pl1++) - (float)(*pl2++);
            break;
          case MRI_FLOAT:
            v = (float)(*pf1++) - (float)(*pf2++);
            break;
          }

          switch (mri_dst->type)
          {
          case MRI_UCHAR:
            (*pdst++)  = (BUFTYPE) nint(v);
            break;
          case MRI_SHORT:
            (*psdst++) = (short)   nint(v);
            break;
          case MRI_INT:
            (*pidst++) = (int)     nint(v);
            break;
          case MRI_LONG:
            (*pldst++) = (long)    nint(v);
            break;
          case MRI_FLOAT:
            (*pfdst++) = (float)    v;
            break;
          }

        }
      }
      s++;
    }
  }

  return(mri_dst) ;
}
#if 0
/*------------------------------------------------------
  MRIsubtract(mri1,mri2,mridiff) - computes mri1-mri2.
  ------------------------------------------------------*/
MRI *MRIsubtract(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     nframes, width, height, depth, x, y, z, f ;
  float v1, v2, vdiff;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;
  nframes = mri1->nframes ;
  if (nframes == 0) nframes = 1;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width, height, depth, mri1->type,nframes) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for (f = 0 ; f < nframes ; f++)
        {
          v1 = MRIgetVoxVal(mri1,x,y,z,f);
          v2 = MRIgetVoxVal(mri2,x,y,z,f);
          vdiff = v1-v2;
          MRIsetVoxVal(mri_dst,x,y,z,f,vdiff);
        }
      }
    }
  }
  return(mri_dst) ;
}
#endif

MRI * MRIabsdiff(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst, v1, v2 ;
  float   f1, f2 ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  if (mri1->type == MRI_UCHAR && mri2->type == MRI_UCHAR)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        p1 = mri1->slices[z][y] ;
        p2 = mri2->slices[z][y] ;
        pdst = mri_dst->slices[z][y] ;
        for (x = 0 ; x < width ; x++)
        {
          v1 = *p1++ ;
          v2 = *p2++ ;
          if (v1 > v2)
            *pdst++ = v1 - v2 ;
          else
            *pdst++ = v2 - v1 ;
        }
      }
    }
  }
  else
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          f1 = MRIgetVoxVal(mri1, x, y, z, 0)  ;
          f2 = MRIgetVoxVal(mri2, x, y, z, 0)  ;
          MRIsetVoxVal(mri_dst, x, y, z, 0, fabs(f1 - f2)) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*!
  \fn MRI *MRIabs(MRI *mri_src, MRI *mri_dst)
  \brief Computes the abs() of each voxel
*/
MRI *MRIabs(MRI *mri_src, MRI *mri_dst)
{
  int   width, height, depth, nframes, x, y, z,f ;
  float val;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  nframes = mri_src->nframes;

  if (!mri_dst)  {
    mri_dst = MRIallocSequence(width, height, depth, mri_src->type, nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)  {
    for (y = 0 ; y < height ; y++)    {
      for (x = 0 ; x < width ; x++)      {
        for (f = 0 ; f < nframes ; f++)        {
          val = fabs(MRIgetVoxVal(mri_src,x,y,z,f));
          MRIsetVoxVal(mri_dst,x,y,z,f,val);
        }
      }
    }
  }
  return(mri_dst) ;
}
/*!
  \fn MRI *MRIpos(MRI *mri_src, MRI *mri_dst)
  \brief If a voxel is negative, sets it to 0.
*/
MRI *MRIpos(MRI *mri_src, MRI *mri_dst)
{
  int   width, height, depth, nframes, x, y, z,f ;
  float val;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  nframes = mri_src->nframes;

  if (!mri_dst){
    mri_dst = MRIallocSequence(width, height, depth, mri_src->type, nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)  {
    for (y = 0 ; y < height ; y++)    {
      for (x = 0 ; x < width ; x++)      {
        for (f = 0 ; f < nframes ; f++) {
          val = MRIgetVoxVal(mri_src,x,y,z,f);
          if(val < 0.0) val = 0.0;
          MRIsetVoxVal(mri_dst,x,y,z,f,val);
        }
      }
    }
  }
  return(mri_dst) ;
}
/*!
  \fn MRI *MRIpos(MRI *mri_src, MRI *mri_dst)
  \brief If a voxel is postive, sets it to 0.
*/
MRI *MRIneg(MRI *mri_src, MRI *mri_dst)
{
  int   width, height, depth, nframes, x, y, z,f ;
  float val;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  nframes = mri_src->nframes;

  if (!mri_dst){
    mri_dst = MRIallocSequence(width, height, depth, mri_src->type, nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)  {
    for (y = 0 ; y < height ; y++)    {
      for (x = 0 ; x < width ; x++)      {
        for (f = 0 ; f < nframes ; f++) {
          val = MRIgetVoxVal(mri_src,x,y,z,f);
          if(val > 0.0) val = 0.0;
          MRIsetVoxVal(mri_dst,x,y,z,f,val);
        }
      }
    }
  }
  return(mri_dst) ;
}


/*-----------------------------------------------------*/
MRI * MRIadd(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     nframes, width, height, depth, x, y, z, f, s ;
  float v1, v2, v=0.0;
  BUFTYPE *p1=NULL, *p2=NULL, *pdst=NULL ;
  short *ps1=NULL, *ps2=NULL, *psdst=NULL ;
  int   *pi1=NULL, *pi2=NULL, *pidst=NULL ;
  long  *pl1=NULL, *pl2=NULL, *pldst=NULL ;
  float *pf1=NULL, *pf2=NULL, *pfdst=NULL ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;
  nframes = mri1->nframes ;
  if (nframes == 0) nframes = 1;

  if (!mri_dst)
  {
    mri_dst = MRIallocSequence(width, height, depth, mri1->type,nframes) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  if (mri1->type == MRI_UCHAR || (mri1->type != mri2->type))
  {
    /* Generic but slow */
    for (f = 0 ; f < nframes ; f++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            v1 = MRIgetVoxVal(mri1,x,y,z,f);
            v2 = MRIgetVoxVal(mri2,x,y,z,f);
            v = v1+v2;
            if (mri_dst->type == MRI_UCHAR && v > 255)
              v = 255 ;
            if (mri_dst->type == MRI_UCHAR && v < 0)
              v = 0 ;
            MRIsetVoxVal(mri_dst,x,y,z,f,v);
          }
        }
      }
    }
    return(mri_dst) ;
  }


  s = 0;
  for (f = 0 ; f < nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {

        switch (mri_dst->type)
        {
        case MRI_UCHAR:
          pdst  = mri_dst->slices[s][y] ;
          break;
        case MRI_SHORT:
          psdst = (short *) mri_dst->slices[s][y] ;
          break;
        case MRI_INT:
          pidst = (int *)   mri_dst->slices[s][y] ;
          break;
        case MRI_LONG:
          pldst = (long *)  mri_dst->slices[s][y] ;
          break;
        case MRI_FLOAT:
          pfdst = (float *) mri_dst->slices[s][y] ;
          break;
        }

        switch (mri1->type)
        {
        case MRI_UCHAR:
          p1 = mri1->slices[s][y] ;
          p2 = mri2->slices[s][y] ;
          break;
        case MRI_SHORT:
          ps1 = (short *) mri1->slices[s][y] ;
          ps2 = (short *) mri2->slices[s][y] ;
          break;
        case MRI_INT:
          pi1 = (int *) mri1->slices[s][y] ;
          pi2 = (int *) mri2->slices[s][y] ;
          break;
        case MRI_LONG:
          pl1 = (long *) mri1->slices[s][y] ;
          pl2 = (long *) mri2->slices[s][y] ;
          break;
        case MRI_FLOAT:
          pf1 = (float *) mri1->slices[s][y] ;
          pf2 = (float *) mri2->slices[s][y] ;
          break;
        }

        for (x = 0 ; x < width ; x++)
        {

          /* note, the code below does not check for over/underflow! */
          switch (mri_dst->type)
          {
          case MRI_UCHAR:
            switch (mri1->type)
            {
            case MRI_UCHAR:
              (*pdst++) = (BUFTYPE) ((*p1++)+(*p2++));
              break;
            case MRI_SHORT:
              (*pdst++) = (BUFTYPE) ((*ps1++)+(*ps2++));
              break;
            case MRI_INT:
              (*pdst++) = (BUFTYPE) ((*pi1++)+(*pi2++));
              break;
            case MRI_LONG:
              (*pdst++) = (BUFTYPE) ((*pl1++)+(*pl2++));
              break;
            case MRI_FLOAT:
              (*pdst++) = (BUFTYPE) nint((*pf1++)+(*pf2++));
              break;
            }
            break;
          case MRI_SHORT:
            switch (mri1->type)
            {
            case MRI_UCHAR:
              (*psdst++) = ((short)(*p1++)+(*p2++));
              break;
            case MRI_SHORT:
              (*psdst++) = (short) ((*ps1++)+(*ps2++));
              break;
            case MRI_INT:
              (*psdst++) = (short) ((*pi1++)+(*pi2++));
              break;
            case MRI_LONG:
              (*psdst++) = (short) ((*pl1++)+(*pl2++));
              break;
            case MRI_FLOAT:
              (*psdst++) = (short) nint((*pf1++)+(*pf2++));
              break;
            }
            break;
          case MRI_INT:
            switch (mri1->type)
            {
            case MRI_UCHAR:
              (*pidst++) = ((int)(*p1++)+(*p2++));
              break;
            case MRI_SHORT:
              (*pidst++) = ((int)(*ps1++)+(*ps2++));
              break;
            case MRI_INT:
              (*pidst++) = (int) ((*pi1++)+(*pi2++));
              break;
            case MRI_LONG:
              (*pidst++) = (int) ((*pl1++)+(*pl2++));
              break;
            case MRI_FLOAT:
              (*pidst++) = (int) nint((*pf1++)+(*pf2++));
              break;
            }
            break;
          case MRI_LONG:
            switch (mri1->type)
            {
            case MRI_UCHAR:
              (*pldst++) = ((long)(*p1++)+(*p2++));
              break;
            case MRI_SHORT:
              (*pldst++) = ((long)(*ps1++)+(*ps2++));
              break;
            case MRI_INT:
              (*pldst++) = ((long)(*pi1++)+(*pi2++));
              break;
            case MRI_LONG:
              (*pldst++) = (long) ((*pl1++)+(*pl2++));
              break;
            case MRI_FLOAT:
              (*pldst++) = (long) nint((*pf1++)+(*pf2++));
              break;
            }
            break;
          case MRI_FLOAT:
            switch (mri1->type)
            {
            case MRI_UCHAR:
              (*pfdst++) = ((float)(*p1++)+(*p2++));
              break;
            case MRI_SHORT:
              (*pfdst++) = ((float)(*ps1++)+(*ps2++));
              break;
            case MRI_INT:
              (*pfdst++) = ((float)(*pi1++)+(*pi2++));
              break;
            case MRI_LONG:
              (*pfdst++) = ((float)(*pl1++)+(*pl2++));
              break;
            case MRI_FLOAT:
              (*pfdst++) =         (*pf1++)+(*pf2++);
              break;
            }
            break;
          }
        }
      }
      s++;
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------------
  MRIaverage() - computes average of source and destination.
  ------------------------------------------------------*/
MRI *
MRIaverage(MRI *mri_src, int dof, MRI *mri_dst)
{
  int     width, height, depth, x, y, z, f ;
  double    src, dst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  if (!MRIcheckSize(mri_src, mri_dst,0,0,0))
    ErrorReturn(NULL,
                (ERROR_BADPARM,"MRIaverage: incompatible volume dimensions"));
#if 0
  if ((mri_src->type != MRI_UCHAR) || (mri_dst->type != MRI_UCHAR))
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRISaverage: unsupported voxel format %d",mri_src->type));
#endif

  //  for (f = 0 ; f < mri_src->nframes ; f++)
  for (f = 0 ; f < mri_dst->nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          MRIsampleVolumeFrameType
          (mri_src, x,  y, z, f, SAMPLE_NEAREST, &src) ;
          MRIsampleVolumeFrameType
          (mri_dst, x,  y, z, f, SAMPLE_NEAREST, &dst) ;
          MRIsetVoxVal
          (mri_dst, x, y, z, f, (dst*dof+src)/(double)(dof+1))  ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  ------------------------------------------------------*/
MRI *
MRImultiply(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  float   f1, f2 ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        f1 = MRIgetVoxVal(mri1, x, y, z, 0) ;
        f2 = MRIgetVoxVal(mri2, x, y, z, 0) ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, f1*f2) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  ------------------------------------------------------*/
MRI *
MRIscaleAndMultiply(MRI *mri1, float scale, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;
  float   out_val ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      p1 = mri1->slices[z][y] ;
      p2 = mri2->slices[z][y] ;
      pdst = mri_dst->slices[z][y] ;
      for (x = 0 ; x < width ; x++)
      {
        out_val = *p1++ * (*p2++/scale) ;
        if (out_val > 255)
          out_val = 255 ;
        else if (out_val < 0)
          out_val = 0 ;
        *pdst++ = (BUFTYPE)nint(out_val) ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *MRIdivide(MRI *mri1, MRI *mri2, MRI *mri_dst)
{
  int     width, height, depth, x, y, z ;
  BUFTYPE *p1, *p2, *pdst ;

  width = mri1->width ;
  height = mri1->height ;
  depth = mri1->depth ;

  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, mri1->type) ;
    MRIcopyHeader(mri1, mri_dst) ;
  }

  if (mri1->type != MRI_UCHAR || mri2->type != MRI_UCHAR)
  {
    double val1, val2, dst ;

    for (z = 0 ; z < depth ; z++)    {
      for (y = 0 ; y < height ; y++)      {
        for (x = 0 ; x < width ; x++)        {
          if (x == Gx && y == Gy && z==Gz)
            DiagBreak() ;
          val1 = MRIgetVoxVal(mri1, x, y, z, 0) ;
          val2 = MRIgetVoxVal(mri2, x, y, z, 0) ;
          if(FZERO(val2)) dst = 0.0;
          else            dst = val1 / val2 ;
          if(abs(dst) > 1000)
            DiagBreak() ;
          MRIsetVoxVal(mri_dst, x, y, z, 0, dst) ;
        }
      }
    }
  }
  else   /* both UCHAR volumes */
  {
    for (z = 0 ; z < depth ; z++)    {
      for (y = 0 ; y < height ; y++)      {
        p1 = mri1->slices[z][y] ;
        p2 = mri2->slices[z][y] ;
        pdst = mri_dst->slices[z][y] ;
        for (x = 0 ; x < width ; x++)        {
          if  (!*p2)          {
            *pdst = FZERO(*p1) ? 0 : 255 ;
            p2++ ;
          }
          else {
            *pdst++ = *p1++ / *p2++ ;
	    //printf("%d / %d = %d\n", *pdst, *p1, *p2);
	  }
        }
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  MRIclone() - create a copy of an mri struct. Copies
  header info and allocs the pixel space (but does not
  copy pixel data).
  ------------------------------------------------------*/
MRI *MRIclone( const MRI *mri_src, MRI *mri_dst)
{
  if (!mri_dst)
    mri_dst =
      MRIallocSequence(mri_src->width, mri_src->height,mri_src->depth,
                       mri_src->type, mri_src->nframes);

  MRIcopyHeader(mri_src, mri_dst) ;
  return(mri_dst) ;
}
/*!
  \fn MRI *MRIcloneBySpace(MRI *mri_src, int type, int nframes)
  \brief Copies mri struct, allocs but not copy pixels (type,nframes)
  \param mri_src - source MRI struct
  \param type - clone will be of this type (-1 to use source)
  \param nframes - clone will have this many frames  (-1 to use source)
  \return New MRI struct.

  Does not copy pixel data.
*/
MRI *MRIcloneBySpace(MRI *mri_src, int type, int nframes)
{
  MRI *mri_dst;
  if (type < 0) type = mri_src->type;
  if (nframes < 0) nframes = mri_src->nframes;
  mri_dst = MRIallocSequence(mri_src->width, mri_src->height,
                             mri_src->depth, type, nframes);
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_dst->nframes = nframes;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Description
  Copy one MRI into another (including header info)
  ------------------------------------------------------*/
MRI *
MRIcloneRoi(MRI *mri_src, MRI *mri_dst)
{
  int  w, h, d ;

  w = mri_src->width - mri_src->roi.x ;
  h = mri_src->height - mri_src->roi.y ;
  d = mri_src->depth - mri_src->roi.z ;
  mri_dst = MRIallocSequence(w, h, d, MRI_FLOAT, mri_src->nframes) ;
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_dst->xstart = mri_src->xstart + mri_src->roi.x * mri_src->xsize ;
  mri_dst->ystart = mri_src->ystart + mri_src->roi.y * mri_src->ysize ;
  mri_dst->zstart = mri_src->zstart + mri_src->roi.z * mri_src->zsize ;
  mri_dst->xend = mri_src->xstart + w * mri_src->xsize ;
  mri_dst->yend = mri_src->ystart + h * mri_src->ysize ;
  mri_dst->zend = mri_src->zstart + d * mri_src->zsize ;
  return(mri_dst) ;
}
/*----------------------------------------------------------*/
/*!
  \fn int MRIchunk(MRI **pmri)
  \brief Change input MRI memory allocation to be "chunked".
  \return 0 on success, 1 if error
  Has no effect if input is already chuncked
*/
int MRIchunk(MRI **pmri)
{
  MRI *mritmp;
  if ((*pmri)->ischunked) return(0);
  printf("Chunking\n");
  mritmp = MRIallocChunk((*pmri)->width, (*pmri)->height, (*pmri)->depth,
                         (*pmri)->type, (*pmri)->nframes);
  if (mritmp == NULL) return(1);
  MRIcopy(*pmri, mritmp);
  MRIfree(pmri);
  *pmri = mritmp;
  return(0);
}


/*----------------------------------------------------------
  Copy one MRI into another (including header info and data)
  -----------------------------------------------------------*/
MRI *
MRIcopy(MRI *mri_src, MRI *mri_dst)
{
  int     width, height, depth, bytes, x, y, z, frame, val ;
  float   *fdst, *fsrc ;
  BUFTYPE *csrc, *cdst ;
  int     dest_ptype, *isrc;
  short   *ssrc, *sdst ;

  if (mri_src == mri_dst)
    return(mri_dst) ;
  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
  {
    if (mri_src->slices)
      mri_dst = MRIallocSequence(width, height, depth, mri_src->type,
                                 mri_src->nframes) ;
    else{
      mri_dst = MRIallocHeader(width, height, depth, mri_src->type, 1);
      mri_dst->nframes = mri_src->nframes ;
    }
  }
  dest_ptype = mri_dst->ptype;
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_dst->ptype = dest_ptype;

  if (!mri_src->slices)
    return(mri_dst);

  if (mri_src->type == mri_dst->type)
  {
    bytes = width ;
    switch (mri_src->type)
    {
    case MRI_UCHAR:
      bytes *= sizeof(BUFTYPE) ;
      break ;
    case MRI_SHORT:
      bytes *= sizeof(short);
      break;
    case MRI_FLOAT:
      bytes *= sizeof(float) ;
      break ;
    case MRI_INT:
      bytes *= sizeof(int) ;
      break ;
    case MRI_LONG:
      bytes *= sizeof(long) ;
      break ;
    }

    for (frame = 0 ; frame < mri_src->nframes ; frame++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          memmove(mri_dst->slices[z+frame*depth][y],
                 mri_src->slices[z+frame*depth][y], bytes) ;
        }
      }
    }
  }
  else
  {
    switch (mri_src->type)
    {
    case MRI_FLOAT:
      switch (mri_dst->type)
      {
      case MRI_SHORT: /* float --> short */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              sdst = &MRISseq_vox(mri_dst, 0, y, z, frame) ;
              fsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                val = nint(*fsrc++) ;
                *sdst++ = (short)val ;
              }
            }
          }
        }
        break ;
      case MRI_UCHAR:  /* float --> unsigned char */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              cdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
              fsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                val = nint(*fsrc++) ;
                if (val > 255)
                  val = 255 ;
                *cdst++ = (BUFTYPE)val ;
              }
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    case MRI_UCHAR:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* unsigned char --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              csrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
                *fdst++ = (float)*csrc++ ;
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    case MRI_SHORT:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* short --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              ssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                if (z == 113 && y == 143 && x == 161)
                  DiagBreak() ;
                *fdst++ = (float)*ssrc++ ;
              }
            }
          }
        }
        break ;
      case MRI_UCHAR:
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              cdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
              ssrc = &MRISseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
              {
                *cdst++ = (float)*ssrc++ ;
              }
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    case MRI_INT:
      switch (mri_dst->type)
      {
      case MRI_FLOAT:   /* unsigned char --> float */
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          for (z = 0 ; z < depth ; z++)
          {
            for (y = 0 ; y < height ; y++)
            {
              fdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
              isrc = &MRIIseq_vox(mri_src, 0, y, z, frame) ;
              for (x = 0 ; x < width ; x++)
                *fdst++ = (float)*isrc++ ;
            }
          }
        }
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIcopy: src type %d & dst type %d unsupported",
                     mri_src->type, mri_dst->type)) ;
        break ;
      }
      break ;
    default:
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "MRIcopy: src type %d & dst type %d unsupported",
                   mri_src->type, mri_dst->type)) ;
      break ;  /* in case someone removes the errorreturn */
    }
  }
  strcpy(mri_dst->fname,mri_src->fname);
  return(mri_dst) ;
}
/*
  make MAX_INDEX way larger than it has to be. This will give
  some headroom for bad (e.g. poorly registered) images without
  sacrificing too much space.
*/
#define MAX_INDEX    500
/*-----------------------------------------------------
  Allocate a lookup table that allows indices which
  are outside the image region.
  ------------------------------------------------------*/
int
MRIallocIndices(MRI *mri)
{
  int width, height, depth, i ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  mri->xi = (int *)calloc(width+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->xi)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIallocIndices: could not allocate %d elt index array",
              width+2*MAX_INDEX) ;
  mri->yi = (int *)calloc(height+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->yi)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIallocIndices: could not allocate %d elt index array",
              height+2*MAX_INDEX) ;
  mri->zi = (int *)calloc(depth+2*MAX_INDEX, sizeof(int)) ;
  if (!mri->zi)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIallocIndices: could not allocate %d elt index array",
              depth+2*MAX_INDEX) ;

  /*
    indexing into these arrays returns valid pixel indices from
    -MAX_INDEX to width+MAX_INDEX
  */
  mri->xi += MAX_INDEX ;
  mri->yi += MAX_INDEX ;
  mri->zi += MAX_INDEX ;
  for (i = -MAX_INDEX ; i < width+MAX_INDEX ; i++)
  {
    if (i <= 0)           mri->xi[i] = 0 ;
    else if (i >= width)  mri->xi[i] = width-1 ;
    else                  mri->xi[i] = i ;
  }
  for (i = -MAX_INDEX ; i < height+MAX_INDEX ; i++)
  {
    if (i <= 0)            mri->yi[i] = 0 ;
    else if (i >= height)  mri->yi[i] = height-1 ;
    else                    mri->yi[i] = i ;
  }
  for (i = -MAX_INDEX ; i < depth+MAX_INDEX ; i++)
  {
    if (i <= 0)           mri->zi[i] = 0 ;
    else if (i >= depth)  mri->zi[i] = depth-1 ;
    else                  mri->zi[i] = i ;
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
/*!
\fn MRI *MRIallocChunk(int width, int height, int depth, int type, int nframes)
\brief Alloc pixel data in MRI struct as one big buffer.
*/
MRI *MRIallocChunk(int width, int height, int depth, int type, int nframes)
{
  MRI *mri ;
  int  slice, row;
  void *p;

  if( sizeof(mri->bytes_total) != sizeof(size_t) ) {
    fprintf( stderr, "%s: WARNING\nbytes_total is not a size_t\n",
             __FUNCTION__ );
  }

  mris_alloced++ ;

  if ((width <= 0) || (height <= 0) || (depth <= 0))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIallocChunk(%d, %d, %d): bad parm",
                 width, height, depth)) ;
  mri = MRIallocHeader(width, height, depth, type, 1) ;
  mri->nframes = nframes ;
  MRIinitHeader(mri) ;

  // Allocate a big chunk of memory
  mri->ischunked = 1;
  mri->bytes_per_row   = mri->bytes_per_vox   * mri->width;
  mri->bytes_per_slice = mri->bytes_per_row   * mri->height;
  mri->bytes_per_vol   = mri->bytes_per_slice * mri->depth;
  mri->bytes_total     = mri->bytes_per_vol   * mri->nframes;
  mri->chunk = calloc(mri->bytes_total,1);
  if (mri->chunk == NULL)
  {
    printf("ERROR: MRIallocChunk(): could not alloc %lu\n",
	   (unsigned long)mri->bytes_total);
    return(NULL);
  }
  //printf("Allocing MRI with Chunk\n");

  MRIallocIndices(mri) ;  // not sure what this does
  mri->outside_val = 0 ;
  mri->slices = (BUFTYPE ***)calloc(depth*nframes, sizeof(BUFTYPE **)) ;
  if (!mri->slices)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIalloc: could not allocate %d slices\n", mri->depth) ;

  p = mri->chunk;
  for (slice = 0 ; slice < depth*nframes ; slice++)
  {
    /* allocate pointer to array of rows */
    mri->slices[slice] = (BUFTYPE **)calloc(mri->height, sizeof(BUFTYPE *)) ;
    if (!mri->slices[slice])
      ErrorExit
      (ERROR_NO_MEMORY,
       "MRIallocChunk(%d, %d, %d): could not allocate "
       "%d bytes for %dth slice\n",
       height, width, depth, mri->height*sizeof(BUFTYPE *), slice) ;
    /* Instead of allocating each row, just point to the
       correct location in the chunk. */
    for (row = 0 ; row < mri->height ; row++)
    {
      mri->slices[slice][row] = p;
      p += mri->bytes_per_row;
    }
  }
  return(mri) ;
}
/*-------------------------------------------------------------*/
/*!
  \fn MRI *MRIallocSequence(int width, int height, int depth, int type, int nframes)
  \brief Allocate header and buffer for MRI struct. Maybe chunked or not
   depending on the FS_USE_MRI_CHUNK environment variable.
*/
/*-------------------------------------------------------------*/
MRI *MRIallocSequence(int width, int height, int depth, int type, int nframes)
{
  MRI     *mri ;
  int     slice, row, bpp;
  BUFTYPE *buf ;

  if (getenv("FS_USE_MRI_CHUNK") != NULL)
  {
    mri = MRIallocChunk(width, height, depth, type, nframes);
    return(mri);
  }

  mris_alloced++ ;

  if ((width <= 0) || (height <= 0) || (depth <= 0))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MRIalloc(%d, %d, %d): bad parm",
                 width, height, depth)) ;
#if 1
  mri = MRIallocHeader(width, height, depth, type, nframes) ;
  MRIinitHeader(mri) ;
#else
  mri = (MRI *)calloc(1, sizeof(MRI)) ;
  if (!mri)
    ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI\n") ;

  mri->scale = 1 ;
  mri->height = height ;
  mri->width = width ;
  mri->yinvert = 1 ;
  mri->depth = depth ;
  mri->type = type ;
  mri->frames = (MRI_FRAME *)calloc(nframes, sizeof(MRI_FRAME)) ;
  if (!mri->frames)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIalloc: could not allocate %d frames\n", nframes) ;
  for (int i = 0 ; i < nframes ; i++)
    mri->frames[i].m_ras2vox = MatrixAlloc(4,4, MATRIX_REAL) ;
#endif
  mri->nframes = nframes ;
  MRIallocIndices(mri) ;
  mri->outside_val = 0 ;
  mri->slices = (BUFTYPE ***)calloc(depth*nframes, sizeof(BUFTYPE **)) ;
  if (!mri->slices)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIalloc: could not allocate %d slices\n", mri->depth) ;

  for (slice = 0 ; slice < depth*nframes ; slice++)
  {
    /* allocate pointer to array of rows */
    mri->slices[slice] = (BUFTYPE **)calloc(mri->height, sizeof(BUFTYPE *)) ;
    if (!mri->slices[slice])
      ErrorExit
      (ERROR_NO_MEMORY,
       "MRIalloc(%d, %d, %d): could not allocate "
       "%d bytes for %dth slice\n",
       height, width, depth, mri->height*sizeof(BUFTYPE *), slice) ;

#if USE_ELECTRIC_FENCE
    switch (mri->type)
    {
    case MRI_BITMAP:
      bpp = 1 ;
      break ;
    case MRI_UCHAR:
      bpp = 8 ;
      break ;
    case MRI_TENSOR:
    case MRI_FLOAT:
      bpp = sizeof(float) * 8 ;
      break ;
    case MRI_INT:
      bpp = sizeof(int) * 8 ;
      break ;
    case MRI_SHORT:
      bpp = sizeof(short) * 8 ;
      break ;
    case MRI_LONG:
      bpp = sizeof(long) * 8 ;
      break ;
    default:
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "MRIalloc(%d, %d, %d, %d): unknown type",
                   width, height, depth, mri->type)) ;
      break ;
    }
    bpp /= 8 ;
    buf = (BUFTYPE *)calloc((mri->width*mri->height*bpp), 1) ;
    if (buf == NULL)
      ErrorExit
      (ERROR_NO_MEMORY,
       "MRIalloc(%d, %d, %d): could not allocate "
       "%d bytes for %dth slice\n",
       height, width, depth, (mri->width*mri->height*bpp), slice) ;
    for (row = 0 ; row < mri->height ; row++)
    {
      mri->slices[slice][row] = buf+(row*mri->width*bpp) ;
    }
#else
    /* allocate each row */
    for (row = 0 ; row < mri->height ; row++)
    {
      switch (mri->type)
      {
      case MRI_BITMAP:
        mri->slices[slice][row] =
          (BUFTYPE*)calloc(mri->width/8,sizeof(BUFTYPE));
        break ;
      case MRI_UCHAR:
        mri->slices[slice][row] =
          (BUFTYPE*)calloc(mri->width,sizeof(BUFTYPE));
        break ;
      case MRI_TENSOR:
      case MRI_FLOAT:
        mri->slices[slice][row] =
          (BUFTYPE *)calloc(mri->width, sizeof(float));
        break ;
      case MRI_INT:
        mri->slices[slice][row] =
          (BUFTYPE *)calloc(mri->width, sizeof(int)) ;
        break ;
      case MRI_SHORT:
        mri->slices[slice][row] =
          (BUFTYPE *)calloc(mri->width, sizeof(short));
        break ;
      case MRI_LONG:
        mri->slices[slice][row] =
          (BUFTYPE *)calloc(mri->width, sizeof(long)) ;
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM,
                     "MRIalloc(%d, %d, %d, %d): unknown type",
                     width, height, depth, mri->type)) ;
        break ;
      }

      if (!mri->slices[slice][row])
        ErrorExit
        (ERROR_NO_MEMORY,
         "MRIalloc(%d,%d,%d): could not allocate "
         "%dth row in %dth slice\n",
         width,height,depth, slice, row) ;
    }
#endif
  }

  return(mri) ;
}
/*-----------------------------------------------------*/
/*!
  \fn MRI *MRIallocHeader(int width, int height, int depth, int type, int nframes)
  \brief allocate an MRI data structure but not space for  the image data
*/
MRI *MRIallocHeader(int width, int height, int depth, int type, int nframes)
{
  MRI  *mri ;
  int  i ;

  mri = (MRI *)calloc(1, sizeof(MRI)) ;
  if (!mri)
    ErrorExit(ERROR_NO_MEMORY, "MRIalloc: could not allocate MRI\n") ;

  // Note: changes here may need to be reflected in MRISeqchangeType()
  mri->frames = (MRI_FRAME *)calloc(nframes, sizeof(MRI_FRAME)) ;
  if (!mri->frames)
    ErrorExit(ERROR_NO_MEMORY,
              "MRIalloc: could not allocate %d frame\n", nframes) ;
  for (i = 0 ; i < mri->nframes ; i++)
    mri->frames[i].m_ras2vox = MatrixAlloc(4,4, MATRIX_REAL) ;

  mri->imnr0 = 1 ;
  mri->imnr1 = depth;
  mri->fov = width ;
  mri->thick = 1.0 ;
  mri->scale = 1 ;
  mri->roi.dx = mri->width = width ;
  mri->roi.dy = mri->height = height ;
  mri->roi.dz = mri->depth = depth ;
  mri->yinvert = 1 ;
  mri->xsize = mri->ysize = mri->zsize = 1 ;
  mri->type = type ;
  mri->nframes = nframes ;
  mri->xi = mri->yi = mri->zi = NULL ;
  mri->slices = NULL ;
  mri->ps = 1 ;
  mri->xstart = -mri->width/2.0 ;
  mri->xend = mri->width/2.0 ;
  mri->ystart = -mri->height/2.0 ;
  mri->yend = mri->height/2.0 ;
  mri->zstart = -mri->depth/2.0 ;
  mri->zend = mri->depth/2 ;
  //
  mri->x_r = -1;
  mri->x_a = 0.;
  mri->x_s = 0.;
  //
  mri->y_r = 0.;
  mri->y_a = 0.;
  mri->y_s = -1;
  //
  mri->z_r = 0.;
  mri->z_a = 1.;
  mri->z_s = 0.;
  //
  mri->c_r = mri->c_a = mri->c_s = 0.0;
  mri->ras_good_flag = 0;
  mri->brightness = 1;
  mri->register_mat = NULL;
  mri->subject_name[0] = '\0';
  mri->path_to_t1[0] = '\0';
  mri->fname_format[0] = '\0';
  mri->gdf_image_stem[0] = '\0';
  mri->tag_data = NULL;
  mri->tag_data_size = 0;
  mri->transform_fname[0] = '\0';
  mri->pedir = NULL;

  MATRIX *tmp;
  tmp = extract_i_to_r(mri);
  AffineMatrixAlloc( &(mri->i_to_r__ ) );
  SetAffineMatrix( mri->i_to_r__, tmp );
  MatrixFree( &tmp );

  mri->r_to_i__ = extract_r_to_i(mri);
  mri->AutoAlign = NULL;

  // Chunking memory management
  mri->ischunked = 0;
  mri->chunk = NULL;
  mri->bytes_per_vox   = MRIsizeof(type);

  // These things are explicitly set to 0 here because we
  // do not yet know the true number of frames, and they
  // are set in MRIallocChunk().
  mri->bytes_per_row   = 0;
  mri->bytes_per_slice = 0;
  mri->bytes_per_vol   = 0;
  mri->bytes_total     = 0;

  mri->bvals = NULL; // For DWI
  mri->bvecs = NULL;

  return(mri) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  allocate an MRI data structure as well as space for
  the image data
  ------------------------------------------------------*/
MRI *
MRIalloc(int width, int height, int depth, int type)
{
  return(MRIallocSequence(width, height, depth, type, 1)) ;
}
/*-----------------------------------------------------
  Free and MRI data structure and all its attached memory
  ------------------------------------------------------*/
int
MRIfree(MRI **pmri)
{
  MRI *mri ;
  int slice, i ;
#if !USE_ELECTRIC_FENCE
  int  row ;
#endif

  mris_alloced-- ;
  mri = *pmri ;
  if (!mri)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIfree: null pointer\n")) ;

  if (mri->frames)
  {
    int i ;
    for (i = 0 ; i < mri->nframes ; i++)
      if (mri->frames[i].m_ras2vox)
        MatrixFree(&mri->frames[i].m_ras2vox) ;
    free(mri->frames) ;
  }
  if (mri->xi)    free(mri->xi-MAX_INDEX) ;
  if (mri->yi)    free(mri->yi-MAX_INDEX) ;
  if (mri->zi)    free(mri->zi-MAX_INDEX) ;

  if (!mri->ischunked)
  {
    if (mri->slices)
    {
      for (slice = 0 ; slice < mri->depth*mri->nframes ; slice++)
      {
        if (mri->slices[slice])
        {
#if USE_ELECTRIC_FENCE
          free(mri->slices[slice][0]) ;
#else
          for (row = 0 ; row < mri->height ; row++)
            free(mri->slices[slice][row]) ;
#endif
          free(mri->slices[slice]) ;
        }
      }
      free(mri->slices) ;
    }
  }
  else
  {
    //printf("Freeing MRI Chunk\n");
    free(mri->chunk);
    mri->chunk = NULL;
    free(mri->slices) ;
  }

  if (mri->free_transform)
    delete_general_transform(&mri->transform) ;

  if (mri->register_mat != NULL)
    MatrixFree(&(mri->register_mat));

  if (mri->i_to_r__) { 
    AffineMatrixFree( &mri->i_to_r__ );
  }

  if (mri->r_to_i__)    MatrixFree(&mri->r_to_i__);
  if (mri->AutoAlign)   MatrixFree(&mri->AutoAlign);

  for (i = 0 ; i < mri->ncmds ; i++) free(mri->cmdlines[i]) ;

  if(mri->bvals) MatrixFree(&mri->bvals);
  if(mri->bvecs) MatrixFree(&mri->bvecs);
  
  if(mri->pedir) free(mri->pedir);

  free(mri) ;
  *pmri = NULL ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
/*!
  \fn MRIfreeFrames(MRI *mri, int start_frame)
  \brief Frees frames of the mri struct starting at start_frame.
  Cannot be run on mri struct with chunked mem management.
*/
int MRIfreeFrames(MRI *mri, int start_frame)
{
  int slice, row, end_frame ;

  if (mri->ischunked)
  {
    printf("WARNING: MRIfreeFrames(): cannot free frames "
           "from chunked memory\n");
    return(1);
  }

  end_frame = mri->nframes-1 ;
  if (!mri)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIfree: null pointer\n")) ;

  if (mri->slices)
  {
    for (slice = start_frame*mri->depth ;
         slice < mri->depth*end_frame ; slice++)
    {
      if (mri->slices[slice])
      {
        for (row = 0 ; row < mri->height ; row++)
          free(mri->slices[slice][row]) ;
        free(mri->slices[slice]) ;
      }
    }
  }

  mri->nframes -= (end_frame-start_frame+1) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
/*!
  \fn MRIdump(MRI *mri, FILE *fp)
  \brief Dump the MRI header to a file
*/
int MRIdump(MRI *mri, FILE *fp)
{
  fprintf(fp, "%6.6s = %s\n", "fname", mri->fname);
  fprintf(fp, "%6.6s = %d\n", "height", mri->height);
  fprintf(fp, "%6.6s = %d\n", "width", mri->width);
  fprintf(fp, "%6.6s = %d\n", "depth", mri->depth);
  fprintf(fp, "%6.6s = %d\n", "nframes", mri->nframes);
  fprintf(fp, "%6.6s = %d\n", "imnr0", mri->imnr0);
  fprintf(fp, "%6.6s = %d\n", "imnr1", mri->imnr1);
  fprintf(fp, "%6.6s = %d\n", "xnum", mri->width);
  fprintf(fp, "%6.6s = %d\n", "ynum", mri->height);
  fprintf(fp, "%6.6s = %f\n", "fov", mri->fov);
  fprintf(fp, "%6.6s = %f\n", "thick", mri->thick);
  fprintf(fp, "%6.6s = %f\n", "xstart", mri->xstart); /* strtx */
  fprintf(fp, "%6.6s = %f\n", "xend", mri->xend); /* endx */
  fprintf(fp, "%6.6s = %f\n", "ystart", mri->ystart); /* strty */
  fprintf(fp, "%6.6s = %f\n", "yend", mri->yend); /* endy */
  fprintf(fp, "%6.6s = %f\n", "zstart", mri->zstart); /* strtz */
  fprintf(fp, "%6.6s = %f\n", "zend", mri->zend); /* endz */
  fprintf(fp, "%6.6s = %d\n", "type", mri->type);
  fprintf(fp, "%6.6s = %f\n", "xsize", mri->xsize);
  fprintf(fp, "%6.6s = %f\n", "ysize", mri->ysize);
  fprintf(fp, "%6.6s = %f\n", "zsize", mri->zsize);
  fprintf(fp, "%6.6s = %f %f %f\n", "x ras", mri->x_r, mri->x_a, mri->x_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "y ras", mri->y_r, mri->y_a, mri->y_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "z ras", mri->z_r, mri->z_a, mri->z_s);
  fprintf(fp, "%6.6s = %f %f %f\n", "c ras", mri->c_r, mri->c_a, mri->c_s);
  fprintf(fp, "%s = %f\n", "det(xyz_ras)", MRIvolumeDeterminant(mri));
  fprintf(fp, "%s = %d\n", "ras_good_flag", mri->ras_good_flag);
  fprintf(fp, "%s = %d\n", "brightness", mri->brightness);
  fprintf(fp, "%s = %s\n", "subject_name", mri->subject_name);
  fprintf(fp, "%s = %s\n", "path_to_t1", mri->path_to_t1);
  fprintf(fp, "%s = %s\n", "fname_format", mri->fname_format);
  if (mri->register_mat != NULL)
  {
    fprintf(fp, "%s = \n", "register_mat");
    MatrixPrint(fp, mri->register_mat);
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
/*!
  \fn int MRIdumpBuffer(MRI *mri, FILE *fp)
  \brief Dump the non-zero elements of an MRI buffer to a file
*/
int MRIdumpBuffer(MRI *mri, FILE *fp)
{
  int  x, y, z ;

  for (z = 0 ; z < mri->depth ; z++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (x = 0 ; x < mri->width ; x++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          if (!FZERO(mri->slices[z][y][x]))
            fprintf
            (fp,
             "[%d][%d][%d]: %d\n",
             x,y,z,(int)mri->slices[z][y][x]);
          break ;
        case MRI_FLOAT:
          /*          if (!FZERO(MRIFvox(mri,x,y,z)))*/
          fprintf
          (fp,
           "[%d][%d][%d]: %2.3f\n",
           x,y,z, MRIFvox(mri,x,y,z));
          break ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------*/
/*!
  \fn int MRIpeak(MRI *mri, int *px, int *py, int *pz)
  \brief Find the peak intensity in an MRI image.
*/
int MRIpeak(MRI *mri, int *px, int *py, int *pz)
{
  int      max_row, max_col, max_slice, row, col, slice, width, height, depth ;
  BUFTYPE  val, max_val, *im ;
  long     lval, lmax_val, *lim ;
  int      ival, imax_val, *iim ;
  float    fval, fmax_val, *fim ;

  max_val = 0 ;
  lmax_val = 0L ;
  imax_val = 0 ;
  fmax_val = 0.0f ;
  max_row = max_col = max_slice = -1 ;   /* to prevent compiler warning */
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  for (slice = 0 ; slice < depth ; slice++)
  {
    for (row = 0 ; row < height ; row++)
    {
      switch (mri->type)
      {
      case MRI_UCHAR:
        im = mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          val = *im++ ;
          if (val > max_val)
          {
            max_val = val ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_LONG:
        lim = (long *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          lval = *lim++ ;
          if (lval > lmax_val)
          {
            lmax_val = lval ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_FLOAT:
        fim = (float *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          fval = *fim++ ;
          if (fval > fmax_val)
          {
            fmax_val = fval ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      case MRI_INT:
        iim = (int *)mri->slices[slice][row] ;
        for (col = 0 ; col < width ; col++)
        {
          ival = *iim++ ;
          if (ival > imax_val)
          {
            imax_val = ival ;
            max_row = row ;
            max_col = col ;
            max_slice = slice ;
          }
        }
        break ;
      }
    }
  }

  *px = max_col ;
  *py = max_row ;
  *pz = max_slice ;
  return(NO_ERROR) ;
}
/*
  compare two headers to see if they are the same voxel and ras coords
*/
int
MRIcompareHeaders(MRI *mri1, MRI *mri2)
{
  if (mri1 == NULL || mri2 == NULL)
    return(1) ; // not the same
  if (mri1->xsize != mri2->xsize)
    return(1) ;
  if (mri1->ysize != mri2->ysize)
    return(1) ;
  if (mri1->zsize != mri2->zsize)
    return(1) ;
  if (mri1->ptype != mri2->ptype)
    return(1) ;
  if (mri1->fov != mri2->fov)
    return(1) ;
  if (mri1->thick != mri2->thick)
    return(1) ;
  if (mri1->ps != mri2->ps)
    return(1) ;
  if (mri1->location != mri2->location)
    return(1) ;
  if (mri1->xstart != mri2->xstart)
    return(1) ;
  if (mri1->xend != mri2->xend)
    return(1) ;
  if (mri1->ystart != mri2->ystart)
    return(1) ;
  if (mri1->yend != mri2->yend)
    return(1) ;
  if (mri1->zstart != mri2->zstart)
    return(1) ;
  if (mri1->zend != mri2->zend)
    return(1) ;
  if (mri1->flip_angle != mri2->flip_angle)
    return(1) ;
  if (mri1->tr != mri2->tr)
    return(1) ;
  if (mri1->te != mri2->te)
    return(1) ;
  if (mri1->ti != mri2->ti)
    return(1) ;
  if (mri1->x_r != mri2->x_r)
    return(1) ;
  if (mri1->x_a != mri2->x_a)
    return(1) ;
  if (mri1->x_s != mri2->x_s)
    return(1) ;
  if (mri1->y_r != mri2->y_r)
    return(1) ;
  if (mri1->y_a != mri2->y_a)
    return(1) ;
  if (mri1->y_s != mri2->y_s)
    return(1) ;
  if (mri1->z_r != mri2->z_r)
    return(1) ;
  if (mri1->z_a != mri2->z_a)
    return(1) ;
  if (mri1->z_s != mri2->z_s)
    return(1) ;
  if (mri1->c_r != mri2->c_r)
    return(1) ;
  if (mri1->c_a != mri2->c_a)
    return(1) ;
  if (mri1->c_s != mri2->c_s)
    return(1) ;
  if (mri1->ras_good_flag != mri2->ras_good_flag)
    return(1) ;
  return(0) ;   // they are the same
}

/*--------------------------------------------------------------
  Description: Copy the header information from one MRI into another.
  Does not copy the dimension lengths, only the geometry, pulse seq,
  etc. Does not copy ischunked or chunk pointer.
  ------------------------------------------------------*/
MRI *
MRIcopyHeader( const MRI *mri_src, MRI *mri_dst )
{
  int i ;

  if (mri_dst == NULL)
    mri_dst = MRIallocHeader(mri_src->width, mri_src->height, mri_src->depth, mri_src->type, mri_src->nframes) ;

  mri_dst->dof = mri_src->dof ;
  mri_dst->mean = mri_src->mean ;
  mri_dst->xsize = mri_src->xsize ;
  mri_dst->ysize = mri_src->ysize ;
  mri_dst->zsize = mri_src->zsize ;
  
  if (mri_dst->free_transform)
    delete_general_transform(&mri_dst->transform) ;
  if (mri_src->linear_transform)
  {
    copy_general_transform(&(((MRI*)mri_src)->transform),
			   &mri_dst->transform) ;
    mri_dst->linear_transform = mri_src->linear_transform ;
    mri_dst->inverse_linear_transform = mri_src->inverse_linear_transform ;
    mri_dst->linear_transform =
      get_linear_transform_ptr(&mri_dst->transform) ;
    mri_dst->inverse_linear_transform =
      get_inverse_linear_transform_ptr(&mri_dst->transform) ;
    mri_dst->free_transform = 1 ;
  }
  strcpy(mri_dst->transform_fname, mri_src->transform_fname) ;
  if (mri_dst->depth == mri_src->depth)
  {
    mri_dst->imnr0 = mri_src->imnr0 ;
    mri_dst->imnr1 = mri_src->imnr1 ;
  }
  mri_dst->ptype = mri_src->ptype ;
  mri_dst->fov = mri_src->fov ;
  mri_dst->thick = mri_src->thick ;
  mri_dst->ps = mri_src->ps ;
  mri_dst->location = mri_src->location ;
  mri_dst->xstart = mri_src->xstart ;
  mri_dst->xend = mri_src->xend ;
  mri_dst->ystart = mri_src->ystart ;
  mri_dst->yend = mri_src->yend ;
  mri_dst->zstart = mri_src->zstart ;
  mri_dst->zend = mri_src->zend ;
  mri_dst->flip_angle = mri_src->flip_angle ;
  mri_dst->tr = mri_src->tr ;
  mri_dst->te = mri_src->te ;
  mri_dst->ti = mri_src->ti ;
  strcpy(mri_dst->fname, mri_src->fname) ;
  mri_dst->x_r = mri_src->x_r;
  mri_dst->x_a = mri_src->x_a;
  mri_dst->x_s = mri_src->x_s;
  mri_dst->y_r = mri_src->y_r;
  mri_dst->y_a = mri_src->y_a;
  mri_dst->y_s = mri_src->y_s;
  mri_dst->z_r = mri_src->z_r;
  mri_dst->z_a = mri_src->z_a;
  mri_dst->z_s = mri_src->z_s;
  mri_dst->c_r = mri_src->c_r;
  mri_dst->c_a = mri_src->c_a;
  mri_dst->c_s = mri_src->c_s;
  mri_dst->ras_good_flag = mri_src->ras_good_flag;

  mri_dst->bytes_per_vox   = MRIsizeof(mri_dst->type);
  mri_dst->bytes_per_row   = mri_dst->bytes_per_vox   * mri_dst->width;
  mri_dst->bytes_per_slice = mri_dst->bytes_per_row   * mri_dst->height;
  mri_dst->bytes_per_vol   = mri_dst->bytes_per_slice * mri_dst->depth;
  mri_dst->bytes_total     = mri_dst->bytes_per_vol   * mri_dst->nframes;

  mri_dst->brightness = mri_src->brightness;
  if (mri_src->register_mat != NULL)
    MatrixCopy(mri_src->register_mat, mri_dst->register_mat);
  else
    mri_dst->register_mat = NULL;
  strcpy(mri_dst->subject_name, mri_src->subject_name);
  strcpy(mri_dst->path_to_t1, mri_src->path_to_t1);
  strcpy(mri_dst->fname_format, mri_src->fname_format);

  strcpy(mri_dst->gdf_image_stem, mri_src->gdf_image_stem);

  mri_dst->i_to_r__ = AffineMatrixCopy( mri_src->i_to_r__,
					mri_dst->i_to_r__ );

  mri_dst->r_to_i__ = MatrixCopy(mri_src->r_to_i__, mri_dst->r_to_i__);
  if(mri_src->AutoAlign != NULL){
    mri_dst->AutoAlign = MatrixCopy(mri_src->AutoAlign,NULL);
  }


  for (i = 0 ; i < mri_dst->ncmds ; i++) free(mri_dst->cmdlines[i]) ;
  for (i = 0 ; i < mri_src->ncmds ; i++)
  {
    mri_dst->cmdlines[i] =
      (char *)calloc(strlen(mri_src->cmdlines[i])+1, sizeof(char)) ;
    strcpy(mri_dst->cmdlines[i], mri_src->cmdlines[i]) ;
  }
  mri_dst->ncmds = mri_src->ncmds ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Translate the MRI image mri_src by the vector
  dx, dy, dz  and store the result in mri_dst.
  ------------------------------------------------------*/
MRI *
MRItranslate(MRI *mri_src, MRI *mri_dst, double dx, double dy, double dz)
{
  int      y1, y2, y3, width, height, depth ;
  BUFTYPE *pdst ;
  double     x1, x2, x3, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  for (y3 = 0 ; y3 < depth ; y3++)
  {
    x3 = (double)y3 - dz ;
    if (x3 < 0 || x3 >= depth)
      continue ;

    for (y2 = 0 ; y2 < height ; y2++)
    {
      x2 = (double)y2 - dy ;
      if (x2 < 0 || x2 >= height)
        continue ;

      pdst = &MRIvox(mri_dst, 0, y2, y3) ;
      for (y1 = 0 ; y1 < width ; y1++, pdst++)
      {
        x1 = (double)y1 - dx ;
        if (x1 >= 0 && x1 < width)
        {
          MRIsampleVolume(mri_src, x1, x2, x3, &val) ;
          *pdst = (BUFTYPE)nint(val) ;
        }
      }
    }
  }

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Y axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateX(MRI *mri_src, MRI *mri_dst, float x_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, x_angle, X_ROTATION) ;

  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = (double)mri_src->width / 2.0 ;
  mO->rptr[2][1] = (double)mri_src->height / 2.0 ;
  mO->rptr[3][1] = (double)mri_src->depth / 2.0 ;

  /* build rotation matrix */
  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Y axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateY(MRI *mri_src, MRI *mri_dst, float y_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  /* origin of coordinate system */
  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = (double)mri_src->width / 2.0 ;
  mO->rptr[2][1] = (double)mri_src->height / 2.0 ;
  mO->rptr[3][1] = (double)mri_src->depth / 2.0 ;

  m = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Z axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateZ(MRI *mri_src, MRI *mri_dst, float z_angle)
{
  int    width, height, depth ;
  MATRIX *m, *mO ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
  mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
  mO->rptr[1][1] = (double)mri_src->width / 2.0 ;
  mO->rptr[2][1] = (double)mri_src->height / 2.0 ;
  mO->rptr[3][1] = (double)mri_src->depth / 2.0 ;

  mri_dst = MRIrotate(mri_src, NULL, m, mO) ;
  MatrixFree(&m) ;
  MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Y axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateX_I(MRI *mri_src, MRI *mri_dst, float x_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, x_angle, X_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Y axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateY_I(MRI *mri_src, MRI *mri_dst, float y_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, y_angle, Y_ROTATION) ;

  /* build rotation matrix */
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate mri_src around the Z axis and return the
  result in mri_dst
  ------------------------------------------------------*/
MRI *
MRIrotateZ_I(MRI *mri_src, MRI *mri_dst, float z_angle)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  m = MatrixAllocRotation(3, z_angle, Z_ROTATION) ;
  mri_dst = MRIrotate_I(mri_src, NULL, m, NULL) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale the MRI image mri_src by sx,sy,sz in the
  x, y, and z directions respectively.
  ------------------------------------------------------*/
MRI *
MRIscale(MRI *mri_src, MRI *mri_dst, float sx, float sy, float sz)
{
  int    width, height, depth ;
  MATRIX *m ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  m = MatrixAlloc(4, 4, MATRIX_REAL) ;

  /* build rotation matrix */
  m->rptr[1][1] = sx ;
  m->rptr[2][2] = sy ;
  m->rptr[3][3] = sz ;
  m->rptr[4][4] = 1.0 ;
  mri_dst = MRIlinearTransform(mri_src, NULL, m) ;
  MatrixFree(&m) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate about the point mO
  ------------------------------------------------------*/
MRI *
MRIrotate(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO)
{
  int    x1, x2, x3, y1, y2, y3, width, height, depth, y1o, y2o, y3o, freeit ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */
  MATRIX *mRinv ;   /* inverse of R */

  mRinv = MatrixInverse(mR, NULL) ;
  if (!mRinv)
    ErrorReturn(NULL,(ERROR_BADPARM,"MRIrotate: rotation matrix is singular"));

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  if (!mO)
  {
    mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
    mO->rptr[1][1] = (double)mri_src->width / 2.0 ;
    mO->rptr[2][1] = (double)mri_src->height / 2.0 ;
    mO->rptr[3][1] = (double)mri_src->depth / 2.0 ;
    freeit = 1 ;
  }
  else
    freeit = 0 ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */


  y1o = mO->rptr[1][1] ;
  y2o = mO->rptr[2][1] ;
  y3o = mO->rptr[3][1] ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    mY->rptr[3][1] = y3 - y3o ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      mY->rptr[2][1] = y2 - y2o ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        mY->rptr[1][1] = y1 - y1o ;
        MatrixMultiply(mRinv, mY, mX) ;
        MatrixAdd(mX, mO, mX) ;

        /* should do bilinear interpolation here */
        x1 = nint(mX->rptr[1][1]) ;
        x2 = nint(mX->rptr[2][1]) ;
        x3 = nint(mX->rptr[3][1]) ;
        if (x1 >= 0 && x1 < width &&
            x2 >= 0 && x2 < height &&
            x3 >= 0 && x3 < depth)
          mri_dst->slices[y3][y2][y1] = mri_src->slices[x3][x2][x1] ;
      }
    }
  }

  MatrixFree(&mX) ;
  MatrixFree(&mRinv) ;
  MatrixFree(&mY) ;
  if (freeit)
    MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Rotate about the point mO using trilinear interpolation
  ------------------------------------------------------*/
MRI *
MRIrotate_I(MRI *mri_src, MRI *mri_dst, MATRIX *mR, MATRIX *mO)
{
  int    y1, y2, y3, width, height, depth, y1o, y2o, y3o, freeit ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */
  MATRIX *mRinv ;   /* inverse of R */
  float  x1, x2, x3 ;
  double   val ;

  mRinv = MatrixInverse(mR, NULL) ;
  if (!mRinv)
    ErrorReturn(NULL,(ERROR_BADPARM,
                      "MRIrotate_I: rotation matrix is singular"));

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  if (!mO)
  {
    mO = MatrixAlloc(3, 1, MATRIX_REAL) ;
    mO->rptr[1][1] = (double)mri_src->width / 2.0 ;
    mO->rptr[2][1] = (double)mri_src->height / 2.0 ;
    mO->rptr[3][1] = (double)mri_src->depth / 2.0 ;
    freeit = 1 ;
  }
  else
    freeit = 0 ;

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */


  y1o = mO->rptr[1][1] ;
  y2o = mO->rptr[2][1] ;
  y3o = mO->rptr[3][1] ;
  if (Gdiag == 99)
    MatrixPrint(stdout, mRinv) ;
  if (Gdiag == 99)
    MatrixPrint(stdout, mO) ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    mY->rptr[3][1] = y3 - y3o ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      mY->rptr[2][1] = y2 - y2o ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        mY->rptr[1][1] = y1 - y1o ;
        MatrixMultiply(mRinv, mY, mX) ;
        MatrixAdd(mX, mO, mX) ;

        /* do trilinear interpolation here */
        x1 = mX->rptr[1][1] ;
        x2 = mX->rptr[2][1] ;
        x3 = mX->rptr[3][1] ;

        MRIsampleVolume(mri_src, x1, x2, x3, &val) ;
        mri_dst->slices[y3][y2][y1] = (BUFTYPE)nint(val) ;
      }
    }
  }

  MatrixFree(&mX) ;
  MatrixFree(&mRinv) ;
  MatrixFree(&mY) ;
  if (freeit)
    MatrixFree(&mO) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Perform an affine coordinate transformation x' = Ax + B on
  the MRI image mri_src into mri_dst
  ------------------------------------------------------*/
MRI *
MRIaffine(MRI *mri_src, MRI *mri_dst, MATRIX *mA, MATRIX *mB)
{
  int    x1, x2, x3, y1, y2, y3, width, height, depth ;
  MATRIX *mX, *mY ;   /* original and transformed coordinate systems */

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  mX = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* input coordinates */
  mY = MatrixAlloc(3, 1, MATRIX_REAL) ;  /* transformed coordinates */

  for (x3 = 0 ; x3 < depth ; x3++)
  {
    mX->rptr[3][1] = x3 ;
    for (x2 = 0 ; x2 < height ; x2++)
    {
      mX->rptr[2][1] = x2 ;
      for (x1 = 0 ; x1 < width ; x1++)
      {
        mX->rptr[1][1] = x1 ;
        if (mA)
          MatrixMultiply(mA, mX, mY) ;
        else
          MatrixCopy(mX, mY) ;
        if (mB)
          MatrixAdd(mY, mB, mY) ;
        y1 = nint(mY->rptr[1][1]) ;
        y2 = nint(mY->rptr[2][1]) ;
        y3 = nint(mY->rptr[3][1]) ;
        if (y1 >= 0 && y1 < width &&
            y2 >= 0 && y2 < height &&
            y3 >= 0 && y3 < depth)
          mri_dst->slices[y3][y2][y1] = mri_src->slices[x3][x2][x1] ;
      }
    }
  }


  MatrixFree(&mX) ;
  MatrixFree(&mY) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Convert a slice of an MRI data structure into a HIPS image.
  ------------------------------------------------------*/
IMAGE *
MRItoImageView(MRI *mri, IMAGE *I, int slice, int view, int frame)
{
  int      width, height, depth, x, y, yp, w, h, d,
  xm, ym, zm, format ;
  float    fmin, fmax, frac ;
  double     val ;
  int src_slice_direction;
  int xsign, ysign;

  xsign=ysign=1; // let compiler be satisfied

  src_slice_direction = getSliceDirection(mri);
  // only strict slice direction is supported
  if (src_slice_direction == MRI_UNDEFINED)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRItoImageView(%d, %d): unsupported view/slice direction %d",
                 slice, view, src_slice_direction)) ;

  // generic routines to get the right axes
  //
  // should be able to do this in generic terms.
  // don't have time to do that
  //
  // Define the view to be the following:
  // coronal    = (-R, -S)
  // sagittal   = ( A, -S)
  // horizontal = (-R,  A)
  //
  // translate these direction into (x,y)
  d = w = h = xm = ym = zm = 0 ;  /* silly compiler warnings */
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  frac = .6 ;
  switch (src_slice_direction)
  {
  case MRI_CORONAL: // x direction can be -R or R,
    // y direction can be -S or S, z direction can be -A or A
    //     +/-1  0  0
    //      0    0  +/-1
    //      0  +/-1 0
    switch (view)
    {
    case MRI_CORONAL:
      frac = .4 ;
      w = width;
      h = height;
      d = depth;
      xsign = (mri->x_r > 0) ? -1 : 1;
      ysign = (mri->y_s > 0) ? -1 : 1;
      break;
    case MRI_SAGITTAL:
      w = depth;
      h = height;
      d = width;
      xsign = (mri->z_a > 0) ?  1 : -1;
      ysign = (mri->y_s > 0) ? -1 : 1;
      break;
    case MRI_HORIZONTAL:
      w = width;
      h = depth;
      d = height;
      xsign = (mri->x_r > 0) ? -1 : 1;
      ysign = (mri->z_a > 0) ? 1 : -1;
      break;
    }
    break;
  case MRI_SAGITTAL: // x direction can be -A or A,
    // y direction can be -S or S, z direction can be -R or R
    //    0    0   +/-1
    //  +/-1   0    0
    //    0  +/-1   0
    switch (view)
    {
    case MRI_CORONAL:
      w = depth;
      h = height;
      d = width;
      xsign =
        (mri->z_r > 0) ? -1 : 1;
      ysign = (mri->y_s > 0) ? -1 : 1;
      break;
    case MRI_SAGITTAL:
      w = width;
      h = height;
      frac = .4 ;
      d = depth;
      xsign = (mri->x_a > 0) ? 1 : -1;
      ysign = (mri->y_s > 0) ? -1 : 1;
      break;
    case MRI_HORIZONTAL:
      w = depth;
      h = width;
      d = height;
      xsign = (mri->z_r > 0) ? -1 : 1;
      ysign = (mri->x_a > 0) ? 1 : -1;
      break;
    }
    break;
  case MRI_HORIZONTAL: // x direction can be -R or R,
    // y direction can be -A or A, z direction can be -S or S
    //   +/-1  0   0
    //    0   +/-1 0
    //    0    0  +/-1
    switch (view)
    {
    case MRI_CORONAL:
      w = width;
      h = depth;
      d = height;
      xsign = (mri->x_r > 0) ? -1 : 1;
      ysign = (mri->z_s > 0) ? -1 : 1;
      break;
    case MRI_SAGITTAL:
      w = height;
      h = depth;
      d = width;
      xsign = (mri->y_a > 0) ? 1 : -1;
      ysign = (mri->z_s > 0) ? -1 : 1;
      break;
    case MRI_HORIZONTAL:
      w = width;
      h = height;
      frac = .4 ;
      d = depth;
      xsign = (mri->x_r > 0) ? -1 : 1;
      ysign = (mri->y_a > 0) ?  1 : -1;
      break;
    }
    break;
  default:
    ErrorReturn
    (NULL,
     (ERROR_UNSUPPORTED,
      "MRItoImageView(%d, %d): unsupported view/slice direction %d",
      slice, view, src_slice_direction)) ;
  }
  if (slice < 0)
    slice = nint(frac*d) ;
  else if (slice >= d)
    ErrorReturn(NULL, (ERROR_BADPARM, "MRItoImageView: bad slice %d\n",slice));
    

#if 0
  format = (mri->type == MRI_UCHAR) ? PFBYTE :
           mri->type == MRI_INT   ? PFINT :
           mri->type == MRI_FLOAT ? PFFLOAT :
           mri->type == MRI_SHORT ? PFFLOAT : PFBYTE ;
#else
  format = (mri->type == MRI_UCHAR) ? PFBYTE : PFFLOAT ;
  format = PFBYTE ;
#endif

  if (I && ((I->rows != h) || (I->cols != w) || (I->pixel_format != format)))
  {
    ImageFree(&I);
    I = NULL ;  /* must allocate a new one */
  }
  if (!I)
    I = ImageAlloc(h, w, format, 1) ;

  fmin = 10000000 ;
  fmax = -fmin ;

  // Image values assigned from MRI
  // we need to gert min and max in IMAGE
  for (y = 0 ; y < h ; y++)
  {
    for (x = 0 ; x < w ; x++)
    {
      switch (src_slice_direction)
      {
      case MRI_CORONAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = (xsign> 0) ? x : (w - 1 -x);
          ym = (ysign> 0) ? y : (h - 1 -y);
          zm = slice ;
          break ;
        case MRI_SAGITTAL:
          xm = slice ;
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        case MRI_HORIZONTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = slice ;
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        }
        break;
      case MRI_SAGITTAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = slice;
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        case MRI_SAGITTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = slice ;
          break ;
        case MRI_HORIZONTAL:
          xm = (ysign>0) ? y : (h - 1 - y);
          ym = slice;
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        }
        break;
      case MRI_HORIZONTAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = slice ;
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        case MRI_SAGITTAL:
          xm = slice ;
          ym = (xsign>0) ? x : (w - 1 - x);
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        case MRI_HORIZONTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = slice ;
          break ;
        }
        break;
      }
      MRIsampleVolumeFrame(mri, xm, ym, zm, frame, &val) ;
      // check min max
      if (val > fmax)
        fmax = val ;
      if (val < fmin)
        fmin = val ;
    }
  }
  if (FZERO(fmax-fmin))
    ErrorReturn(I, (ERROR_BADPARM, "MRItoImageView: constant image")) ;

  // after all these calculation, we are going to do it again?
  for (y = 0 ; y < h ; y++)
  {

    for (x = 0 ; x < w ; x++)
    {
      switch (src_slice_direction)
      {
      case MRI_CORONAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = (xsign> 0) ? x : (w - 1 -x);
          ym = (ysign> 0) ? y : (h - 1 -y);
          zm = slice ;
          break ;
        case MRI_SAGITTAL:
          xm = slice ;
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        case MRI_HORIZONTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = slice ;
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        }
        break;
      case MRI_SAGITTAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = slice;
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        case MRI_SAGITTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = slice ;
          break ;
        case MRI_HORIZONTAL:
          xm = (ysign>0) ? y : (h - 1 - y);
          ym = slice;
          zm = (xsign>0) ? x : (w - 1 - x);
          break ;
        }
        break;
      case MRI_HORIZONTAL:
        switch (view)   /* calculate coordinates in MR structure */
        {
        case MRI_CORONAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = slice ;
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        case MRI_SAGITTAL:
          xm = slice ;
          ym = (xsign>0) ? x : (w - 1 - x);
          zm = (ysign>0) ? y : (h - 1 - y);
          break ;
        case MRI_HORIZONTAL:
          xm = (xsign>0) ? x : (w - 1 - x);
          ym = (ysign>0) ? y : (h - 1 - y);
          zm = slice ;
          break ;
        }
        break;
      }
      MRIsampleVolumeFrame(mri, xm, ym, zm, frame, &val) ;
      yp = h - (y+1) ;   /* hips coordinate system is inverted */
      if (format == PFBYTE)
        *IMAGEpix(I, x, yp) = (byte)(255.0 * (val - fmin) / (fmax - fmin)) ;
      else
        *IMAGEFpix(I, x, yp) = (byte)(255.0 * (val - fmin) / (fmax - fmin)) ;
    }
  }

  return(I) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Convert a slice of an MRI data structure into a HIPS image.
  ------------------------------------------------------*/
IMAGE *
MRItoImage(MRI *mri, IMAGE *I, int slice)
{
  int           width, height, y, yp ;

  width = mri->width ;
  height = mri->height ;
  if (slice < 0 || slice >= mri->depth)
    ErrorReturn(NULL, (ERROR_BADPARM, "MRItoImage: bad slice %d\n", slice)) ;

  if (!I)
  {
    I = ImageAlloc(height, width,
                   mri->type == MRI_UCHAR ? PFBYTE :
                   mri->type == MRI_INT   ? PFINT :
                   mri->type == MRI_FLOAT ? PFFLOAT : PFBYTE, 1) ;
  }

  for (y = 0 ; y < height ; y++)
  {
    yp = height - (y+1) ;

#if 0
    if (mri->type == MRI_UCHAR)
      image_to_buffer(I->image, mri, slice) ;
    else
#endif
      switch (mri->type)
      {
      case MRI_INT:
        memmove
        (IMAGEIpix(I, 0, yp),mri->slices[slice][y],width*sizeof(int)) ;
        break ;
      case MRI_FLOAT:
        memmove
        (IMAGEFpix(I, 0, yp),mri->slices[slice][y],width*sizeof(float));
        break ;
      case MRI_UCHAR:
        memmove(IMAGEpix(I, 0, yp), mri->slices[slice][y],
               width*sizeof(unsigned char)) ;
        break ;
      default:
      case MRI_LONG:
        ImageFree(&I) ;
        ErrorReturn
        (NULL,
         (ERROR_UNSUPPORTED,
          "MRItoImage: unsupported type %d",
          mri->type)) ;
        break ;
      }
  }

  return(I) ;
}

MRI *
ImageToMRI(IMAGE *I)
{
  MRI *mri;
  int  width, height, depth, type, nframes, y, yp, x ;
  int frames;

  type = MRI_UCHAR; // to make compiler happy
  width = I->ocols;
  height = I->orows;
  depth = 1;
  nframes = I->num_frame;

  switch (I->pixel_format)
  {
  case PFBYTE:
    type = MRI_UCHAR;
    break;
  case PFSHORT:
    type = MRI_INT;
    break;
  case PFINT:
    type = MRI_INT;
    break;
  case PFFLOAT:
    type = MRI_FLOAT;
    break;
  case PFDOUBLE:
  case PFCOMPLEX:
  default:
    ErrorExit
    (ERROR_BADPARM, "IMAGE type = %d not supported\n", I->pixel_format);
    break;
  }
  // allocate memory
  if (nframes <= 1)
    mri = MRIalloc(width, height, depth, type);
  else
    mri = MRIallocSequence(width, height, depth, type, nframes);

  // just fake the size
  mri->nframes = nframes;
  mri->imnr0 = 1;
  mri->imnr1 = 1;
  mri->thick = mri->ps = 1.0;
  mri->xsize = mri->ysize = mri->zsize = 1.0;
  mri->xend = mri->width  * mri->xsize / 2.0;
  mri->xstart = -mri->xend;
  mri->yend = mri->height * mri->ysize / 2.0;
  mri->ystart = -mri->yend;
  mri->zend = mri->depth  * mri->zsize / 2.0;
  mri->zstart = -mri->zend;
  mri->fov = ((mri->xend - mri->xstart) > (mri->yend - mri->ystart) ?
              (mri->xend - mri->xstart) : (mri->yend - mri->ystart));
  // set orientation to be coronal
  mri->x_r = -1;
  mri->y_r =  0;
  mri->z_r =  0;
  mri->c_r = 0;
  mri->x_a =  0;
  mri->y_a =  0;
  mri->z_a =  1;
  mri->c_a = 0;
  mri->x_s =  0;
  mri->y_s = -1;
  mri->z_s =  0;
  mri->c_s = 0;

  // hips coordinate system is inverted
  for (frames = 0; frames < nframes; ++frames)
  {

    for (y = 0 ; y < height ; y++)
    {
      yp = height - (y+1) ;

      for (x = 0; x < width; ++x)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          MRIseq_vox(mri, x, y, 0, frames) = *IMAGEpix(I, x, yp);
          break ;
        case MRI_SHORT:
          MRISseq_vox(mri, x, y, 0, frames) = *IMAGESpix(I, x, yp);
          break;
        case MRI_INT:
          {
            int val ;

            if (*IMAGESpix(I, x, yp) < 0)
              DiagBreak() ;
            if (I->pixel_format == PFSHORT)
              MRIIseq_vox(mri, x, y, 0, frames) = 
                (int)((unsigned short)(*IMAGESpix(I, x, yp)));
            else
              MRIIseq_vox(mri, x, y, 0, frames) = *IMAGEIpix(I, x, yp);
            val = MRIIseq_vox(mri, x, y, 0, frames);
            break ;
          }
        case MRI_FLOAT:
          MRIFseq_vox(mri, x, y, 0, frames) = *IMAGEFpix(I, x, yp);
          break ;
        default:
          ErrorReturn
          (NULL,
           (ERROR_UNSUPPORTED,
            "MRItoImage: unsupported type %d",
            mri->type)) ;
          break ;
        }
      }
    }
  } // frames
  return(mri) ;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Build an MRI from all values s.t. min_val <= val <= max_val
  ------------------------------------------------------*/
MRI *
MRIextractValues(MRI *mri_src, MRI *mri_dst, float min_val, float max_val)
{
  BUFTYPE  *psrc, *pdst, val ;
  float    *pfsrc, *pfdst, fval ;
  int      frame, x, y, z, width, height, depth ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (frame = 0 ; frame < mri_src->nframes ; frame++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        switch (mri_src->type)
        {
        case MRI_UCHAR:
          psrc = &MRIseq_vox(mri_src, 0, y, z, frame) ;
          pdst = &MRIseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            val = *psrc++ ;
            if ((val < min_val) || (val > max_val))
              val = 0 ;
            *pdst++ = val ;
          }
          break ;
        case MRI_FLOAT:
          pfsrc = &MRIFseq_vox(mri_src, 0, y, z, frame) ;
          pfdst = &MRIFseq_vox(mri_dst, 0, y, z, frame) ;
          for (x = 0 ; x < width ; x++)
          {
            fval = *pfsrc++ ;
            if ((fval < min_val) || (fval > max_val))
              fval = 0.0f ;
            *pfdst++ = fval ;
          }
          break ;
        default:
          ErrorReturn
          (NULL,
           (ERROR_UNSUPPORTED,
            "MRIextractValues: unsupported type %d",mri_src->type));
        }
      }
    }
  }

  return(mri_dst) ;
}


/*-----------------------------------------------------
  Wrapper around MRIupsampleN for N=2
  ------------------------------------------------------*/
MRI *
MRIupsample2(MRI *mri_src, MRI *mri_dst)
{
  return(MRIupsampleN(mri_src, mri_dst, 2)) ;
}
/*-----------------------------------------------------
  MRI *MRIupsampleN(MRI *mri_src, MRI *mri_dst, int N)
  Upsample volume by integer factor. No error checking, upsample
  factor must be valid.
  ------------------------------------------------------*/
MRI *MRIupsampleN(MRI *mri_src, MRI *mri_dst, int N)
{
  int     width, depth, height, x, y, z, f ;
  double val;
  MATRIX *Vox2RAS,*CRS0,*RAS0;

  width  = N*mri_src->width ;
  height = N*mri_src->height ;
  depth  = N*mri_src->depth ;

  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, mri_src->type, mri_src->nframes) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  // Recompute geometry for finer resolution
  // Only the xsize and cras change
  mri_dst->xsize = mri_src->xsize/N;
  mri_dst->ysize = mri_src->ysize/N;
  mri_dst->zsize = mri_src->zsize/N;
  mri_dst->x_r = mri_src->x_r;
  mri_dst->x_a = mri_src->x_a;
  mri_dst->x_s = mri_src->x_s;
  mri_dst->y_r = mri_src->y_r;
  mri_dst->y_a = mri_src->y_a;
  mri_dst->y_s = mri_src->y_s;
  mri_dst->z_r = mri_src->z_r;
  mri_dst->z_a = mri_src->z_a;
  mri_dst->z_s = mri_src->z_s;
  mri_dst->xstart = mri_src->xstart;
  mri_dst->ystart = mri_src->ystart;
  mri_dst->zstart = mri_src->zstart;
  mri_dst->xend = mri_src->xend;
  mri_dst->yend = mri_src->yend;
  mri_dst->zend = mri_src->zend;
  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;

  // Computes CRAS based on location of the 1st voxel in upsampled space
  // The new location is 1/(2*Nth) of a voxel from the corner. The RAS of
  // a voxel is at the center of the voxel (unfortunately)
  Vox2RAS = MRIxfmCRS2XYZ(mri_src,0); // scanner vox2ras of source mri
  CRS0 = MatrixZero(4,1,NULL);
  CRS0->rptr[1][1] = -1.0/(2*N);
  CRS0->rptr[2][1] = -1.0/(2*N);
  CRS0->rptr[3][1] = -1.0/(2*N);
  CRS0->rptr[4][1] =  1.0;
  RAS0 = MatrixMultiply(Vox2RAS,CRS0,NULL);
  MRIp0ToCRAS(mri_dst, RAS0->rptr[1][1],RAS0->rptr[2][1],RAS0->rptr[3][1]);
  MatrixFree(&Vox2RAS);
  MatrixFree(&CRS0);
  MatrixFree(&RAS0);

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
	for(f = 0; f < mri_src->nframes; f++){
	  val = MRIgetVoxVal(mri_src, x/N, y/N, z/N, f);
	  MRIsetVoxVal(mri_dst,x,y,z,f,val);
	}
      }
    }
  }
  MRIreInitCache(mri_dst) ;

  return(mri_dst) ;
}


MRI *MRIdownsample2LabeledVolume(MRI *mri_src, MRI *mri_dst)
{
  int     width, depth, height, x, y, z, x1, y1, z1, counts[256], label,
  max_count, out_label ;
  BUFTYPE *psrc ;

  if (mri_src->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,
                 "MRIdownsample2LabeledVolume: source must be UCHAR"));

  width = mri_src->width/2 ;
  height = mri_src->height/2 ;
  depth = mri_src->depth/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  MRIclear(mri_dst) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        memset(counts, 0, sizeof(counts)) ;
        if (x == 96 && y == 66 && z == 56)
          DiagBreak() ;
        for (z1 = 2*z ; z1 <= 2*z+1 ; z1++)
        {
          for (y1 = 2*y ; y1 <= 2*y+1 ; y1++)
          {
            psrc = &MRIvox(mri_src, 2*x, y1, z1) ;
            for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
            {
              label = *psrc++ ;
              counts[label]++ ;
            }
          }
        }
        for (out_label = label = 0, max_count = counts[0] ;
             label <= 255 ;
             label++)
        {
          if (counts[label] > max_count)
          {
            out_label = label ;
            max_count = counts[label] ;
          }
        }
        MRIvox(mri_dst, x, y, z) = out_label ;
      }
    }
  }

  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;
  mri_dst->xsize = mri_src->xsize*2 ;
  mri_dst->ysize = mri_src->ysize*2 ;
  mri_dst->zsize = mri_src->zsize*2 ;
  mri_dst->thick = mri_src->thick*2 ;
  mri_dst->ps    = mri_src->ps*2 ;

  VECTOR* C = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(C,1) = mri_src->width/2+0.5;
  VECTOR_ELT(C,2) = mri_src->height/2+0.5;
  VECTOR_ELT(C,3) = mri_src->depth/2+0.5;
  VECTOR_ELT(C,4) = 1.0;
  MATRIX* V2R     = extract_i_to_r(mri_src);
  MATRIX* P       = MatrixMultiply(V2R,C,NULL);
  mri_dst->c_r    = P->rptr[1][1];
  mri_dst->c_a    = P->rptr[2][1];
  mri_dst->c_s    = P->rptr[3][1];
  MatrixFree(&P);
  MatrixFree(&V2R);
  VectorFree(&C);
  
  MRIreInitCache(mri_dst) ;

  //mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIdownsample2(MRI *mri_src, MRI *mri_dst)
{
  int     width, depth, height, x, y, z, x1, y1, z1 ;
  BUFTYPE *psrc ;
  short   *pssrc ;
  float   *pfsrc ;
  float   val ;

  if (mri_dst && mri_src->type != mri_dst->type)
    ErrorReturn
    (NULL,
     (ERROR_UNSUPPORTED,
      "MRIdownsample2: source and dst must be same type"));

  width  = mri_src->width/2 ;
  height = mri_src->height/2 ;
  depth  = mri_src->depth/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, depth, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  MRIclear(mri_dst) ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        for (val = 0.0f, z1 = 2*z ; z1 <= 2*z+1 ; z1++)
        {
          for (y1 = 2*y ; y1 <= 2*y+1 ; y1++)
          {
            switch (mri_src->type)
            {
            case MRI_UCHAR:
              psrc = &MRIvox(mri_src, 2*x, y1, z1) ;
              for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
                val += *psrc++ ;
              break ;
            case MRI_SHORT:
              pssrc = &MRISvox(mri_src, 2*x, y1, z1) ;
              for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
                val += *pssrc++ ;
              break ;
            case MRI_FLOAT:
              pfsrc = &MRIFvox(mri_src, 2*x, y1, z1) ;
              for (x1 = 2*x ; x1 <= 2*x+1 ; x1++)
                val += *pfsrc++ ;
              break ;
            default:
              ErrorReturn
              (NULL,
               (ERROR_UNSUPPORTED,
                "MRIdownsample2: unsupported input type %d",
                mri_src->type));
            }
          }
        }
        switch (mri_src->type)
        {
        case MRI_UCHAR:
          MRIvox(mri_dst, x, y, z) = (BUFTYPE)nint(val/8.0f) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_dst, x, y, z) = val/8.0f ;
          break ;
        case MRI_SHORT:
          MRISvox(mri_dst, x, y, z) = (short)nint(val/8.0f) ;
          break ;
        }
      }
    }
  }

  mri_dst->imnr0 = mri_src->imnr0 ;
  mri_dst->imnr1 = mri_src->imnr0 + mri_dst->depth - 1 ;
  mri_dst->xsize = mri_src->xsize*2 ;
  mri_dst->ysize = mri_src->ysize*2 ;
  mri_dst->zsize = mri_src->zsize*2 ;
  mri_dst->thick = mri_src->thick*2 ;
  mri_dst->ps    = mri_src->ps*2 ;
  
  // adjust cras
  //printf("COMPUTING new CRAS\n") ;
  VECTOR* C = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(C,1) = mri_src->width/2+0.5;
  VECTOR_ELT(C,2) = mri_src->height/2+0.5;
  VECTOR_ELT(C,3) = mri_src->depth/2+0.5;
  VECTOR_ELT(C,4) = 1.0;
  MATRIX* V2R     = extract_i_to_r(mri_src);
  MATRIX* P       = MatrixMultiply(V2R,C,NULL);
  mri_dst->c_r    = P->rptr[1][1];
  mri_dst->c_a    = P->rptr[2][1];
  mri_dst->c_s    = P->rptr[3][1];
  MatrixFree(&P);
  MatrixFree(&V2R);
  VectorFree(&C);
  
  MRIreInitCache(mri_dst) ;
  //printf("CRAS new: %2.3f %2.3f %2.3f\n",mri_dst->c_r,mri_dst->c_a,mri_dst->c_s);

  //  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-------------------------------------------------------------
  MRI MRIvalueFill(MRI *mri, float value) -- fills an mri volume
  with the given value. Type of mri can be anything.
  -------------------------------------------------------------*/
int MRIvalueFill(MRI *mri, float value)
{
  int c,r,s,f;

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        for (f=0; f < mri->nframes; f++)
        {
          switch (mri->type)
          {
          case MRI_UCHAR:
            MRIseq_vox(mri,c,r,s,f) = (unsigned char) (nint(value));
            break;
          case MRI_SHORT:
            MRISseq_vox(mri,c,r,s,f) = (short) (nint(value));
            break;
          case MRI_INT:
            MRIIseq_vox(mri,c,r,s,f) = (int) (nint(value));
            break;
          case MRI_FLOAT:
            MRIFseq_vox(mri,c,r,s,f) = value;
            break;
          }
        }
      }
    }
  }

  return(0);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Iteratively set all voxels in mri_dst that neighbor
  a voxel that has already been filled (starting with the seed),
  and for which the corresponding voxel in mri_src is
  below threshold.
  ------------------------------------------------------*/
MRI *
MRIfill(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
        int seed_z, int threshold, int fill_val, int max_count)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
  zmin, zmax, on, x0, x1, y0, y1, z0, z1 ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = height-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = depth-1 ;

  xmin = MAX(seed_x-1,0) ;
  xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ;
  ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ;
  zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  do
  {
    z0 = zmin ;
    z1 = zmax ;
    y0 = ymin ;
    y1 = ymax ;
    x0 = xmin ;
    x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          val = *psrc ;
          if ((val > threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ;
    zmax = z1 ;
    ymin = y0 ;
    ymax = y1 ;
    xmin = x0 ;
    xmax = x1 ;
    /*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
    if ((max_count > 0) && (nfilled > max_count))
      break ;
  }
  while (nfilled > 0) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIfillFG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
          int seed_z, int threshold, int fill_val, int *npix)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
  zmin, zmax, on, x0, x1, y0, y1, z0, z1, total_pix ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = width-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = width-1 ;

  xmin = MAX(seed_x-1,0) ;
  xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ;
  ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ;
  zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  total_pix = 1 ;  /* include the seed point */
  do
  {
    z0 = zmin ;
    z1 = zmax ;
    y0 = ymin ;
    y1 = ymax ;
    x0 = xmin ;
    x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          val = *psrc ;
          if ((val < threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ;
    zmax = z1 ;
    ymin = y0 ;
    ymax = y1 ;
    xmin = x0 ;
    xmax = x1 ;
    total_pix += nfilled ;
    /*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
  }
  while (nfilled > 0) ;

  if (npix)
    *npix = total_pix ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIfillBG(MRI *mri_src, MRI *mri_dst, int seed_x, int seed_y,
          int seed_z, int threshold, int fill_val, int *npix)
{
  int     width, height, depth, x, y, z, nfilled, xmin, xmax, ymin, ymax,
  zmin, zmax, on, x0, x1, y0, y1, z0, z1, total_pix ;
  BUFTYPE *psrc, *pdst, val ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  if (seed_x < 0)
    seed_x = 0 ;
  else if (seed_x >= width)
    seed_x = width-1 ;
  if (seed_y < 0)
    seed_y = 0 ;
  else if (seed_y >= height)
    seed_y = width-1 ;
  if (seed_z < 0)
    seed_z = 0 ;
  else if (seed_z >= depth)
    seed_z = width-1 ;

  xmin = MAX(seed_x-1,0) ;
  xmax = MIN(seed_x+1, width-1) ;
  ymin = MAX(seed_y-1,0) ;
  ymax = MIN(seed_y+1, height-1) ;
  zmin = MAX(seed_z-1,0) ;
  zmax = MIN(seed_z+1, depth-1) ;

  /* replace all occurrences of fill_val with fill_val-1 */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++, pdst++)
      {
        val = *pdst ;
        if (val == fill_val)
          *pdst = val-1 ;
      }
    }
  }

  MRIvox(mri_dst, seed_x, seed_y, seed_z) = fill_val ;
  total_pix = 0 ;
  do
  {
    z0 = zmin ;
    z1 = zmax ;
    y0 = ymin ;
    y1 = ymax ;
    x0 = xmin ;
    x1 = xmax ;
    nfilled = 0 ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        psrc = &MRIvox(mri_src, xmin, y, z) ;
        pdst = &MRIvox(mri_dst, xmin, y, z) ;
        for (x = xmin ; x <= xmax ; x++, psrc++, pdst++)
        {
          if (z == 130 && (y == 145 || y == 146) && x == 142)
            DiagBreak() ;
          val = *psrc ;
          if ((val > threshold) || (*pdst == fill_val))
            continue ;

          on = 0 ;
          if (x > 0)
            on = (MRIvox(mri_dst, x-1, y, z) == fill_val) ;
          if (!on && (x < width-1))
            on = (MRIvox(mri_dst, x+1, y, z) == fill_val) ;
          if (!on && (y > 0))
            on = (MRIvox(mri_dst, x, y-1, z) == fill_val) ;
          if (!on && (y < height-1))
            on = (MRIvox(mri_dst, x, y+1, z) == fill_val) ;
          if (!on && (z > 0))
            on = (MRIvox(mri_dst, x, y, z-1) == fill_val) ;
          if (!on && (z < depth-1))
            on = (MRIvox(mri_dst, x, y, z+1) == fill_val) ;
          if (on)
          {
            if (z == 130 && (y == 145 || y == 146) && x == 142)
              DiagBreak() ;
            if (x <= x0)
              x0 = MAX(x-1,0) ;
            if (x >= x1)
              x1 = MIN(x+1,width-1) ;
            if (y <= y0)
              y0 = MAX(y-1,0) ;
            if (y >= y1)
              y1 = MIN(y+1,height-1) ;
            if (z <= z0)
              z0 = MAX(z-1,0) ;
            if (z >= z1)
              z1 = MIN(z+1,depth-1) ;
            nfilled++ ;
            *pdst = fill_val ;
          }
        }
      }
    }
    zmin = z0 ;
    zmax = z1 ;
    ymin = y0 ;
    ymax = y1 ;
    xmin = x0 ;
    xmax = x1 ;
    total_pix += nfilled ;
    /*    fprintf(stderr, "# filled = %d\n", nfilled) ;*/
  }
  while (nfilled > 0) ;

  if (npix)
    *npix = total_pix ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIsetResolution(MRI *mri, float xres, float yres, float zres)
{
  mri->ps = (xres+yres+zres) / 3.0f ;
  mri->xsize = xres ;
  mri->ysize = yres ;
  mri->zsize = zres ;
  mri->xstart *= xres ;
  mri->xend *= xres ;
  mri->ystart *= yres ;
  mri->yend *= yres ;
  mri->zstart *= zres ;
  mri->zend *= zres ;
  mri->thick = MAX(MAX(xres,yres),zres) ;
#if 0
  mri->x_r = mri->x_r * xres ;
  mri->x_a = mri->x_a * xres ;
  mri->x_s = mri->x_s * xres ;
  mri->y_r = mri->y_r * yres ;
  mri->y_a = mri->y_a * yres ;
  mri->y_s = mri->y_s * yres ;
  mri->z_r = mri->z_r * zres ;
  mri->z_a = mri->z_a * zres ;
  mri->z_s = mri->z_s * zres ;
#endif
  MRIreInitCache(mri) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIsetTransform(MRI *mri,   General_transform *transform)
{
  mri->transform = *transform ;
  mri->linear_transform = get_linear_transform_ptr(transform) ;
  mri->inverse_linear_transform = get_inverse_linear_transform_ptr(transform) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIextractTalairachPlane(MRI *mri_src, MRI *mri_dst, int orientation,
                         int x, int y, int z, int wsize)
{
  double     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;
  double     ex, ey, ez, len, x0, y0, z0 ;

  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }

  MRIvoxelToTalairachVoxel(mri_src, x, y, z, &x0, &y0, &z0) ;
  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (double)x0+1 ;
    ey = (double)y0 ;
    ez = (double)z0 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x0 ;
    ey = (double)y0+1 ;
    ez = (double)z0 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (double)x0+1 ;
    ey = (double)y0 ;
    ez = (double)z0 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x0 ;
    ey = (double)y0 ;
    ez = (double)z0+1 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (double)x0 ;
    ey = (double)y0 ;
    ez = (double)z0+1.0 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x0 ;
    ey = (double)y0+1.0 ;
    ez = (double)z0 ;
    MRItalairachVoxelToVoxel(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  }

  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
  /*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
  /*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIextractArbitraryPlane(MRI *mri_src, MRI *mri_dst,
                         double e1_x, double e1_y, double e1_z,
                         double e2_x, double e2_y, double e2_z,
                         int x, int y, int z, int wsize)
{
  double     xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;

  whalf = (wsize-1)/2 ;

  if (!mri_dst)
  {
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->xstart = x-whalf*mri_dst->xsize ;
    mri_dst->ystart = y-whalf*mri_dst->ysize ;
    mri_dst->zstart = z-whalf*mri_dst->zsize ;
    mri_dst->xend = mri_dst->xstart + wsize*mri_dst->xsize ;
    mri_dst->yend = mri_dst->ystart + wsize*mri_dst->ysize ;
    mri_dst->zend = mri_dst->zstart + wsize*mri_dst->zsize ;
    mri_dst->imnr0 = z + mri_src->imnr0 ;
    mri_dst->imnr1 = mri_dst->imnr0 ;
  }

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk*e1_x)] ;
      yi = mri_src->yi[nint(ybase + xk*e1_y)] ;
      zi = mri_src->zi[nint(zbase + xk*e1_z)] ;
      MRIvox(mri_dst, xk+whalf,yk+whalf,0) = MRIvox(mri_src, xi, yi, zi) ;
    }
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIeraseTalairachPlaneNew(MRI *mri, MRI *mri_mask, int orientation, int x,
                          int y, int z, int wsize, int fill_val)
{
  double     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi, xki, yki, x0, y0 ;
  double     ex, ey, ez, len, xt0, yt0, zt0 ;

  whalf = (wsize-1)/2 ;

  x0 = mri_mask->width/2 ;
  y0 = mri_mask->height/2 ;
  MRIvoxelToTalairachVoxel(mri, x, y, z, &xt0, &yt0, &zt0) ;
  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (double)xt0+1 ;
    ey = (double)yt0 ;
    ez = (double)zt0 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)xt0 ;
    ey = (double)yt0+1 ;
    ez = (double)zt0 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (double)xt0+1 ;
    ey = (double)yt0 ;
    ez = (double)zt0 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)xt0 ;
    ey = (double)yt0 ;
    ez = (double)zt0+1 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (double)xt0 ;
    ey = (double)yt0 ;
    ez = (double)zt0+1.0 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)xt0 ;
    ey = (double)yt0+1.0 ;
    ez = (double)zt0 ;
    MRItalairachVoxelToVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  }

  /*
    don't want to normalize basis - they are orthonormal in magnet space,
    not necessarily Talairach space.
  */
  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
  /*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
  /*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri->xi[nint(xbase + xk*e1_x)] ;
      yi = mri->yi[nint(ybase + xk*e1_y)] ;
      zi = mri->zi[nint(zbase + xk*e1_z)] ;
      xki = mri_mask->xi[xk+x0] ;
      yki = mri_mask->yi[yk+y0] ;
      if (MRIvox(mri_mask, xki, yki,0))
        MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIeraseTalairachPlane(MRI *mri, MRI *mri_mask, int orientation, int x, int y,
                       int z, int wsize, int fill_val)
{
  double     e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase ;
  int      whalf, xk, yk, xi, yi, zi ;
  double     ex, ey, ez, len ;

  whalf = (wsize-1)/2 ;

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    /* the 'x' basis vector */
    ex = (double)x+1 ;
    ey = (double)y ;
    ez = (double)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x ;
    ey = (double)y+1 ;
    ez = (double)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    /* the 'x' basis vector */
    ex = (double)x+1 ;
    ey = (double)y ;
    ez = (double)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x ;
    ey = (double)y ;
    ez = (double)z+1 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    /* the 'x' basis vector */
    ex = (double)x ;
    ey = (double)y ;
    ez = (double)z+1.0 ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z) ;
    e1_x -= (double)x ;
    e1_y -= (double)y ;
    e1_z -= (double)z ;

    /* the 'y' basis vector */
    ex = (double)x ;
    ey = (double)y+1.0 ;
    ez = (double)z ;
    MRIvoxelToTalairachVoxel(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z) ;
    e2_x -= (double)x ;
    e2_y -= (double)y ;
    e2_z -= (double)z ;
    break ;
  }

  /*
    don't want to normalize basis - they are orthonormal in magnet space,
    not necessarily Talairach space.
  */
  len = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
  /*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  len = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
  /*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf ; yk <= whalf ; yk++)
  {
    xbase = (float)x + (float)yk * e2_x ;
    ybase = (float)y + (float)yk * e2_y ;
    zbase = (float)z + (float)yk * e2_z ;
    for (xk = -whalf ; xk <= whalf ; xk++)
    {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri->xi[nint(xbase + xk*e1_x)] ;
      yi = mri->yi[nint(ybase + xk*e1_y)] ;
      zi = mri->zi[nint(zbase + xk*e1_z)] ;
      if (MRIvox(mri_mask, xk+whalf,yk+whalf,0))
        MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIextractPlane(MRI *mri_src, MRI *mri_dst, int orientation, int where)
{
  int      x, y, z, width, height ;

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    width = mri_src->width ;
    height = mri_src->height ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    width = mri_src->width ;
    height = mri_src->depth ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    width = mri_src->depth ;
    height = mri_src->height ;
    break ;
  }

  if (!mri_dst)
  {
    mri_dst = MRIalloc(width, height, 1, MRI_UCHAR) ;
    MRIcopyHeader(mri_src, mri_dst) ;
    mri_dst->zstart = where ;
    mri_dst->zend = where+1 ;
    mri_dst->imnr0 = 0 ;
    mri_dst->imnr1 = 1 ;
  }
  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
        MRIvox(mri_dst, x, y, 0) = MRIvox(mri_src, x, y, where) ;
    }
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
        MRIvox(mri_dst, x, z, 0) = MRIvox(mri_src, x, where, z) ;
    }
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    for (z = 0 ; z < mri_src->depth ; z++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
        MRIvox(mri_dst, z, y, 0) = MRIvox(mri_src, where, y, z) ;
    }
    break ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI *
MRIfillPlane
(MRI *mri_mask, MRI *mri_dst, int orientation, int where, int fillval)
{
  int      x, y, z, width, height ;

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    width = mri_mask->width ;
    height = mri_mask->height ;
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    width = mri_mask->width ;
    height = mri_mask->depth ;
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    width = mri_mask->depth ;
    height = mri_mask->height ;
    break ;
  }

  switch (orientation)
  {
  default:
  case MRI_CORONAL:   /* basis vectors in x-y plane */
    for (x = 0 ; x < mri_dst->width ; x++)
    {
      for (y = 0 ; y < mri_dst->height ; y++)
        if (MRIvox(mri_mask, x, y, 0))
          MRIvox(mri_dst, x, y, where) = fillval ;
    }
    break ;
  case MRI_HORIZONTAL:  /* basis vectors in x-z plane */
    for (x = 0 ; x < mri_dst->width ; x++)
    {
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (MRIvox(mri_mask, x, z, 0))
          MRIvox(mri_dst, x, where, z) = fillval ;
      }
    }
    break ;
  case MRI_SAGITTAL:    /* basis vectors in y-z plane */
    for (z = 0 ; z < mri_dst->depth ; z++)
    {
      for (y = 0 ; y < mri_dst->height ; y++)
        if (MRIvox(mri_mask, z, y, 0))
          MRIvox(mri_dst, where, y, z) = fillval ;
    }
    break ;
  }

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIerasePlane(MRI *mri, float x0, float y0, float z0,
              float dx, float dy, float dz, int fill_val)
{
  int      *pxi, *pyi, *pzi, xi, yi, zi, x, y, z ;
  float    e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, maxlen, l1, l2 ;

  maxlen = MAX(mri->width, mri->height) ;
  maxlen = MAX(maxlen, mri->depth) ;

  /* don't care about sign (right-hand rule) */
  e1_x = dz*dz - dx*dy ;
  e1_y = dx*dx - dy*dz ;
  e1_z = dy*dy - dz*dx ;
  l1 = sqrt(e1_x*e1_x + e1_y*e1_y + e1_z*e1_z) ;
  e1_x /= l1 ;
  e1_y /= l1 ;
  e1_z /= l1 ;

  e2_x = e1_y*dz - e1_z*dy ;
  e2_y = e1_x*dz - e1_z*dx ;
  e2_z = e1_y*dx - e1_x*dy ;
  l2 = sqrt(e2_x*e2_x + e2_y*e2_y + e2_z*e2_z) ;
  e2_x /= l2 ;
  e2_y /= l2 ;
  e2_z /= l2 ;

  pxi = mri->xi ;
  pyi = mri->yi ;
  pzi = mri->zi ;
  maxlen *= 1.5 ;  /* make sure to get entire extent */
  for (l1 = -maxlen/2 ; l1 <= maxlen/2 ; l1 += 0.5f)
  {
    for (l2 = -maxlen/2 ; l2 <= maxlen/2 ; l2 += 0.5f)
    {
      x = nint(x0 + l1 * e1_x + l2 * e2_x) ;
      xi = pxi[x] ;
      y = nint(y0 + l1 * e1_y + l2 * e2_y) ;
      yi = pyi[y] ;
      z = nint(z0 + l1 * e1_z + l2 * e2_z) ;
      zi = pzi[z] ;
      MRIvox(mri, xi, yi, zi) = fill_val ;
    }
  }

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume to cubic voxels.
  ------------------------------------------------------*/
MRI *
MRIinterpolate(MRI *mri_src, MRI *mri_dst)
{
  int     xs, ys, zs, xd, yd, zd, max_dim, xorig, yorig, zorig, dorig ;
  float   sx, sy, sz, psize ;
  int     width, height, depth, i ;
  float   xsmd, ysmd, zsmd, xspd, yspd, zspd, weights[8], fout ;
  int     xsp, xsm, ysp, ysm, zsp, zsm ;  /* surrounding coordinates */
  float vals[8], outval ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth;
  if (width > height)
  {
    max_dim = width > depth ? width : depth ;
    psize = width > depth ? mri_src->xsize : mri_src->zsize ;
  }
  else
  {
    max_dim = height > depth ? height : depth ;
    psize = height > depth ? mri_src->ysize : mri_src->zsize ;
  }

  if (!mri_dst)
  {
    mri_dst = MRIalloc(max_dim, max_dim, max_dim, mri_src->type) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  mri_dst->xsize = mri_dst->ysize = mri_dst->zsize = psize ;
  sx = (float)width / (float)max_dim  ;
  sy = (float)height / (float)max_dim  ;
  sz = (float)depth / (float)max_dim  ;

  mri_dst->imnr0 = 1;
  mri_dst->imnr1 = mri_dst->depth;
  mri_dst->xstart = -mri_dst->xsize * (double)mri_dst->width / 2.0;
  mri_dst->xend = -mri_dst->xstart;
  mri_dst->ystart = -mri_dst->ysize * (double)mri_dst->height / 2.0;
  mri_dst->yend = -mri_dst->ystart;
  mri_dst->zstart = -mri_dst->zsize * (double)mri_dst->depth / 2.0;
  mri_dst->zend = -mri_dst->zstart;

  xorig = (width-1)/2 ;
  yorig = (height-1)/2 ;
  zorig = (depth-1)/2 ;
  dorig = (max_dim-1)/2 ;

  /*
    for each output voxel, find the 8 nearest input voxels and interpolate
    the output voxel from them.
  */
  for (zd = 0 ; zd < max_dim ; zd++)
  {
    printf("interpolate: %d/%d\n",zd+1,max_dim);
    for (yd = 0 ; yd < max_dim ; yd++)
    {
      for (xd = 0 ; xd < max_dim ; xd++)
      {
        /* do trilinear interpolation here */
        xs = sx*(xd-dorig) + xorig ;
        ys = sy*(yd-dorig) + yorig ;
        zs = sz*(zd-dorig) + zorig ;

        /*
          these boundary conditions will cause
          reflection across the border
          for the 1st negative pixel.
        */
        if (xs > -1 && xs < width &&
            ys > -1 && ys < height &&
            zs > -1 && zs < depth)
        {
          xsm = (int)xs ;
          xsp = MIN(width-1, xsm+1) ;
          ysm = (int)ys ;
          ysp = MIN(height-1, ysm+1) ;
          zsm = (int)zs ;
          zsp = MIN(depth-1, zsm+1) ;
          xsmd = xs - (float)xsm ;
          ysmd = ys - (float)ysm ;
          zsmd = zs - (float)zsm ;
          xspd = (1.0f - xsmd) ;
          yspd = (1.0f - ysmd) ;
          zspd = (1.0f - zsmd) ;

          /*          vals[0] = mri_src->slices[zsm][ysm][xsm] ;
                      vals[1] = mri_src->slices[zsm][ysm][xsp] ;
                      vals[2] = mri_src->slices[zsm][ysp][xsm] ;
                      vals[3] = mri_src->slices[zsm][ysp][xsp] ;
                      vals[4] = mri_src->slices[zsp][ysm][xsm] ;
                      vals[5] = mri_src->slices[zsp][ysm][xsp] ;
                      vals[6] = mri_src->slices[zsp][ysp][xsm] ;
                      vals[7] = mri_src->slices[zsp][ysp][xsp] ;
          */
          /* different vox types... */
          vals[0] = MRIFvox(mri_src, xsm, ysm, zsm);
          vals[1] = MRIFvox(mri_src, xsp, ysm, zsm);
          vals[2] = MRIFvox(mri_src, xsm, ysp, zsm);
          vals[3] = MRIFvox(mri_src, xsp, ysp, zsm);
          vals[4] = MRIFvox(mri_src, xsm, ysm, zsp);
          vals[5] = MRIFvox(mri_src, xsp, ysm, zsp);
          vals[6] = MRIFvox(mri_src, xsm, ysp, zsp);
          vals[7] = MRIFvox(mri_src, xsp, ysp, zsp);

          weights[0] = zsmd * ysmd * xsmd ;
          weights[1] = zsmd * ysmd * xspd ;
          weights[2] = zsmd * yspd * xsmd ;
          weights[3] = zsmd * yspd * xspd ;
          weights[4] = zspd * ysmd * xsmd ;
          weights[5] = zspd * ysmd * xspd ;
          weights[6] = zspd * yspd * xsmd ;
          weights[7] = zspd * yspd * xspd ;
          /*
            for(i = 0;i < 8;i++)
            printf("%d, %f, %f\n",i,vals[i], weights[i]);
          */
          for (fout = 0.0f, i = 0 ; i < 8 ; i++)
            fout += (float)vals[i] * weights[i] ;
          outval = (float)nint(fout) ;
          MRIvox(mri_dst, xd, yd, zd) = (BUFTYPE)nint(fout) ;

        }
      }
    }
  }

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIsampleVolumeFrame( const MRI *mri,
		      double x, double y, double z,
		      const int frame,
		      double *pval )
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  double val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    return(MRIsampleVolumeFrameType
           (mri, x, y, z, frame, SAMPLE_NEAREST, pval)) ;

  if (frame >= mri->nframes)
  {
    *pval = mri->outside_val ;
    return(NO_ERROR) ;
  }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguoulsy out of bounds */
    *pval = mri->outside_val ;
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = val =
              xpd * ypd * zpd * (double)MRIseq_vox(mri, xm, ym, zm, frame) +
              xpd * ypd * zmd * (double)MRIseq_vox(mri, xm, ym, zp, frame) +
              xpd * ymd * zpd * (double)MRIseq_vox(mri, xm, yp, zm, frame) +
              xpd * ymd * zmd * (double)MRIseq_vox(mri, xm, yp, zp, frame) +
              xmd * ypd * zpd * (double)MRIseq_vox(mri, xp, ym, zm, frame) +
              xmd * ypd * zmd * (double)MRIseq_vox(mri, xp, ym, zp, frame) +
              xmd * ymd * zpd * (double)MRIseq_vox(mri, xp, yp, zm, frame) +
              xmd * ymd * zmd * (double)MRIseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_FLOAT:
    *pval = val =
              xpd * ypd * zpd * (double)MRIFseq_vox(mri, xm, ym, zm, frame) +
              xpd * ypd * zmd * (double)MRIFseq_vox(mri, xm, ym, zp, frame) +
              xpd * ymd * zpd * (double)MRIFseq_vox(mri, xm, yp, zm, frame) +
              xpd * ymd * zmd * (double)MRIFseq_vox(mri, xm, yp, zp, frame) +
              xmd * ypd * zpd * (double)MRIFseq_vox(mri, xp, ym, zm, frame) +
              xmd * ypd * zmd * (double)MRIFseq_vox(mri, xp, ym, zp, frame) +
              xmd * ymd * zpd * (double)MRIFseq_vox(mri, xp, yp, zm, frame) +
              xmd * ymd * zmd * (double)MRIFseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_SHORT:
    *pval = val =
              xpd * ypd * zpd * (double)MRISseq_vox(mri, xm, ym, zm, frame) +
              xpd * ypd * zmd * (double)MRISseq_vox(mri, xm, ym, zp, frame) +
              xpd * ymd * zpd * (double)MRISseq_vox(mri, xm, yp, zm, frame) +
              xpd * ymd * zmd * (double)MRISseq_vox(mri, xm, yp, zp, frame) +
              xmd * ypd * zpd * (double)MRISseq_vox(mri, xp, ym, zm, frame) +
              xmd * ypd * zmd * (double)MRISseq_vox(mri, xp, ym, zp, frame) +
              xmd * ymd * zpd * (double)MRISseq_vox(mri, xp, yp, zm, frame) +
              xmd * ymd * zmd * (double)MRISseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_INT:
    *pval = val =
              xpd * ypd * zpd * (double)MRIIseq_vox(mri, xm, ym, zm, frame) +
              xpd * ypd * zmd * (double)MRIIseq_vox(mri, xm, ym, zp, frame) +
              xpd * ymd * zpd * (double)MRIIseq_vox(mri, xm, yp, zm, frame) +
              xpd * ymd * zmd * (double)MRIIseq_vox(mri, xm, yp, zp, frame) +
              xmd * ypd * zpd * (double)MRIIseq_vox(mri, xp, ym, zm, frame) +
              xmd * ypd * zmd * (double)MRIIseq_vox(mri, xp, ym, zp, frame) +
              xmd * ymd * zpd * (double)MRIIseq_vox(mri, xp, yp, zm, frame) +
              xmd * ymd * zmd * (double)MRIIseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  case MRI_LONG:
    *pval = val =
              xpd * ypd * zpd * (double)MRILseq_vox(mri, xm, ym, zm, frame) +
              xpd * ypd * zmd * (double)MRILseq_vox(mri, xm, ym, zp, frame) +
              xpd * ymd * zpd * (double)MRILseq_vox(mri, xm, yp, zm, frame) +
              xpd * ymd * zmd * (double)MRILseq_vox(mri, xm, yp, zp, frame) +
              xmd * ypd * zpd * (double)MRILseq_vox(mri, xp, ym, zm, frame) +
              xmd * ypd * zmd * (double)MRILseq_vox(mri, xp, ym, zp, frame) +
              xmd * ymd * zpd * (double)MRILseq_vox(mri, xp, yp, zm, frame) +
              xmd * ymd * zmd * (double)MRILseq_vox(mri, xp, yp, zp, frame) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIsampleVolumeFrame: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}



/*---------------------------------------------------------------------------
  Purpose: to return the approximate fraction of a voxel centered
  at the given point
  is labeled with a given label by the labeled volume: mri

  Input: mri  is the labeled volume. Ea voxel contains the uchar label index
  x,y,z is the floating point location of the center of a voxel
  whose labeling is to be determined. The voxel is examined to see
  how much of it is labeled with the label, ucharLabel

  Output: pval is the fraction which the given voxel location is labeled
  by ucharLabel
  returns NO_ERROR or ERROR_UNSUPPORTED if an unsupported
  (non-uchar) mri labeled volume is passed in
  AAM: 7/26/00
  --------------------------------------------------------------------------*/

#ifndef uchar
#define uchar unsigned char
#endif
int
MRIsampleLabeledVolume(MRI *mri,
                       double x, double y, double z,
                       double *pval, unsigned char ucharLabel)
{
  /* m's are the mri grid locations less than x (or y or z)
     (i.e. floor(x), p's are essentially rounding up  to the next
     grid location  greater than x  */
  int  xm, xp, ym, yp, zm, zp;
  int width, height, depth ;
  double xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  uchar ucharDmmm;
  uchar ucharDmmp;
  uchar ucharDmpm;
  uchar ucharDmpp;
  uchar ucharDpmm;
  uchar ucharDpmp;
  uchar ucharDppm;
  uchar ucharDppp;

  *pval=0.0;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  /*  if (x >= width)
      x = width - 1.0 ;
      if (y >= height)
      y = height - 1.0 ;
      if (z >= depth)
      z = depth - 1.0 ;
      if (x < 0.0)
      x = 0.0 ;
      if (y < 0.0)
      y = 0.0 ;
      if (z < 0.0)
      z = 0.0 ;
  */

  /* if the x,y,z point is out of range
     then return that none of the given voxel was labeled by ucharLabel */
  if (x >= width) return(NO_ERROR) ;
  if (y >= height) return(NO_ERROR) ;
  if (z >= depth) return(NO_ERROR) ;
  if (x < 0.0) return(NO_ERROR) ;
  if (y < 0.0) return(NO_ERROR) ;
  if (z < 0.0) return(NO_ERROR) ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  ucharDmmm = MRIvox(mri, xm, ym, zm) == ucharLabel ? 1 : 0;
  ucharDmmp = MRIvox(mri, xm, ym, zp) == ucharLabel ? 1 : 0;
  ucharDmpm = MRIvox(mri, xm, yp, zm) == ucharLabel ? 1 : 0;
  ucharDmpp = MRIvox(mri, xm, yp, zp) == ucharLabel ? 1 : 0;
  ucharDpmm = MRIvox(mri, xp, ym, zm) == ucharLabel ? 1 : 0;
  ucharDpmp = MRIvox(mri, xp, ym, zp) == ucharLabel ? 1 : 0;
  ucharDppm = MRIvox(mri, xp, yp, zm) == ucharLabel ? 1 : 0;
  ucharDppp = MRIvox(mri, xp, yp, zp) == ucharLabel ? 1 : 0;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval =
      xpd * ypd * zpd * (double)ucharDmmm +
      xpd * ypd * zmd * (double)ucharDmmp +
      xpd * ymd * zpd * (double)ucharDmpm +
      xpd * ymd * zmd * (double)ucharDmpp +
      xmd * ypd * zpd * (double)ucharDpmm +
      xmd * ypd * zmd * (double)ucharDpmp +
      xmd * ymd * zpd * (double)ucharDppm +
      xmd * ymd * zmd * (double)ucharDppp ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIsampleVolume: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIsampleVolumeType( const MRI *mri,
                     double x, double y, double z,
                     double *pval, int type)
{
  int   xv, yv, zv ;
  int OutOfBounds;

  switch (type)
  {
  default:
  case SAMPLE_NEAREST:
    break ;
  case SAMPLE_TRILINEAR:
    return(MRIsampleVolume(mri, x, y, z, pval)) ;
  case SAMPLE_CUBIC:
    return(MRIcubicSampleVolume(mri, x, y, z, pval)) ;
  case SAMPLE_SINC:
    return(MRIsincSampleVolume(mri, x, y, z, 5, pval)) ;
  }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val;
    return(NO_ERROR) ;
  }

  xv = nint(x) ;
  yv = nint(y) ;
  zv = nint(z) ;
  if (xv < 0)
    xv = 0 ;
  if (xv >= mri->width)
    xv = mri->width-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= mri->height)
    yv = mri->height-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= mri->depth)
    zv = mri->depth-1 ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = (float)MRIvox(mri, xv, yv, zv) ;
    break ;
  case MRI_SHORT:
    *pval = (float)MRISvox(mri, xv, yv, zv) ;
    break ;
  case MRI_INT:
    *pval = (float)MRIIvox(mri, xv, yv, zv) ;
    break ;
  case MRI_FLOAT:
    *pval = MRIFvox(mri, xv, yv, zv) ;
    break ;
  default:
    *pval = 0 ;
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeType: unsupported volume type %d",
      mri->type)) ;
  }
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIsampleVolumeFrameType
( const MRI *mri,
  const double x, const double y, const double z,
  const int frame, int type,
  double *pval )
{
  int   xv, yv, zv ;
  int OutOfBounds;

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    type = SAMPLE_NEAREST ;

  switch (type)
  {
  case SAMPLE_NEAREST:
    break ;
  case SAMPLE_TRILINEAR:
    return(MRIsampleVolumeFrame(mri, x, y, z, frame, pval)) ;
  case SAMPLE_CUBIC_BSPLINE:
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeFrameType(%d): First create coeff image and then use it to sample (see mriBSpline).",
      type));    
  case SAMPLE_SINC:
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeFrameType(%d): unsupported interpolation type",
      type));
  default: break;
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeFrameType(%d): unsupported interpolation type",
      type));
    /*E* add SAMPLE_CUBIC here? */
    /*    return(MRIsincSampleVolume(mri, x, y, z, 5, pval)) ;*/
  }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val;
    return(NO_ERROR) ;
  }

  xv = nint(x) ;
  yv = nint(y) ;
  zv = nint(z) ;
  if (xv < 0)
    xv = 0 ;
  if (xv >= mri->width)
    xv = mri->width-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= mri->height)
    yv = mri->height-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= mri->depth)
    zv = mri->depth-1 ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = (float)MRIseq_vox(mri, xv, yv, zv, frame) ;
    break ;
  case MRI_SHORT:
    *pval = (float)MRISseq_vox(mri, xv, yv, zv, frame) ;
    break ;
  case MRI_INT:
    *pval = (float)MRIIseq_vox(mri, xv, yv, zv, frame) ;
    break ;
  case MRI_FLOAT:
    *pval = MRIFseq_vox(mri, xv, yv, zv, frame) ;
    break ;
  default:
    *pval = 0 ;
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeFrameType: unsupported volume type %d",
      mri->type)) ;
  }
  return(NO_ERROR) ;
}

int
MRIinterpolateIntoVolume(MRI *mri, double x, double y, double z, double val)
{
  return(MRIinterpolateIntoVolumeFrame(mri, x, y, z,  0, val)) ;
}
int
MRIinterpolateIntoVolumeFrame(MRI *mri, double x, double y, double z, int frame, double val)
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  double xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;
  //  printf("MRIinterpolateIntoVolume: (xpd, ypd, zpd)%f, %f, %f\n", xpd, ypd, zpd) ; //LZ
  //  printf("MRIinterpolateIntoVolume: (xmd, ymd, zmd)%f, %f, %f\n", xmd, ymd, zmd) ; //LZ
  switch (mri->type)
  {
  case MRI_UCHAR:
    MRIseq_vox(mri, xm, ym, zm, frame) += nint(xpd * ypd * zpd * val) ;
    MRIseq_vox(mri, xm, ym, zp, frame) += nint(xpd * ypd * zmd * val) ;
    MRIseq_vox(mri, xm, yp, zm, frame) += nint(xpd * ymd * zpd * val) ;
    MRIseq_vox(mri, xm, yp, zp, frame) += nint(xpd * ymd * zmd * val) ;
    MRIseq_vox(mri, xp, ym, zm, frame) += nint(xmd * ypd * zpd * val) ;
    MRIseq_vox(mri, xp, ym, zp, frame) += nint(xmd * ypd * zmd * val) ;
    MRIseq_vox(mri, xp, yp, zm, frame) += nint(xmd * ymd * zpd * val) ;
    MRIseq_vox(mri, xp, yp, zp, frame) += nint(xmd * ymd * zmd * val) ;
    break ;
  case MRI_FLOAT:
    MRIFseq_vox(mri, xm, ym, zm, frame) += (xpd * ypd * zpd * val) ;
    MRIFseq_vox(mri, xm, ym, zp, frame) += (xpd * ypd * zmd * val) ;
    MRIFseq_vox(mri, xm, yp, zm, frame) += (xpd * ymd * zpd * val) ;
    MRIFseq_vox(mri, xm, yp, zp, frame) += (xpd * ymd * zmd * val) ;
    MRIFseq_vox(mri, xp, ym, zm, frame) += (xmd * ypd * zpd * val) ;
    MRIFseq_vox(mri, xp, ym, zp, frame) += (xmd * ypd * zmd * val) ;
    MRIFseq_vox(mri, xp, yp, zm, frame) += (xmd * ymd * zpd * val) ;
    MRIFseq_vox(mri, xp, yp, zp, frame) += (xmd * ymd * zmd * val) ;
    break ;
  case MRI_SHORT:
    MRISseq_vox(mri, xm, ym, zm, frame) += nint(xpd * ypd * zpd * val) ;
    MRISseq_vox(mri, xm, ym, zp, frame) += nint(xpd * ypd * zmd * val) ;
    MRISseq_vox(mri, xm, yp, zm, frame) += nint(xpd * ymd * zpd * val) ;
    MRISseq_vox(mri, xm, yp, zp, frame) += nint(xpd * ymd * zmd * val) ;
    MRISseq_vox(mri, xp, ym, zm, frame) += nint(xmd * ypd * zpd * val) ;
    MRISseq_vox(mri, xp, ym, zp, frame) += nint(xmd * ypd * zmd * val) ;
    MRISseq_vox(mri, xp, yp, zm, frame) += nint(xmd * ymd * zpd * val) ;
    MRISseq_vox(mri, xp, yp, zp, frame) += nint(xmd * ymd * zmd * val) ;
    break ;
  case MRI_INT:
    MRIIseq_vox(mri, xm, ym, zm, frame) += nint(xpd * ypd * zpd * val) ;
    MRIIseq_vox(mri, xm, ym, zp, frame) += nint(xpd * ypd * zmd * val) ;
    MRIIseq_vox(mri, xm, yp, zm, frame) += nint(xpd * ymd * zpd * val) ;
    MRIIseq_vox(mri, xm, yp, zp, frame) += nint(xpd * ymd * zmd * val) ;
    MRIIseq_vox(mri, xp, ym, zm, frame) += nint(xmd * ypd * zpd * val) ;
    MRIIseq_vox(mri, xp, ym, zp, frame) += nint(xmd * ypd * zmd * val) ;
    MRIIseq_vox(mri, xp, yp, zm, frame) += nint(xmd * ymd * zpd * val) ;
    MRIIseq_vox(mri, xp, yp, zp, frame) += nint(xmd * ymd * zmd * val) ;
    break ;
  case MRI_LONG:
    MRILseq_vox(mri, xm, ym, zm, frame) += nint(xpd * ypd * zpd * val) ;
    MRILseq_vox(mri, xm, ym, zp, frame) += nint(xpd * ypd * zmd * val) ;
    MRILseq_vox(mri, xm, yp, zm, frame) += nint(xpd * ymd * zpd * val) ;
    MRILseq_vox(mri, xm, yp, zp, frame) += nint(xpd * ymd * zmd * val) ;
    MRILseq_vox(mri, xp, ym, zm, frame) += nint(xmd * ypd * zpd * val) ;
    MRILseq_vox(mri, xp, ym, zp, frame) += nint(xmd * ypd * zmd * val) ;
    MRILseq_vox(mri, xp, yp, zm, frame) += nint(xmd * ymd * zpd * val) ;
    MRILseq_vox(mri, xp, yp, zp, frame) += nint(xmd * ymd * zmd * val) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIsampleVolume: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}

/*-------------------------------------------------------------------
  MRIsampleVolume() - performs trilinear interpolation on a
  single-frame volume. See MRIsampleSeqVolume() for sampling
  multi-frame.
  -------------------------------------------------------------------*/
int
MRIsampleVolume( const MRI *mri, double x, double y, double z, double *pval)
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  double val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    return(MRIsampleVolumeType(mri, x, y, z, pval, SAMPLE_NEAREST)) ;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val ;
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  switch (mri->type)
  {
  case MRI_UCHAR:
    *pval = val =
              xpd * ypd * zpd * (double)MRIvox(mri, xm, ym, zm) +
              xpd * ypd * zmd * (double)MRIvox(mri, xm, ym, zp) +
              xpd * ymd * zpd * (double)MRIvox(mri, xm, yp, zm) +
              xpd * ymd * zmd * (double)MRIvox(mri, xm, yp, zp) +
              xmd * ypd * zpd * (double)MRIvox(mri, xp, ym, zm) +
              xmd * ypd * zmd * (double)MRIvox(mri, xp, ym, zp) +
              xmd * ymd * zpd * (double)MRIvox(mri, xp, yp, zm) +
              xmd * ymd * zmd * (double)MRIvox(mri, xp, yp, zp) ;
    break ;
  case MRI_FLOAT:
    *pval = val =
              xpd * ypd * zpd * (double)MRIFvox(mri, xm, ym, zm) +
              xpd * ypd * zmd * (double)MRIFvox(mri, xm, ym, zp) +
              xpd * ymd * zpd * (double)MRIFvox(mri, xm, yp, zm) +
              xpd * ymd * zmd * (double)MRIFvox(mri, xm, yp, zp) +
              xmd * ypd * zpd * (double)MRIFvox(mri, xp, ym, zm) +
              xmd * ypd * zmd * (double)MRIFvox(mri, xp, ym, zp) +
              xmd * ymd * zpd * (double)MRIFvox(mri, xp, yp, zm) +
              xmd * ymd * zmd * (double)MRIFvox(mri, xp, yp, zp) ;
    break ;
  case MRI_SHORT:
    *pval = val =
              xpd * ypd * zpd * (double)MRISvox(mri, xm, ym, zm) +
              xpd * ypd * zmd * (double)MRISvox(mri, xm, ym, zp) +
              xpd * ymd * zpd * (double)MRISvox(mri, xm, yp, zm) +
              xpd * ymd * zmd * (double)MRISvox(mri, xm, yp, zp) +
              xmd * ypd * zpd * (double)MRISvox(mri, xp, ym, zm) +
              xmd * ypd * zmd * (double)MRISvox(mri, xp, ym, zp) +
              xmd * ymd * zpd * (double)MRISvox(mri, xp, yp, zm) +
              xmd * ymd * zmd * (double)MRISvox(mri, xp, yp, zp) ;
    break ;
  case MRI_INT:
    *pval = val =
              xpd * ypd * zpd * (double)MRIIvox(mri, xm, ym, zm) +
              xpd * ypd * zmd * (double)MRIIvox(mri, xm, ym, zp) +
              xpd * ymd * zpd * (double)MRIIvox(mri, xm, yp, zm) +
              xpd * ymd * zmd * (double)MRIIvox(mri, xm, yp, zp) +
              xmd * ypd * zpd * (double)MRIIvox(mri, xp, ym, zm) +
              xmd * ypd * zmd * (double)MRIIvox(mri, xp, ym, zp) +
              xmd * ymd * zpd * (double)MRIIvox(mri, xp, yp, zm) +
              xmd * ymd * zmd * (double)MRIIvox(mri, xp, yp, zp) ;
    break ;
  case MRI_LONG:
    *pval = val =
              xpd * ypd * zpd * (double)MRILvox(mri, xm, ym, zm) +
              xpd * ypd * zmd * (double)MRILvox(mri, xm, ym, zp) +
              xpd * ymd * zpd * (double)MRILvox(mri, xm, yp, zm) +
              xpd * ymd * zmd * (double)MRILvox(mri, xm, yp, zp) +
              xmd * ypd * zpd * (double)MRILvox(mri, xp, ym, zm) +
              xmd * ypd * zmd * (double)MRILvox(mri, xp, ym, zp) +
              xmd * ymd * zpd * (double)MRILvox(mri, xp, yp, zm) +
              xmd * ymd * zmd * (double)MRILvox(mri, xp, yp, zp) ;
    break ;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "MRIsampleVolume: unsupported type %d", mri->type)) ;
    break ;
  }
  return(NO_ERROR) ;
}

/*------------------------------------------------------------------
  MRIsampleSeqVolume() - performs trilinear interpolation on a
  multi-frame volume. valvect is a vector of length nframes. No error
  checking for first and last frame. The caller is able to specify
  first frame and last frame so that all frames do not have to be
  sampled at the same time (this can be important in time-sensitive
  applications).
  -------------------------------------------------------------------*/
int MRIsampleSeqVolume( const MRI *mri, double x, double y, double z, float *valvect,
                        int firstframe, int lastframe )
{
  int  OutOfBounds;
  int  f,xm, xp, ym, yp, zm, zp, width, height, depth ;
  double xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    for (f=firstframe; f <= lastframe; f++) valvect[f] = mri->outside_val;
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  for (f = firstframe; f <= lastframe; f++)
  {
    switch (mri->type)
    {
    case MRI_UCHAR:
      valvect[f] =
        xpd * ypd * zpd * (double)MRIseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRIseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRIseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRIseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRIseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRIseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRIseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRIseq_vox(mri, xp, yp, zp, f) ;
      break ;
    case MRI_FLOAT:
      valvect[f] =
        xpd * ypd * zpd * (double)MRIFseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRIFseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRIFseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRIFseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRIFseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRIFseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRIFseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRIFseq_vox(mri, xp, yp, zp, f) ;
      break ;
    case MRI_SHORT:
      valvect[f] =
        xpd * ypd * zpd * (double)MRISseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRISseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRISseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRISseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRISseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRISseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRISseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRISseq_vox(mri, xp, yp, zp, f) ;
      break ;
    case MRI_INT:
      valvect[f] =
        xpd * ypd * zpd * (double)MRIIseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRIIseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRIIseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRIIseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRIIseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRIIseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRIIseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRIIseq_vox(mri, xp, yp, zp, f) ;
      break ;
    case MRI_LONG:
      valvect[f] =
        xpd * ypd * zpd * (double)MRILseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRILseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRILseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRILseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRILseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRILseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRILseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRILseq_vox(mri, xp, yp, zp, f) ;
      break ;
    default:
      ErrorReturn(ERROR_UNSUPPORTED,
                  (ERROR_UNSUPPORTED,
                   "MRIsampleSeqVolume: unsupported type %d", mri->type)) ;
      break ;
    }
  }/* end loop over frames */

  return(NO_ERROR) ;
}

//testing - LZ
int MRIsampleSeqVolumeType(MRI *mri, double x, double y, double z, float *valvect,
			   int firstframe, int lastframe, int type)
{
  int  OutOfBounds;
  int  /*f,xm, xp, ym, yp, zm, zp, */ f, width, height, depth ;
  //double xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  int   xv, yv, zv ;

  switch (type)
    {
    default:
    case SAMPLE_NEAREST:
      break ;
    case SAMPLE_TRILINEAR:
      return(MRIsampleSeqVolume(mri, x, y, z, valvect, firstframe, lastframe)) ;
    case SAMPLE_CUBIC:
      {
	printf("Cubic interpolation is not implemented yet on multi-frame data. Going ahead with tri-linear") ;
	return(MRIsampleSeqVolume(mri, x, y, z, valvect, firstframe, lastframe)) ;
      }
    case SAMPLE_SINC:
      {
	printf("Sinc interpolation is not implemented yet on multi-frame data. Going ahead with tri-linear") ;
	return(MRIsampleSeqVolume(mri, x, y, z, valvect, firstframe, lastframe)) ;
      }
    case SAMPLE_CUBIC_BSPLINE:
      {
	printf("Cubic Bspline interpolation not implemented here. Needs to be done by calling function!") ;
	exit(1);
      }
    }

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    for (f=firstframe; f <= lastframe; f++) 
      valvect[f] = mri->outside_val;
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  //
  xv = nint(x) ;
  yv = nint(y) ;
  zv = nint(z) ;
  if (xv < 0)
    xv = 0 ;
  if (xv >= width)
    xv = width-1 ;
  if (yv < 0)
    yv = 0 ;
  if (yv >= height)
    yv = height-1 ;
  if (zv < 0)
    zv = 0 ;
  if (zv >= depth)
    zv = depth-1 ;
  //

  //if (x >= width)    x = width - 1.0 ;
  //if (y >= height)   y = height - 1.0 ;
  //if (z >= depth)    z = depth - 1.0 ;
  //if (x < 0.0)       x = 0.0 ;
  //if (y < 0.0)       y = 0.0 ;
  //if (z < 0.0)       z = 0.0 ;

  /*  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;*/

  for (f = firstframe; f <= lastframe; f++)
  {
    switch (mri->type)
    {
    case MRI_UCHAR:
      valvect[f] = (double)MRIseq_vox(mri, xv, yv, zv, f) ;
      /*        xpd * ypd * zpd * (double)MRIseq_vox(mri, xm, ym, zm, f) +
		xpd * ypd * zmd * (double)MRIseq_vox(mri, xm, ym, zp, f) +
		xpd * ymd * zpd * (double)MRIseq_vox(mri, xm, yp, zm, f) +
		xpd * ymd * zmd * (double)MRIseq_vox(mri, xm, yp, zp, f) +
		xmd * ypd * zpd * (double)MRIseq_vox(mri, xp, ym, zm, f) +
		xmd * ypd * zmd * (double)MRIseq_vox(mri, xp, ym, zp, f) +
		xmd * ymd * zpd * (double)MRIseq_vox(mri, xp, yp, zm, f) +
		xmd * ymd * zmd * (double)MRIseq_vox(mri, xp, yp, zp, f) ;*/
      break ;
    case MRI_FLOAT:
      valvect[f] = (double)MRIFseq_vox(mri, xv, yv, zv, f);
      /*xpd * ypd * zpd * (double)MRIFseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRIFseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRIFseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRIFseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRIFseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRIFseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRIFseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRIFseq_vox(mri, xp, yp, zp, f) ;*/
      break ;
    case MRI_SHORT:
      valvect[f] = (double)MRISseq_vox(mri, xv, yv, zv, f);
      /*xpd * ypd * zpd * (double)MRISseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRISseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRISseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRISseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRISseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRISseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRISseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRISseq_vox(mri, xp, yp, zp, f) ;*/
      break ;
    case MRI_INT:
      valvect[f] = (double)MRIIseq_vox(mri, xv, yv, zv, f);
      /*xpd * ypd * zpd * (double)MRIIseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRIIseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRIIseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRIIseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRIIseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRIIseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRIIseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRIIseq_vox(mri, xp, yp, zp, f) ;*/
      break ;
    case MRI_LONG:
      valvect[f] = (double)MRILseq_vox(mri, xv, yv, zv, f);
      /*xpd * ypd * zpd * (double)MRILseq_vox(mri, xm, ym, zm, f) +
        xpd * ypd * zmd * (double)MRILseq_vox(mri, xm, ym, zp, f) +
        xpd * ymd * zpd * (double)MRILseq_vox(mri, xm, yp, zm, f) +
        xpd * ymd * zmd * (double)MRILseq_vox(mri, xm, yp, zp, f) +
        xmd * ypd * zpd * (double)MRILseq_vox(mri, xp, ym, zm, f) +
        xmd * ypd * zmd * (double)MRILseq_vox(mri, xp, ym, zp, f) +
        xmd * ymd * zpd * (double)MRILseq_vox(mri, xp, yp, zm, f) +
        xmd * ymd * zmd * (double)MRILseq_vox(mri, xp, yp, zp, f) ;*/
      break ;
    default:
      valvect[f] = 0;
      ErrorReturn(ERROR_UNSUPPORTED,
                  (ERROR_UNSUPPORTED,
                   "MRIsampleSeqVolumeType: unsupported type %d", mri->type)) ;
      break ;
    }
  }/* end loop over frames */

  return(NO_ERROR) ;
}


/*---------------------------------------------------------------------*/
/*!
  \fn double *MRItrilinKernel(MRI *mri, double c, double r, double s, double *kernel)
  \brief Computes the kernel used for trilinear interpolation
  \param mri - used to get volume dimensions
  \param c - column (continuous)
  \param r - row (continuous)
  \param s - slice (continuous)
  \param kernel - array of length 8. Can be NULL.
  \return kernel - array of length 8.
  \description Computes the kernel used for trilinear interpolation. See
  also MRIsampleSeqVolume().
 */
double *MRItrilinKernel(MRI *mri, double c, double r, double s, double *kernel)
{
  int  OutOfBounds;
  int  f,xm, xp, ym, yp, zm, zp, width, height, depth ;
  double xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */

  if (kernel == NULL) kernel = (double *) calloc(8,sizeof(double));

  OutOfBounds = MRIindexNotInVolume(mri, c, r, s);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    for (f=0; f < 8; f++) kernel[f] = 0;
    return(kernel) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (c >= width)   c = width - 1.0 ;
  if (r >= height)  r = height - 1.0 ;
  if (s >= depth)   s = depth - 1.0 ;
  if (c < 0.0)      c = 0.0 ;
  if (r < 0.0)      r = 0.0 ;
  if (s < 0.0)      s = 0.0 ;

  xm = MAX((int)c, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)r, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)s, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = c - (double)xm ;
  ymd = r - (double)ym ;
  zmd = s - (double)zm ;
  xpd = (1.0 - xmd) ;
  ypd = (1.0 - ymd) ;
  zpd = (1.0 - zmd) ;

  kernel[0] = xpd * ypd * zpd; // xm ym zm
  kernel[1] = xpd * ypd * zmd; // xm ym zp
  kernel[2] = xpd * ymd * zpd; // xm yp zm
  kernel[3] = xpd * ymd * zmd; // xm yp zp
  kernel[4] = xmd * ypd * zpd; // xp ym zm
  kernel[5] = xmd * ypd * zmd; // xp ym zp
  kernel[6] = xmd * ymd * zpd; // xp yp zm
  kernel[7] = xmd * ymd * zmd; // xp yp zp

  return(kernel) ;
}

/*-----------------------------------------------------
  used by MRIcubicSampleVolume
  ------------------------------------------------------*/

double
localeval(double x, int iter)
{
  double p;
  switch (iter)
  {
  case 0:
    p = ((2-x)*x-1)*x;
    break;
  case 1:
    p = (3*x-5)*x*x+2;
    break;
  case 2:
    p = ((4-3*x)*x+1)*x;
    break;
  case 3:
    p = (x-1)*x*x;
    break;
  default:
    ErrorReturn(ERROR_UNSUPPORTED,
                (ERROR_UNSUPPORTED,
                 "localeval: called wrong by MRIcubicSampleVolume!")) ;
  }
  return(p);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  by analogy with
  /usr/pubsw/common/matlab/6.5/toolbox/matlab/polyfun/interp3.m

  uses localeval above

  ------------------------------------------------------*/
int
MRIcubicSampleVolume( const MRI *mri,
                      double x, double y, double z, double *pval )
{
  int  OutOfBounds;
  int  width, height, depth ;
  int ix_low,iy_low,iz_low,ix,iy,iz;
  double val,xx,yy,zz,fx,fy,fz,vv[4][4][4];

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    return(MRIsampleVolumeType(mri, x, y, z, pval, SAMPLE_NEAREST)) ;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = val = mri->outside_val;
    return(NO_ERROR) ;
  }
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  /*E* I suppose these are for "ambiguously out of bounds" - within .5vox */

  /*E* I think this needs an edit - x is double, whatever that is, not
    int, so any x>= width-1 should be set to width-1.
    if (x >= width)    x = width - 1.0 ;
    if (y >= height)   y = height - 1.0 ;
    if (z >= depth)    z = depth - 1.0 ;
  */

  if (x > width-1.0)    x = width - 1.0 ;
  if (y > height-1.0)   y = height - 1.0 ;
  if (z > depth-1.0)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  ix_low = floor((double)x);
  if ((ix_low = floor((double)x)) < width-1)
    xx = x - ix_low;
  else
  {
    ix_low--;
    xx = 1;
  }
  iy_low = floor((double)y);
  if ((iy_low = floor((double)y)) < height-1)
    yy = y - iy_low;
  else
  {
    iy_low--;
    yy = 1;
  }
  iz_low = floor((double)z);
  if ((iz_low = floor((double)z)) < depth-1)
    zz = z - iz_low;
  else
  {
    iz_low--;
    zz = 1;
  }

  /*E* build a little box of the local points plus boundary stuff -
    for this rev accept zeroes for the border expansion */

  for (iz= MAX(0,1-iz_low); iz<MIN(4,depth+1-iz_low); iz++)
  {
    for (iy= MAX(0,1-iy_low); iy<MIN(4,height+1-iy_low); iy++)
    {
      for (ix= MAX(0,1-ix_low); ix<MIN(4,width+1-ix_low); ix++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          vv[ix][iy][iz] =
            (double)MRIvox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz);
          break;
        case MRI_FLOAT:
          vv[ix][iy][iz] =
            (double)MRIFvox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz);
          break;
        case MRI_SHORT:
          vv[ix][iy][iz] =
            (double)MRISvox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz);
          break;
        case MRI_INT:
          vv[ix][iy][iz] =
            (double)MRIIvox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz);
          break;
        case MRI_LONG:
          vv[ix][iy][iz] =
            (double)MRILvox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz);
          break;
        default:
          ErrorReturn(ERROR_UNSUPPORTED,
                      (ERROR_UNSUPPORTED,
                       "MRIcubicSampleVolume: unsupported type %d",
                       mri->type)) ;
          break ;
        }
      }
    }
  }

  val = 0;

  for (iz=0; iz<=3; iz++)
  {
    fz = localeval(zz,iz);
    for (iy=0; iy<=3; iy++)
    {
      fy = localeval(yy,iy);
      for (ix=0; ix<=3; ix++)
      {
        fx = localeval(xx,ix);
        val += (double)(vv[ix][iy][iz]*fx*fy*fz);
      }
    }
  }

  *pval = val/8.;

  return(NO_ERROR);
}

int
MRIcubicSampleVolumeFrame(MRI *mri, double x, double y, double z, int frame, double *pval)
{
  int  OutOfBounds;
  int  width, height, depth ;
  int ix_low,iy_low,iz_low,ix,iy,iz;
  double val,xx,yy,zz,fx,fy,fz,vv[4][4][4];

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    return(MRIsampleVolumeFrameType(mri, x, y, z, frame, SAMPLE_NEAREST, pval)) ;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = val = mri->outside_val;
    return(NO_ERROR) ;
  }
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  /*E* I suppose these are for "ambiguously out of bounds" - within .5vox */

  /*E* I think this needs an edit - x is double, whatever that is, not
    int, so any x>= width-1 should be set to width-1.
    if (x >= width)    x = width - 1.0 ;
    if (y >= height)   y = height - 1.0 ;
    if (z >= depth)    z = depth - 1.0 ;
  */

  if (x > width-1.0)    x = width - 1.0 ;
  if (y > height-1.0)   y = height - 1.0 ;
  if (z > depth-1.0)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  ix_low = floor((double)x);
  if ((ix_low = floor((double)x)) < width-1)
    xx = x - ix_low;
  else
  {
    ix_low--;
    xx = 1;
  }
  iy_low = floor((double)y);
  if ((iy_low = floor((double)y)) < height-1)
    yy = y - iy_low;
  else
  {
    iy_low--;
    yy = 1;
  }
  iz_low = floor((double)z);
  if ((iz_low = floor((double)z)) < depth-1)
    zz = z - iz_low;
  else
  {
    iz_low--;
    zz = 1;
  }

  /*E* build a little box of the local points plus boundary stuff -
    for this rev accept zeroes for the border expansion */

  for (iz= MAX(0,1-iz_low); iz<MIN(4,depth+1-iz_low); iz++)
  {
    for (iy= MAX(0,1-iy_low); iy<MIN(4,height+1-iy_low); iy++)
    {
      for (ix= MAX(0,1-ix_low); ix<MIN(4,width+1-ix_low); ix++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          vv[ix][iy][iz] =
            (double)MRIseq_vox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz,frame);
          break;
        case MRI_FLOAT:
          vv[ix][iy][iz] =
            (double)MRIFseq_vox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz,frame);
          break;
        case MRI_SHORT:
          vv[ix][iy][iz] =
            (double)MRISseq_vox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz,frame);
          break;
        case MRI_INT:
          vv[ix][iy][iz] =
            (double)MRIIseq_vox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz,frame);
          break;
        case MRI_LONG:
          vv[ix][iy][iz] =
            (double)MRILseq_vox(mri,ix_low-1+ix,iy_low-1+iy,iz_low-1+iz,frame);
          break;
        default:
          ErrorReturn(ERROR_UNSUPPORTED,
                      (ERROR_UNSUPPORTED,
                       "MRIcubicSampleVolumeFrame: unsupported type %d",
                       mri->type)) ;
          break ;
        }
      }
    }
  }

  val = 0;

  for (iz=0; iz<=3; iz++)
  {
    fz = localeval(zz,iz);
    for (iy=0; iy<=3; iy++)
    {
      fy = localeval(yy,iy);
      for (ix=0; ix<=3; ix++)
      {
        fx = localeval(xx,ix);
        val += (double)(vv[ix][iy][iz]*fx*fy*fz);
      }
    }
  }

  *pval = val/8.;

  return(NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define IMIN(a,b) (a < b ? a : b)
#define IMAX(a,b) (a > b ? a : b)
double ham_sinc(double x,double fullwidth)
{
  double ham;
  if ( fabs(x) < 1.0e-5)
    ham = 1.0;
  else
  {
    ham = sin(PI*x)/(PI*x);
    ham *= 0.54 + 0.46 * cos(2.0*PI*x/fullwidth);
  }
  return ham;
}

/*-------------------------------------------------------------------------*/
int
MRIsincSampleVolume( const MRI *mri,
                     double x, double y, double z,
                     int hw, double *pval )
{
  int  OutOfBounds;
  int  width, height, depth ;
  int nwidth;
  int ix_low,ix_high,iy_low,iy_high,iz_low,iz_high;
  int jx1,jy1,jz1,jx_rel,jy_rel,jz_rel;
  double coeff_x[128],coeff_y[128],coeff_z[128];
  double coeff_x_sum,coeff_y_sum,coeff_z_sum;
  double sum_x,sum_y,sum_z;
  double xsize,ysize,zsize;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val ;
    return(NO_ERROR) ;
  }

  xsize = mri->xsize;
  ysize=mri->ysize;
  zsize=mri->zsize;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  nwidth = hw;
  ix_low = floor((double)x);
  ix_high = ceil((double)x);
  iy_low = floor((double)y);
  iy_high = ceil((double)y);
  iz_low = floor((double)z);
  iz_high = ceil((double)z);

  coeff_x_sum = coeff_y_sum = coeff_z_sum = 0;
  if (iz_low>=0 && iz_high < depth)
  {
    for (jx1=IMAX(ix_high-nwidth,0), jx_rel=0;
         jx1<IMIN(ix_low+nwidth,width-1);
         jx1++,jx_rel++)
    {
      coeff_x[jx_rel] = ham_sinc((double)(x-jx1),2*nwidth);
      coeff_x_sum += coeff_x[jx_rel];
    }
    for (jy1=IMAX(iy_high-nwidth,0), jy_rel=0;
         jy1<IMIN(iy_low+nwidth,height-1);
         jy1++,jy_rel++)
    {
      coeff_y[jy_rel] = ham_sinc((double)(y-jy1),2*nwidth);
      coeff_y_sum += coeff_y[jy_rel];
    }
    for (jz1=IMAX(iz_high-nwidth,0), jz_rel=0;
         jz1<IMIN(iz_low+nwidth,depth-1);
         jz1++,jz_rel++)
    {
      coeff_z[jz_rel] = ham_sinc((double)(z-jz1),2*nwidth);
      coeff_z_sum += coeff_z[jz_rel];
    }

    for (sum_z=0., jz1=IMAX(iz_high-nwidth,0), jz_rel = 0;
         jz1 < IMIN(iz_low+nwidth,depth-1);
         jz1++, jz_rel++)
    {

      for (sum_y=0., jy1=IMAX(iy_high-nwidth,0), jy_rel = 0;
           jy1 < IMIN(iy_low+nwidth,height-1);
           jy1++, jy_rel++)
      {
        for (sum_x=0., jx1=IMAX(ix_high-nwidth,0), jx_rel = 0;
             jx1 < IMIN(ix_low+nwidth,width-1);
             jx1++, jx_rel++)
        {

          switch (mri->type)
          {
          case MRI_UCHAR:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
                     * (double)MRIvox(mri,jx1,jy1,jz1);
            break;
          case MRI_SHORT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
                     * (double)MRISvox(mri,jx1,jy1,jz1);
            break;
          case MRI_INT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
                     * (double)MRIIvox(mri,jx1,jy1,jz1);
            break;
          case MRI_LONG:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
                     * (double)MRILvox(mri,jx1,jy1,jz1);
            break;
          case MRI_FLOAT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
                     * (double)MRIFvox(mri,jx1,jy1,jz1);
            break;
          default:
            ErrorReturn(ERROR_UNSUPPORTED,
                        (ERROR_UNSUPPORTED,
                         "MRIsincSampleVolume: unsupported type %d",
                         mri->type)) ;
            break;
          }
        }
        sum_y += sum_x * (coeff_y[jy_rel]/coeff_y_sum);
      }
      sum_z += sum_y * (coeff_z[jz_rel]/coeff_z_sum);
    }
    if ((mri->type == MRI_UCHAR || mri->type == MRI_SHORT) && sum_z<0.0)
      *pval = 0.0;
    else if (mri->type == MRI_UCHAR && sum_z >255.0)
      *pval = 255.0;
    else if (mri->type == MRI_SHORT && sum_z > 65535.0)
      *pval = 65535.0;
    else
      *pval = sum_z;
  }
  else
    *pval = 0.0;

  return(NO_ERROR);
}

int
MRIsincSampleVolumeFrame(MRI *mri,
		       double x, double y, double z, int frame, 
		       int hw, double *pval)
{
  int  OutOfBounds;
  int  width, height, depth ;
  int nwidth;
  int ix_low,ix_high,iy_low,iy_high,iz_low,iz_high;
  int jx1,jy1,jz1,jx_rel,jy_rel,jz_rel;
  double coeff_x[128],coeff_y[128],coeff_z[128];
  double coeff_x_sum,coeff_y_sum,coeff_z_sum;
  double sum_x,sum_y,sum_z;
  double xsize,ysize,zsize;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val ;
    return(NO_ERROR) ;
  }

  xsize = mri->xsize;
  ysize=mri->ysize;
  zsize=mri->zsize;
  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  nwidth = hw;
  ix_low = floor((double)x);
  ix_high = ceil((double)x);
  iy_low = floor((double)y);
  iy_high = ceil((double)y);
  iz_low = floor((double)z);
  iz_high = ceil((double)z);

  coeff_x_sum = coeff_y_sum = coeff_z_sum = 0;
  if (iz_low>=0 && iz_high < depth)
  {
    for (jx1=IMAX(ix_high-nwidth,0), jx_rel=0;
         jx1<IMIN(ix_low+nwidth,width-1);
         jx1++,jx_rel++)
    {
      coeff_x[jx_rel] = ham_sinc((double)(x-jx1),2*nwidth);
      coeff_x_sum += coeff_x[jx_rel];
    }
    for (jy1=IMAX(iy_high-nwidth,0), jy_rel=0;
         jy1<IMIN(iy_low+nwidth,height-1);
         jy1++,jy_rel++)
    {
      coeff_y[jy_rel] = ham_sinc((double)(y-jy1),2*nwidth);
      coeff_y_sum += coeff_y[jy_rel];
    }
    for (jz1=IMAX(iz_high-nwidth,0), jz_rel=0;
         jz1<IMIN(iz_low+nwidth,depth-1);
         jz1++,jz_rel++)
    {
      coeff_z[jz_rel] = ham_sinc((double)(z-jz1),2*nwidth);
      coeff_z_sum += coeff_z[jz_rel];
    }

    for (sum_z=0., jz1=IMAX(iz_high-nwidth,0), jz_rel = 0;
         jz1 < IMIN(iz_low+nwidth,depth-1);
         jz1++, jz_rel++)
    {

      for (sum_y=0., jy1=IMAX(iy_high-nwidth,0), jy_rel = 0;
           jy1 < IMIN(iy_low+nwidth,height-1);
           jy1++, jy_rel++)
      {
        for (sum_x=0., jx1=IMAX(ix_high-nwidth,0), jx_rel = 0;
             jx1 < IMIN(ix_low+nwidth,width-1);
             jx1++, jx_rel++)
        {

          switch (mri->type)
          {
          case MRI_UCHAR:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
	      * (double)MRIseq_vox(mri,jx1,jy1,jz1,frame);
            break;
          case MRI_SHORT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
	      * (double)MRISseq_vox(mri,jx1,jy1,jz1,frame);
            break;
          case MRI_INT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
	      * (double)MRIIseq_vox(mri,jx1,jy1,jz1,frame);
            break;
          case MRI_LONG:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
	      * (double)MRILseq_vox(mri,jx1,jy1,jz1,frame);
            break;
          case MRI_FLOAT:
            sum_x += (coeff_x[jx_rel]/coeff_x_sum)
	      * (double)MRIFseq_vox(mri,jx1,jy1,jz1,frame);
            break;
          default:
            ErrorReturn(ERROR_UNSUPPORTED,
                        (ERROR_UNSUPPORTED,
                         "MRIsincSampleVolumeFrame: unsupported type %d",
                         mri->type)) ;
            break;
          }
        }
        sum_y += sum_x * (coeff_y[jy_rel]/coeff_y_sum);
      }
      sum_z += sum_y * (coeff_z[jz_rel]/coeff_z_sum);
    }
    if ((mri->type == MRI_UCHAR || mri->type == MRI_SHORT) && sum_z<0.0)
      *pval = 0.0;
    else if (mri->type == MRI_UCHAR && sum_z >255.0)
      *pval = 255.0;
    else if (mri->type == MRI_SHORT && sum_z > 65535.0)
      *pval = 65535.0;
    else
      *pval = sum_z;
  }
  else
    *pval = 0.0;

  return(NO_ERROR);
}
/*-----------------------------------------------------------------
  MRIindexNotInVolume() - determines whether a col, row, slice point is
  in the mri volume. If it is unambiguously in the volume, then 0
  is returned. If it is within 0.5 of the edge of the volume, -1
  is returned. Otherwise 1 is returned. Flagging the case where
  the point is within 0.5 of the edge can be used for assigning
  a nearest neighbor when the point is outside but close to the
  volume. In this case the index of the nearest neighbor can safely
  be computed as the nearest integer to col, row, and slice.
  -----------------------------------------------------------------*/
int MRIindexNotInVolume( const MRI *mri,
			 const double col, const double row, const double slice )
{
  float nicol, nirow, nislice;

  /* unambiguously in the volume */
  if (col   >= 0 && col   <= mri->width-1  &&
      row   >= 0 && row   <= mri->height-1 &&
      slice >= 0 && slice <= mri->depth-1  )
    return(0);

  /* within 0.5 of the edge of the volume */
  nicol   = rint(col);
  nirow   = rint(row);
  nislice = rint(slice);
  if (nicol   >= 0 && nicol   < mri->width  &&
      nirow   >= 0 && nirow   < mri->height &&
      nislice >= 0 && nislice < mri->depth  )
    return(-1);

  /* unambiguously NOT in the volume */
  return(1);

}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume directional derivative using
  trilinear interpolation.
  ------------------------------------------------------*/
float
MRIsampleCardinalDerivative(MRI *mri, int x, int y, int z,
                            int xk, int yk, int zk)
{
  float d ;

  if (xk)
    d = MRIsampleXDerivative(mri, x, y, z, xk) ;
  else if (yk)
    d = MRIsampleYDerivative(mri, x, y, z, yk) ;
  else
    d = MRIsampleZDerivative(mri, x, y, z, zk) ;
  return(d) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume directional derivative using
  trilinear interpolation.
  ------------------------------------------------------*/
float
MRIsampleXDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dx ;
  int   yk, zk, xi, yi, zi, nvox ;

  dx = 0.0 ;

  xi = mri->xi[x+dir] ;
  for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dx += dir*MRIgetVoxVal(mri, xi, yi, zi,0) ;  /* x+dir */
      dx -= dir*MRIgetVoxVal(mri, x, yi, zi, 0) ;   /* - x */
      nvox += 2 ;
    }
  }
  dx /= ((float)nvox*mri->xsize) ;
  return(dx) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume directional derivative using
  trilinear interpolation.
  ------------------------------------------------------*/
float
MRIsampleYDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dy ;
  int   xk, zk, xi, yi, zi, nvox ;

  dy = 0.0 ;

  yi = mri->yi[y+dir] ;
  for (nvox = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (xk = -1 ; xk <= 1 ; xk++)
    {
      xi = mri->xi[x+xk] ;
      dy += dir*MRIgetVoxVal(mri, xi, yi, zi, 0) ;  /* x+dir */
      dy -= dir*MRIgetVoxVal(mri, x, yi, zi, 0) ;   /* - x */
      nvox += 2 ;
    }
  }
  dy /= ((float)nvox*mri->ysize) ;
  return(dy) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume directional derivative using
  trilinear interpolation.
  ------------------------------------------------------*/
float
MRIsampleZDerivative(MRI *mri, int x, int y, int z, int dir)
{
  float dz ;
  int   xk, yk, xi, yi, zi, nvox ;

  dz = 0.0 ;

  zi = mri->zi[z+dir] ;
  for (nvox = 0, xk = -1 ; xk <= 1 ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dz += dir*MRIgetVoxVal(mri, xi, yi, zi, 0) ;  /* x+dir */
      dz -= dir*MRIgetVoxVal(mri, x, yi, zi, 0) ;   /* - x */
      nvox += 2 ;
    }
  }
  dz /= ((float)nvox*mri->zsize) ;
  return(dz) ;
}
int
MRIsampleVolumeDirectionScale(MRI *mri, double x, double y, double z,
                              double dx, double dy, double dz, double *pmag,
                              double sigma)
{
  int  width, height, depth ;
  double xp1, xm1, yp1, ym1, zp1, zm1, len ;
  double dist, val, k, ktotal, step_size, total_val ;
  int  n ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  step_size = MAX(.5,sigma/2) ;
  for (total_val = ktotal = 0.0,n = 0, len = 0.0, dist = step_size ;
       dist <= MAX(2*sigma,step_size);
       dist += step_size, n++)
  {
    if (FZERO(sigma))
      k = 1.0 ;
    else
      k = exp(-dist*dist/(2*sigma*sigma)) ;
    ktotal += k ;
    len += dist ;
    xp1 = x + dist*dx ;
    yp1 = y + dist*dy ;
    zp1 = z + dist*dz ;
    MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
    total_val += k*val ;

    xm1 = x - dist*dx ;
    ym1 = y - dist*dy ;
    zm1 = z - dist*dz ;
    MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
    total_val += k*val ;
    if (FZERO(step_size))
      break ;
  }
  total_val /= (double)2.0*ktotal ;

  *pmag = total_val ;
  return(NO_ERROR) ;
}

int
MRIsampleVolumeDerivativeScale(MRI *mri,
                               double x, double y, double z,
                               double dx, double dy, double dz,
                               double *pmag, double sigma)
{
  int  width, height, depth ;
  double xp1, xm1, yp1, ym1, zp1, zm1, vp1, vm1, len ;
  double dist, val, k, ktotal, step_size ;
  int  n ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;

  step_size = MAX(.25,sigma/5.0) ;
  for (ktotal = 0.0,n = 0, len = vp1 = vm1 = 0.0, dist = step_size ;
       dist <= MAX(2*sigma,step_size);
       dist += step_size, n++)
  {
    if (FZERO(sigma))
      k = 1.0 ;
    else
      k = exp(-dist*dist/(2*sigma*sigma)) ;
    ktotal += k ;
    len += dist ;
    xp1 = x + dist*dx ;
    yp1 = y + dist*dy ;
    zp1 = z + dist*dz ;
    MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
    vp1 += k*val ;

    xm1 = x - dist*dx ;
    ym1 = y - dist*dy ;
    zm1 = z - dist*dz ;
    MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
    vm1 += k*val ;
    if (FZERO(step_size))
      break ;
  }
  vm1 /= (double)ktotal ;
  vp1 /= (double)ktotal ;
  len /= (double)ktotal ;

  *pmag = (vp1-vm1) / (2.0*len) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume directional derivative using
  trilinear interpolation.
  ------------------------------------------------------*/
int
MRIsampleVolumeDerivative(MRI *mri, double x, double y, double z,
                          double dx, double dy, double dz, double *pmag)
{
  int  width, height, depth ;
  double xp1, xm1, yp1, ym1, zp1, zm1, vp1, vm1, len ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
#if 1
  {
    double dist, val ;
    int  n ;

    for (n = 0, len = vp1 = vm1 = 0.0, dist = .5 ; dist <= 2 ;
         dist += 0.5, n++)
    {
      len += dist ;
      xp1 = x + dist*dx ;
      yp1 = y + dist*dy ;
      zp1 = z + dist*dz ;
      MRIsampleVolume(mri, xp1, yp1, zp1, &val) ;
      vp1 += val ;

      xm1 = x - dist*dx ;
      ym1 = y - dist*dy ;
      zm1 = z - dist*dz ;
      MRIsampleVolume(mri, xm1, ym1, zm1, &val) ;
      vm1 += val ;
    }
    vm1 /= (double)n ;
    vp1 /= (double)n ;
    len /= (double)n ;
  }
#else
  xp1 = x + dx ;
  xm1 = x - dx ;
  yp1 = y + dy ;
  ym1 = y - dy ;
  zp1 = z + dz ;
  zm1 = z - dz ;
  len = sqrt(dx*dx+dy*dy+dz*dz) ;
  MRIsampleVolume(mri, xp1, yp1, zp1, &vp1) ;
  MRIsampleVolume(mri, xm1, ym1, zm1, &vm1) ;
#endif

  *pmag = (vp1-vm1) / (2.0*len) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume gradient to cubic voxels.
  ------------------------------------------------------*/
int
MRIsampleVolumeGradient(MRI *mri, double x, double y, double z,
                        double *pdx, double *pdy, double *pdz)
{
  int  width, height, depth ;
  double xp1, xm1, yp1, ym1, zp1, zm1 ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
#if 0
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
#endif
  if (MRIindexNotInVolume(mri, x+1, y, z))
    MRIsampleVolume(mri, x, y, z, &xp1) ;
  else
    MRIsampleVolume(mri, x+1.0, y, z, &xp1) ;
  if (MRIindexNotInVolume(mri, x-1, y, z))
    MRIsampleVolume(mri, x, y, z, &xm1) ;
  else
    MRIsampleVolume(mri, x-1.0, y, z, &xm1) ;

  if (MRIindexNotInVolume(mri, x, y+1, z))
    MRIsampleVolume(mri, x, y, z, &yp1) ;
  else
    MRIsampleVolume(mri, x, y+1.0, z, &yp1) ;
  if (MRIindexNotInVolume(mri, x, y-1.0, z))
    MRIsampleVolume(mri, x, y, z, &ym1) ;
  else
    MRIsampleVolume(mri, x, y-1.0, z, &ym1) ;

  if (MRIindexNotInVolume(mri, x, y, z+1.0))
    MRIsampleVolume(mri, x, y, z, &zp1) ;
  else
    MRIsampleVolume(mri, x, y, z+1.0, &zp1) ;
  if (MRIindexNotInVolume(mri, x, y, z-1.0))
    MRIsampleVolume(mri, x, y, z, &zm1) ;
  else
    MRIsampleVolume(mri, x, y, z-1.0, &zm1) ;


  *pdx = (xp1-xm1)/(2.0*mri->xsize) ;
  *pdy = (yp1-ym1)/(2.0*mri->ysize) ;
  *pdz = (zp1-zm1)/(2.0*mri->zsize) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Interpolate the volume gradient to cubic voxels.
  ------------------------------------------------------*/
int
MRIsampleVolumeGradientFrame( const MRI *mri,
			      double x, double y, double z,
			      double *pdx, double *pdy, double *pdz,
			      int frame )
{
  int  width, height, depth ;
  double xp1, xm1, yp1, ym1, zp1, zm1 ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
#if 0
  if (x >= width)
    x = width - 1.0 ;
  if (y >= height)
    y = height - 1.0 ;
  if (z >= depth)
    z = depth - 1.0 ;
  if (x < 0.0)
    x = 0.0 ;
  if (y < 0.0)
    y = 0.0 ;
  if (z < 0.0)
    z = 0.0 ;
  if (frame >= mri->nframes)
    frame = mri->nframes-1 ;
  if (frame < 0)
    frame = 0 ;
#endif

  MRIsampleVolumeFrame(mri, x+1.0, y, z, frame, &xp1) ;
  MRIsampleVolumeFrame(mri, x-1.0, y, z, frame, &xm1) ;

  MRIsampleVolumeFrame(mri, x, y+1.0, z, frame, &yp1) ;
  MRIsampleVolumeFrame(mri, x, y-1.0, z, frame, &ym1) ;

  MRIsampleVolumeFrame(mri, x, y, z+1.0, frame, &zp1) ;
  MRIsampleVolumeFrame(mri, x, y, z-1.0, frame, &zm1) ;

  *pdx = (xp1-xm1)/(2.0*mri->xsize) ;
  *pdy = (yp1-ym1)/(2.0*mri->ysize) ;
  *pdz = (zp1-zm1)/(2.0*mri->zsize) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighborsOn(MRI *mri, int x0, int y0, int z0, int min_val)
{
  int   nbrs = 0 ;

  if (MRIgetVoxVal(mri,mri->xi[x0-1],y0,z0, 0) >= min_val)
    nbrs++ ;
  if (MRIgetVoxVal(mri,mri->xi[x0+1],y0,z0, 0) >= min_val)
    nbrs++ ;
  if (MRIgetVoxVal(mri,x0,mri->yi[y0+1],z0, 0) >= min_val)
    nbrs++ ;
  if (MRIgetVoxVal(mri,x0,mri->yi[y0-1],z0, 0) >= min_val)
    nbrs++ ;
  if (MRIgetVoxVal(mri,x0,y0,mri->zi[z0+1], 0) >= min_val)
    nbrs++ ;
  if (MRIgetVoxVal(mri,x0,y0,mri->zi[z0-1], 0) >= min_val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighborsOn3x3(MRI *mri, int x, int y, int z, int min_val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs ;

  for (nbrs = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIgetVoxVal(mri, xi, yi, zi, 0) > min_val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighborsInWindow(MRI *mri, int x, int y, int z, int wsize, int val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs, whalf ;

  whalf = (wsize-1)/2 ;

  for (nbrs = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIgetVoxVal(mri, xi, yi, zi, 0) == val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighbors(MRI *mri, int x0, int y0, int z0, int val)
{
  int   nbrs = 0 ;

  if (nint(MRIgetVoxVal(mri,mri->xi[x0-1],y0,z0,0)) == val)
    nbrs++ ;
  if (nint(MRIgetVoxVal(mri,mri->xi[x0+1],y0,z0,0)) == val)
    nbrs++ ;
  if (nint(MRIgetVoxVal(mri,x0,mri->yi[y0+1],z0,0)) == val)
    nbrs++ ;
  if (nint(MRIgetVoxVal(mri,x0,mri->yi[y0-1],z0,0)) == val)
    nbrs++ ;
  if (nint(MRIgetVoxVal(mri,x0,y0,mri->zi[z0+1],0)) == val)
    nbrs++ ;
  if (nint(MRIgetVoxVal(mri,x0,y0,mri->zi[z0-1],0)) == val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighborsInRange(MRI *mri, 
                    int x0, int y0, int z0,
                    int frame,
                    float low_val,
                    float hi_val)
{
  int   nbrs = 0 ;
  float val ;

  val = MRIgetVoxVal(mri,mri->xi[x0-1],y0,z0,frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  val = MRIgetVoxVal(mri,mri->xi[x0+1],y0,z0,frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  val = MRIgetVoxVal(mri,x0,mri->yi[y0+1],z0,frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  val = MRIgetVoxVal(mri,x0,mri->yi[y0-1],z0,frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  val = MRIgetVoxVal(mri,x0,y0,mri->zi[z0+1],frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  val = MRIgetVoxVal(mri,x0,y0,mri->zi[z0-1],frame) ;
  if (val >= low_val && val <= hi_val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int
MRIneighbors3x3(MRI *mri, int x, int y, int z, int val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs ;

  for (nbrs = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIvox(mri, xi, yi, zi) == val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIneighborsOff(MRI *mri, int x0, int y0, int z0, int min_val)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,x0-1,y0,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0+1,y0,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0+1,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0-1,z0) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0+1) < min_val)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0-1) < min_val)
    nbrs++ ;
  return(nbrs) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
int
MRIneighborsOff3x3(MRI *mri, int x, int y, int z, int min_val)
{
  int   xk, yk, zk, xi, yi, zi, nbrs ;

  for (nbrs = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (!zk && !yk && !xk)
          continue ;
        if (MRIvox(mri, xi, yi, zi) < min_val)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
/*-----------------------------------------------------
  Perform an linear coordinate transformation x' = Ax on
  the MRI image mri_src into mri_dst
  ------------------------------------------------------*/
MRI *
MRIinverseLinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA)
{
  MATRIX  *m_inv ;

  m_inv = MatrixInverse(mA, NULL) ;
  if   (!m_inv)
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MRIinverseLinearTransform: xform is singular!")) ;
  //fprintf
  //(stderr,
 //  "applying the vox-to-vox linear transform (calculated inverse)\n");
 // MatrixPrint(stderr, m_inv);
  mri_dst = MRIlinearTransform(mri_src, mri_dst, m_inv) ;
  MatrixFree(&m_inv) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Convert a transform from RAS to voxel coordinates, then apply
  it to an MRI.
  ------------------------------------------------------*/
MRI *
MRIapplyRASlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *m_ras_xform)
{
  MATRIX   *m_voxel_xform ;

  m_voxel_xform = MRIrasXformToVoxelXform(mri_src, mri_dst, m_ras_xform, NULL);
  //fprintf(stderr, "applying the vox to vox linear transform\n");
  //MatrixPrint(stderr, m_voxel_xform);
  mri_dst = MRIlinearTransform(mri_src, mri_dst, m_voxel_xform) ;
  MatrixFree(&m_voxel_xform) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Convert a transform from RAS to voxel coordinates, then apply
  it to an MRI.
  ------------------------------------------------------*/
MRI *
MRIapplyRASlinearTransformInterp
(MRI *mri_src, MRI *mri_dst, MATRIX *m_ras_xform, int interp)
{
  MATRIX   *m_voxel_xform ;

  m_voxel_xform = MRIrasXformToVoxelXform(mri_src, mri_dst, m_ras_xform, NULL);
  //fprintf(stderr, "applying the vox to vox linear transform\n");
  //MatrixPrint(stderr, m_voxel_xform);
  mri_dst = MRIlinearTransformInterp(mri_src, mri_dst, m_voxel_xform, interp) ;
  MatrixFree(&m_voxel_xform) ;
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Convert a transform from RAS to voxel coordinates, then apply
  it to an MRI.
  ------------------------------------------------------*/
MRI *
MRIapplyRASinverseLinearTransform(MRI *mri_src, MRI *mri_dst,
                                  MATRIX *m_ras_xform)
{
  MATRIX   *m_voxel_xform ;

  m_voxel_xform = MRIrasXformToVoxelXform(mri_src, mri_dst, m_ras_xform, NULL);
  mri_dst = MRIinverseLinearTransform(mri_src, mri_dst, m_voxel_xform) ;
  MatrixFree(&m_voxel_xform) ;
  return(mri_dst) ;
}

/*-----------------------------------------------------
  Perform an linear coordinate transformation x' = Ax on
  the MRI image mri_src into mri_dst using sinc interp.
  ------------------------------------------------------*/
MRI *
MRIsincTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA, int hw)
{
  int    y1, y2, y3, width, height, depth ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  MATRIX *mAinv ;     /* inverse of mA */
  double   val, x1, x2, x3 ;

  mAinv = MatrixInverse(mA, NULL) ;      /* will sample from dst back to src */
  if (!mAinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIsincTransform: xform is singular")) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;

        x1 = V3_X(v_X) ;
        x2 = V3_Y(v_X) ;
        x3 = V3_Z(v_X) ;

        if (nint(y1) == 13 && nint(y2) == 10 && nint(y3) == 7)
          DiagBreak() ;
        if (nint(x1) == 13 && nint(x2) == 10 && nint(x3) == 7)
        {
#if 0
          fprintf
          (stderr,
           "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
           (float)x1, (float)x2, (float)x3,
           (float)y1, (float)y2, (float)y3) ;
#endif
          DiagBreak() ;
        }

        if (x1 > -1 && x1 < width &&
            x2 > -1 && x2 < height &&
            x3 > -1 && x3 < depth)
        {
          MRIsincSampleVolume(mri_src, x1, x2, x3, hw, &val);
          MRIvox(mri_dst,y1,y2,y3) = (BUFTYPE)nint(val) ;
        }
      }
    }
  }

  MatrixFree(&v_X) ;
  MatrixFree(&mAinv) ;
  MatrixFree(&v_Y) ;

  mri_dst->ras_good_flag = 0;

  return(mri_dst) ;
}
/*-----------------------------------------------------------------
  MRIlinearTransform() - for historical reasons, this uses trilinear
  interpolation. This the operations under this function name can
  now (2/20/02) be found under MRIlinearTransformInterp().
  -----------------------------------------------------------------*/
MRI *
MRIlinearTransform(MRI *mri_src, MRI *mri_dst, MATRIX *mA)
{
  mri_dst = MRIlinearTransformInterp(mri_src, mri_dst, mA, SAMPLE_TRILINEAR);
  return(mri_dst);
}
/*-------------------------------------------------------------------
  MRIlinearTransformInterp() Perform linear coordinate transformation
  x' = Ax on the MRI image mri_src into mri_dst using the specified
  interpolation method. A is a voxel-to-voxel transform.
  ------------------------------------------------------------------*/
MRI *
MRIlinearTransformInterp(MRI *mri_src, MRI *mri_dst, MATRIX *mA,
                         int InterpMethod)
{
  int    y1, y2, y3, width, height, depth,frame ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  MATRIX *mAinv ;     /* inverse of mA */
  double   val, x1, x2, x3 ;

  if (InterpMethod != SAMPLE_NEAREST &&
      InterpMethod != SAMPLE_TRILINEAR &&
      InterpMethod != SAMPLE_CUBIC_BSPLINE )
  {
    printf("ERROR: MRIlinearTransformInterp: unrecognized or unsupported interpolation "
           "method %d\n",InterpMethod);
  }

  mAinv = MatrixInverse(mA, NULL) ;      /* will sample from dst back to src */
  if (!mAinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIlinearTransform: xform is singular")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;
  else
    MRIclear(mri_dst) ;

  MRI_BSPLINE * bspline = NULL;
  if (InterpMethod == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(mri_src,NULL,3);

  width  = mri_dst->width ;
  height = mri_dst->height ;
  depth  = mri_dst->depth ;
  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;

        x1 = V3_X(v_X) ;
        x2 = V3_Y(v_X) ;
        x3 = V3_Z(v_X) ;

        if (nint(y1) == Gx && nint(y2) == Gy && nint(y3) == Gz)
          DiagBreak() ;
        if (nint(x1) == Gx && nint(x2) == Gy && nint(x3) == Gz)
        {
#if 0
          fprintf
          (stderr,
           "(%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
           (float)x1, (float)x2, (float)x3,
           (float)y1, (float)y2, (float)y3) ;
#endif
          DiagBreak() ;
        }

        //MRIsampleVolume(mri_src, x1, x2, x3, &val);
        for (frame = 0 ; frame < mri_src->nframes ; frame++)
        {
          if (InterpMethod == SAMPLE_CUBIC_BSPLINE)
            // recommended to externally call this and keep mri_coeff
            // if image is resampled often (e.g. in registration algo)
            MRIsampleBSpline(bspline, x1, x2, x3, frame, &val);
          else
            MRIsampleVolumeFrameType(mri_src, x1, x2, x3, 
                                     frame, InterpMethod, &val);

          // will clip the val according to mri_dst type:
          MRIsetVoxVal(mri_dst, y1, y2, y3, frame, val) ;

#if 0
          // if this gets ever enabled, don't forget to clip val to type...
          switch (mri_dst->type)
          {
          case MRI_UCHAR:
            MRIvox(mri_dst,y1,y2,y3) = (BUFTYPE)nint(val) ;
          break ;
          case MRI_SHORT:
            MRISvox(mri_dst,y1,y2,y3) = (short)nint(val) ;
            break ;
          case MRI_FLOAT:
            MRIFvox(mri_dst,y1,y2,y3) = (float)(val) ;
            break ;
          case MRI_INT:
            MRIIvox(mri_dst,y1,y2,y3) = nint(val) ;
            break ;
          default:
            ErrorReturn(NULL,
                        (ERROR_UNSUPPORTED,
                         "MRIlinearTransform: unsupported dst type %d",
                         mri_dst->type)) ;
            break ;
          }
#endif
        }
      }
    }
  }
  if (bspline) MRIfreeBSpline(&bspline);
  MatrixFree(&v_X) ;
  MatrixFree(&mAinv) ;
  MatrixFree(&v_Y) ;

  mri_dst->ras_good_flag = 1;

  return(mri_dst) ;
}
MRI *
MRIconcatenateFrames(MRI *mri_frame1, MRI *mri_frame2, MRI *mri_dst)
{
  int       width, height, depth, x, y, z ;
  BUFTYPE   *pf1, *pf2, *pdst1, *pdst2 ;

  if (mri_frame1->type != MRI_UCHAR || mri_frame1->type != MRI_UCHAR)
    ErrorReturn(NULL,
                (ERROR_UNSUPPORTED,"MRIconcatenateFrames: src must be UCHAR"));

  width = mri_frame1->width ;
  height = mri_frame1->height ;
  depth = mri_frame1->depth ;
  if (mri_dst == NULL)
  {
    mri_dst = MRIallocSequence(width, height, depth, mri_frame1->type, 2) ;
    MRIcopyHeader(mri_frame1, mri_dst) ;
  }
  if (!mri_dst)
    ErrorExit(ERROR_NOMEMORY, "MRIconcatenateFrames: could not alloc dst") ;


  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst1 = &MRIvox(mri_dst, 0, y, z) ;
      pdst2 = &MRIseq_vox(mri_dst, 0, y, z, 1) ;
      pf1 = &MRIvox(mri_frame1, 0, y, z) ;
      pf2 = &MRIvox(mri_frame2, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        *pdst1++ = *pf1++ ;
        *pdst2++ = *pf2++ ;
      }
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  ------------------------------------------------------*/
MRI *
MRIcopyFrame(MRI *mri_src, MRI *mri_dst, int src_frame, int dst_frame)
{
  int       width, height, depth, y, z, bytes ;
  BUFTYPE   *psrc, *pdst ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  if (!mri_dst)
    mri_dst =
      MRIallocSequence(width, height, depth, mri_src->type, dst_frame+1) ;
  if (!mri_dst)
    ErrorExit(ERROR_NOMEMORY, "MRIcopyFrame: could not alloc dst") ;

  if (mri_src->type != mri_dst->type)
    ErrorReturn(NULL,(ERROR_UNSUPPORTED,
                      "MRIcopyFrame: src and dst must be same type"));
  if (dst_frame >= mri_dst->nframes)
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MRIcopyFrame: dst frame #%d out of range (nframes=%d)\n",
      dst_frame, mri_dst->nframes)) ;
  MRIcopyHeader(mri_src, mri_dst) ;

  switch (mri_src->type)
  {
  case MRI_UCHAR:
    bytes = sizeof(unsigned char) ;
    break ;
  case MRI_SHORT:
    bytes = sizeof(short) ;
    break ;
  case MRI_INT:
    bytes = sizeof(int) ;
    break ;
  case MRI_FLOAT:
    bytes = sizeof(float) ;
    break ;
  default:
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIcopyFrame: unsupported src format %d",
                       mri_src->type));
    break ;
  }
  bytes *= width ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      switch (mri_src->type)
      {
      default:  /* already handled above */
      case MRI_UCHAR:
        psrc = &MRIseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = &MRIseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_SHORT:
        psrc = (BUFTYPE *)&MRISseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRISseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_INT:
        psrc = (BUFTYPE *)&MRIIseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRIIseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      case MRI_FLOAT:
        psrc = (BUFTYPE *)&MRIFseq_vox(mri_src, 0, y, z, src_frame) ;
        pdst = (BUFTYPE *)&MRIFseq_vox(mri_dst, 0, y, z, dst_frame) ;
        break ;
      }
      memmove(pdst, psrc, bytes) ;
    }
  }
  return(mri_dst) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
   compute the mean in a frame of all values
  ------------------------------------------------------*/
double
MRImeanFrame(MRI *mri, int frame)
{
  int       width, height, depth, x, y, z ;
  double    mean ;
  double      val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (mri->nframes <= frame)
    ErrorReturn(0.0,(ERROR_BADPARM,
                     "MRImeanFrame: frame %d out of bounds (%d)",
                     frame, mri->nframes));

  if (frame >= 0)
  {
    for (mean = 0.0, z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          val = MRIgetVoxVal(mri, x, y, z, frame) ;
          mean += (double)val ;
        }
      }
    }
    mean /= (double)(width*height*depth) ;
  }
  else
  {
    for (mean = 0.0, frame = 0 ; frame < mri->nframes ; frame++)
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            val = MRIgetVoxVal(mri, x, y, z, frame) ;
            mean += (double)val ;
          }
        }
      }
    mean /= (double)(width*height*depth*mri->nframes) ;
  }

  return(mean) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
   compute the mean in a frame of all values
  ------------------------------------------------------*/
double
MRImeanFrameNonzeroMask(MRI *mri, int frame, MRI *mri_mask)
{
  int       width, height, depth, x, y, z, nvox ;
  double    mean ;
  double      val ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  if (mri->nframes <= frame)
    ErrorReturn(0.0,(ERROR_BADPARM,
                     "MRImeanFrameNonzeroMask: frame %d out of bounds (%d)",
                     frame, mri->nframes));

  nvox = 0 ;
  if (frame >= 0)
  {
    for (mean = 0.0, z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (MRIgetVoxVal(mri_mask, x, y, z, frame) == 0)
            continue ;
          val = MRIgetVoxVal(mri, x, y, z, frame) ;
          mean += (double)val ;
          nvox++ ;
        }
      }
    }
    mean /= nvox ;
  }
  else
  {
    for (mean = 0.0, frame = 0 ; frame < mri->nframes ; frame++)
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            if (MRIgetVoxVal(mri_mask, x, y, z, frame) == 0)
              continue ;
            val = MRIgetVoxVal(mri, x, y, z, frame) ;
            mean += (double)val ;
          }
        }
      }
    mean /= nvox ;
  }

  return(mean) ;
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
   compute the mean in a frame of all values above thresh
  ------------------------------------------------------*/
double
MRImeanFrameThresh(MRI *mri, int frame, float thresh)
{
  int       width, height, depth, x, y, z, num ;
  double    mean ;
  double      val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;
  if (mri->nframes <= frame)
    ErrorReturn(0.0,(ERROR_BADPARM,
                     "MRImeanFrame: frame %d out of bounds (%d)",
                     frame, mri->nframes));

  for (mean = 0.0, num = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        MRIsampleVolumeType(mri, x, y, z, &val, SAMPLE_NEAREST) ;
        if (val < thresh)
          continue ;
        mean += (double)val ;
        num++ ;
      }
    }
  }
  mean /= (double)(num) ;
  return(mean) ;
}


/*--------------------------------------------------------------
  MRISeqchangeType() - changes the data type for a 3D or 4D volume.
  This simply changes the volume dimensions so that it appears to be a
  3D volume, then calls MRIchangeType(), and then resets the
  dimensions to their original values. The values of the volume can be
  rescaled between f_low and f_high.
  ------------------------------------------------------------*/
MRI *MRISeqchangeType(MRI *vol, int dest_type, float f_low,
                      float f_high, int no_scale_option_flag)
{
  int nslices, nframes, i;
  MRI *mri;

  /* Change vol dimensions to make it look like a single frame */
  // This can cause problems with MRI_FRAME operations, see below
  nslices = vol->depth;
  nframes = vol->nframes;
  vol->depth = nslices*nframes;
  vol->nframes = 1;

  /* Change the type */
  mri = MRIchangeType(vol,dest_type,f_low,f_high,no_scale_option_flag);

  /* Change vol dimensions back to original */
  vol->depth = nslices;
  vol->nframes = nframes;

  /* Check for error */
  if (mri == NULL)
  {
    fprintf(stderr,"ERROR: MRISeqchangeType: MRIchangeType\n");
    return(NULL);
  }

  /* Change mri dimensions back to original */
  mri->depth = nslices;
  mri->nframes = nframes;

  // Alloc MRI_FRAME. This needs to be updated when MRI_FRAME items are added
  mri->frames = (MRI_FRAME *)calloc(mri->nframes, sizeof(MRI_FRAME)) ;
  for (i = 0 ; i < mri->nframes ; i++)
    mri->frames[i].m_ras2vox = MatrixAlloc(4,4, MATRIX_REAL) ;

  return(mri);
}
/*-----------------------------------------------------------
  MRIchangeType() - changes the data type of a 3D MRI volume,
  with optional rescaling. Use MRISeqchangeType() for 3D or
  4D volumes.
  ---------------------------------------------------------*/
MRI *MRIchangeType(MRI *src, int dest_type, float f_low,
                   float f_high, int no_scale_option_flag)
{

  MRI *dest = NULL;
  int i, j, k;
  float val;
  int no_scale_flag = FALSE;
  float scale, dest_min, dest_max; /* new = scale * (val - min) */
  float src_min=0.0, src_max=0.0;
  int hist_bins[N_HIST_BINS];
  float bin_size;
  int bin;
  int nth, n_passed;

  /* ----- shut the compiler up ----- */
  val = 0.0;
  dest_min = dest_max = 0.0;

  if (src->type == dest_type)
  {
    dest = MRIcopy(src, NULL);
    return(dest);
  }

  if (src->type == MRI_UCHAR &&
      (dest_type == MRI_SHORT ||
       dest_type == MRI_INT ||
       dest_type == MRI_LONG ||
       dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if (src->type == MRI_SHORT &&
           (dest_type == MRI_INT ||
            dest_type == MRI_LONG ||
            dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if (src->type == MRI_LONG &&
           (dest_type == MRI_INT ||
            dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else if (src->type == MRI_INT &&
           (dest_type == MRI_LONG ||
            dest_type == MRI_FLOAT))
    no_scale_flag = TRUE;
  else
  {
    MRIlimits(src, &src_min, &src_max);

    if (no_scale_option_flag > 1)
      no_scale_flag = TRUE;
    else if (no_scale_option_flag == 1)
    {
      if (dest_type == MRI_UCHAR &&
          src_min >= UCHAR_MIN &&
          src_max <= UCHAR_MAX)
        no_scale_flag = TRUE;
      if (dest_type == MRI_SHORT &&
          src_min >= SHORT_MIN &&
          src_max <= SHORT_MAX)
        no_scale_flag = TRUE;
      if (dest_type == MRI_INT &&
          src_min >= INT_MIN &&
          src_max <= INT_MAX)
        no_scale_flag = TRUE;
      if (dest_type == MRI_LONG &&
          src_min >= LONG_MIN &&
          src_max <= LONG_MAX)
        no_scale_flag = TRUE;
      if (no_scale_flag == FALSE)
        printf("******* WARNING - forcing scaling of values "
               "to prevent cropping of input [%2.0f, %2.0f]\n",
               src_min, src_max) ;
    }
  }

  if (no_scale_flag)
  {
    dest = MRIalloc(src->width, src->height, src->depth, dest_type);
    MRIcopyHeader(src, dest);
    dest->type = dest_type;

    for (k = 0;k < src->depth;k++)
      for (j = 0;j < src->height;j++)
        for (i = 0;i < src->width;i++)
        {

          if (src->type == MRI_UCHAR)
            val = (float)MRIvox(src, i, j, k);
          if (src->type == MRI_SHORT)
            val = (float)MRISvox(src, i, j, k);
          if (src->type == MRI_INT)
            val = (float)MRIIvox(src, i, j, k);
          if (src->type == MRI_LONG)
            val = (float)MRILvox(src, i, j, k);
          if (src->type == MRI_FLOAT)
            val = (float)MRIFvox(src, i, j, k);

          if (dest->type == MRI_UCHAR)
          {
            if (val < UCHAR_MIN)
              val = UCHAR_MIN;
            if (val > UCHAR_MAX)
              val = UCHAR_MAX;
            MRIvox(dest, i, j, k) = (unsigned char)nint(val);
          }
          if (dest->type == MRI_SHORT)
          {
            if (val < SHORT_MIN)
              val = SHORT_MIN;
            if (val > SHORT_MAX)
              val = SHORT_MAX;
            MRISvox(dest, i, j, k) = (short)nint(val);
          }
          if (dest->type == MRI_INT)
          {
            if (val < INT_MIN)
              val = INT_MIN;
            if (val > INT_MAX)
              val = INT_MAX;
            MRIIvox(dest, i, j, k) = (int)nint(val);
          }
          if (dest->type == MRI_LONG)
          {
            if (val < LONG_MIN)
              val = LONG_MIN;
            if (val > LONG_MAX)
              val = LONG_MAX;
            MRILvox(dest, i, j, k) = (long)nint(val);
          }
          
          if (dest_type == MRI_FLOAT)
            MRIFvox(dest, i, j, k) = (float)val;
        }
  }
  else
  {
    long nonzero = 0 ;

    /* ----- build a histogram ----- */
    printf("MRIchangeType: Building histogram \n");
    bin_size = (src_max - src_min) / (float)N_HIST_BINS;

    for (i = 0;i < N_HIST_BINS;i++)
      hist_bins[i] = 0;

    for (i = 0;i < src->width;i++)
      for (j = 0;j < src->height;j++)
        for (k = 0;k < src->depth;k++)
        {

          if (src->type == MRI_UCHAR)
            val = (float)MRIvox(src, i, j, k);
          if (src->type == MRI_SHORT)
            val = (float)MRISvox(src, i, j, k);
          if (src->type == MRI_INT)
            val = (float)MRIIvox(src, i, j, k);
          if (src->type == MRI_LONG)
            val = (float)MRILvox(src, i, j, k);
          if (src->type == MRI_FLOAT)
            val = (float)MRIFvox(src, i, j, k);

          if (!DZERO(val))
            nonzero++ ;
          bin = (int)((val - src_min) / bin_size);

          if (bin < 0)
            bin = 0;
          if (bin >= N_HIST_BINS)
            bin = N_HIST_BINS-1;

          hist_bins[bin]++;

        }

    nth = (int)(f_low * src->width * src->height * src->depth);
    for (n_passed = 0,bin = 0;n_passed < nth && bin < N_HIST_BINS;bin++)
      n_passed += hist_bins[bin];
    src_min = (float)bin * bin_size + src_min;

#if 1  // handle mostly empty volumes
    nth = (int)((1.0-f_high) * nonzero);
#else
    nth = (int)((1.0-f_high) * src->width * src->height * src->depth);
#endif
    for (n_passed = 0,bin = N_HIST_BINS-1;n_passed < nth && bin > 0;bin--)
      n_passed += hist_bins[bin];
    src_max = (float)bin * bin_size + src_min;

    //if (src_min >= src_max)
    if (src_min > src_max)
    {
      ErrorReturn
      (NULL,
       (ERROR_BADPARM,
        "MRIchangeType(): after hist: src_min = %g, "
        "src_max = %g (f_low = %g, f_high = %g)",
        src_min, src_max, f_low, f_high));
    }

    /* ----- continue ----- */

    if (dest_type == MRI_UCHAR)
    {
      dest_min = UCHAR_MIN;
      dest_max = UCHAR_MAX;
    }
    if (dest_type == MRI_SHORT)
    {
      dest_min = SHORT_MIN;
      dest_max = SHORT_MAX;
    }
    if (dest_type == MRI_INT)
    {
      dest_min = INT_MIN;
      dest_max = INT_MAX;
    }
    if (dest_type == MRI_LONG)
    {
      dest_min = LONG_MIN;
      dest_max = LONG_MAX;
    }

    if (src_max == src_min)
      scale = 1.0;
    else
      scale = (dest_max - dest_min) / (src_max - src_min);

    dest = MRIalloc(src->width, src->height, src->depth, dest_type);
    MRIcopyHeader(src, dest);
    dest->type = dest_type;

    for (i = 0;i < src->width;i++)
      for (j = 0;j < src->height;j++)
        for (k = 0;k < src->depth;k++)
        {

          if (src->type == MRI_SHORT)
            val = MRISvox(src, i, j, k);
          if (src->type == MRI_INT)
            val = MRIIvox(src, i, j, k);
          if (src->type == MRI_LONG)
            val = MRILvox(src, i, j, k);
          if (src->type == MRI_FLOAT)
            val = MRIFvox(src, i, j, k);

          val = dest_min + scale * (val - src_min);

          if (dest->type == MRI_UCHAR)
          {
            if (val < UCHAR_MIN)
              val = UCHAR_MIN;
            if (val > UCHAR_MAX)
              val = UCHAR_MAX;
            MRIvox(dest, i, j, k) = (unsigned char)nint(val);
          }
          if (dest->type == MRI_SHORT)
          {
            if (val < SHORT_MIN)
              val = SHORT_MIN;
            if (val > SHORT_MAX)
              val = SHORT_MAX;
            MRISvox(dest, i, j, k) = (short)nint(val);
          }
          if (dest->type == MRI_INT)
          {
            if (val < INT_MIN)
              val = INT_MIN;
            if (val > INT_MAX)
              val = INT_MAX;
            MRIIvox(dest, i, j, k) = (int)nint(val);
          }
          if (dest->type == MRI_LONG)
          {
            if (val < LONG_MIN)
              val = LONG_MIN;
            if (val > LONG_MAX)
              val = LONG_MAX;
            MRILvox(dest, i, j, k) = (long)nint(val);
          }

        }

  }

  return(dest);

} /* end MRIchangeType() */

/*-----------------------------------------------------*/
MATRIX *MRIgetResampleMatrix(MRI *src, MRI *template_vol)
{

  MATRIX *src_mat, *dest_mat; /* from i to ras */
  float src_det, dest_det;
  MATRIX *src_inv, *m;

  /* ----- fake the ras values if ras_good_flag is not set ----- */
  int src_slice_direction = getSliceDirection(src);
  if (!src->ras_good_flag
      && (src_slice_direction != MRI_CORONAL)
      && (src_slice_direction != MRI_SAGITTAL)
      && (src_slice_direction != MRI_HORIZONTAL))
  {
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MRIresample(): source volume orientation is unknown"));
  }

  /*

  : solve each row of src_mat * [midx;midy;midz;1] = centerr, centera, centers
  and for dest

  S = M * s_v

  where M = [(x_r y_r z_r)(xsize  0     0  )  s14]
  s_v = (center_x)   S = (c_r)   etc.
  [(x_a y_a z_a)(  0  ysize   0  )  s24]        (center_y)       (c_a)
  [(x_s y_s z_s)(  0    0   zsize)  s34]        (center_z)       (c_s)
  [       0        0       0         1 ]        (    1   )       ( 1 )

  Write M = [  m    s],  s_v = (c_v), S = (c_R), then  c_R = m * c_v + s or
  [ 0 0 0 1]         ( 1 )      ( 1 )

  The translation s is given by   s = c_R - m*c_v

  Note the convention c_(r,a,s) being defined in
  terms of c_v = (width/2., height/2, depth/2.),
  not ((width-1)/2, (height-1)/2, (depth-1)/2).

  */

  src_mat = extract_i_to_r(src); // error when allocation fails
  if (src_mat == NULL)
    return NULL; // did ErrorPrintf in extract_i_to_r()

  dest_mat = extract_i_to_r(template_vol); // error when allocation fails
  if (dest_mat == NULL)
  {
    MatrixFree(&src_mat);
    return NULL; // did ErrorPrintf in extract_i_to_r()
  }

  // error check
  src_det = MatrixDeterminant(src_mat);
  dest_det = MatrixDeterminant(dest_mat);

  if (src_det == 0.0)
  {
    errno = 0;
    ErrorPrintf
    (ERROR_BADPARM,
     "MRIresample(): source matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  if (dest_det == 0.0)
  {
    errno = 0;
    ErrorPrintf
    (ERROR_BADPARM,
     "MRIresample(): destination matrix has zero determinant; matrix is:");
    MatrixPrint(stderr, dest_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  src_inv = MatrixInverse(src_mat, NULL);

  if (src_inv == NULL)
  {
    errno = 0;
    ErrorPrintf
    (ERROR_BADPARM,
     "MRIresample(): error inverting matrix; determinant is "
     "%g, matrix is:", src_det);
    MatrixPrint(stderr, src_mat);
    MatrixFree(&src_mat);
    MatrixFree(&dest_mat);
    return(NULL);
  }

  m = MatrixMultiply(src_inv, dest_mat, NULL);
  if (m == NULL)
    return(NULL);

  MatrixFree(&src_inv);
  MatrixFree(&src_mat);
  MatrixFree(&dest_mat);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("MRIresample() matrix is:\n");
    MatrixPrint(stdout, m);
  }
  return(m) ;

} /* end MRIreslice() */


MRI *
MRIresample(MRI *src, MRI *template_vol, int resample_type)
{
  return(MRIresampleFill(src, template_vol, resample_type, 0)) ;
} /* end MRIresample() */


MRI *MRIresampleFill
(MRI *src, MRI *template_vol, int resample_type, float fill_val)
{

  MRI *dest = NULL;
  MATRIX *m;
  int nframe;
  int di, dj, dk;
  int si, sj, sk;
  float si_f, sj_f, sk_f;
  float si_ff, sj_ff, sk_ff;
  MATRIX *sp, *dp;
  float val, val000, val001, val010, val011, val100, val101, val110, val111;
  float w000, w001, w010, w011, w100, w101, w110, w111;
  float si_f2, sj_f2, sk_f2;
  float ii2, ij2, ik2, isi_f2, isj_f2, isk_f2;
  float w[8];
  int wi[8];
  int mwi;
  double pval;
  int i_good_flag, i1_good_flag;
  int j_good_flag, j1_good_flag;
  int k_good_flag, k1_good_flag;

  /* ----- keep the compiler quiet ----- */
  val = 0.0;
  val000 = val001 = val010 = val011 = val100 = val101 = val110 = val111 = 0.0;

#if 0
  if (src->type != template_vol->type)
  {
    ErrorReturn
    (NULL,
     (ERROR_UNSUPPORTED,
      "MRIresample(): source and destination types must be identical"));
  }
#endif

  /* ----- fake the ras values if ras_good_flag is not set ----- */
  if (!src->ras_good_flag)
    printf
    ("MRIresample(): WARNING: ras_good_flag "
     "is not set, changing orientation\n"
     "to default.\n");

  // get dst voxel -> src voxel transform
  m = MRIgetResampleMatrix(src, template_vol);
  if (m == NULL)
    return(NULL);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("MRIresample() matrix is:\n");
    MatrixPrint(stdout, m);
  }

  dest = MRIallocSequence(template_vol->width,
                          template_vol->height,
                          template_vol->depth,
                          src->type,
                          template_vol->nframes);
  if (dest == NULL)
    return(NULL);
  MRIreplaceValues(dest, dest, 0.0f, fill_val) ;

  MRIcopyHeader(template_vol, dest);
  MRIcopyPulseParameters(src, dest) ;

  sp = MatrixAlloc(4, 1, MATRIX_REAL);
  dp = MatrixAlloc(4, 1, MATRIX_REAL);

  *MATRIX_RELT(dp, 4, 1) = 1.0;

  MRI_BSPLINE * bspline = NULL;
  if (resample_type == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(src,NULL,3);
  

  for (nframe = 0; nframe < template_vol->nframes; nframe++)
  {
    for (di = 0;di < template_vol->width;di++)
    {
      for (dj = 0;dj < template_vol->height;dj++)
      {
        for (dk = 0;dk < template_vol->depth;dk++)
        {
          if (di == Gx && dj == Gy && dk == Gz)
            DiagBreak() ;
          *MATRIX_RELT(dp, 1, 1) = (float)di;
          *MATRIX_RELT(dp, 2, 1) = (float)dj;
          *MATRIX_RELT(dp, 3, 1) = (float)dk;

          MatrixMultiply(m, dp, sp);

          si_ff = *MATRIX_RELT(sp, 1, 1);
          sj_ff = *MATRIX_RELT(sp, 2, 1);
          sk_ff = *MATRIX_RELT(sp, 3, 1);

          si = (int)floor(si_ff);
          sj = (int)floor(sj_ff);
          sk = (int)floor(sk_ff);

          if (si == 147 && sj == 91 && sk == 86)
            DiagBreak() ;
          if (di == 129 && dj == 164 && dk == 147)
            DiagBreak() ;
          si_f = si_ff - si;
          sj_f = sj_ff - sj;
          sk_f = sk_ff - sk;

#if 0
          if (di % 20 == 0 && dj == di && dk == dj)
          {
            printf("MRIresample() sample point: %3d %3d %3d: "
                   "%f (%3d + %f) %f (%3d + %f) %f (%3d + %f)\n",
                   di, dj, dk,
                   si_ff, si, si_f,
                   sj_ff, sj, sj_f,
                   sk_ff, sk, sk_f);
          }
#endif

          if (resample_type == SAMPLE_SINC)
          {
            MRIsincSampleVolumeFrame(src, si_ff, sj_ff, sk_ff, nframe, 5, &pval);
            val = (float)pval;
          }
          else if (resample_type == SAMPLE_CUBIC)
          {
            MRIcubicSampleVolumeFrame(src, si_ff, sj_ff, sk_ff, nframe, &pval);
            val = (float)pval;
          }
          else if (resample_type == SAMPLE_CUBIC_BSPLINE)
          {
            MRIsampleBSpline(bspline, si_ff, sj_ff, sk_ff, nframe, &pval);
            val = (float)pval;
          }
          else
          {
            i_good_flag = (si >= 0 && si < src->width);
            i1_good_flag = (si+1 >= 0 && si+1 < src->width);
            j_good_flag = (sj >= 0 && sj < src->height);
            j1_good_flag = (sj+1 >= 0 && sj+1 < src->height);
            k_good_flag = (sk >= 0 && sk < src->depth);
            k1_good_flag = (sk+1 >= 0 && sk+1 < src->depth);

            if (src->type == MRI_UCHAR)
            {
              val000 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIseq_vox(src, si    , sj    , sk    , nframe));
              val001 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIseq_vox(src, si    , sj    , sk + 1, nframe));
              val010 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIseq_vox(src, si    , sj + 1, sk    , nframe));
              val011 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIseq_vox(src, si    , sj + 1, sk + 1, nframe));
              val100 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIseq_vox(src, si + 1, sj    , sk    , nframe));
              val101 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIseq_vox(src, si + 1, sj    , sk + 1, nframe));
              val110 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIseq_vox(src, si + 1, sj + 1, sk    , nframe));
              val111 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIseq_vox(src, si + 1, sj + 1, sk + 1, nframe));
            }

            if (si == 154 && sj == 134 && sk == 136)
              DiagBreak() ;

            if (src->type == MRI_SHORT)
            {
              val000 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRISseq_vox(src, si    , sj    , sk    , nframe));
              val001 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRISseq_vox(src, si    , sj    , sk + 1, nframe));
              val010 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRISseq_vox(src, si    , sj + 1, sk    , nframe));
              val011 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRISseq_vox(src, si    , sj + 1, sk + 1, nframe));
              val100 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRISseq_vox(src, si + 1, sj    , sk    , nframe));
              val101 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRISseq_vox(src, si + 1, sj    , sk + 1, nframe));
              val110 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRISseq_vox(src, si + 1, sj + 1, sk    , nframe));
              val111 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRISseq_vox(src, si + 1, sj + 1, sk + 1, nframe));
            }

            if (src->type == MRI_INT)
            {
              val000 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIIseq_vox(src, si    , sj    , sk    , nframe));
              val001 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIIseq_vox(src, si    , sj    , sk + 1, nframe));
              val010 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIIseq_vox(src, si    , sj + 1, sk    , nframe));
              val011 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIIseq_vox(src, si    , sj + 1, sk + 1, nframe));
              val100 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIIseq_vox(src, si + 1, sj    , sk    , nframe));
              val101 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIIseq_vox(src, si + 1, sj    , sk + 1, nframe));
              val110 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIIseq_vox(src, si + 1, sj + 1, sk    , nframe));
              val111 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIIseq_vox(src, si + 1, sj + 1, sk + 1, nframe));
            }

            if (src->type == MRI_LONG)
            {
              val000 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRILseq_vox(src, si    , sj    , sk    , nframe));
              val001 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRILseq_vox(src, si    , sj    , sk + 1, nframe));
              val010 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRILseq_vox(src, si    , sj + 1, sk    , nframe));
              val011 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRILseq_vox(src, si    , sj + 1, sk + 1, nframe));
              val100 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRILseq_vox(src, si + 1, sj    , sk    , nframe));
              val101 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRILseq_vox(src, si + 1, sj    , sk + 1, nframe));
              val110 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRILseq_vox(src, si + 1, sj + 1, sk    , nframe));
              val111 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRILseq_vox(src, si + 1, sj + 1, sk + 1, nframe));
            }

            if (src->type == MRI_FLOAT)
            {
              val000 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIFseq_vox(src, si    , sj    , sk    , nframe));
              val001 = ( !i_good_flag ||
                         !j_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIFseq_vox(src, si    , sj    , sk + 1, nframe));
              val010 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k_good_flag ? 0.0 :
                         (float)MRIFseq_vox(src, si    , sj + 1, sk    , nframe));
              val011 = ( !i_good_flag ||
                         !j1_good_flag ||
                         !k1_good_flag ? 0.0 :
                         (float)MRIFseq_vox(src, si    , sj + 1, sk + 1, nframe));
              val100 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIFseq_vox(src, si + 1, sj    , sk    , nframe));
              val101 = (!i1_good_flag ||
                        !j_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIFseq_vox(src, si + 1, sj    , sk + 1, nframe));
              val110 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k_good_flag ? 0.0 :
                        (float)MRIFseq_vox(src, si + 1, sj + 1, sk    , nframe));
              val111 = (!i1_good_flag ||
                        !j1_good_flag ||
                        !k1_good_flag ? 0.0 :
                        (float)MRIFseq_vox(src, si + 1, sj + 1, sk + 1, nframe));
            }

            if (resample_type == SAMPLE_TRILINEAR)
            {
              val = (1.0-si_f) * (1.0-sj_f) * (1.0-sk_f) * val000 +
                (1.0-si_f) * (1.0-sj_f) * (    sk_f) * val001 +
                (1.0-si_f) * (    sj_f) * (1.0-sk_f) * val010 +
                (1.0-si_f) * (    sj_f) * (    sk_f) * val011 +
                (    si_f) * (1.0-sj_f) * (1.0-sk_f) * val100 +
                (    si_f) * (1.0-sj_f) * (    sk_f) * val101 +
                (    si_f) * (    sj_f) * (1.0-sk_f) * val110 +
                (    si_f) * (    sj_f) * (    sk_f) * val111;
            }

            if (resample_type == SAMPLE_NEAREST)
            {
              if (si_f < 0.5)
              {
                if (sj_f < 0.5)
                {
                  if (sk_f < 0.5)
                    val = val000;
                  else
                    val = val001;
                }
                else
                {
                  if (sk_f < 0.5)
                    val = val010;
                  else
                    val = val011;
                }
              }
              else
              {
                if (sj_f < 0.5)
                {
                  if (sk_f < 0.5)
                    val = val100;
                  else
                    val = val101;
                }
                else
                {
                  if (sk_f < 0.5)
                    val = val110;
                  else
                    val = val111;
                }
              }
            }

            if (resample_type == SAMPLE_WEIGHTED)
            {
              /* unfinished */
              si_f2 = si_f * si_f;
              sj_f2 = sj_f * sj_f;
              sk_f2 = sk_f * sk_f;

              ii2 = 1. / (1 - 2*si_f + si_f2);
              ij2 = 1. / (1 - 2*sj_f + sj_f2);
              ik2 = 1. / (1 - 2*sk_f + sk_f2);

              isi_f2 = 1 / si_f2;
              isj_f2 = 1 / sj_f2;
              isk_f2 = 1 / sk_f2;

              w000 =   ii2    + ij2    + ik2;
              w001 =   ii2    + ij2    + isk_f2;
              w010 =   ii2    + isj_f2 + ik2;
              w011 =   ii2    + isj_f2 + isk_f2;
              w100 =   isi_f2 + ij2    + ik2;
              w101 =   isi_f2 + ij2    + isk_f2;
              w110 =   isi_f2 + isj_f2 + ik2;
              w111 =   isi_f2 + isj_f2 + isk_f2;

              w[0] = w[1] = w[2] = w[3] =
                w[4] = w[5] = w[6] = w[7] = 0.0;

              wi[0] = 0;
              wi[1] = 1;
              wi[2] = 2;
              wi[3] = 3;
              wi[4] = 4;
              wi[5] = 5;
              wi[6] = 6;
              wi[7] = 7;

              if (val001 == val000)
                wi[1] = 0;

              if (val010 == val001)
                wi[2] = 1;
              if (val010 == val000)
                wi[2] = 0;

              if (val011 == val010)
                wi[3] = 2;
              if (val011 == val001)
                wi[3] = 1;
              if (val011 == val000)
                wi[3] = 0;

              if (val100 == val011)
                wi[4] = 3;
              if (val100 == val010)
                wi[4] = 2;
              if (val100 == val001)
                wi[4] = 1;
              if (val100 == val000)
                wi[4] = 0;

              if (val101 == val100)
                wi[5] = 4;
              if (val101 == val011)
                wi[5] = 3;
              if (val101 == val010)
                wi[5] = 2;
              if (val101 == val001)
                wi[5] = 1;
              if (val101 == val000)
                wi[5] = 0;

              if (val110 == val101)
                wi[6] = 5;
              if (val110 == val100)
                wi[6] = 4;
              if (val110 == val011)
                wi[6] = 3;
              if (val110 == val010)
                wi[6] = 2;
              if (val110 == val001)
                wi[6] = 1;
              if (val110 == val000)
                wi[6] = 0;

              if (val111 == val110)
                wi[7] = 6;
              if (val111 == val101)
                wi[7] = 5;
              if (val111 == val100)
                wi[7] = 4;
              if (val111 == val011)
                wi[7] = 3;
              if (val111 == val010)
                wi[7] = 2;
              if (val111 == val001)
                wi[7] = 1;
              if (val111 == val000)
                wi[7] = 0;

              w[wi[0]] += w000;
              w[wi[1]] += w001;
              w[wi[2]] += w010;
              w[wi[3]] += w011;
              w[wi[4]] += w100;
              w[wi[5]] += w101;
              w[wi[6]] += w110;
              w[wi[7]] += w111;

              mwi = 0;

              if (w[1] > w[mwi])
                mwi = 1;
              if (w[2] > w[mwi])
                mwi = 2;
              if (w[3] > w[mwi])
                mwi = 3;
              if (w[4] > w[mwi])
                mwi = 4;
              if (w[5] > w[mwi])
                mwi = 5;
              if (w[6] > w[mwi])
                mwi = 6;
              if (w[7] > w[mwi])
                mwi = 7;

              if (mwi == 0)
                val = val000;
              if (mwi == 1)
                val = val001;
              if (mwi == 2)
                val = val010;
              if (mwi == 3)
                val = val011;
              if (mwi == 4)
                val = val100;
              if (mwi == 5)
                val = val101;
              if (mwi == 6)
                val = val110;
              if (mwi == 7)
                val = val111;

            }

          }

          if (dest->type == MRI_UCHAR)
          {
            if (val < UCHAR_MIN) val = UCHAR_MIN;
            if (val > UCHAR_MAX) val = UCHAR_MAX;
            MRIseq_vox(dest, di, dj, dk, nframe) = (unsigned char)nint(val);
          }
          if (dest->type == MRI_SHORT)
          {
            if (val < SHORT_MIN) val = SHORT_MIN;
            if (val > SHORT_MAX) val = SHORT_MAX;
            MRISseq_vox(dest, di, dj, dk, nframe) = (short)nint(val);
          }
          if (dest->type == MRI_INT)
          {
            if (val < INT_MIN) val = INT_MIN;
            if (val > INT_MAX) val = INT_MAX;
            MRIIseq_vox(dest, di, dj, dk, nframe) = (int)nint(val);
          }
          if (dest->type == MRI_LONG)
          {
            if (val < LONG_MIN) val = LONG_MIN;
            if (val > LONG_MAX) val = LONG_MAX;
            MRILseq_vox(dest, di, dj, dk, nframe) = (long)nint(val);
          }
          if (dest->type == MRI_FLOAT)
          {
            MRIFseq_vox(dest, di, dj, dk, nframe) = (float)val;
          }

        }
      }
    }
  }

  {
    MATRIX  *m_old_voxel_to_ras, *m_voxel_to_ras,
    *m_old_ras_to_voxel, *v_ras, *v_vox, *m_new_ras_to_voxel, *v_vox2 ;

    m_old_voxel_to_ras = MRIgetVoxelToRasXform(src) ;

    m_voxel_to_ras = MatrixMultiply(m_old_voxel_to_ras, m, NULL) ;

    dest->x_r = *MATRIX_RELT(m_voxel_to_ras,1,1)/dest->xsize;
    dest->x_a = *MATRIX_RELT(m_voxel_to_ras,2,1)/dest->xsize;
    dest->x_s = *MATRIX_RELT(m_voxel_to_ras,3,1)/dest->xsize;

    dest->y_r = *MATRIX_RELT(m_voxel_to_ras,1,2)/dest->ysize;
    dest->y_a = *MATRIX_RELT(m_voxel_to_ras,2,2)/dest->ysize;
    dest->y_s = *MATRIX_RELT(m_voxel_to_ras,3,2)/dest->ysize;

    dest->z_r = *MATRIX_RELT(m_voxel_to_ras,1,3)/dest->zsize;
    dest->z_a = *MATRIX_RELT(m_voxel_to_ras,2,3)/dest->zsize;
    dest->z_s = *MATRIX_RELT(m_voxel_to_ras,3,3)/dest->zsize;

    /* compute the RAS coordinates of the center of the dest. image
       and put them in c_r, c_a, and c_s.

       C = M * c_v
    */
    v_vox = VectorAlloc(4, MATRIX_REAL) ;

    /* voxel coords of center of dest image */
    VECTOR_ELT(v_vox,4) = 1.0 ;
    VECTOR_ELT(v_vox,1) = (dest->width)/2.0 ;
    VECTOR_ELT(v_vox,2) = (dest->height)/2.0 ;
    VECTOR_ELT(v_vox,3) = (dest->depth)/2.0 ;

    v_vox2 =
      MatrixMultiply(m, v_vox, NULL) ; /* voxel coords in source image */
    v_ras =
      MatrixMultiply(m_old_voxel_to_ras, v_vox2, NULL) ; /* ras cntr of dest */

    dest->c_r = VECTOR_ELT(v_ras, 1);
    dest->c_a = VECTOR_ELT(v_ras, 2);
    dest->c_s = VECTOR_ELT(v_ras, 3);
    dest->ras_good_flag = 1 ;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    {
      m_new_ras_to_voxel = MRIgetRasToVoxelXform(dest) ;
      m_old_ras_to_voxel = MRIgetRasToVoxelXform(src) ;
      V3_X(v_ras) = V3_Y(v_ras) = V3_Z(v_ras) = 0.0 ;
      MatrixMultiply(m_old_ras_to_voxel, v_ras, v_vox) ;
      printf("old RAS (0,0,0) -> (%2.0f, %2.0f, %2.0f)\n",
             VECTOR_ELT(v_vox,1), VECTOR_ELT(v_vox,2), VECTOR_ELT(v_vox,3)) ;
      MatrixMultiply(m_new_ras_to_voxel, v_ras, v_vox) ;
      printf("new RAS (0,0,0) -> (%2.0f, %2.0f, %2.0f)\n",
             VECTOR_ELT(v_vox,1), VECTOR_ELT(v_vox,2), VECTOR_ELT(v_vox,3)) ;
      MatrixFree(&m_new_ras_to_voxel) ;
      MatrixFree(&m_old_ras_to_voxel) ;
    }

    MatrixFree(&v_vox) ;
    MatrixFree(&v_vox2) ;
    MatrixFree(&v_ras) ;
    MatrixFree(&m_voxel_to_ras) ;
    MatrixFree(&m_old_voxel_to_ras) ;
  }

  if (bspline) MRIfreeBSpline(&bspline);
  MatrixFree(&dp);
  MatrixFree(&sp);
  MatrixFree(&m);

  return(dest);

} /* end MRIresample() */

int MRIlimits(MRI *mri, float *min, float *max)
{

  float val;
  int i, j, k, f;

  if (mri == NULL)
    return(NO_ERROR);

  if (mri->slices == NULL)
    return(NO_ERROR);

#if 1
  *min = *max = (float)MRIgetVoxVal(mri, 0, 0, 0, 0);
  for (f = 0 ; f < mri->nframes ; f++)
    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
        {
          val = (float)MRIgetVoxVal(mri, i, j, k, f);
          if (val < *min)
            *min = val;
          if (val > *max)
            *max = val;
        }
#else
  if (mri->type == MRI_UCHAR)
  {
    *min = *max = (float)MRIvox(mri, 0, 0, 0);
    for (f = 0 ; f < mri->nframes ; f++)
      for (i = 0;i < mri->width;i++)
        for (j = 0;j < mri->height;j++)
          for (k = 0;k < mri->depth;k++)
          {
            val = (float)MRIseq_vox(mri, i, j, k, f);
            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;
          }

  }

  if (mri->type == MRI_SHORT)
  {

    *min = *max = (float)MRISvox(mri, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
        {
          val = (float)MRISvox(mri, i, j, k);
          if (val < *min)
            *min = val;
          if (val > *max)
            *max = val;
        }

  }

  if (mri->type == MRI_INT)
  {

    *min = *max = (float)MRILvox(mri, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
        {
          val = (float)MRILvox(mri, i, j, k);
          if (val < *min)
            *min = val;
          if (val > *max)
            *max = val;
        }

  }

  if (mri->type == MRI_LONG)
  {

    *min = *max = (float)MRILvox(mri, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
        {
          val = (float)MRILvox(mri, i, j, k);
          if (val < *min)
            *min = val;
          if (val > *max)
            *max = val;
        }

  }

  if (mri->type == MRI_FLOAT)
  {

    *min = *max = (float)MRIFvox(mri, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
        {
          val = (float)MRIFvox(mri, i, j, k);
          if (val < *min)
            *min = val;
          if (val > *max)
            *max = val;
        }

  }
#endif

  return(NO_ERROR);

} /* end MRIlimits() */

int MRIprintStats(MRI *mri, FILE *stream)
{

  float min, max, mean, std;
  int n;
  double com[3];

  MRIstats(mri, &min, &max, &n, &mean, &std);

  fprintf(stream, "%d values\n", n);
  fprintf(stream, "min = %g\n", min);
  fprintf(stream, "max = %g\n", max);
  fprintf(stream, "mean = %g\n", mean);
  fprintf(stream, "std = %g\n", std);

  MRIcenterOfMass(mri, com, 0);

  fprintf(stream, "com = %g %g %g\n", com[0], com[1], com[2]);

  return(NO_ERROR);

} /* end MRIprintStats() */

int MRIstats
(MRI *mri, float *min, float *max, int *n_voxels, float *mean, float *std)
{

  float val;
  float sum, sq_sum;
  int i, j, k, t;
  float n;

  if (mri == NULL)
    return(NO_ERROR);

  if (mri->slices == NULL)
    return(NO_ERROR);

  sum = sq_sum = 0.0;
  *n_voxels = mri->width * mri->height * mri->depth * mri->nframes;
  n = (float)(*n_voxels);

  if (mri->type == MRI_UCHAR)
  {

    *max = *min = MRIseq_vox(mri, 0, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
          for (t = 0;t < mri->nframes;t++)
          {

            val = (float)MRIseq_vox(mri, i, j, k, t);

            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if (mri->type == MRI_SHORT)
  {

    *min = *max = (float)MRISseq_vox(mri, 0, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
          for (t = 0;t < mri->nframes;t++)
          {

            val = (float)MRISseq_vox(mri, i, j, k, t);

            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if (mri->type == MRI_INT)
  {

    *min = *max = (float)MRILseq_vox(mri, 0, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
          for (t = 0;t < mri->nframes;t++)
          {

            val = (float)MRILseq_vox(mri, i, j, k, t);

            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if (mri->type == MRI_LONG)
  {

    *min = *max = (float)MRILseq_vox(mri, 0, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
          for (t = 0;t < mri->nframes;t++)
          {

            val = (float)MRILseq_vox(mri, i, j, k, t);

            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  if (mri->type == MRI_FLOAT)
  {

    *min = *max = (float)MRIFseq_vox(mri, 0, 0, 0, 0);

    for (i = 0;i < mri->width;i++)
      for (j = 0;j < mri->height;j++)
        for (k = 0;k < mri->depth;k++)
          for (t = 0;t < mri->nframes;t++)
          {

            val = (float)MRIFseq_vox(mri, i, j, k, t);

            if (val < *min)
              *min = val;
            if (val > *max)
              *max = val;

            sum += val;
            sq_sum += val * val;

          }

  }

  *mean = sum / n;
  *std = (sq_sum - n * (*mean) * (*mean)) / (n - 1.0);

  *std = sqrt(*std);

  return(NO_ERROR);

} /* end MRIstats() */

float MRIvolumeDeterminant(MRI *mri)
{

  MATRIX *m;
  float det;

  m = MatrixAlloc(4, 4, MATRIX_REAL);
  if (m == NULL)
    return(0.0);

  stuff_four_by_four(m, mri->x_r, mri->y_r, mri->z_r, 0.0,
                     mri->x_a, mri->y_a, mri->z_a, 0.0,
                     mri->x_s, mri->y_s, mri->z_s, 0.0,
                     0.0,      0.0,      0.0, 1.0);

  det = MatrixDeterminant(m);

  MatrixFree(&m);

  return(det);

} /* end MRIvolumeDeterminant() */

int stuff_four_by_four(MATRIX *m, float m11, float m12, float m13, float m14,
                       float m21, float m22, float m23, float m24,
                       float m31, float m32, float m33, float m34,
                       float m41, float m42, float m43, float m44)
{

  if (m == NULL)
  {
    ErrorReturn
    (ERROR_BADPARM,
     (ERROR_BADPARM,
      "stuff_four_by_four(): matrix is NULL"));
  }

  if (m->rows != 4 || m->cols != 4)
  {
    ErrorReturn
    (ERROR_BADPARM,
     (ERROR_BADPARM,
      "stuff_four_by_four(): matrix is not four-by-four"));
  }

  *MATRIX_RELT(m, 1, 1) = m11;
  *MATRIX_RELT(m, 1, 2) = m12;
  *MATRIX_RELT(m, 1, 3) = m13;
  *MATRIX_RELT(m, 1, 4) = m14;
  *MATRIX_RELT(m, 2, 1) = m21;
  *MATRIX_RELT(m, 2, 2) = m22;
  *MATRIX_RELT(m, 2, 3) = m23;
  *MATRIX_RELT(m, 2, 4) = m24;
  *MATRIX_RELT(m, 3, 1) = m31;
  *MATRIX_RELT(m, 3, 2) = m32;
  *MATRIX_RELT(m, 3, 3) = m33;
  *MATRIX_RELT(m, 3, 4) = m34;
  *MATRIX_RELT(m, 4, 1) = m41;
  *MATRIX_RELT(m, 4, 2) = m42;
  *MATRIX_RELT(m, 4, 3) = m43;
  *MATRIX_RELT(m, 4, 4) = m44;

  return(NO_ERROR);

} /* end stuff_four_by_four() */

int apply_i_to_r(MRI *mri, MATRIX *m)
{

  float x_r, x_a, x_s;
  float y_r, y_a, y_s;
  float z_r, z_a, z_s;
  float mag;
  MATRIX *origin, *c;

  x_r = *MATRIX_RELT(m, 1, 1);
  x_a = *MATRIX_RELT(m, 2, 1);
  x_s = *MATRIX_RELT(m, 3, 1);
  mag = sqrt(x_r*x_r + x_a*x_a + x_s*x_s);
  mri->x_r = x_r / mag;
  mri->x_a = x_a / mag;
  mri->x_s = x_s / mag;
  mri->xsize = mag;

  y_r = *MATRIX_RELT(m, 1, 2);
  y_a = *MATRIX_RELT(m, 2, 2);
  y_s = *MATRIX_RELT(m, 3, 2);
  mag = sqrt(y_r*y_r + y_a*y_a + y_s*y_s);
  mri->y_r = y_r / mag;
  mri->y_a = y_a / mag;
  mri->y_s = y_s / mag;
  mri->ysize = mag;

  z_r = *MATRIX_RELT(m, 1, 3);
  z_a = *MATRIX_RELT(m, 2, 3);
  z_s = *MATRIX_RELT(m, 3, 3);
  mag = sqrt(z_r*z_r + z_a*z_a + z_s*z_s);
  mri->z_r = z_r / mag;
  mri->z_a = z_a / mag;
  mri->z_s = z_s / mag;
  mri->zsize = mag;

  origin = MatrixAlloc(4, 1, MATRIX_REAL);
  if (origin == NULL)
  {
    ErrorReturn
    (ERROR_BADPARM,
     (ERROR_BADPARM,
      "apply_i_to_r(): error allocating matrix"));
  }
  *MATRIX_RELT(origin, 1, 1) = (double)mri->width / 2.0;
  *MATRIX_RELT(origin, 2, 1) = (double)mri->height / 2.0;
  *MATRIX_RELT(origin, 3, 1) = (double)mri->depth / 2.0;
  *MATRIX_RELT(origin, 4, 1) = 1.0;

  c = MatrixMultiply(m, origin, NULL);
  if (c == NULL)
  {
    MatrixFree(&origin);
    ErrorReturn
    (ERROR_BADPARM,
     (ERROR_BADPARM,
      "apply_i_to_r(): error multiplying matrices"));
  }

  mri->c_r = *MATRIX_RELT(c, 1, 1);
  mri->c_a = *MATRIX_RELT(c, 2, 1);
  mri->c_s = *MATRIX_RELT(c, 3, 1);

  MatrixFree(&origin);
  MatrixFree(&c);

  mri->ras_good_flag = 1;

  return(NO_ERROR);

} /* end apply_i_to_r() */

MATRIX *
MRIrasXformToVoxelXform(MRI *mri_src, MRI *mri_dst, MATRIX *m_ras_xform,
                        MATRIX *m_voxel_xform)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras, *m_tmp ;

  if (!mri_dst)
    mri_dst = mri_src ;   /* assume they will be in the same space */

  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_src) ;
  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_dst) ;

  m_tmp = MatrixMultiply(m_ras_xform, m_voxel_to_ras, NULL) ;
  m_voxel_xform = MatrixMultiply(m_ras_to_voxel, m_tmp, m_voxel_xform) ;

  MatrixFree(&m_voxel_to_ras);
  MatrixFree(&m_ras_to_voxel);
  MatrixFree(&m_tmp);

  return(m_voxel_xform) ;
}

MATRIX *
MRIvoxelXformToRasXform(MRI *mri_src, MRI *mri_dst, MATRIX *m_voxel_xform,
                        MATRIX *m_ras_xform)
{
  MATRIX *m_ras_to_voxel, *m_voxel_to_ras, *m_tmp ;

  if (!mri_dst)
    mri_dst = mri_src ;  /* assume they will be in the same space */

  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_src) ;
  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_dst) ;

  m_tmp = MatrixMultiply(m_voxel_xform, m_ras_to_voxel, NULL) ;
  m_ras_xform = MatrixMultiply(m_voxel_to_ras, m_tmp, m_ras_xform) ;

  MatrixFree(&m_voxel_to_ras);
  MatrixFree(&m_ras_to_voxel);
  MatrixFree(&m_tmp);

  return(m_ras_xform) ;
}


int
MRIsetVoxelToRasXform(MRI *mri, MATRIX *m_vox2ras)
{
  float ci, cj, ck;

  mri->x_r = *MATRIX_RELT(m_vox2ras, 1, 1) / mri->xsize ;
  mri->y_r = *MATRIX_RELT(m_vox2ras, 1, 2) / mri->ysize ;
  mri->z_r = *MATRIX_RELT(m_vox2ras, 1, 3) / mri->zsize ;

  mri->x_a = *MATRIX_RELT(m_vox2ras, 2, 1) / mri->xsize ;
  mri->y_a = *MATRIX_RELT(m_vox2ras, 2, 2) / mri->ysize ;
  mri->z_a = *MATRIX_RELT(m_vox2ras, 2, 3) / mri->zsize ;

  mri->x_s = *MATRIX_RELT(m_vox2ras, 3, 1) / mri->xsize ;
  mri->y_s = *MATRIX_RELT(m_vox2ras, 3, 2) / mri->ysize ;
  mri->z_s = *MATRIX_RELT(m_vox2ras, 3, 3) / mri->zsize ;

  ci = (mri->width) / 2.0;
  cj = (mri->height) / 2.0;
  ck = (mri->depth) / 2.0;
  mri->c_r = *MATRIX_RELT(m_vox2ras, 1, 4) +
             (*MATRIX_RELT(m_vox2ras, 1, 1) * ci +
              *MATRIX_RELT(m_vox2ras, 1, 2) * cj +
              *MATRIX_RELT(m_vox2ras, 1, 3) * ck);
  mri->c_a = *MATRIX_RELT(m_vox2ras, 2, 4) +
             (*MATRIX_RELT(m_vox2ras, 2, 1) * ci +
              *MATRIX_RELT(m_vox2ras, 2, 2) * cj +
              *MATRIX_RELT(m_vox2ras, 2, 3) * ck);
  mri->c_s = *MATRIX_RELT(m_vox2ras, 3, 4) +
             (*MATRIX_RELT(m_vox2ras, 3, 1) * ci +
              *MATRIX_RELT(m_vox2ras, 3, 2) * cj +
              *MATRIX_RELT(m_vox2ras, 3, 3) * ck);
  mri->ras_good_flag = 1;
  MRIreInitCache(mri);
  return(NO_ERROR) ;
}


/* eof */
MRI *
MRIscaleMeanIntensities(MRI *mri_src, MRI *mri_ref, MRI *mri_dst)
{
  int    width, height, depth, x, y, z, val ;
  double ref_mean, src_mean, nref_vox, nsrc_vox, scale ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  width = mri_dst->width ;
  height = mri_dst->height ;
  depth = mri_dst->depth;

  nref_vox = nsrc_vox = src_mean = ref_mean = 0.0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_ref,x,y,z) > 10)
        {
          nref_vox++ ;
          ref_mean += (double)MRIvox(mri_ref, x, y, z) ;
        }
        if (MRIvox(mri_src,x,y,z) > 10)
        {
          src_mean += (double)MRIvox(mri_src, x, y, z) ;
          nsrc_vox++ ;
        }
      }
    }
  }

  ref_mean /= nref_vox ;
  src_mean /= nsrc_vox ;
  fprintf(stderr, "mean brightnesses: ref = %2.1f, in = %2.1f\n",
          ref_mean, src_mean) ;
  scale = ref_mean / src_mean ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        val = MRIvox(mri_src, x, y, z) ;
        val = nint(val*scale) ;
        if (val > 255)
          val = 255 ;
        MRIvox(mri_src, x, y, z) = val ;
      }
    }
  }

  return(mri_dst) ;
}

MRI *MRIsmoothParcellation(MRI *mri, int smooth_parcellation_count)
{

  MRI *mri2;
  int i, j, k;
  short vals[26];
  int counts[32768];
  int c;

  if (mri->type != MRI_SHORT)
  {
    ErrorReturn
    (NULL,
     (ERROR_UNSUPPORTED,
      "MRIsmoothParcellation(): only supported for shorts data"));
  }

  mri2 = MRIcopy(mri, NULL);
  if (mri2 == NULL)
  {
    ErrorReturn
    (NULL,
     (ERROR_NOMEMORY,
      "MRIsmoothParcellation(): error copying structre"));
  }

  for (i = 1;i < mri->width-1;i++)
  {
    for (j = 1;j < mri->height-1;j++)
    {
      for (k = 1;k < mri->depth-1;k++)
      {

        memset(counts, 0x00, 32768 * sizeof(int));

        vals[ 0] = MRISvox(mri, i+1, j+1, k+1);
        vals[ 1] = MRISvox(mri, i+1, j+1, k  );
        vals[ 2] = MRISvox(mri, i+1, j+1, k-1);
        vals[ 3] = MRISvox(mri, i+1, j  , k+1);
        vals[ 4] = MRISvox(mri, i+1, j  , k  );
        vals[ 5] = MRISvox(mri, i+1, j  , k-1);
        vals[ 6] = MRISvox(mri, i+1, j-1, k+1);
        vals[ 7] = MRISvox(mri, i+1, j-1, k  );
        vals[ 8] = MRISvox(mri, i+1, j-1, k-1);

        vals[ 9] = MRISvox(mri, i  , j+1, k+1);
        vals[10] = MRISvox(mri, i  , j+1, k  );
        vals[11] = MRISvox(mri, i  , j+1, k-1);
        vals[12] = MRISvox(mri, i  , j  , k+1);
        /* --- ignore the voxel itself --- */
        vals[13] = MRISvox(mri, i  , j  , k-1);
        vals[14] = MRISvox(mri, i  , j-1, k+1);
        vals[15] = MRISvox(mri, i  , j-1, k  );
        vals[16] = MRISvox(mri, i  , j-1, k-1);

        vals[17] = MRISvox(mri, i-1, j+1, k+1);
        vals[18] = MRISvox(mri, i-1, j+1, k  );
        vals[19] = MRISvox(mri, i-1, j+1, k-1);
        vals[20] = MRISvox(mri, i-1, j  , k+1);
        vals[21] = MRISvox(mri, i-1, j  , k  );
        vals[22] = MRISvox(mri, i-1, j  , k-1);
        vals[23] = MRISvox(mri, i-1, j-1, k+1);
        vals[24] = MRISvox(mri, i-1, j-1, k  );
        vals[25] = MRISvox(mri, i-1, j-1, k-1);

        for (c = 0;c < 26;c++)
          counts[vals[c]]++;

        for (c = 0;c < 26;c++)
          if (counts[vals[c]] >= smooth_parcellation_count)
            MRISvox(mri2, i, j, k) = vals[c];

      }
    }
  }

  return(mri2);

} /* end MRIsmoothParcellation() */

int
MRIeraseBorderPlanes(MRI *mri, int border_size)
{
  int  x, y, z, i, f ;

  for (f = 0 ; f < mri->nframes ; f++)
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
      {
        for (i = 0 ; i < border_size ; i++)
        {
          MRIsetVoxVal(mri, x, y, i, f, 0)  ;
          MRIsetVoxVal(mri, x, y, mri->depth-i-1, f, 0) ;
      }
    }

  for (f = 0 ; f < mri->nframes ; f++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        for (i = 0 ; i < border_size ; i++)
        {
          MRIsetVoxVal(mri, i, y, z, f, 0) ;
          MRIsetVoxVal(mri, mri->width-i-1, y, z, f,0) ;
        }
      }

  for (f = 0 ; f < mri->nframes ; f++)
    for (x = 0 ; x < mri->width ; x++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        for (i = 0 ; i < border_size ; i++)
        {
          MRIsetVoxVal(mri, x, i, z, f, 0) ;
          MRIsetVoxVal(mri, x, mri->height-i-1, z, f, 0) ;
        }
      }

  return(NO_ERROR) ;
}
int
MRIcopyPulseParameters(MRI *mri_src, MRI *mri_dst)
{
  mri_dst->flip_angle = mri_src->flip_angle ;
  mri_dst->tr = mri_src->tr ;
  mri_dst->te = mri_src->te ;
  mri_dst->ti = mri_src->ti ;
  return(NO_ERROR) ;
}
float
MRIfindNearestNonzero(MRI *mri,
                      int wsize,
                      double xr, double yr, double zr,
                      float max_dist)
{
  int   xk, yk, zk, xi, yi, zi, whalf, x, y, z ;
  float dist, min_dist, min_val, dx, dy, dz ;
  dist = 0.0;

  x = mri->xi[nint(xr)] ;
  y = mri->yi[nint(yr)] ;
  z = mri->zi[nint(zr)] ;
  min_val = MRIgetVoxVal(mri, x, y, z, 0) ;
  if (min_val > 0)
    return(min_val) ;

  min_dist = 100000 ;
  min_val = 0 ;
  whalf = (wsize-1)/2 ;
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    dz = zi-zr ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dy = yi-yr ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        dx = xi-xr ;
        if (MRIgetVoxVal(mri, xi, yi, zi,0) > 0)
        {
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            min_val = MRIgetVoxVal(mri, xi, yi, zi,0) ;
          }
        }
      }
    }
  }
  if (max_dist > 0 && dist > max_dist)
    return(0) ;
  return(min_val) ;
}
int
MRIcountNonzeroInNbhd(MRI *mri, int wsize, int x, int y, int z)
{
  int   xk, yk, zk, xi, yi, zi, whalf, total ;

  whalf = (wsize-1)/2 ;
  for (total = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi,0) > 0)
          total++ ;
      }
    }
  }
  return(total) ;
}
int
MRIareNonzeroInNbhd(MRI *mri, int wsize, int x, int y, int z)
{
  int   xk, yk, zk, xi, yi, zi, whalf ;

  whalf = (wsize-1)/2 ;
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi,0) > 0)
          return(1) ;
      }
    }
  }
  return(0) ;
}
float
MRIfindNearestNonzeroLocation(MRI *mri, int wsize,
                              double xr, double yr, double zr,
                              int *pxv, int *pyv, int *pzv)
{
  int   xk, yk, zk, xi, yi, zi, whalf, x, y, z ;
  float dist, min_dist, min_val, dx, dy, dz ;

  x = nint(xr) ;
  y = nint(yr) ;
  z = nint(zr) ;
  if (MRIvox(mri, x, y, z) > 0)
    return((float)MRIvox(mri, x, y, z)) ;

  min_dist = 100000 ;
  min_val = 0 ;
  whalf = (wsize-1)/2 ;
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    dz = zi-zr ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      dy = yi-yr ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        dx = xi-xr ;
        if (MRIvox(mri, xi, yi, zi) > 0)
        {
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < min_dist)
          {
            if (pxv)
            {
              *pxv = xi ;
              *pyv = yi ;
              *pzv = zi ;
            }
            min_dist = dist ;
            min_val = MRIvox(mri, xi, yi, zi) ;
          }
        }
      }
    }
  }
  return(min_val) ;
}

MRI *
MRIfromTalairach(MRI *mri_src, MRI *mri_dst)
{
  int  x, y, z, xv, yv, zv ;
  double xt, yt, zt, xn, yn, zn, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (x = 0 ; x < mri_dst->width ; x++)
  {
    xn = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      yn = (double)y ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        zn = (double)z ;
        MRIvoxelToTalairachVoxel(mri_src, xn, yn, zn, &xt, &yt, &zt) ;
        xv = nint(xt) ;
        yv = nint(yt) ;
        zv = nint(zt) ;
        if ((xv >= 0 && xv < mri_src->width) &&
            (yv >= 0 && yv < mri_src->height) &&
            (zv >= 0 && zv < mri_src->depth))
        {
          MRIsampleVolume(mri_src, xt, yt, zt, &val) ;
          MRIvox(mri_dst, x, y, z) = val ;
        }
        else
          MRIvox(mri_dst, x, y, z) = 0 ;
      }
    }
  }

  return(mri_dst) ;
}
MRI *
MRItoTalairach(MRI *mri_src, MRI *mri_dst)
{
  int  x, y, z, xv, yv, zv ;
  double xt, yt, zt, xn, yn, zn, val ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (x = 0 ; x < mri_dst->width ; x++)
  {
    xt = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      yt = (double)y ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        zt = (double)z ;
        MRItalairachVoxelToVoxel(mri_src, xt, yt, zt, &xn, &yn, &zn) ;
        xv = nint(xn) ;
        yv = nint(yn) ;
        zv = nint(zt) ;
        if ((xv >= 0 && xv < mri_src->width) &&
            (yv >= 0 && yv < mri_src->height) &&
            (zv >= 0 && zv < mri_src->depth))
        {
          MRIsampleVolume(mri_src, xn, yn, zn, &val) ;
          MRIvox(mri_dst, x, y, z) = val ;
        }
        else
          MRIvox(mri_dst, x, y, z) = 0 ;
      }
    }
  }

  return(mri_dst) ;
}
/*-------------------------------------------------------------------
  MRIlog10() - computes the log10 of the abs value at each voxel and
  frame. If a value is zero, the result is set to 10000000000.0. If
  the negflag is set, then -log10 is computed. If a mask is specified,
  voxels outside of the mask are 0'ed. The result is stored
  in outmri. If outmri is NULL, the output MRI is alloced and its
  pointer returned.
  ------------------------------------------------------------------*/
MRI *MRIlog10(MRI *inmri, MRI *mask, MRI *outmri, int negflag)
{
  int c, r, s, f;
  double val,m;

  if (outmri==NULL)
  {
    outmri = MRIallocSequence(inmri->width, inmri->height, inmri->depth,
                              MRI_FLOAT, inmri->nframes);

    if (outmri==NULL)
    {
      printf("ERROR: fMRIlog10: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(inmri,outmri);
  }
  else
  {
    if (inmri->width   != outmri->width  ||
        inmri->height  != outmri->height ||
        inmri->depth   != outmri->depth  ||
        inmri->nframes != outmri->nframes)
    {
      printf("ERROR: MRIlog10: output dimension mismatch\n");
      return(NULL);
    }
    if (outmri->type != MRI_FLOAT)
    {
      printf("ERROR: MRIlog10: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for (c=0; c < inmri->width; c++)
  {
    for (r=0; r < inmri->height; r++)
    {
      for (s=0; s < inmri->depth; s++)
      {
        if (mask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5)
          {
            for (f=0; f < inmri->nframes; f++)
              MRIsetVoxVal(outmri,c,r,s,f,0);
            continue;
          }
        }
        for (f=0; f < inmri->nframes; f++)
        {
          val = MRIgetVoxVal(inmri, c, r, s, f);
          if (val == 0)     MRIFseq_vox(outmri,c,r,s,f) = 10000000000.0;
          else
          {
            if (negflag)
            {
              if (val < 0) MRIFseq_vox(outmri,c,r,s,f) =  log10(fabs(val));
              else        MRIFseq_vox(outmri,c,r,s,f) = -log10(val);
            }
            else
            {
              if (val < 0) MRIFseq_vox(outmri,c,r,s,f) = -log10(fabs(val));
              else        MRIFseq_vox(outmri,c,r,s,f) =  log10(val);
            }
          }
        }
      }
    }
  }

  return(outmri);
}
/*-------------------------------------------------------
  MRIlog() - computes natural log: out = a*log(abs(b*in)).
  If a=0, then it is ignored (ie, set to 1).
  If b=0, then it is ignored (ie, set to 1).
  If in=0, then EPSILON is used.
  If mask is non-NULL, then sets vox=0 where mask < 0.5
  Input can be any data type.
  Output is float.
  Can be run in-place.
  Note: 3 to 8 times faster than using MRIgetVoxVal(), which
  is important because this is run on raw data (dti).
  -------------------------------------------------------*/
#define EPSILON 0.25
MRI *MRIlog(MRI *in, MRI *mask, double a, double b, MRI *out)
{
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  float m;
  float   *pout=NULL;
  void    *pin=NULL;
  double  v, aa, bb;
  int sz;

  if (a == 0) aa = 1;
  else       aa = a;
  if (b == 0) bb = 1;
  else       bb = b;
  v = 0.0;

  ncols = in->width;
  nrows = in->height;
  nslices = in->depth;
  nframes = in->nframes;

  if (out==NULL)
  {
    out = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (out==NULL)
    {
      printf("ERROR: MRIlog: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(in,out); // ordinarily would need to change nframes
  }
  if (out->type != MRI_FLOAT)
  {
    printf("ERROR: MRIlog: output must be of type float\n");
    return(NULL);
  }
  if (out->width != ncols   || out->height != nrows ||
      out->depth != nslices || out->nframes != nframes)
  {
    printf("ERROR: MRIlog: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the input data type
  sz = 0;
  switch (in->type)
  {
  case MRI_UCHAR:
    sz = sizeof(char);
    break;
  case MRI_SHORT:
    sz = sizeof(short);
    break;
  case MRI_INT:
    sz = sizeof(int);
    break;
  case MRI_LONG:
    sz = sizeof(long);
    break;
  case MRI_FLOAT:
    sz = sizeof(float);
    break;
  }

  n = 0;
  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        // Pointers to the start of the column
        pin  = (void *)  in->slices[n][r];
        pout = (float *) out->slices[n][r];
        for (c=0; c<ncols; c++)
        {
          if (mask)
          {
            m = MRIgetVoxVal(mask,c,r,s,0);
            if (m < 0.5)
            {
              // must increment pointers
              *pout++ = 0.0;
              pin += sz;
              continue;
            }
          }
          // Get input value
          switch (in->type)
          {
          case MRI_UCHAR:
            v = (double)(*((char *)pin));
            break;
          case MRI_SHORT:
            v = (double)(*((short*)pin));
            break;
          case MRI_INT:
            v = (double)(*((int  *)pin));
            break;
          case MRI_LONG:
            v = (double)(*((long *)pin));
            break;
          case MRI_FLOAT:
            v = (double)(*((float*)pin));
            break;
          }
          if (v==0) v = EPSILON;
          // Compute and set output
          *pout = aa*log(fabs(bb*v));
          // Increment
          pout ++;
          pin += sz;
        } // cols
      } // rows
      n++;
    } // slices
  } // frames
  return(out);
}

/*---------------------------------------------------------------------
  MRIrandn() - fills an MRI structure with values sampled from a
  normal distribution with mean avg and standard devation stddev.
  --------------------------------------------------------*/
MRI *MRIrandn(int ncols, int nrows, int nslices, int nframes,
              float avg, float stddev, MRI *mri)
{
  int c, r, s, f;

  if (mri==NULL)
  {
    mri = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (mri==NULL)
    {
      printf("ERROR: MRIrandn: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mri->width != ncols   || mri->height != nrows ||
        mri->depth != nslices || mri->nframes != nframes)
    {
      printf("ERROR: MRIrandn: dimension mismatch\n");
      return(NULL);
    }
    if (mri->type != MRI_FLOAT)
    {
      printf("ERROR: MRIrandn: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        for (c=0; c<ncols; c++)
        {
          MRIFseq_vox(mri,c,r,s,f) =
            stddev*PDFgaussian() + avg;
        }
      }
    }
  }

  return(mri);
}
/*---------------------------------------------------------------------
  MRIrande() - fills an MRI structure with values sampled from an
  Erlang distribution with mean avg and order order. The variance
  will be (avg^2)/order. Theoretical distribution is
  r = order
  pdf(x) = r*((r*(x-avg+1))^(r-1)) * exp(-r*(x-avg+1)) / (r-1)!
  when order=1, this generates an exponential distribution.
  --------------------------------------------------------*/
MRI *MRIrande(int ncols, int nrows, int nslices, int nframes,
              float avg, int order, MRI *mri)
{
  int c, r, s, f;

  if (mri==NULL)
  {
    mri = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (mri==NULL)
    {
      printf("ERROR: MRIrande: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mri->width != ncols   || mri->height != nrows ||
        mri->depth != nslices || mri->nframes != nframes)
    {
      printf("ERROR: MRIrande: dimension mismatch\n");
      return(NULL);
    }
    if (mri->type != MRI_FLOAT)
    {
      printf("ERROR: MRIrande: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  avg = avg - 1; // PDFrande() already has average of 1
  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        for (c=0; c<ncols; c++)
        {
          MRIFseq_vox(mri,c,r,s,f) =
            PDFerlang(order) + avg;
        }
      }
    }
  }

  return(mri);
}
/*---------------------------------------------------------------------
  MRIdrand48() - fills an MRI structure with values sampled from a
  the drand48 uniform random number generator. The user must specify
  the min and max. Note: the stddev = (max-min)*0.289. If mri is NULL,
  it will alloc a MRI_FLOAT volume, otherwise, it will use the type
  as specified in mri.
  --------------------------------------------------------*/
MRI *MRIdrand48(int ncols, int nrows, int nslices, int nframes,
                float min, float max, MRI *mri)
{
  int c, r, s, f, n;
  float range, v;
  BUFTYPE *pmri=NULL;
  short   *psmri=NULL;
  int     *pimri=NULL;
  long    *plmri=NULL;
  float   *pfmri=NULL;

  if (mri==NULL)
  {
    mri = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (mri==NULL)
    {
      printf("ERROR: MRIdrand48: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mri->width != ncols   || mri->height != nrows ||
        mri->depth != nslices || mri->nframes != nframes)
    {
      printf("ERROR: MRIdrand48: dimension mismatch\n");
      return(NULL);
    }
  }

  range = max-min;
  n = 0;
  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          pmri  =           mri->slices[n][r];
          break;
        case MRI_SHORT:
          psmri = (short *) mri->slices[n][r];
          break;
        case MRI_INT:
          pimri = (int *)   mri->slices[n][r];
          break;
        case MRI_LONG:
          plmri = (long *)  mri->slices[n][r];
          break;
        case MRI_FLOAT:
          pfmri = (float *) mri->slices[n][r];
          break;
        }
        for (c=0; c<ncols; c++)
        {
          v = range*drand48() + min;
          switch (mri->type)
          {
          case MRI_UCHAR:
            *pmri++  = (BUFTYPE) nint(v);
            break;
          case MRI_SHORT:
            *psmri++ = (short) nint(v);
            break;
          case MRI_INT:
            *pimri++ = (int)   nint(v);
            break;
          case MRI_LONG:
            *plmri++ = (long)  nint(v);
            break;
          case MRI_FLOAT:
            *pfmri++ = (float) v;
            break;
          }
        }
      }
      n++;
    }
  }

  return(mri);
}
/*---------------------------------------------------------------------
  MRIsampleCDF() - fills an MRI structure with values sampled from a
  the given CDF. See PDFsampleCDF(). CDF[n] is the probability that
  the random number is <= xCDF[n].
  --------------------------------------------------------------------*/
MRI *MRIsampleCDF(int ncols, int nrows, int nslices, int nframes,
                  double *xCDF, double *CDF, int nCDF, MRI *mri)
{
  int c, r, s, f;

  if (mri==NULL)
  {
    mri = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (mri==NULL)
    {
      printf("ERROR: MRIsampleCDF: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mri->width != ncols   || mri->height != nrows ||
        mri->depth != nslices || mri->nframes != nframes)
    {
      printf("ERROR: MRIsampleCDF: dimension mismatch\n");
      return(NULL);
    }
    if (mri->type != MRI_FLOAT)
    {
      printf("ERROR: MRIsampleCDF: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        for (c=0; c<ncols; c++)
        {
          MRIFseq_vox(mri,c,r,s,f) = PDFsampleCDF(xCDF,CDF,nCDF);
        }
      }
    }
  }

  return(mri);
}
/*---------------------------------------------------------------------
  MRIconst() - fills an MRI structure with the given value. If mri is
  NULL, it will alloc a MRI_FLOAT volume, otherwise, it will use the type
  as specified in mri.
  --------------------------------------------------------*/
MRI *MRIconst(int ncols, int nrows, int nslices, int nframes,
              float val, MRI *mri)
{
  int c, r, s, f, n;
  BUFTYPE *pmri=NULL;
  short   *psmri=NULL;
  int     *pimri=NULL;
  long    *plmri=NULL;
  float   *pfmri=NULL;

  if (mri==NULL)
  {
    mri = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (mri==NULL)
    {
      printf("ERROR: MRIdconst: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mri->width != ncols   || mri->height != nrows ||
        mri->depth != nslices || mri->nframes != nframes)
    {
      printf("ERROR: MRIconst: dimension mismatch\n");
      return(NULL);
    }
  }

  n = 0;
  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        switch (mri->type)
        {
        case MRI_UCHAR:
          pmri  =           mri->slices[n][r];
          break;
        case MRI_SHORT:
          psmri = (short *) mri->slices[n][r];
          break;
        case MRI_INT:
          pimri = (int *)   mri->slices[n][r];
          break;
        case MRI_LONG:
          plmri = (long *)  mri->slices[n][r];
          break;
        case MRI_FLOAT:
          pfmri = (float *) mri->slices[n][r];
          break;
        }
        for (c=0; c<ncols; c++)
        {
          switch (mri->type)
          {
          case MRI_UCHAR:
            *pmri++  = (BUFTYPE) nint(val);
            break;
          case MRI_SHORT:
            *psmri++ = (short) nint(val);
            break;
          case MRI_INT:
            *pimri++ = (int)   nint(val);
            break;
          case MRI_LONG:
            *plmri++ = (long)  nint(val);
            break;
          case MRI_FLOAT:
            *pfmri++ = (float) val;
            break;
          }
        }
      }
      n++;
    }
  }

  return(mri);
}
/*--------------------------------------------------------------*/

int
MRInormalizeSequence(MRI *mri, float target)
{
  int    x, y, z, frame ;
  double norm ;
  double   val ;

  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        for (frame = 0, norm = 0 ; frame < mri->nframes ; frame++)
        {
          MRIsampleVolumeFrame(mri, x, y, z, frame, &val) ;
          norm += (val*val) ;
        }
        norm = sqrt(norm) / target ;
        if (FZERO(norm))
          norm = 1 ;
        for (frame = 0 ; frame < mri->nframes ; frame++)
        {
          switch (mri->type)
          {
          default:
            ErrorReturn
            (ERROR_UNSUPPORTED,
             (ERROR_UNSUPPORTED,
              "MRInormalizeSequence: unsupported input type %d",
              mri->type)) ;
            break ;
          case MRI_SHORT:
            MRISseq_vox(mri, x, y, z, frame) =
              MRISseq_vox(mri, x, y, z, frame) / norm ;
            break ;
          case MRI_FLOAT:
            MRIFseq_vox(mri, x, y, z, frame) /= norm ;
            break ;
          case MRI_UCHAR:
            MRIseq_vox(mri, x, y, z, frame) /= norm ;
            break ;
          }
        }
      }
    }
  }

  return(NO_ERROR) ;
}

double 
MRIstdInLabel(MRI *mri_src, MRI *mri_labeled, MRI *mri_mean, int label)
{
  int  x, y, z, nvox, l ;
  double mean, var = 0.0 ;
  float  val ;

  nvox = 0 ;
  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        l = nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) ;
        if (l == label)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          mean = MRIgetVoxVal(mri_mean, x, y, z, 0) ;
          val -= mean ;
          var += (val*val) ;
          nvox++ ;
        }
      }
    }
  }
  if (!nvox)
    nvox = 1 ;
  return(sqrt(var/nvox)) ;
}

double
MRImeanInLabel(MRI *mri_src, MRI *mri_labeled, int label)
{
  int  x, y, z, nvox, l ;
  double mean = 0.0 ;
  float  val ;

  nvox = 0 ;
  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        l = nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) ;
        if (l == label)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          mean += val ;
          nvox++ ;
        }
      }
    }
  }
  if (!nvox)
    nvox = 1 ;
  return(mean/nvox) ;
}

double
MRImeanInLabelInRegion(MRI *mri_src, MRI *mri_labeled,
                       int label, int x0, int y0, int z0, int whalf)
{
  int  x, y, z, nvox, l ;
  double mean = 0.0 ;
  float  val ;

  nvox = 0 ;
  for (x = x0-whalf ; x <= x0+whalf ; x++)
  {
    if (x < 0 || x >= mri_src->width)
      continue ;
    for (y = y0-whalf ; y <= y0+whalf ; y++)
    {
      if (y < 0 || y >= mri_src->height)
        continue ;
      for (z = z0-whalf ; z <= z0+whalf ; z++)
      {
        if (z < 0 || z >= mri_src->depth)
          continue ;
        l = nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) ;
        if (l == label)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          mean += val ;
          nvox++ ;
        }
      }
    }
  }
  if (!nvox)
    nvox = 1 ;
  return(mean/nvox) ;
}
double
MRImaxInLabelInRegion(MRI *mri_src, MRI *mri_labeled,
                      int label, int x0, int y0, int z0, int whalf)
{
  int  x, y, z;
  int l;
  float  val ;
  int nvox = 0;
  double max = 0.0 ;

  for (x = x0-whalf ; x <= x0+whalf ; x++)
  {
    if (x < 0 || x >= mri_src->width)
      continue ;
    for (y = y0-whalf ; y <= y0+whalf ; y++)
    {
      if (y < 0 || y >= mri_src->height)
        continue ;
      for (z = z0-whalf ; z <= z0+whalf ; z++)
      {
        if (z < 0 || z >= mri_src->depth)
          continue ;
        l = nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) ;
        if (l == label)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
          if (val > max)
            max = val ;
          nvox++ ;
        }
      }
    }
  }
  if (!nvox)
    nvox = 1 ;
  return(max) ;
}

MRI *
MRImakePositive(MRI *mri_src, MRI *mri_dst)
{
  float   fmin, fmax, val ;
  int     x, y, z,f  ;

  MRIvalRange(mri_src, &fmin, &fmax) ;
  mri_dst = MRIcopy(mri_src, mri_dst) ;
  if (fmin >= 0)
    return(mri_dst) ;

  for (f = 0 ; f < mri_dst->nframes ; f++)
  {
    for (x = 0 ; x < mri_dst->width ; x++)
    {
      for (y = 0 ; y < mri_dst->height ; y++)
      {
        for (z = 0 ; z < mri_dst->depth ; z++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          val -= fmin ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val) ;
        }
      }
    }
  }

  return(mri_dst) ;
}
MRI *
MRIeraseNegative(MRI *mri_src, MRI *mri_dst)
{
  int  x, y, z ;
  float  val ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        if (val < 0)
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
      }
    }
  }
  return(mri_dst) ;
}
int
MRIsampleVolumeSlice
(MRI *mri, double x, double y, double z, double *pval, int slice_direction)
{
  int  OutOfBounds;
  int  xm, xp, ym, yp, zm, zp, width, height, depth ;
  double val, xmd, ymd, xpd, ypd ;  /* d's are distances */
  double val11, val12, val21, val22 ;

  if (FEQUAL((int)x,x) && FEQUAL((int)y,y) && FEQUAL((int)z, z))
    return(MRIsampleVolumeType(mri, x, y, z, pval, SAMPLE_NEAREST)) ;

  OutOfBounds = MRIindexNotInVolume(mri, x, y, z);
  if (OutOfBounds == 1)
  {
    /* unambiguously out of bounds */
    *pval = mri->outside_val;
    return(NO_ERROR) ;
  }

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  if (x >= width)    x = width - 1.0 ;
  if (y >= height)   y = height - 1.0 ;
  if (z >= depth)    z = depth - 1.0 ;
  if (x < 0.0)       x = 0.0 ;
  if (y < 0.0)       y = 0.0 ;
  if (z < 0.0)       z = 0.0 ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;

  switch (slice_direction)
  {
  case MRI_CORONAL:
    zp = nint(z) ;
    val11 = MRIgetVoxVal(mri, xm, ym, zp, 0) ;
    val12 = MRIgetVoxVal(mri, xm, yp, zp, 0) ;
    val21 = MRIgetVoxVal(mri, xp, ym, zp, 0) ;
    val22 = MRIgetVoxVal(mri, xp, yp, zp, 0) ;
    *pval = val = xpd * ypd * val11 + xpd * ymd * val12 + xmd * ypd * val21 +
                  xmd * ymd * val22 ;
    break ;
  default:
    ErrorReturn
    (ERROR_UNSUPPORTED,
     (ERROR_UNSUPPORTED,
      "MRIsampleVolumeSlice: unsupported direction %d", slice_direction)) ;
    break ;
  }
  return(NO_ERROR) ;
}


#define MRI_VOX_LABEL_PARTIAL_VOLUME_OUTPUT 0

float MRIvoxelsInLabelWithPartialVolumeEffects( const MRI *mri,
						const MRI *mri_vals, 
						const int label, 
						MRI *mri_mixing_coef, 
						MRI *mri_nbr_labels ) {
  enum { maxlabels = 20000 };
  float volume;
#if (defined(FS_CUDA) && defined(GCAMORPH_ON_GPU))
  volume = MRIvoxelsInLabelWithPartialVolumeEffectsGPU( mri,
							mri_vals,
							label,
							mri_mixing_coef,
							mri_nbr_labels );
  
#else
  int     x, y, z;
  MRI     *mri_border ;
  // DNG 6/7/07 : had to use maxlabels instead of MAX_CMA_LABELS here
  // so that segmentations with values > MAX_CMA_LABELS can be
  // accessed. This includes the cortical segmentations as well as
  // white matter segs. Currently, the max seg no is 4181, but this
  // could easily change.
  // NJS 2/17/10 : Indeed, it did change... the Destrieux a2009s atlas has
  // label values up to about 15000.

  if(label >= maxlabels){
    printf("ERROR: MRIvoxelsInLabelWithPartialVolumeEffects()\n");
    printf(" label %d exceeds maximum label number %d\n",label,maxlabels);
    return(-100000);
  }

#if MRI_VOX_LABEL_PARTIAL_VOLUME_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( ( nCalls%outputFreq ) == 0 ) {
    char fname[STRLEN];
    const unsigned int nOut = nCalls / outputFreq;

    snprintf( fname, STRLEN-1, "MRIvoxelsInLabelWithPartialVolumeEffectsInput%04u.mgz", nOut );
    fname[STRLEN-1] = '\0';
    MRIwrite( (MRI*)mri, fname );

    snprintf( fname, STRLEN-1, "MRIvoxelsInLabelWithPartialVolumeEffectsVals%04u.mgz", nOut );
    fname[STRLEN-1] = '\0';
    MRIwrite( (MRI*)mri_vals, fname );
  }
  nCalls++; 
#endif

  const float vox_vol = mri->xsize*mri->ysize*mri->zsize ;

  /* first find border voxels */
  mri_border = MRImarkLabelBorderVoxels(mri, NULL, label, 1, 1);

  if( DIAG_VERBOSE_ON && (Gdiag & DIAG_WRITE) ) {
    MRIwrite(mri_border, "b.mgz");
  }

  volume = 0 ;
  for( x = 0 ; x < mri->width ; x++ ) {
    for( y = 0 ; y < mri->height ; y++ ) {
      for( z = 0 ; z < mri->depth ; z++ ) {

        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
	}

        const int vox_label = MRIgetVoxVal(mri, x, y, z, 0);
        const int border = MRIgetVoxVal(mri_border, x, y, z, 0);

	/*
	  Note that these are all zeroed at the start
	  of MRIcomputeLabelNbhd
	*/
	int nbr_label_counts[maxlabels], label_counts[maxlabels];
	float label_means[maxlabels];

        if( (vox_label != label) && (border == 0) ) {
          continue;
	}

        if( border == 0 ) {
	  volume += vox_vol;
        } else { /* compute partial volume */
          MRIcomputeLabelNbhd( mri, mri_vals, x, y, z,
			       nbr_label_counts, label_means, 1, maxlabels);

          MRIcomputeLabelNbhd( mri, mri_vals, x, y, z,
			       label_counts, label_means, 7, maxlabels );

	  // Compute partial volume based on intensity
          const float val = MRIgetVoxVal( mri_vals, x, y, z, 0 );

	  

          float mean_label = label_means[vox_label];
          int nbr_label = -1 ;
          int max_count = 0 ;
	  float pv, mean_nbr;

          /* look for a label that is a nbr and is
             on the other side of val from the label mean */
	  int this_label;

          for( this_label = 0 ; this_label < maxlabels;  this_label++ ) {

            if( this_label == vox_label ) {
              continue ;
	    }

            if( nbr_label_counts[this_label] == 0 ) {
              continue ; /* not a nbr */
	    }

            if( (label_counts[this_label] > max_count) &&
                ((label_means[this_label] - val) *
                 (mean_label - val) < 0) ) {
              //max_count = label_means[this_label] ; // bug
              max_count = label_counts[this_label] ;
              nbr_label = this_label ;
            }
          }

          if( vox_label != label && nbr_label != label ) {
            continue; // this struct not in voxel 
	  }

          if( max_count == 0 ) {
            volume += vox_vol ; // couldn't find an appropriate label

            if (mri_nbr_labels) {
	      // find max nbr label anyway for caller
              for( this_label = 0; this_label < maxlabels;  this_label++ ) {

                if( this_label == vox_label ) {
                  continue;
		}

                if( nbr_label_counts[this_label] == 0 ) {
                  continue ; /* not a nbr */
		}
                
                if( label_counts[this_label] > max_count ) {
                  //max_count = label_means[this_label] ; //bug
                  max_count = label_counts[this_label] ;
                  nbr_label = this_label ;
                }
              }

              MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, nbr_label) ;
              if( mri_mixing_coef ) {
                MRIsetVoxVal(mri_mixing_coef, x, y, z, 0, 1.0);
	      }

            }

          } else {
	    // compute partial volume pct 
            mean_nbr = label_means[nbr_label] ;
            pv = (val - mean_nbr) / (mean_label - mean_nbr) ;

            if (pv < 0 || pv > 1) {
              DiagBreak() ;
	    }

            if (pv > 1) {
              pv = 1 ;
	    }

            if (pv < 0) {
              continue ;  // shouldn't happen
	    }

            if( vox_label != label ) {
              pv = 1-pv ;
	    }

            volume += vox_vol * pv ;

            if( mri_mixing_coef ) {
              MRIsetVoxVal(mri_mixing_coef, x, y, z, 0, pv);
	    }

            if (mri_nbr_labels) {
	      // return nbr label to caller
              if (vox_label != label) {
                MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, vox_label);
              } else {
                MRIsetVoxVal(mri_nbr_labels, x, y, z, 0, nbr_label);
	      }
            }
          }
        }
      }
    }
  }
  
  MRIfree(&mri_border) ;
#endif

  return(volume) ;
}

MRI *
MRImakeDensityMap(MRI *mri,
                  MRI *mri_vals,
                  int label,
                  MRI *mri_dst,
                  float orig_res)
{
  float   vox_vol, volume, current_res ;
  int     x, y, z, nbr_label_counts[MAX_CMA_LABELS], ndilates;
  int     label_counts[MAX_CMA_LABELS], this_label, border;
  int     nbr_label, max_count, vox_label ;
  MRI     *mri_border ;
  float   label_means[MAX_CMA_LABELS], pv, mean_label, mean_nbr, val ;

  memset(nbr_label_counts,0,sizeof(nbr_label_counts));

  if (mri_dst == NULL)
  {
    mri_dst = MRIalloc(mri->width, mri->height, mri->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri, mri_dst) ;
  }

  /* first find border voxels */
  mri_border = MRImarkLabelBorderVoxels(mri, NULL, label, 1, 1) ;
  current_res = mri->xsize ;
  ndilates = 0 ;
  while (current_res < orig_res)
  {
    //    MRIdilateLabel(mri_border, mri_border, 1, 1) ;
    MRIdilate(mri_border, mri_border) ;
    //    MRImask(mri_border, mri, mri_border, 0, 0) ; // turn off exterior
    current_res += mri->xsize ;
    ndilates++ ;
  }
  if (DIAG_VERBOSE_ON && (Gdiag & DIAG_WRITE))
    MRIwrite(mri_border, "b.mgz") ;
  vox_vol = mri->xsize*mri->ysize*mri->zsize ;
  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        vox_label = MRIgetVoxVal(mri, x, y, z, 0) ;
        border = MRIgetVoxVal(mri_border, x, y, z, 0) ;
        if ((vox_label != label) && (border == 0))
          continue ;

        volume = vox_vol ;
        if (border == 0)
          volume = vox_vol ;
        else  /* compute partial volume */
        {
#if 0
          MRIcomputeLabelNbhd
          (mri, mri_vals, x, y, z,
           nbr_label_counts, label_means, 1, MAX_CMA_LABELS) ;
#endif
          MRIcomputeLabelNbhd
          (mri, mri_vals, x, y, z,
           label_counts, label_means, 7+(ndilates), MAX_CMA_LABELS) ;
          val =
            MRIgetVoxVal(mri_vals, x, y, z, 0) ;  /* compute partial
                                                     volume based on
                                                     intensity */
          mean_label = label_means[vox_label] ;
          nbr_label = -1 ;
          max_count = 0 ;
          /* look for a label that is a nbr and is
             on the other side of val from the label mean */
          for (this_label = 0 ;
               this_label < MAX_CMA_LABELS ;
               this_label++)
          {
            if (this_label == vox_label)
              continue ;
            if (nbr_label_counts[this_label] == 0)   /* not a nbr */
              continue ;

            if ((label_counts[this_label] > max_count) &&
                ((label_means[this_label] - val) *
                 (mean_label - val) < 0))
            {
              max_count = label_means[this_label] ;
              nbr_label = this_label ;
            }
          }
          if (vox_label != label &&
              nbr_label != label)  /* this struct not in voxel */
            continue ;

          if (max_count > 0)  /* compute partial volume pct */
          {
            mean_nbr = label_means[nbr_label] ;
            pv = (val - mean_nbr) / (mean_label - mean_nbr) ;
            if (vox_label == label)
              volume = pv*vox_vol ;
            else
              volume = vox_vol * (1-pv) ;
            if (pv < 0 || pv > 1)
              DiagBreak() ;
          }
        }
        MRIsetVoxVal(mri_dst, x, y, z, 0, volume) ;
      }
    }
  }

  MRIfree(&mri_border) ;
  return(mri_dst) ;
}



#define MRI_MARK_LABEL_BORDER_VOXELS_OUTPUT 0

MRI *
MRImarkLabelBorderVoxels( const MRI *mri_src,
			  MRI *mri_dst,
			  int label,
			  int mark,
			  int six_connected ) {

#if MRI_MARK_LABEL_BORDER_VOXELS_OUTPUT
  const unsigned int outputFreq = 10;
  static unsigned int nCalls = 0;
  if( ( nCalls%outputFreq ) == 0 ) {
    char fname[STRLEN];
    const unsigned int nOut = nCalls / outputFreq;
    snprintf( fname, STRLEN-1, "mriMarkLabelBorderVoxelsInput%04u.mgz", nOut );
    fname[STRLEN-1] = '\0';
    MRIwrite( (MRI*)mri_src, fname );
  }
  nCalls++; 
#endif
  if (mri_dst == NULL) {
    mri_dst = MRIclone(mri_src, NULL) ;
  }

#if (defined(FS_CUDA) && defined(GCAMORPH_ON_GPU))
  MRImarkLabelBorderVoxelsGPU( mri_src, mri_dst, label, mark, six_connected );
#else
  int  x, y, z, xk, yk, zk, xi, yi, zi;
  int this_label, that_label, border;
  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        this_label = MRIgetVoxVal(mri_src, x, y, z, 0) ;
        border = 0 ;
        for (xk = -1 ; xk <= 1 && !border ; xk++)
        {
          xi = mri_src->xi[x+xk] ;
          for (yk = -1 ; yk <= 1 && !border ; yk++)
          {
            yi = mri_src->yi[y+yk] ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              if (six_connected && (abs(xk)+abs(yk)+abs(zk) != 1))
                continue ;
              zi = mri_src->zi[z+zk] ;
              that_label = MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
              if (((this_label == label) &&
                   (that_label != label)) ||
                  ((this_label != label) &&
                   (that_label == label)))
              {
                border = 1 ;
                break ;
              }
            }
          }
        }
        if (border)
          MRIsetVoxVal(mri_dst, x, y, z, 0, mark) ;
      }
    }
  }
#endif

  return(mri_dst) ;
}


int
MRIcomputeLabelNbhd( const MRI *mri_labels,
		     const MRI *mri_vals,
		     const int x, const int y, const int z,
		     int *label_counts, float *label_means,
		     const int whalf, const int max_labels ) {
  int    xi, yi, zi, xk, yk, zk, label ;
  float  val ;

  // Zero the input arrays
  memset( label_counts, 0, sizeof(label_counts[0])*max_labels );
  if( label_means ) {
    memset( label_means, 0, sizeof(label_means[0])*max_labels );
  }

  for (xk = -whalf ; xk <= whalf ; xk++) {
    xi = mri_labels->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++) {
      yi = mri_labels->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++) {
        zi = mri_labels->zi[z+zk] ;
        label = MRIgetVoxVal(mri_labels, xi, yi, zi, 0) ;
        label_counts[label]++ ;

        if( mri_vals ) {
          val = MRIgetVoxVal(mri_vals, xi, yi, zi, 0) ;
          label_means[label] += val ;
        }
      }
    }
  }

  if( mri_vals ) {
    for (label = 0 ; label < max_labels ; label++) {
      if (label_counts[label] > 0) {
        label_means[label] /= label_counts[label] ;
      }
    }
  }
  return(NO_ERROR) ;
}

/**
 * void MRIcalCRASforSampledVolume
 *
 * @param src MRI* volume
 * @param dst MRI* sampled volume
 * @param pr  output ptr to c_r
 * @param pa  output ptr to c_a
 * @param ps  output ptr to c_s
 */
void MRIcalcCRASforSampledVolume
(MRI *src, MRI *dst, double *pr, double *pa, double *ps)
{
  // get the voxel position of the "center" voxel of the dst in the src volume
  // i.e. sample is 2, then get the voxel position 64 in the src volume
  //      thus  it is 128.5 (in the src) for 64 (in the dst)
  double sx, sy, sz;
  int dx, dy, dz;
  int samplex, sampley, samplez;

#if 0
  // error check first
  if ((src->width)%(dst->width) != 0)
    ErrorExit
    (ERROR_BADPARM, "src width must be integer multiple of dst width");
  if ((src->height)%(dst->height) != 0)
    ErrorExit
    (ERROR_BADPARM, "src height must be integer multiple of dst height");
  if ((src->depth)%(dst->depth) != 0)
    ErrorExit
    (ERROR_BADPARM, "src depth must be integer multiple of dst depth");
#endif

  samplex = src->width/dst->width;
  sampley = src->height/dst->height;
  samplez = src->depth/dst->depth;

  // "center" voxel position in dst
  dx = dst->width/2; // if the length is odd,
  // then it does the right thing (truncation)
  dy = dst->height/2;// i.e.  0 1 2 then 3/2 = 1,  0 1 2 3 then 4/2=2
  dz = dst->depth/2;
  // corresponding position in src
  sx = dx*samplex + (samplex-1.)/2.; // ((dx*samplex - 0.5) +
  // ((dx+1)*samplex -0.5))/2.
  sy = dy*sampley + (sampley-1.)/2.;
  sz = dz*samplez + (samplez-1.)/2.;
  //
  // Example
  // | 0 1 2 | 3 4 5 | 6 7 8 | -(sample by 3).
  // |   0   |   1   |   2   |
  //                                       1 (3/2) in dst corresponds
  //                                       to 1*3+(3-1)/2 = 4! in src
  //
  // | 0 1 2 3 | 4 5 6 7 |     -(sample by 4)->0 1
  // |    0    |    1    |                        1 (2/2) in dst corresponds
  //                                              to 1*4+(4-1)/2 = 5.5 in src

  // get ras of the "center" voxel position in dst
  if (!src->i_to_r__)
  {
    MATRIX *tmp = extract_i_to_r( src );
    AffineMatrixAlloc( &(src->i_to_r__ ) );
    SetAffineMatrix( src->i_to_r__, tmp );
    MatrixFree( &tmp );

    src->r_to_i__ = extract_r_to_i(src);
  }

  AffineVector s, p;
  SetAffineVector( &s, sx, sy, sz );
  
  AffineMV( &p, src->i_to_r__, &s );

  float pxf, pyf, pzf;
  GetAffineVector( &p, &pxf, &pyf, &pzf );
  *pr = pxf;
  *pa = pyf;
  *ps = pzf;



  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "c_ras for sample volume is (%f, %f, %f) "
            "compared with the src (%f, %f, %f)\n",
            *pr, *pa, *ps, src->c_r, src->c_a, src->c_s);
}

/**
 * MRIcalcCRASfroExtractedVolume
 *
 * @param src  MRI* src volume
 * @param x0   src start position of the extraction region
 * @param y0
 * @param z0
 * @param x1   target start position of the extracted region
 * @param y1
 * @param z1
 * @param pr   output double*  c_r
 * @param pa                 c_a
 * @param ps                 c_s
 */
void MRIcalcCRASforExtractedVolume
(MRI *src, MRI *dst, int x0, int y0, int z0, int x1, int y1, int z1,
 double *pr, double *pa, double *ps)
{
  double cx, cy, cz;
  // The "center" voxel position of the
  // extracted volume in the original voxel position
  //        x1 of dst corresponds to x0 of src
  // Thus, the "center" of dst corresponds to that of the src is
  cx = x0+(double)dst->width/2.0 - x1;  // integer divide cutoff extra
  cy = y0+(double)dst->height/2.0 - y1;
  cz = z0+(double)dst->depth/2.0 - z1;

  if (!src->i_to_r__)
  {
    MATRIX *tmp;
    tmp = extract_i_to_r(src);
    AffineMatrixAlloc( &(src->i_to_r__) );
    SetAffineMatrix( src->i_to_r__, tmp );
    MatrixFree( &tmp );

    src->r_to_i__ = extract_r_to_i(src);
  }

  AffineVector c, p;
  float pxf, pyf,pzf;

  SetAffineVector( &c, cx, cy,cz );
  AffineMV( &p, src->i_to_r__, &c );

  GetAffineVector( &p, &pxf, &pyf, &pzf );
  *pr = pxf;
  *pa = pyf;
  *ps = pzf;

  // got where the RAS position of the new volume position
  // we have to translate so that we can get the same value
  // under the new volume

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "c_ras for sample volume is (%f, %f, %f) "
            "compared with the src (%f, %f, %f)\n",
            *pr, *pa, *ps, src->c_r, src->c_a, src->c_s);
}

// transform is the hires to lowres vox-to-vox transform
void MRIcalcCRASforHiresVolume(MRI *hires, MRI *lowres,
                               MATRIX *vox_xform,
                               double *pr, double *pa, double *ps)
{
  // get where the center of hires volume goes to in the lowres volume
  double cx, cy, cz;
  double dx, dy, dz;
  cx = (double)hires->width/2.0;
  cy = (double)hires->height/2.0;
  cz = (double)hires->depth/2.0;
  TransformWithMatrix(vox_xform, cx, cy, cz, &dx, &dy, &dz);
  // get the c_ras values for this position
  AffineVector d, p;
  float prf, paf, psf;
  SetAffineVector( &d, dx, dy, dz );
  AffineMV( &p, lowres->i_to_r__, &d );
  GetAffineVector( &p, &prf, &paf, &psf );
  *pr = prf;
  *pa = paf;
  *ps = psf;

  // if we use this c_ras value for the transformed hires volume, then
  // the volume will be containing the original points
}

// transform is src to dst vox to vox transform
// return rotated src whose center is at the right place
// Just rotating the original volume make the non-zero voxels go
// outside of the rotated volume.  This routine will keep the
// center of the rotated volume at the right location.
MRI *
MRIsrcTransformedCentered
(MRI *src, MRI *dst, MATRIX *stod_voxtovox, int interp_method)
{
  double cr, ca, cs;
  MRI *rotated;
  MATRIX *stosrotVox;
  MATRIX *tmp;

  // get the c_ras value for the rotated volume
  MRIcalcCRASforHiresVolume(src, dst, stod_voxtovox, &cr, &ca, &cs);
  rotated = MRIclone(src, NULL);
  // reset the rotated volume
  rotated->c_r = cr;
  rotated->c_a = ca;
  rotated->c_s = cs;
  MRIreInitCache(rotated); // recalculate the transform
  // the map is the following
  //
  //     src -->    RAS
  //       |         |
  //       |given    |
  //       V         V
  //     dst -->    RAS
  //       |         |
  //       |       (identity)
  //       V         |
  //     src'  -->  RAS
  // new rotated src volume whose center is at the right location
  //
  MATRIX *tmp2 = MatrixAlloc( 4, 4, MATRIX_REAL );
  GetAffineMatrix( tmp2, dst->i_to_r__ );
  tmp = MatrixMultiply(rotated->r_to_i__, tmp2, NULL);
  MatrixFree( &tmp2 );

  stosrotVox = MatrixMultiply(tmp, stod_voxtovox, NULL);
  MRIlinearTransformInterp(src, rotated, stosrotVox, interp_method);
  return rotated;
}

////////////////////////////////////////////////////////////////////////////
// No sample method ////////////////////////////////////////////////////////
// using the src and orig_dst to modify
// the direction cosines and c_(ras) value
// so that it will be rotated in the RAS space but no pixel is sampled
MRI *MRITransformedCenteredMatrix(MRI *src, MRI *orig_dst, MATRIX *m_L)
{
  LTA *lta ;
  MRI *mri_dst ;

  lta = LTAalloc(1, NULL) ;
  MatrixCopy(m_L, lta->xforms[0].m_L) ;
  lta->type = LINEAR_VOX_TO_VOX ;
  mri_dst = MRITransformedCentered(src, orig_dst, lta) ;
  LTAfree(&lta) ;
  return(mri_dst) ;
}
MRI *MRITransformedCentered(MRI *src, MRI *orig_dst, LTA *lta)
{
  MRI *dst = 0;
  double cx, cy, cz;
  double cr, ca, cs;
  MATRIX *dstToRas = 0;
  MATRIX *SI = 0;
  MATRIX *D = 0;
  LT *tran = &lta->xforms[0];

  dst = MRIcopy(src, NULL);

  if (lta->num_xforms > 1)
    ErrorExit(ERROR_BADPARM, "The LTA contains more than one transforms.");

  // first verify the consistency of the transform stored geometry vs. argument
  if (tran->dst.valid == 1)
  {
    // compare the one with orig_dst to the stored one
    VG dstvg, storedvg;
    getVolGeom(orig_dst, &dstvg);
    storedvg.valid = tran->dst.valid;
    storedvg.width = tran->dst.width;
    storedvg.height = tran->dst.height;
    storedvg.depth = tran->dst.depth;
    storedvg.xsize = tran->dst.xsize;
    storedvg.ysize = tran->dst.ysize;
    storedvg.zsize = tran->dst.zsize;
    storedvg.x_r = tran->dst.x_r;
    storedvg.x_a = tran->dst.x_a;
    storedvg.x_s = tran->dst.x_s;
    storedvg.y_r = tran->dst.y_r;
    storedvg.y_a = tran->dst.y_a;
    storedvg.y_s = tran->dst.y_s;
    storedvg.z_r = tran->dst.z_r;
    storedvg.z_a = tran->dst.z_a;
    storedvg.z_s = tran->dst.z_s;
    storedvg.c_r = tran->dst.c_r;
    storedvg.c_a = tran->dst.c_a;
    storedvg.c_s = tran->dst.c_s;
    if (!vg_isEqual(&dstvg, &storedvg))
    {
      fprintf
      (stderr,
       "WARNING: stored destination volume for the "
       "transform differs from the argument.");
    }
  }
  if (tran->src.valid == 1)
  {
    // compare the one with orig_dst to the stored one
    VG srcvg, storedvg;
    getVolGeom(src, &srcvg);
    storedvg.valid = tran->src.valid;
    storedvg.width = tran->src.width;
    storedvg.height = tran->src.height;
    storedvg.depth = tran->src.depth;
    storedvg.xsize = tran->src.xsize;
    storedvg.ysize = tran->src.ysize;
    storedvg.zsize = tran->src.zsize;
    storedvg.x_r = tran->src.x_r;
    storedvg.x_a = tran->src.x_a;
    storedvg.x_s = tran->src.x_s;
    storedvg.y_r = tran->src.y_r;
    storedvg.y_a = tran->src.y_a;
    storedvg.y_s = tran->src.y_s;
    storedvg.z_r = tran->src.z_r;
    storedvg.z_a = tran->src.z_a;
    storedvg.z_s = tran->src.z_s;
    storedvg.c_r = tran->src.c_r;
    storedvg.c_a = tran->src.c_a;
    storedvg.c_s = tran->src.c_s;
    if (!vg_isEqual(&srcvg, &storedvg))
    {
      fprintf(stderr, "WARNING: stored destination volume for the "
              "transform differes from the argument.");
    }
  }

  if (lta->type == LINEAR_RAS_TO_RAS)
  {
    MATRIX *tmp;
    //
    //      src   -->   RAS
    //       |           |
    //       | M         | Y
    //       V           V
    //    orig_dst -->  RAS
    //
    // transform it to the vox-to-vox
    // now calculate M
    MATRIX *tmp2 = MatrixAlloc( 4, 4, MATRIX_REAL );
    GetAffineMatrix( tmp2, src->i_to_r__ );
    tmp = MatrixMultiply(tran->m_L, tmp2, NULL);
    MatrixFree( &tmp2 );

    MatrixFree(&tran->m_L);
    tran->m_L = MatrixMultiply(orig_dst->r_to_i__, tmp, NULL);
    MatrixFree(&tmp);
    lta->type = LINEAR_VOX_TO_VOX;
  }

  if (lta->type != LINEAR_VOX_TO_VOX)
    ErrorExit
    (ERROR_BADPARM,
     "The LTA does not contain LINEASR_RAS_TO_RAS nor LINEAR_VOX_TO_VOX.");
  //
  //      src   -->   RAS
  //       |           |
  //       | M         | Y
  //       V           V
  //    orig_dst -->  RAS
  //       |           |
  //       |       (identity)
  //       V           |
  //      dst   -->   RAS
  // new rotated src volume whose center is at the right location
  //            ?
  // we try src->dst is identity (no sampling)
  //
  //    Y = i_to_r(orig_dst) * M * r_to_i(src)
  //
  // Thus we have simpler picture
  //
  //      src   -----> RAS
  //       |            |
  //    (identity)     (Y)
  //       |            |
  //       V            V
  //      dst   -----> RAS
  //            (???)
  //
  // Thus
  //      (???) = Y * i_to_r(src)
  //            = i_to_r(orig_dst) * M * r_to_i(src) * i_to_r(src)
  //
  //      (???) = i_to_r(orig_dst) * M           (A)
  //
  // Thus we derive direction cosines and c_(ras) from ???
  //
  //      (???)  =  ( X   T ) ( S  0 )    where S is the voxel size
  //                ( 0   1 ) ( 0  1 )
  // or
  //      ( X  T ) = (???) * (S^(-1) 0 )         (B)
  //      ( 0  1 )           ( 0     1 )
  //
  //
  MATRIX *tmp2 = MatrixAlloc( 4, 4, MATRIX_REAL );
  GetAffineMatrix( tmp2, orig_dst->i_to_r__ );
  dstToRas = MatrixMultiply( tmp2, tran->m_L, NULL);
  MatrixFree( &tmp2 );

  SI = MatrixAlloc(4, 4, MATRIX_REAL);
  *MATRIX_RELT(SI, 1, 1) = 1./dst->xsize ;
  *MATRIX_RELT(SI, 2, 2) = 1./dst->ysize ;
  *MATRIX_RELT(SI, 3, 3) = 1./dst->zsize ;
  *MATRIX_RELT(SI, 4, 4) = 1. ;

  D = MatrixMultiply(dstToRas, SI, NULL);

  dst->x_r = *MATRIX_RELT(D, 1, 1);
  dst->x_a = *MATRIX_RELT(D, 2, 1);
  dst->x_s = *MATRIX_RELT(D, 3, 1);
  dst->y_r = *MATRIX_RELT(D, 1, 2);
  dst->y_a = *MATRIX_RELT(D, 2, 2);
  dst->y_s = *MATRIX_RELT(D, 3, 2);
  dst->z_r = *MATRIX_RELT(D, 1, 3);
  dst->z_a = *MATRIX_RELT(D, 2, 3);
  dst->z_s = *MATRIX_RELT(D, 3, 3);

  MatrixFree(&D);
  MatrixFree(&SI);

  // c_ras is calculated by
  //
  //      ( X  T )(S 0)(dstwidth/2 )  = (c_r)
  // or  i_to_r(orig_dst) * M * (dstwidth/2 ) = C'
  //      ( 0  1 )(0 1)(dstheight/2)    (c_a)
  //                               (dstheight/2)
  //                   (dstdepth/2 )    (c_s)
  //                                                  (dstdepth/2 )
  //                   (    1      )    ( 1 )               (    1      )
  //
  // This is the same as the original center point
  // is mapped to orig_dst and then obtaining
  // its ras position! (because src and dst has the
  // same volume size).  The only difference
  // is its orientation with respect to orig_dst RAS.  Good
  //
  // get where the center of hires volume goes to in the lowres volume
  cx = (double)dst->width/2.0;
  cy = (double)dst->height/2.0;
  cz = (double)dst->depth/2.0;
  // get the c_ras values for this position
  TransformWithMatrix(dstToRas, cx, cy, cz, &cr, &ca, &cs);
  // if we use this c_ras value for the transformed hires volume, then
  // the volume will be containing the original points
  dst->c_r = cr;
  dst->c_a = ca;
  dst->c_s = cs;

  // when you change the direction cosine, you have to do this
  MRIreInitCache(dst);

  MatrixFree(&dstToRas);

  return dst;
}
int
MRIcropBoundingBox(MRI *mri,    MRI_REGION  *box)
{
  box->x = MAX(0, box->x) ;
  box->y = MAX(0, box->y) ;
  box->z = MAX(0, box->z) ;
  box->dx = MIN(mri->width-box->x, box->dx) ;
  box->dy = MIN(mri->height-box->y, box->dy) ;
  box->dz = MIN(mri->depth-box->z, box->dz) ;
  return(NO_ERROR) ;
}

MATRIX *
MRIgetVoxelToVoxelXform(MRI *mri_src, MRI *mri_dst)
{
  MATRIX *m_ras2vox_dst, *m_vox2ras_src, *m_vox2vox ;

  m_vox2ras_src = MRIgetVoxelToRasXform(mri_src) ;
  m_ras2vox_dst = MRIgetRasToVoxelXform(mri_dst) ;
  m_vox2vox = MatrixMultiply(m_ras2vox_dst, m_vox2ras_src, NULL) ;
  MatrixFree(&m_vox2ras_src) ;
  MatrixFree(&m_ras2vox_dst) ;
  return(m_vox2vox) ;
}

/*--------------------------------------------------------------
  MRIfovCol(mri) - computes the edge-to-edge FOV in the column
  direction. fov is in mm.
  -------------------------------------------------------------*/
float MRIfovCol(MRI *mri)
{
  MATRIX *M,*v,*a,*b,*d;
  float fov;

  M = MRIgetVoxelToRasXform(mri) ;
  v = MatrixAlloc(4,1,MATRIX_REAL);
  v->rptr[1][4] = 1;
  v->rptr[1][1] = mri->width-1+0.5; // edge of last column
  a = MatrixMultiply(M,v,NULL);     // xyz of last column
  v->rptr[1][1] = -0.5;             // edge of first column
  b = MatrixMultiply(M,v,NULL);     // xyz of first column
  d = MatrixSubtract(a,b,NULL); // xyz difference
  fov = VectorLen(d);           // fov is in mm
  MatrixFree(&M);
  MatrixFree(&v);
  MatrixFree(&a);
  MatrixFree(&b);
  MatrixFree(&d);

  //printf("MRIfovCol() %g\n",fov);
  return(fov);
}
/* ---------------------------------------------------------------------
   MRIorientationStringToDircos() - sets the direction cosines of to
   be that dictated by the Orientation String. This is helpful for
   setting the direction cosines when the information is not present
   in a header (eg, with FSL analyze format).  The Orientation String
   is a three character string indicating the primary direction of
   each axis in the 3d matrix. The characters can be L,R,A,P,I,S.  The
   string must be valid (see MRIcheckOrientationString). If the string
   is not valid, the errors are printed and a 1 is returned. Eg, of
   valid strings are RAS, LPS, LAI. Invalid are ABC (B and C are not
   valid), RAP (the AP axis is represented twice, IS axis not at all).
   There are 48 possible valid strings. This should only be used to
   get the direction cosines "about right".
   --------------------------------------------------------------------*/

int MRIorientationStringToDircos(MRI *mri, char *ostr)
{
  int c,r=0;
  double Mdc[3][3], v=0;
  char *errstr;

  errstr = MRIcheckOrientationString(ostr);
  if (errstr != NULL)
  {
    printf("ERROR: in orientation string %s\n",ostr);
    printf("%s",errstr);
    free(errstr);
    return(1);
  }

  // Initialize the matrix of direction cosines (Mdc)
  for (c=0; c<3; c++) for (r=0; r<3; r++) Mdc[r][c] = 0;

  // Each column of Mdc corresponds to a different axis which corresonds to
  // a different letter in the orientation string. The value of the letter
  // determine which row of the Mdc will be affected.
  for (c=0; c<3; c++)
  {
    switch (toupper(ostr[c]))
    {
    case 'L':
      r=0;
      v=-1;
      break;
    case 'R':
      r=0;
      v=+1;
      break;
    case 'P':
      r=1;
      v=-1;
      break;
    case 'A':
      r=1;
      v=+1;
      break;
    case 'I':
      r=2;
      v=-1;
      break;
    case 'S':
      r=2;
      v=+1;
      break;
    }
    Mdc[r][c] = v;
  }
  mri->x_r = Mdc[0][0];
  mri->x_a = Mdc[1][0];
  mri->x_s = Mdc[2][0];

  mri->y_r = Mdc[0][1];
  mri->y_a = Mdc[1][1];
  mri->y_s = Mdc[2][1];

  mri->z_r = Mdc[0][2];
  mri->z_a = Mdc[1][2];
  mri->z_s = Mdc[2][2];

  return(0);
}
/* ---------------------------------------------------------------
   MRIcheckOrientationString() - this checks the orientation string
   to make sure that it is valid. "Valid" means that all axes are
   represented exactly once and no invalid characters are present
   in the string. Case-insensitive. Returns NULL if everything is
   ok, otherwise is returns a string that lists all the errors
   it encountered.
   ---------------------------------------------------------------*/
char *MRIcheckOrientationString(char *ostr)
{
  int c, nsag=0, ncor=0, nax=0, err;
  char errstr[1000], *errstrret=NULL;
  char tmpstr[1000];

  errstr[0] = '\0';
  tmpstr[0] = '\0';
  err = 0;

  for (c=0; c<3; c++)
  {
    switch (toupper(ostr[c]))
    {
    case 'L':
      nsag++;
      break;
    case 'R':
      nsag++;
      break;
    case 'P':
      ncor++;
      break;
    case 'A':
      ncor++;
      break;
    case 'I':
      nax++;
      break;
    case 'S':
      nax++;
      break;
    default:
      sprintf(tmpstr,"%s  Character %c in position %d is invalid.\n",
              errstr,ostr[c],c+1);
      strcpy(errstr,tmpstr);
      err = 1;
      break;
    }
  }

  if (nsag > 1)
  {
    strcat(errstr,"  LR axis represented multiple times.\n");
    err = 1;
  }
  if (ncor > 1)
  {
    strcat(errstr,"  PA axis represented multiple times.\n");
    err = 1;
  }
  if (nax > 1)
  {
    strcat(errstr,"  IS axis represented multiple times.\n");
    err = 1;
  }
  if (nsag == 0)
  {
    strcat(errstr,"  LR axis not represented.\n");
    err = 1;
  }
  if (ncor == 0)
  {
    strcat(errstr,"  PA axis not represented.\n");
    err = 1;
  }
  if (nax == 0)
  {
    strcat(errstr,"  IS axis not represented.\n");
    err = 1;
  }

  if (err)
  {
    errstrret = (char *) calloc(sizeof(char),strlen(errstr)+1);
    memmove(errstrret,errstr,strlen(errstr));
  }

  return(errstrret);
}
/*------------------------------------------------------------------
  MRIdircosToOrientationString() - examines the direction cosines and
  creates an Orientation String. The Orientation String is a three
  character string indicating the primary direction of each axis
  in the 3d matrix. The characters can be L,R,A,P,I,S. Case is not
  important, but upper case is used here. If ras_good_flag == 0,
  then ostr = ??? and 1 is returned.
  ------------------------------------------------------------------*/
int MRIdircosToOrientationString(MRI *mri, char *ostr)
{
  int c;
  float Mdc[3][3], sag, cor, ax;

  ostr[3] = 0 ;
  if (! mri->ras_good_flag)
  {
    ostr[0] = '?';
    ostr[1] = '?';
    ostr[2] = '?';
    return(1);
  }

  Mdc[0][0] = mri->x_r;
  Mdc[1][0] = mri->x_a;
  Mdc[2][0] = mri->x_s;

  Mdc[0][1] = mri->y_r;
  Mdc[1][1] = mri->y_a;
  Mdc[2][1] = mri->y_s;

  Mdc[0][2] = mri->z_r;
  Mdc[1][2] = mri->z_a;
  Mdc[2][2] = mri->z_s;

  for (c=0; c<3; c++) ostr[c] = '\0';

  for (c=0; c<3; c++)
  {
    sag = Mdc[0][c]; // LR axis
    cor = Mdc[1][c]; // PA axis
    ax  = Mdc[2][c]; // IS axis
    //printf("c = %d, sag = %g, cor = %g, ax = %g\n",c,sag,cor,ax);
    if (fabs(sag) > fabs(cor) && fabs(sag) > fabs(ax))
    {
      if (sag > 0) ostr[c] = 'R';
      else        ostr[c] = 'L';
      continue;
    }
    if (fabs(cor) > fabs(ax))
    {
      if (cor > 0) ostr[c] = 'A';
      else        ostr[c] = 'P';
      continue;
    }
    if (ax > 0) ostr[c] = 'S';
    else       ostr[c] = 'I';
  }
  return(0);

}
/*-------------------------------------------------------------------
  MRIsliceDirectionName() - returns the name of the primary slice
  orientation base on the direction cosine. If mri->ras_good_flag=0,
  then "unknown" is returned.
  -------------------------------------------------------------------*/
char *MRIsliceDirectionName(MRI *mri)
{
  char ostr[4];
  char *slicedir = NULL;
  char *rtstr;
  int len;

  ostr[3] = '\0';

  MRIdircosToOrientationString(mri, ostr);
  if (toupper(ostr[2]) == 'L'|| toupper(ostr[2]) == 'R')
    slicedir = "sagittal";
  if (toupper(ostr[2]) == 'P'|| toupper(ostr[2]) == 'A')
    slicedir = "coronal";
  if (toupper(ostr[2]) == 'I'|| toupper(ostr[2]) == 'S')
    slicedir = "axial";
  if (! mri->ras_good_flag) slicedir = "unknown";

  len = strlen(slicedir);
  rtstr = (char *) calloc(sizeof(char),len+1);
  memmove(rtstr,slicedir,len);

  return(rtstr);
}

/** 
 * This is deprecated.  Please use MRIextractDistanceMap in fastmarching.h 
 * instead 
 **/
MRI *
MRIdistanceTransform(MRI *mri_src,
                     MRI *mri_dist,
                     int label,
                     float max_dist,
                     int mode,
                     MRI *mri_mask)
{

  const int width = mri_src->width;
  const int height = mri_src->height;
  const int depth = mri_src->depth;

  if (mri_dist == NULL)
  {
    mri_dist = MRIalloc(width, height, depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dist) ;
  }
  else
    MRIclear(mri_dist) ;
    
  // these are the modes in fastmarching...
  const int outside = 1;
  // this one isn't used in this function
  // const int inside = 2;
  const int both_signed = 3;
  const int both_unsigned = 4;
  
  if( mode == DTRANS_MODE_SIGNED ) {
    // DTRANS_MODE_SIGNED is negative inside and positive outside.  This 
    // corresponds to both.
    mode = both_signed;
  } else if(mode == DTRANS_MODE_UNSIGNED ) {
    // DTRANS_MODE_UNSIGNED is positive inside and positive outside
    mode = both_unsigned;
  } else if( mode == DTRANS_MODE_OUTSIDE ) {
    // DTRANS_MODE_OUTSIDE is zero inside and positive outside
    mode = outside;
  }

  // Not to get an error within MRIextractDistanceMap
  if (mri_src->type != MRI_FLOAT)
    {
      MRI *mri_tmp = MRIchangeType(mri_src, MRI_FLOAT, 0, 1, 1) ;
      mri_dist = MRIextractDistanceMap
        ( mri_tmp, mri_dist, label, max_dist, mode, mri_mask);
      MRIfree(&mri_tmp); 
    }
  else
    mri_dist = MRIextractDistanceMap
      ( mri_src, mri_dist, label, max_dist, mode, mri_mask);

  mri_dist->outside_val = max_dist ;  
  MRIscalarMul(mri_dist, mri_dist, mri_src->xsize) ;
  return mri_dist;  
}

/*-------------------------------------------------------------------
  MRIreverseSliceOrder() - reverses the order of the slices WITHOUT
  changing the gemoetry information. This is specially desgined to
  undo the reversal that Siemens sometimes makes to volumes. If
  using FreeSurfer unpacking (ie, DICOMRead.c), it should only need
  to be done to mosaics. Note: cannot be done in-place!
  -------------------------------------------------------------------*/
MRI *MRIreverseSliceOrder(MRI *invol, MRI *outvol)
{
  int c,r,s1,s2,f;
  double val;

  if (invol == outvol)
  {
    printf("ERROR: MRIreverseSliceOrder: cannot be done in-place\n");
    return(NULL);
  }

  outvol = MRIclone(invol,outvol);

  s2 = invol->depth;
  for (s1=0; s1 < invol->depth; s1++)
  {
    s2--;
    for (c=0; c < invol->width; c++)
    {
      for (r=0; r < invol->height; r++)
      {
        for (f=0; f < invol->nframes; f++)
        {
          val = MRIgetVoxVal(invol,c,r,s1,f);
          MRIsetVoxVal(outvol,c,r,s2,f,val);
        }
      }
    }
  }

  return(outvol);
}

#define NDIRS 3
static int dirs[NDIRS][3] =
  {
    {
      0, 0, 1
    },
    { 0, 1, 0},
    { 1, 0, 0}
  } ;

MRI *
MRInonMaxSuppress(MRI *mri_src, MRI *mri_sup, float thresh, int thresh_dir)
{
  int  x, y, z, i, max_i ;
  double val, dx, dy, dz, Ix, Iy, Iz;
  double dot, max_dot, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
  double p3_x, p3_y, p3_z, nx, ny, nz, u, xf, yf, zf, oval, numer, denom ;

  mri_sup = MRIcopy(mri_src, mri_sup) ;

  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri_src, x, y, z,0) ;
        val = val*thresh_dir ;
        if (val <= thresh)
        {
          MRIsetVoxVal(mri_sup, x, y, z, 0, 0.0) ;
          continue ;
        }
        MRIsampleVolumeGradient(mri_src, x, y, z,  &Ix, &Iy, &Iz) ;
        dx = dirs[0][0] ;
        dy = dirs[0][1] ;
        dz = dirs[0][2] ;
        max_i = 0 ;
        max_dot = fabs(dx*Ix+dy*Iy+dz*Iz) ;
        for (i = 1 ; i < NDIRS ; i++)
        {
          // compute intersection of gradient
          // direction with coordinate planes
          dx = dirs[i][0] ;
          dy = dirs[i][1] ;
          dz = dirs[i][2] ;
          dot = fabs(dx*Ix+dy*Iy+dz*Iz) ;
          if (dot > max_dot)
          {
            max_i = i ;
            max_dot = dot ;
          }
        }
        // normal to coordinate plane -
        // compute intersections of gradient dir with plane
        nx = dirs[max_i][0] ;
        ny = dirs[max_i][1] ;
        nz = dirs[max_i][2] ;
        p1_x = x ;
        p1_y = y ;
        p1_z = z ;    // point on line
        p2_x = x+Ix ;
        p2_y = y+Iy ;
        p2_z = z+Iz ; // 2nd point on
        //                                           gradient line
        p3_x = x+nx ;
        p3_y = y+ny ;
        p3_z = z+nz ; // point on plane
        numer =
          nx * (p3_x - p1_x) +
          ny * (p3_y - p1_y) +
          nz * (p3_z - p1_z) ;
        denom =
          nx * (p2_x - p1_x) +
          ny * (p2_y - p1_y) +
          nz * (p2_z - p1_z) ;
        if (DZERO(denom))
        {
          DiagBreak() ;
          continue ;
        }
        u = numer / denom ;
        xf = p1_x + u * Ix ;
        yf = p1_y + u * Iy ;
        zf = p1_z + u * Iz ;
        MRIsampleVolume(mri_src, xf, yf, zf, &oval) ;
        oval = oval*thresh_dir ;
        if (oval > val)   // not at a maximum
        {
          MRIsetVoxVal(mri_sup, x, y, z, 0, 0) ;
          continue ;
        }
        xf = p1_x - u * Ix ;
        yf = p1_y - u * Iy ;
        zf = p1_z - u * Iz ;
        MRIsampleVolume(mri_src, xf, yf, zf, &oval) ;
        oval = oval*thresh_dir ;
        if (oval > val)   // not at a maximum
        {
          MRIsetVoxVal(mri_sup, x, y, z, 0, 0) ;
          continue ;
        }
      }
    }
  }

  return(mri_sup) ;
}
MRI *
MRIextractRegionAndPad(MRI *mri_src, MRI *mri_dst, MRI_REGION *region, int pad)
{
  MRI *mri_tmp ;
  MRI_REGION box ;
  MATRIX *m_src_vox2ras, *m_trans, *m_dst_vox2ras ;

  if (region == NULL)
  {
    region = & box ;
    box.x = box.y = box.z = 0 ;
    box.dx = mri_src->width ;
    box.dy = mri_src->height ;
    box.dz = mri_src->depth ;
  }
  mri_dst =
    MRIallocSequence
    (region->dx+2*pad, region->dy+2*pad, region->dz+2*pad, mri_src->type,mri_src->nframes) ;
  MRIcopyHeader(mri_src, mri_dst) ;
  mri_tmp = MRIextractInto(mri_src, NULL, region->x, region->y, region->z,
                           region->dx, region->dy, region->dz, 0, 0, 0) ;
  MRIextractInto(mri_tmp, mri_dst, 0, 0, 0, region->dx, region->dy, region->dz, pad, pad, pad);
  MRIfree(&mri_tmp) ;
  m_src_vox2ras = MRIgetVoxelToRasXform(mri_src) ;
  m_trans = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_trans, 1, 4) = -pad ;
  *MATRIX_RELT(m_trans, 2, 4) = -pad ;
  *MATRIX_RELT(m_trans, 3, 4) = -pad ;
  m_dst_vox2ras = MatrixMultiply(m_src_vox2ras, m_trans, NULL) ;
  MRIsetVox2RASFromMatrix(mri_dst, m_dst_vox2ras) ;

  MatrixFree(&m_src_vox2ras) ; MatrixFree(&m_trans) ; MatrixFree(&m_dst_vox2ras) ;
  
  return(mri_dst) ;
}

MRI *
MRIsetValuesOutsideRegion(MRI *mri_src,
                          MRI_REGION *region,
                          MRI *mri_dst,
                          float val)
{
  int x, y, z ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (x = 0 ; x < mri_dst->width ; x++)
  {
    if (x >= region->x && x < region->x+region->dx)
      continue ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      if (y >= region->y && y < region->y+region->dy)
        continue ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (z >= region->z && z < region->z+region->dz)
          continue ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, val) ;
      }
    }
  }

  return(mri_dst) ;
}

int
MRIlabelsInNbhd6(MRI *mri, int x, int y, int z, int label) 
{
  int xi, yi, zi, xk, yk, zk, count ;

  for (count = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
	if (fabs(xk) + fabs(yk) + fabs(zk) > 1)
	  continue ;
        if (nint(MRIgetVoxVal(mri, xi, yi, zi,0)) == label)
          count++;
      }
    }
  }
  return(count) ;
}

int
MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label)
{
  int xi, yi, zi, xk, yk, zk, count ;

  for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (nint(MRIgetVoxVal(mri, xi, yi, zi,0)) == label)
          count++;
      }
    }
  }
  return(count) ;
}

int
MRIlabelsInPlanarNbhd(MRI *mri,
                      int x, int y, int z,
                      int whalf,
                      int label,
                      int which)
{
  int xi, yi, zi, xk, yk, zk, count ;

  switch (which)
  {
  case MRI_HORIZONTAL:
    for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
    {
      zi = mri->zi[z+zk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (nint(MRIgetVoxVal(mri, xi, y, zi,0)) == label)
          count++;
      }
    }
    break ;
  case MRI_SAGITTAL:
    for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
    {
      zi = mri->zi[z+zk] ;
      for (yk = -whalf ; yk <= whalf ; yk++)
      {
        yi = mri->yi[y+yk] ;
        if (nint(MRIgetVoxVal(mri, x, yi, zi,0)) == label)
          count++;
      }
    }
    break ;
  case MRI_CORONAL:
    for (count = 0, xk = -whalf ; xk <= whalf ; xk++)
    {
      xi = mri->xi[x+xk] ;
      for (yk = -whalf ; yk <= whalf ; yk++)
      {
        yi = mri->yi[y+yk] ;
        if (nint(MRIgetVoxVal(mri, xi, yi, z,0)) == label)
          count++;
      }
    }
    break ;
  default:
    ErrorReturn
      (0, 
       (ERROR_UNSUPPORTED, 
        "MRIlabelsInPlanarNbhd: which=%d not supported yet", 
        which)) ;
    break ;
  }
  return(count) ;
}

MRI *
MRImatchMeanIntensity(MRI *mri_source, MRI *mri_target, MRI *mri_source_scaled)
{
  float   mean_source, mean_target, scale ;

  mri_source_scaled = MRIcopy(mri_source, mri_source_scaled) ;

  mean_source = MRImeanFrameThresh(mri_source, 0, 1) ;
  mean_target = MRImeanFrameThresh(mri_target, 0, 1) ;
  scale = mean_target / mean_source ;
  printf("source mean = %2.1f, target mean = %2.1f, scaling by %2.3f\n",
         mean_source, mean_target, scale) ;
  MRIscalarMul(mri_source, mri_source_scaled, scale) ;
  return(mri_source_scaled) ;
}

/*-----------------------------------------------------
  Parameters:    LTA or XFM filename of atlas transform, and scale factor
                 optional ptrs to determinant

  Returns value: estimated total intracranial volume, per Buckner et al. 2004
                 also, if supplied with ptrs, the scale factor and
                 determinant can be returned (for debug)

  Description:   https://surfer.nmr.mgh.harvard.edu/fswiki/eTIV
  ------------------------------------------------------*/
// Buckner version: eTIV_scale_factor = 1755;
// old version with talairach_with_skull.lta: eTIV_scale_factor = 2889.2;
// NJS calculated scale factor, with newest talairach_with_skull.lta: 2150
// NJS: scale factor must be passed by caller, since it depends on which
// transform is being passed-in
double MRIestimateTIV(char* theLtaFile,
                      double theScaleFactor,
                      double* theAtlasDet)
{
  LTA *atlas_lta ;
  double atlas_det=0, atlas_icv=0 ;

  atlas_lta = LTAreadEx(theLtaFile);

  if (atlas_lta == NULL)
    ErrorExit
    (ERROR_NOFILE,
     "%s: could not open atlas transform file %s",
     Progname, theLtaFile) ;

  atlas_det = MatrixDeterminant(atlas_lta->xforms[0].m_L) ;

  LTAfree(&atlas_lta) ;

  atlas_icv = (theScaleFactor*(10*10*10)) / atlas_det ;

  // if given ptr, return determinant, for debug purposes
  if (theAtlasDet) *theAtlasDet = atlas_det;

  return atlas_icv;
}

int MRInormalizeFrames(MRI *mri)
{
  int    c,r,s,f;
  double mean, var, val, std ;

  if (NULL==mri)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"MRInormalize: no MRI\n"));

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        for (mean = var = 0.0, f = 0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, c, r, s, f);
          mean += val;
          var += val*val;
        }
        mean /= mri->nframes ;
        var = var / mri->nframes - mean*mean ;
        std = sqrt(var) ;

        for (f = 0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, c, r, s, f);
          val = (val - mean) / std ;
          MRIsetVoxVal(mri, c, r, s, f, val);
        }
      }
    }
  }

  return(ERROR_NONE);
}

// divides by mean (per voxel across frames)
int MRInormalizeFramesMean(MRI *mri)
{
  int    c,r,s,f;
  double mean, val ;

  if (NULL==mri)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"MRInormalize: no MRI\n"));

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        for (mean = 0.0, f = 0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, c, r, s, f);
          mean += val;
        }
        mean /= mri->nframes ;
				if (fabs(mean) < 0.000001) mean = 1;
        for (f = 0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, c, r, s, f);
          val = val / mean ;
          MRIsetVoxVal(mri, c, r, s, f, val);
        }
      }
    }
  }

  return(ERROR_NONE);
}

// divides by all frames by first
int MRInormalizeFramesFirst(MRI *mri)
{
  int    c,r,s,f;
  double  val,val1 ;

  if (NULL==mri)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"MRInormalize: no MRI\n"));

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {

        val1 = MRIgetVoxVal(mri, c, r, s, 0);

        for (f = 0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, c, r, s, f);
          val = val / val1 ;
          MRIsetVoxVal(mri, c, r, s, f, val);
        }
      }
    }
  }

  return(ERROR_NONE);
}

MRI *
MRIzeroMean(MRI *mri_src, MRI *mri_dst)
{
  int     x, y, z, f ;
  double  mean ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;

  for (mean = 0.0, f = 0 ; f < mri_src->nframes ; f++)
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
      {
        for (z = 0 ; z < mri_src->height ; z++)
        {
          mean += MRIgetVoxVal(mri_src, x, y, z, f) ;
        }
      }
    }
  mean /= (mri_src->width*mri_src->height*mri_src->depth*mri_src->nframes) ;
  for (f = 0 ; f < mri_src->nframes ; f++)
    for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
      {
        for (z = 0 ; z < mri_src->height ; z++)
        {
          MRIsetVoxVal(mri_dst, x, y, z, f, 
                       MRIgetVoxVal(mri_src, x, y, z, f)-mean) ;
        }
      }
    }
  return(mri_dst) ;
}
/*
  remove the mean from the timecourse at each voxel.
*/
MRI *
MRIzeroMeanTimecourse(MRI *mri_src, MRI *mri_dst)
{
  int     x, y, z, f ;
  double  mean, val ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;

  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->height ; z++)
      {
        for (mean = 0.0, f = 0 ; f < mri_src->nframes ; f++)
          mean += MRIgetVoxVal(mri_src, x, y, z, f) ;
        mean /= mri_src->nframes ;
        for (f = 0 ; f < mri_src->nframes ; f++)
        {
          val = MRIgetVoxVal(mri_src, x, y, z, f) ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val-mean) ;
        }
      }
    }
  }
  return(mri_dst) ;
}
/*
  compute the mean from the timecourse at each voxel.
*/
MRI *
MRImeanTimecourse(MRI *mri_src, MRI *mri_dst)
{
  int     x, y, z, f ;
  double  mean ;

  if (mri_dst == NULL)
  {
    mri_dst = 
      MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_src, mri_dst) ;
  }

  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->height ; z++)
      {
        for (mean = 0.0, f = 0 ; f < mri_src->nframes ; f++)
          mean += MRIgetVoxVal(mri_src, x, y, z, f) ;
        mean /= mri_src->nframes ;
        MRIsetVoxVal(mri_dst, x, y, z, 0, mean) ;
      }
    }
  }
  return(mri_dst) ;
}

/*---------------------------------------------------------------
  MRIsort() - sorts the frames of each voxel in ascending order.
  ---------------------------------------------------------------*/
MRI *MRIsort(MRI *in, MRI *mask, MRI *sorted)
{
  int c, r, s, f, ncols, nrows, nslices,nframes;
  float m;
  double  *vf;

  ncols   = in->width;
  nrows   = in->height;
  nslices = in->depth;
  nframes = in->nframes;

  if (sorted==NULL)
  {
    sorted = MRIallocSequence(ncols, nrows, nslices, in->type, nframes);
    if (sorted==NULL)
    {
      printf("ERROR: MRIsort: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(in,sorted); // ordinarily would need to change nframes
  }
  if (sorted->width != ncols   || sorted->height  != nrows ||
      sorted->depth != nslices || sorted->nframes != nframes)
  {
    printf("ERROR: MRIsort: dimension mismatch\n");
    return(NULL);
  }

  vf = (double *) calloc(in->nframes,sizeof(double));
  if (vf==NULL)
  {
    printf("ERROR: MRIsort: could not alloc vf of size %d\n",in->nframes);
    return(NULL);
  }

  for (s=0; s<nslices; s++)
  {
    for (r=0; r<nrows; r++)
    {
      for (c=0; c<ncols; c++)
      {
        if (mask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5) continue;
        }
        for (f=0; f<nframes; f++)
          vf[f] = MRIgetVoxVal(in,c,r,s,f);
        qsort(vf,nframes,sizeof(double),CompareDoubles);
        for (f=0; f<nframes; f++)
          MRIsetVoxVal(sorted,c,r,s,f,vf[f]);
      } // cols
    } // rows
  } // slices

  return(sorted);
}


/*------------------------------------------------
  CompareDoubles() - simply compares two doubles in a
  way compatible with qsort.
  ------------------------------------------------*/
int CompareDoubles(const void *a, const void *b)
{
  double ad, bd;

  ad = *((double *)a);
  bd = *((double *)b);

  if (ad < bd) return(-1);
  if (ad > bd) return(+1);
  return(0);
}
MRI *
MRIaddScalar(MRI *mri_src, MRI *mri_dst, float scalar)
{
  float  val ;
  int    x, y, z, f ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (f = 0 ; f < mri_dst->nframes ; f++)
    for (x = 0 ; x < mri_dst->width ; x++)
      for (y = 0 ; y < mri_dst->height ; y++)
        for (z = 0 ; z < mri_dst->depth ; z++)
        {
          val = MRIgetVoxVal(mri_dst, x, y, z, f) ;
          MRIsetVoxVal(mri_dst, x, y, z, f, val+scalar) ;
        }
  return(mri_dst) ;
}
int
MRIgeometryMatched(MRI *mri1, MRI *mri2)
{
  if (mri1->width != mri2->width ||
      mri1->height != mri2->height ||
      mri1->depth != mri2->depth)
    return(0) ;
  
  if (!FEQUAL(mri1->xsize, mri2->xsize) || 
      !FEQUAL(mri1->ysize, mri2->ysize) || 
      !FEQUAL(mri1->zsize, mri2->zsize))
    return(0) ;

  if (!FEQUAL(mri1->c_r, mri2->c_r) || 
      !FEQUAL(mri1->c_a, mri2->c_a) || 
      !FEQUAL(mri1->c_s, mri2->c_s))
    return(0) ;
  return(1) ;
}

      
MRI *
MRIfillHoles(MRI *mri_src, MRI *mri_fill, int thresh) 
{
  int  nfilled, ntotal = 0, cnt, cntmax= thresh ;
  int im0,x0,i0,z,i,x;
  int v,vmax, ylim0, ylim1, xlim0, xlim1;

  mri_fill = MRIcopy(mri_src, mri_fill) ;

  ylim0 = 0 ; ylim1 = mri_src->height-1 ;
  xlim0 = 0 ; xlim1 = mri_src->width-1 ;
  do 
  {
    nfilled = 0;
    for (z=1;z!=mri_fill->depth-1;z++)
      for (i=1;i!=mri_fill->height-1;i++)
        for (x=1;x!=mri_fill->width-1;x++)
          if (MRIgetVoxVal(mri_fill, x, i, z, 0) == 0 &&
              i>ylim0-10 && i<ylim1+10 && x>xlim0-10 && x<xlim1+10) 
          {
            cnt = 0;
            vmax = 0;
            for (im0= -1;im0<=1;im0++)
              for (i0= -1;i0<=1;i0++)
                for (x0= -1;x0<=1;x0++) {
                  v = MRIgetVoxVal(mri_fill, x+x0, i+i0, z+im0, 0) ;
                  if (v>vmax) vmax = v;
                  if (v == 0) cnt++;/* count # of nbrs which are off */
                  if (cnt>cntmax) im0=i0=x0=1;  /* break out
                                                   of all 3 loops */
                }
            if (cnt<=cntmax)   /* toggle pixel (off to on, or on to off) */
            {
              MRIsetVoxVal(mri_fill, x, i, z, 0, vmax);
              nfilled++;
              ntotal++;
            }
          }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%d holes filled\n",nfilled);
  } while (nfilled > 0) ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "total of %d holes filled\n",ntotal);
  return(mri_fill) ;
}
int *
MRIhistogramLabels(MRI *mri, int *counts, int max_label)
{
  int   x, y, z, label ;

  if (counts == NULL)
    counts = (int *)calloc(max_label+1, sizeof(*counts)) ;
  if (counts == NULL)
    ErrorExit(ERROR_NOMEMORY, 
              "MRIhistogramLabels: could not allocate %d counts",
              max_label+1) ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        label = (int)MRIgetVoxVal(mri, x, y, z, 0) ;
        if (label <= max_label)
          counts[label]++ ;
      }
  return(counts) ;
}

int
MRIlabelBoundingBox(MRI *mri, int label, MRI_REGION *region)
{
  int    x, y, z, l, x1, y1, z1 ;

  region->x = mri->width ;
  region->y = mri->height ;
  region->z = mri->depth ;
  x1 = y1 = z1 = -1 ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        l = (int)MRIgetVoxVal(mri, x, y, z, 0) ;
        if (l == label)
        {
          region->x = MIN(x, region->x) ;
          region->y = MIN(y, region->y) ;
          region->z = MIN(z, region->z) ;
          x1 = MAX(x, x1) ; y1 = MAX(y, y1) ; z1 = MAX(z, z1) ;
        }
      }

  region->dx = x1-region->x+1 ;
  region->dy = y1-region->y+1 ;
  region->dz = z1-region->z+1 ;
  return(NO_ERROR) ;
}

/*-------------------------------------------------------------------
  MRIcopyVolGeomToMRI - copies the volume geometry passed in into the MRI
  structure.
  -------------------------------------------------------------------*/
int
MRIcopyVolGeomToMRI(MRI *mri, const VOL_GEOM *vg)
{
  mri->xsize = vg->xsize ;
  mri->ysize = vg->ysize ;
  mri->zsize = vg->zsize ;
  mri->x_r = vg->x_r;
  mri->y_r = vg->y_r;
  mri->z_r = vg->z_r;
  mri->c_r = vg->c_r;
  mri->x_a = vg->x_a;
  mri->y_a = vg->y_a;
  mri->z_a = vg->z_a;
  mri->c_a = vg->c_a;
  mri->x_s = vg->x_s;
  mri->y_s = vg->y_s;
  mri->z_s = vg->z_s;
  mri->c_s = vg->c_s;
  return(NO_ERROR) ;
}
/*-------------------------------------------------------------------
  MRIcopyVolGeomToMRI - copies the volume geometry passed in into the MRI
  structure.
  -------------------------------------------------------------------*/
int
MRIcopyVolGeomFromMRI( const MRI *mri, VOL_GEOM *vg)
{
  vg->xsize = mri->xsize ;
  vg->ysize = mri->ysize ;
  vg->zsize = mri->zsize ;
  vg->x_r = mri->x_r ;
  vg->y_r = mri->y_r ;
  vg->z_r = mri->z_r ;
  vg->c_r = mri->c_r ;
  vg->x_a = mri->x_a ;
  vg->y_a = mri->y_a ;
  vg->z_a = mri->z_a ;
  vg->c_a = mri->c_a ;
  vg->x_s = mri->x_s ;
  vg->y_s = mri->y_s ;
  vg->z_s = mri->z_s ;
  vg->c_s = mri->c_s ;
  vg->valid = mri->ras_good_flag ;
  strcpy(vg->fname, mri->fname) ;
  return(NO_ERROR) ;
}
int 
MRIfillRegion(MRI *mri, int x,int y,int z,float fill_val,int whalf)
{
  int   xi, xk, yi, yk, zi, zk ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri->zi[z+zk] ;
        MRIsetVoxVal(mri, xi, yi, zi, 0, fill_val) ;
      }
    }
  }
  return(NO_ERROR) ;
}
MRI  *
MRImatchIntensityRatio(MRI *mri_source, MRI *mri_target, MRI *mri_matched, 
                       double min_scale, double max_scale, double low_thresh,
                       double high_thresh)
{
  HISTOGRAM *h ;
  int       x, y, z, bin, peak ;
  float     val_source, val_target, ratio ;
  
#define NBINS  256
  h = HISTOalloc(NBINS) ;
  h->bin_size = (max_scale - min_scale) / ((float)NBINS-1) ;
  for (bin = 0 ; bin < h->nbins ; bin++)
    h->bins[bin] = min_scale + h->bin_size*(float)bin ;
  
  for (x = 0 ; x < mri_source->width ; x++)
    for (y = 0 ; y < mri_source->height ; y++)
      for (z = 0 ; z < mri_source->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val_source = MRIgetVoxVal(mri_source, x, y, z, 0); 
        val_target = MRIgetVoxVal(mri_target, x, y, z, 0); 
        if (FZERO(val_source) || 
            val_source < low_thresh || val_target < low_thresh ||
            val_source > high_thresh || val_target > high_thresh)
          continue ;
        ratio = val_target / val_source ;
        bin = nint((ratio - min_scale) / h->bin_size) ;
        if (bin < 0 || bin >= NBINS)
          continue ;
        h->counts[bin]++ ;
        if (bin == Gdiag_no)
          DiagBreak() ;
      }

  if (Gdiag & DIAG_WRITE)
  {
    printf("saving h.plt for ratio histogram\n") ;
    HISTOplot(h, "h.plt") ;
  }

  peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins) ;
  if (peak < 0)
  {
    printf("warning: MRImatchIntensityRatio could not find any peaks\n");
    ratio = 1.0 ;
  }
  ratio = h->bins[peak] ;
  printf("peak ratio at %2.3f\n",ratio) ;
  if (FZERO(ratio))
  {
    printf("warning: MRImatchIntensityRatio peak at 0 - disabling\n");
    ratio = 1.0 ;
  }
  mri_matched = MRIscalarMul(mri_source, mri_matched, ratio) ;
  HISTOfree(&h) ;
  return(mri_matched) ;
}

int
MRIcountThreshInNbhd(MRI *mri, int wsize, int x, int y, int z, float thresh)
{
  int   xk, yk, zk, xi, yi, zi, whalf, total ;

  whalf = (wsize-1)/2 ;
  for (total = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIgetVoxVal(mri, xi, yi, zi,0) > thresh)
          total++ ;
      }
    }
  }
  return(total) ;
}
int
MRIfillBox(MRI *mri, MRI_REGION *box, float fillval)
{
  int   x, y, z, xmin, xmax, ymin, ymax, zmin, zmax ;

  xmin = MAX(0, box->x) ;
  xmax = MIN(mri->width-1, box->x+box->dx-1) ;
  ymin = MAX(0, box->y) ;
  ymax = MIN(mri->height-1, box->y+box->dy-1) ;
  zmin = MAX(0, box->z) ;
  zmax = MIN(mri->depth-1, box->z+box->dz-1) ;
  for (x = xmin ; x <= xmax ; x++)
    for (y = ymin ; y <= ymax ; y++)
      for (z = zmin ; z <= zmax ; z++)
        MRIsetVoxVal(mri, x, y, z, 0, fillval) ;
  return(NO_ERROR) ;
}

MRI *
MRIcloneDifferentType(MRI *mri_src, int type)
{
  MRI *mri_dst ;

  if (type == mri_src->type)
    return(MRIclone(mri_src, NULL)) ;

  mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, 
                             type, mri_src->nframes) ;
  MRIcopyHeader(mri_src, mri_dst) ;
  return(mri_dst) ;
}

double
MRImaxNorm(MRI *mri)
{
  double max_norm, norm, val ;
  int    x, y, z, f ;

  max_norm = 0 ;
  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
        for (f = 0, norm = 0.0 ; f < mri->nframes ; f++)
        {
          val = MRIgetVoxVal(mri, x, y, z, f) ;
          norm += (val*val) ;
        }
        if (norm > max_norm)
          max_norm = norm ;
      }
  return(sqrt(max_norm)) ;
}

static int compare_sort_mri(const void *plp1, const void *plp2);
typedef struct
{
  unsigned char x, y, z, val ;
}
SORT_VOXEL ;


static int
compare_sort_mri(const void *psv1, const void *psv2)
{
  SORT_VOXEL  *sv1, *sv2 ;

  sv1 = (SORT_VOXEL *)psv1 ;
  sv2 = (SORT_VOXEL *)psv2 ;

  if (sv1->val > sv2->val)
    return(1) ;
  else if (sv1->val < sv2->val)
    return(-1) ;

  return(0) ;
}
int
MRIorderIndices(MRI *mri, short *x_indices, short *y_indices, short *z_indices)
{
  int         width, height, depth, nindices, index, x, y, z ;
  SORT_VOXEL  *sort_voxels ;

  width = mri->width, height = mri->height ;
  depth = mri->depth ;
  nindices = width*height*depth ;

  sort_voxels = (SORT_VOXEL *)calloc(nindices, sizeof(SORT_VOXEL)) ;
  if (!sort_voxels)
    ErrorExit(ERROR_NOMEMORY,"MRIorderIndices: could not allocate sort table");

  for (index = x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++, index++)
      {
        sort_voxels[index].x = x ;
        sort_voxels[index].y = y ;
        sort_voxels[index].z = z ;
        sort_voxels[index].val = MRIgetVoxVal(mri, x, y, z,0) ;
      }
    }
  }
  qsort(sort_voxels, nindices, sizeof(SORT_VOXEL), compare_sort_mri) ;

  for (index = 0 ; index < nindices ; index++)
  {
    x_indices[index] = sort_voxels[index].x ;
    y_indices[index] = sort_voxels[index].y ;
    z_indices[index] = sort_voxels[index].z ;
  }

  free(sort_voxels) ;
  return(NO_ERROR) ;
}
int
MRIcomputeVoxelPermutation(MRI *mri, short *x_indices, short *y_indices,
                           short *z_indices)
{
  int width, height, depth, tmp, nindices, i, index ;

  width = mri->width, height = mri->height ;
  depth = mri->depth ;
  nindices = width*height*depth ;

  for (i = 0 ; i < nindices ; i++)
  {
    x_indices[i] = i % width ;
    y_indices[i] = (i/width) % height ;
    z_indices[i] = (i / (width*height)) % depth ;
  }
  for (i = 0 ; i < nindices ; i++)
  {
    index = (int)randomNumber(0.0, (double)(nindices-0.0001)) ;

    tmp = x_indices[index] ;
    x_indices[index] = x_indices[i] ;
    x_indices[i] = tmp ;

    tmp = y_indices[index] ;
    y_indices[index] = y_indices[i] ;
    y_indices[i] = tmp ;

    tmp = z_indices[index] ;
    z_indices[index] = z_indices[i] ;
    z_indices[i] = tmp ;
  }
  return(NO_ERROR) ;
}
MRI *
MRImaskZero(MRI *mri_src, MRI *mri_mask, MRI *mri_dst)
{
  int   x, y, z ;

  if (mri_dst == NULL)
    mri_dst = MRIclone(mri_src, NULL) ;

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (MRIgetVoxVal(mri_mask, x, y, z, 0) > 0)
          MRIsetVoxVal(mri_dst, x, y, z, 0, 
                       MRIgetVoxVal(mri_src, x, y, z, 0)) ;
        else
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
      }
  return(mri_dst) ;
}


const char* MRItype2str(int type)
{
  switch (type)
  {
  case MRI_UCHAR:
    return("uchar");
  case MRI_SHORT:
    return("short");
  case MRI_INT:
    return("int");
  case MRI_LONG:
    return("long");
  case MRI_FLOAT:
    return("float");
  case MRI_BITMAP:
    return("bitmap");
  case MRI_TENSOR:
    return("tensor");
  }
  return("unknown data type");
}
#define MAX_VOX 5000
#define MAX_FOV 300
int
MRIfindSliceWithMostStructure(MRI *mri_aseg, int slice_direction, int label)
{
  int     max_slice, max_vox, cor_vox[MAX_VOX], hor_vox[MAX_VOX], sag_vox[MAX_VOX], x, y, z, i, vox = 0 ;
  float   r, a, s, vsize ;
  VECTOR  *v_ras, *v_vox ;
  MATRIX  *m_ras2vox ;

  v_ras = VectorAlloc(4, MATRIX_REAL) ;
  v_vox = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_vox, 4) = 1.0 ; 
  VECTOR_ELT(v_ras, 4) = 1.0 ; 
  m_ras2vox = MRIgetRasToVoxelXform(mri_aseg) ;
  memset(cor_vox, 0, sizeof(cor_vox)) ; memset(sag_vox, 0, sizeof(sag_vox)) ; memset(hor_vox, 0, sizeof(hor_vox)) ;
  
  vsize = MIN(mri_aseg->xsize, MIN(mri_aseg->ysize, mri_aseg->zsize)) ;
  for (r = -MAX_FOV ; r <= MAX_FOV ; r += vsize)
  {
    V3_X(v_ras) = r ;
    for (a = -MAX_FOV ; a <= MAX_FOV ; a += vsize)
    {
      V3_Y(v_ras) = a ;
      for (s = -MAX_FOV ; s <= MAX_FOV ; s += vsize)
      {
        V3_Z(v_ras) = s ;
        MatrixMultiply(m_ras2vox, v_ras, v_vox) ;
        x = nint(V3_X(v_vox)) ; y = nint(V3_Y(v_vox)) ; z = nint(V3_Z(v_vox)) ;
        if (x < 0 || x >= mri_aseg->width ||
            y < 0 || y >= mri_aseg->height ||
            z < 0 || z >= mri_aseg->depth)
          continue ;
        if (MRIgetVoxVal(mri_aseg, x, y, z, 0) != label)
          continue ;
        cor_vox[z]++ ;
        sag_vox[x]++ ;
        hor_vox[y]++ ;
      }
    }
  }
  VectorFree(&v_ras) ; VectorFree(&v_vox) ; MatrixFree(&m_ras2vox) ;
  
  for (max_vox = max_slice =i = 0 ; i < MAX_VOX ; i++)
  {
    switch (slice_direction)
    {
    default:
      ErrorExit(ERROR_UNSUPPORTED, "MRIfindSliceWithMostStructure: unknown slice direction %d", slice_direction);
    case MRI_CORONAL:    vox = cor_vox[i] ; break ; 
    case MRI_SAGITTAL:   vox = sag_vox[i] ; break ; 
    case MRI_HORIZONTAL: vox = hor_vox[i] ; break ; 
    }
    if (vox > max_vox)
    {
      max_vox = vox ;
      max_slice = i ;
    }
  }
  return(max_slice) ;
}
double
MRIrmsDiff(MRI *mri1, MRI *mri2)
{
  double  rms, val1, val2 ;
  int     x, y, z, nvox ;

  for (rms = 0.0, nvox = x = 0 ;  x < mri1->width ; x++)
    for (y = 0 ; y < mri1->height ; y++)
      for (z = 0 ; z < mri1->depth ; z++)
      {
        val1 = MRIgetVoxVal(mri1, x, y, z, 0) ;
        val2 = MRIgetVoxVal(mri2, x, y, z, 0) ;
        if (!FZERO(val1) || !FZERO(val2))
        {
          nvox++ ;
          rms += (val1-val2)*(val1-val2) ;
        }
      }
  if (nvox > 0)
    rms = sqrt(rms/nvox) ;
  return(rms) ;
}


// compute root mean square of 'in', which is assumed to be multi-framed
// writes to 'out'
void MRIrms(MRI *in, MRI *out)
{
  int f,z,y,x;
  int width = in->width ;
  int height = in->height ;
  int depth = in->depth ;
  int nframes = in->nframes ;
  if (nframes == 0) nframes = 1;

  // square and sum each frame of input
  // then divide by nframes and take sqrt
  for (f = 0 ; f < nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          double vin = MRIgetVoxVal(in,x,y,z,f);
          double vout = 0; // output summation
          if (f != 0)
          {
            // after first frame, output gets summed
            vout = MRIgetVoxVal(out,x,y,z,0);
          }
          double v = (vin*vin) + vout; // square and sum
          if (f == (nframes - 1)) // if last frame, div and sqrt
          {
            v /= nframes; // divide
            v = sqrt(v); // square root
          }
          MRIsetVoxVal(out,x,y,z,0,v);
        }
      }
    }
  }
}

int
MRImaskLabel(MRI *mri_src, MRI *mri_dst, MRI *mri_labeled, int label_to_mask, float out_val) 
{
  int    x, y, z, label ;

  if (mri_dst != mri_src)
    mri_dst = MRIcopy(mri_src, mri_dst) ;
  for (x = 0 ; x < mri_src->width; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth; z++)
      {
	label = MRIgetVoxVal(mri_labeled, x, y, z, 0) ;
	if (label == label_to_mask)
	  MRIsetVoxVal(mri_dst, x, y, z, 0, out_val) ;
      }
  return(NO_ERROR) ;
}
int
MRIcomputeVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox, 
                                 MRI *mri_seg, MRI *mri_fractions)
{
  int    x, y, z, xs, ys, zs, label ;
  VECTOR *v1, *v2 ;
  MRI    *mri_counts ;
  float  val, count ;
  MATRIX *m_inv ;

  m_inv = MatrixInverse(m_vox2vox, NULL) ;
  if (m_inv == NULL)
  {
    MatrixPrint(stdout, m_vox2vox) ;
    ErrorExit(ERROR_BADPARM, "MRIcomputePartialVolumeFractions: non-invertible vox2vox matrix");
  }
  mri_counts = MRIcloneDifferentType(mri_src, MRI_INT) ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  for (x = 0 ; x < mri_seg->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_seg->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
        if (xs == Gx && ys == Gy && zs == Gz)
          DiagBreak() ;
        if (xs >= 0 && ys >= 0 && zs >= 0 &&
            xs < mri_src->width && ys < mri_src->height && zs < mri_src->depth)
        {
          val = MRIgetVoxVal(mri_counts, xs, ys, zs, 0) ;
          if (val > 0)
            DiagBreak() ;
          MRIsetVoxVal(mri_counts, xs, ys, zs, 0, val+1) ;

          label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
          if (label >= 0 && label <= mri_fractions->nframes)
          {
            val = MRIgetVoxVal(mri_fractions, xs, ys, zs, label-1) ;  // labels are frame+1
            MRIsetVoxVal(mri_fractions, xs, ys, zs, label-1, val+1) ;
          }
          else
            DiagBreak() ;
        }
      }
    }
  }

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        count = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
        if (count >= 1)
        {
          for (label = 0 ; label < mri_fractions->nframes ; label++)
          {
            if (x == Gx && y == Gy && z == Gz)
              DiagBreak() ;
            val = MRIgetVoxVal(mri_fractions, x, y, z, label) ;
            MRIsetVoxVal(mri_fractions, x, y, z, label, val/count) ;
          }
        }
        else  // sample in other direction
        {
          V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
          MatrixMultiply(m_inv, v1, v2) ;
          MatrixMultiply(m_inv, v1, v2) ;
          xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
          if (xs >= 0 && ys >= 0 && zs >= 0 &&
              xs < mri_seg->width && ys < mri_seg->height && zs < mri_seg->depth)
          {
            label = MRIgetVoxVal(mri_seg, xs, ys, zs, 0) ;
            if (label >= 0 && label < mri_fractions->nframes)
              MRIsetVoxVal(mri_fractions, x, y, z, label-1, 1) ;
            else
              DiagBreak() ;
          }
        }
      }
  VectorFree(&v1) ; VectorFree(&v2) ;
  MRIfree(&mri_counts) ; MatrixFree(&m_inv) ;
  
  return(NO_ERROR) ;
}
