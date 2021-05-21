/**
 * @brief new talairach related routines with Ex
 *
 * takes lta as the talairach transform (use LTAreadEx routine)
 * doesn't rely on COR volume type
 */
/*
 * Original Author: Yasunari Tosa
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

#include "talairachex.h"
#include "diag.h"
#include "error.h"
#include "proto.h"

extern const char *Progname;

////////////////////////////////////////////////////////////////////

#define V4_LOAD(v, x, y, z, r) (VECTOR_ELT(v, 1) = x, VECTOR_ELT(v, 2) = y, VECTOR_ELT(v, 3) = z, VECTOR_ELT(v, 4) = r);

int ModifyTalairachCRAS(MRI *mri_tal, const LTA *lta)
{
  LT *tran = 0;
  // if lta is given
  if (lta != 0) {
    if (lta->num_xforms == 1) {
      if (lta->type == LINEAR_RAS_TO_RAS) {
        tran = &lta->xforms[0];
        if (tran->dst.valid == 1)  // transform dst is valid
        {
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("INFO: Modifying dst c_(r,a,s), using the transform dst\n");
          mri_tal->c_r = tran->dst.c_r;
          mri_tal->c_a = tran->dst.c_a;
          mri_tal->c_s = tran->dst.c_s;
        }
        else if (getenv("NO_AVERAGE305"))  // if this is set
          fprintf(stderr, "INFO: tal c_(r,a,s) not modified\n");
        else {
          // use average_305 value
          if (DIAG_VERBOSE_ON) {
            fprintf(stderr, "INFO: Modifying talairach volume c_(r,a,s) based on average_305.\n");
            fprintf(stderr, "INFO: If not preferred, set environmental varible NO_AVERAGE305 true.\n");
          }
          mri_tal->c_r = -0.095;
          mri_tal->c_a = -16.51;
          mri_tal->c_s = 9.75;
        }
      }
      else
        ErrorExit(ERROR_BADPARM, "%s: xfm passed is not of RAS-to-RAS type", Progname);
    }
    else
      ErrorExit(ERROR_BADPARM, "%s: xfm has more than one xfrm", Progname);
  }  // lta != 0
  else {
    if (getenv("NO_AVERAGE305"))  // if this is set
      fprintf(stderr, "INFO: tal c_(r,a,s) not modified\n");
    else {
      // use average_305 value
      if (DIAG_VERBOSE_ON) fprintf(stderr, "INFO: Modifying talairach volume c_(r,a,s) based on average_305\n");
      mri_tal->c_r = -0.095;
      mri_tal->c_a = -16.51;
      mri_tal->c_s = 9.75;
    }
  }  // lta == 0
  // when you modify c_(ras), you must recalculate i_to_r__ and r_to_i__
  if (mri_tal->i_to_r__) {
    AffineMatrixFree(&mri_tal->i_to_r__);
  }

  if (mri_tal->r_to_i__) {
    MatrixFree(&mri_tal->r_to_i__);
  }

  AffineMatrixAlloc(&(mri_tal->i_to_r__));
  MATRIX *tmp = extract_i_to_r(mri_tal);
  SetAffineMatrix(mri_tal->i_to_r__, tmp);
  MatrixFree(&tmp);

  mri_tal->r_to_i__ = extract_r_to_i(mri_tal);

  return (NO_ERROR);
}
///////////////////////////////////////////////////////////////////////////////
//
//            src  --->  RAS
//             |  \       |
//             |   \      |
//             |    \     |
//             V     V    V
//         talVol  ---> Talairach
//
// voxel -> RAS -> talairach RAS
// needs information only on the src volume
// does not need information on talairach volume
////////////////////////////////////////////////////////////////////
// Matrix routines
///////////////////////////////////////////////////////////////////
MATRIX *MtalairachFromVoxel(MRI *mri_src, const LTA *lta)
{
  MATRIX *RASfromVoxel = 0;
  MATRIX *talRASfromVoxel = 0;
  RASfromVoxel = extract_i_to_r(mri_src);
  if (lta) {
    if (lta->num_xforms == 0) ErrorExit(ERROR_BADPARM, "lta does not have xform");
    if (lta->type != LINEAR_RAS_TO_RAS) ErrorExit(ERROR_BADPARM, "lta must be RAS_TO_RAS transform");

    talRASfromVoxel = MatrixMultiply(lta->xforms[0].m_L, RASfromVoxel, NULL);  // allocate memory
  }
  else  // no transform.  just copy
    talRASfromVoxel = MatrixCopy(RASfromVoxel, NULL);

  MatrixFree(&RASfromVoxel);
  return talRASfromVoxel;
}

MATRIX *MtalVoxelFromVoxel(MRI *mri_src, const LTA *lta)
{
  MATRIX *talRASFromVoxel = 0;
  MATRIX *talVoxelFromTalRAS = 0;
  MATRIX *res;
  MRI *mri_talvol;

  talRASFromVoxel = MtalairachFromVoxel(mri_src, lta);
  // use lta->xform[0].dst to construct the talairach volume matrix
  mri_talvol = MRIallocHeader(mri_src->width, mri_src->height, mri_src->depth, mri_src->type, 1);
  MRIcopyHeader(mri_src, mri_talvol);
  ModifyTalairachCRAS(mri_talvol, lta);

  talVoxelFromTalRAS = extract_r_to_i(mri_talvol);
  res = MatrixMultiply(talVoxelFromTalRAS, talRASFromVoxel, NULL);
  MatrixFree(&talRASFromVoxel);
  MatrixFree(&talVoxelFromTalRAS);
  MRIfree(&mri_talvol);
  return res;
}

MATRIX *MvoxelFromTalairach(MRI *mri_dst, const LTA *lta)
{
  MATRIX *RASFromTalairach = 0;
  MATRIX *voxelFromRAS = 0;
  MATRIX *res = 0;
  if (lta)
    RASFromTalairach = MatrixInverse(lta->xforms[0].m_L, NULL);
  else
    RASFromTalairach = MatrixIdentity(4, NULL);
  voxelFromRAS = extract_r_to_i(mri_dst);
  res = MatrixMultiply(voxelFromRAS, RASFromTalairach, NULL);
  MatrixFree(&RASFromTalairach);
  MatrixFree(&voxelFromRAS);
  return res;
}

MATRIX *MvoxelFromTalVoxel(MRI *mri_dst, const LTA *lta)
{
  MATRIX *talairachFromTalVol = 0;
  MATRIX *voxelFromTalRAS = 0;
  MATRIX *res = 0;
  MRI *mri_talvol = 0;

  mri_talvol = MRIallocHeader(mri_dst->width, mri_dst->height, mri_dst->depth, mri_dst->type, 1);
  MRIcopyHeader(mri_dst, mri_talvol);
  ModifyTalairachCRAS(mri_talvol, lta);

  talairachFromTalVol = extract_i_to_r(mri_talvol);
  voxelFromTalRAS = MvoxelFromTalairach(mri_dst, lta);
  res = MatrixMultiply(voxelFromTalRAS, talairachFromTalVol, NULL);
  MatrixFree(&talairachFromTalVol);
  MatrixFree(&voxelFromTalRAS);
  MRIfree(&mri_talvol);
  return res;
}

MATRIX *MRASFromTalVoxel(MRI *mri, const LTA *lta)
{
  MRI *mriTal = 0;
  MATRIX *talRASfromTalVoxel = 0;
  MATRIX *RASfromTalRAS = 0;
  MATRIX *res = 0;

  mriTal = MRIallocHeader(mri->width, mri->height, mri->depth, mri->type, 1);
  MRIcopyHeader(mri, mriTal);
  ModifyTalairachCRAS(mriTal, lta);

  talRASfromTalVoxel = extract_i_to_r(mriTal);
  RASfromTalRAS = MatrixInverse(lta->xforms[0].m_L, NULL);
  res = MatrixMultiply(RASfromTalRAS, talRASfromTalVoxel, NULL);
  MatrixFree(&talRASfromTalVoxel);
  MatrixFree(&RASfromTalRAS);
  MRIfree(&mriTal);

  return res;
}

void TransformWithMatrix(
    const MATRIX *mat, const double x, const double y, const double z, double *px, double *py, double *pz)
{
#if 0
  // VECTOR *src, *dst;
  static VECTOR *src__ = 0;
  static VECTOR *dst__ = 0;
  if (src__ == 0)
    src__ = VectorAlloc(4, MATRIX_REAL);
  if (dst__ == 0)
    dst__ = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(src__, x, y, z, 1.);
  MatrixMultiply(mat, src__, dst__);
  *px = V3_X(dst__);
  *py = V3_Y(dst__);
  *pz = V3_Z(dst__);
  // VectorFree(&src);
  // VectorFree(&dst);
#else

/*
  Original version uses float internally.
  Set this to zero to recover original behaviour
*/
#define DOUBLE_INTERNAL 0

#if DOUBLE_INTERNAL
  *px = mat->data[0] * x + mat->data[1] * y + mat->data[2] * z + mat->data[3];
  *py = mat->data[4] * x + mat->data[5] * y + mat->data[6] * z + mat->data[7];
  *pz = mat->data[8] * x + mat->data[9] * y + mat->data[10] * z + mat->data[11];
#else
  float xf, yf, zf;

  xf = x;
  yf = y;
  zf = z;

  *px = mat->data[0] * xf + mat->data[1] * yf + mat->data[2] * zf + mat->data[3];
  *py = mat->data[4] * xf + mat->data[5] * yf + mat->data[6] * zf + mat->data[7];
  *pz = mat->data[8] * xf + mat->data[9] * yf + mat->data[10] * zf + mat->data[11];
#endif

#endif
}

//////////////////////////////////////////////////////////////////////////////////
// point to point routines
//////////////////////////////////////////////////////////////////////////////////
int MRIvoxelToTalairachEx(
    MRI *mri_src, double xv, double yv, double zv, double *pxt, double *pyt, double *pzt, const LTA *lta)
{
  MATRIX *talairachFromVoxel = MtalairachFromVoxel(mri_src, lta);
  TransformWithMatrix(talairachFromVoxel, xv, yv, zv, pxt, pyt, pzt);

  MatrixFree(&talairachFromVoxel);
  return (NO_ERROR);
}

// voxel -> RAS -> talairach RAS -> talairachVolume
// mri must be the source volume
int MRIvoxelToTalairachVoxelEx(
    MRI *mri_src, double xv, double yv, double zv, double *pxt, double *pyt, double *pzt, const LTA *lta)
{
  MATRIX *talVoxelFromVoxel = MtalVoxelFromVoxel(mri_src, lta);
  TransformWithMatrix(talVoxelFromVoxel, xv, yv, zv, pxt, pyt, pzt);

  MatrixFree(&talVoxelFromVoxel);
  return (NO_ERROR);
}

// talairachRAS -> RAS -> voxel
// needs the target non-tal volume
int MRItalairachToVoxelEx(
    MRI *mri_dst, double xt, double yt, double zt, double *pxv, double *pyv, double *pzv, const LTA *lta)
{
  MATRIX *voxelFromTalairach = MvoxelFromTalairach(mri_dst, lta);
  TransformWithMatrix(voxelFromTalairach, xt, yt, zt, pxv, pyv, pzv);

  MatrixFree(&voxelFromTalairach);
  return (NO_ERROR);
}

// talairachVolume-> talairachRAS -> RAS
// dst is the non-talairach volume
int MRItalairachVoxelToWorldEx(
    MRI *mri_dst, double xt, double yt, double zt, double *pxw, double *pyw, double *pzw, const LTA *lta)
{
  MATRIX *RASfromTalVoxel = MRASFromTalVoxel(mri_dst, lta);
  TransformWithMatrix(RASfromTalVoxel, xt, yt, zt, pxw, pyw, pzw);

  MatrixFree(&RASfromTalVoxel);
  return (NO_ERROR);
}

// talairachVolume-> talairach RAS -> RAS -> voxel
int MRItalairachVoxelToVoxelEx(
    MRI *mri_dst, double xtv, double ytv, double ztv, double *pxv, double *pyv, double *pzv, const LTA *lta)
{
  MATRIX *voxelFromTalVoxel = MvoxelFromTalVoxel(mri_dst, lta);
  TransformWithMatrix(voxelFromTalVoxel, xtv, ytv, ztv, pxv, pyv, pzv);

  MatrixFree(&voxelFromTalVoxel);
  return (NO_ERROR);
}

////////////////////////////////////////////////////////////////////////////////
// volume to volume routines
////////////////////////////////////////////////////////////////////////////////
MRI *MRItoTalairachExInterp(MRI *mri_src, MRI *mri_tal, const LTA *lta, int interp)
{
  MATRIX *voxToTalvoxel = MtalVoxelFromVoxel(mri_src, lta);
  fprintf(stderr, "voxel to talairach voxel transform\n");
  MatrixPrint(stderr, voxToTalvoxel);
  if (!mri_tal) {
    mri_tal = MRIclone(mri_src, NULL);  // data is not copied
    ModifyTalairachCRAS(mri_tal, lta);
  }
  MRIlinearTransformInterp(mri_src, mri_tal, voxToTalvoxel, interp);
  MatrixFree(&voxToTalvoxel);

  return (mri_tal);
}
// volume -> Talairach volume
MRI *MRItoTalairachEx(MRI *mri_src, MRI *mri_tal, const LTA *lta)
{
  MATRIX *voxToTalvoxel = MtalVoxelFromVoxel(mri_src, lta);
  fprintf(stderr, "voxel to talairach voxel transform\n");
  MatrixPrint(stderr, voxToTalvoxel);
  if (!mri_tal) {
    mri_tal = MRIclone(mri_src, NULL);  // data is not copied
    ModifyTalairachCRAS(mri_tal, lta);
  }
  MRIlinearTransform(mri_src, mri_tal, voxToTalvoxel);
  MatrixFree(&voxToTalvoxel);

  return (mri_tal);
}

// assumes mri contain xform
// transform talairach volume into the dst volume
MRI *MRIfromTalairachEx(MRI *mri_tal, MRI *mri_dst, const LTA *lta)
{
  MATRIX *talVoxelToVoxel = MvoxelFromTalVoxel(mri_dst, lta);
  fprintf(stderr, "talairach voxel to voxel transform\n");
  MatrixPrint(stderr, talVoxelToVoxel);
  if (!mri_dst) {
    ErrorExit(ERROR_BADPARM, "%s: Needs target volume to recover c_(ras)", Progname);
    // mri_dst = MRIclone(mri_src, NULL) ;
  }
  MRIlinearTransform(mri_tal, mri_dst, talVoxelToVoxel);
  MatrixFree(&talVoxelToVoxel);

  return (mri_dst);
}

// extract a talairach plane at point (x,y,z)
MRI *MRIextractTalairachPlaneEx(MRI *mri_src, MRI *mri_dst, int orientation, int x, int y, int z, int wsize, LTA *lta)
{
  double e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int whalf, xk, yk, xi, yi, zi;
  double ex, ey, ez, x0, y0, z0;
  // double len;

  whalf = (wsize - 1) / 2;

  if (!mri_dst) {
    // allocate a plane
    mri_dst = MRIalloc(wsize, wsize, 1, MRI_UCHAR);
    MRIcopyHeader(mri_src, mri_dst);
    mri_dst->xstart = x - whalf * mri_dst->xsize;
    mri_dst->ystart = y - whalf * mri_dst->ysize;
    mri_dst->zstart = z - whalf * mri_dst->zsize;
    mri_dst->xend = mri_dst->xstart + wsize * mri_dst->xsize;
    mri_dst->yend = mri_dst->ystart + wsize * mri_dst->ysize;
    mri_dst->zend = mri_dst->zstart + wsize * mri_dst->zsize;
    mri_dst->imnr0 = z + mri_src->imnr0;
    mri_dst->imnr1 = mri_dst->imnr0;
  }
  // get a poisition in talairach volume
  MRIvoxelToTalairachVoxelEx(mri_src, x, y, z, &x0, &y0, &z0, lta);
  switch (orientation) {
    default:
    case MRI_CORONAL: /* basis vectors in x-y plane */
      /* the 'x' basis vector in talairach space */
      ex = (double)x0 + 1;
      ey = (double)y0;
      ez = (double)z0;
      // get the vector in src volume
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector in talairach space */
      ex = (double)x0;
      ey = (double)y0 + 1;
      ez = (double)z0;
      // get the vector
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
    case MRI_HORIZONTAL: /* basis vectors in x-z plane */
      /* the 'x' basis vector in talairach space */
      ex = (double)x0 + 1;
      ey = (double)y0;
      ez = (double)z0;
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector in talairach space */
      ex = (double)x0;
      ey = (double)y0;
      ez = (double)z0 + 1;
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
    case MRI_SAGITTAL: /* basis vectors in y-z plane */
      /* the 'x' basis vector */
      ex = (double)x0;
      ey = (double)y0;
      ez = (double)z0 + 1.0;
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector */
      ex = (double)x0;
      ey = (double)y0 + 1.0;
      ez = (double)z0;
      MRItalairachVoxelToVoxelEx(mri_src, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
  }
  // calculate the length of the vector in x direction
  // len = sqrt(e1_x * e1_x + e1_y * e1_y + e1_z * e1_z);
  /*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  // calculate the length of the vector in y direction
  // len = sqrt(e2_x * e2_x + e2_y * e2_y + e2_z * e2_z);
  /*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri_src->xi[nint(xbase + xk * e1_x)];
      yi = mri_src->yi[nint(ybase + xk * e1_y)];
      zi = mri_src->zi[nint(zbase + xk * e1_z)];
      MRIvox(mri_dst, xk + whalf, yk + whalf, 0) = MRIvox(mri_src, xi, yi, zi);
    }
  }

  return (mri_dst);
}

int MRIeraseTalairachPlaneNewEx(
    MRI *mri, MRI *mri_mask, int orientation, int x, int y, int z, int wsize, int fill_val, LTA *lta)
{
  double e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, xbase, ybase, zbase;
  int whalf, xk, yk, xi, yi, zi, xki, yki, x0, y0;
  double ex, ey, ez, xt0, yt0, zt0;
  // double len;

  whalf = (wsize - 1) / 2;

  x0 = mri_mask->width / 2;
  y0 = mri_mask->height / 2;
  MRIvoxelToTalairachVoxelEx(mri, x, y, z, &xt0, &yt0, &zt0, lta);
  switch (orientation) {
    default:
    case MRI_CORONAL: /* basis vectors in x-y plane */
      /* the 'x' basis vector */
      ex = (double)xt0 + 1;
      ey = (double)yt0;
      ez = (double)zt0;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector */
      ex = (double)xt0;
      ey = (double)yt0 + 1;
      ez = (double)zt0;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
    case MRI_HORIZONTAL: /* basis vectors in x-z plane */
      /* the 'x' basis vector */
      ex = (double)xt0 + 1;
      ey = (double)yt0;
      ez = (double)zt0;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector */
      ex = (double)xt0;
      ey = (double)yt0;
      ez = (double)zt0 + 1;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
    case MRI_SAGITTAL: /* basis vectors in y-z plane */
      /* the 'x' basis vector */
      ex = (double)xt0;
      ey = (double)yt0;
      ez = (double)zt0 + 1.0;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e1_x, &e1_y, &e1_z, lta);
      e1_x -= (double)x;
      e1_y -= (double)y;
      e1_z -= (double)z;

      /* the 'y' basis vector */
      ex = (double)xt0;
      ey = (double)yt0 + 1.0;
      ez = (double)zt0;
      MRItalairachVoxelToVoxelEx(mri, ex, ey, ez, &e2_x, &e2_y, &e2_z, lta);
      e2_x -= (double)x;
      e2_y -= (double)y;
      e2_z -= (double)z;
      break;
  }

  /*
     don't want to normalize basis - they are orthonormal in magnet space,
     not necessarily Talairach space.
  */
  // len = sqrt(e1_x * e1_x + e1_y * e1_y + e1_z * e1_z);
  /*  e1_x /= len ; e1_y /= len ; e1_z /= len ;*/
  // len = sqrt(e2_x * e2_x + e2_y * e2_y + e2_z * e2_z);
  /*  e2_x /= len ; e2_y /= len ; e2_z /= len ;*/

  for (yk = -whalf; yk <= whalf; yk++) {
    xbase = (float)x + (float)yk * e2_x;
    ybase = (float)y + (float)yk * e2_y;
    zbase = (float)z + (float)yk * e2_z;
    for (xk = -whalf; xk <= whalf; xk++) {
      /* in-plane vect. is linear combination of scaled basis vects */
      xi = mri->xi[nint(xbase + xk * e1_x)];
      yi = mri->yi[nint(ybase + xk * e1_y)];
      zi = mri->zi[nint(zbase + xk * e1_z)];
      xki = mri_mask->xi[xk + x0];
      yki = mri_mask->yi[yk + y0];
      if (MRIvox(mri_mask, xki, yki, 0)) MRIvox(mri, xi, yi, zi) = fill_val;
    }
  }
  return (NO_ERROR);
}
