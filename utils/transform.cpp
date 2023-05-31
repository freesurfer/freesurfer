/**
 * @brief utilities for linear transforms
 *
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define _TRANSFORM_SRC
#include "transform.h"
#undef _TRANSFORM_SRC

#include "diag.h"
#include "error.h"
#include "fio.h"
#include "gcamorph.h"
#include "macros.h"
#include "matrix.h"
#include "mri.h"
#include "mri_circulars.h"
#include "mrinorm.h"
#include "proto.h"
#include "registerio.h"
#include "resample.h"
#include "talairachex.h"
#include "timer.h"
#include "pointset.h"

extern const char *Progname;

#define MAX_TRANSFORMS (1024 * 4)

static LTA *ltaMNIread(const char *fname);
static LTA *ltaFSLread(const char *fname);
static int ltaMNIwrite(const LTA *lta, const char *fname);
static int ltaFSLwrite(const LTA *lta, const char *fname);
static LTA *ltaReadFile(const char *fname);

static LTA *ltaMNIreadEx(const char *fname);
static LTA *ltaReadFileEx(const char *fname);

/*!
\fn MATRIX *LTA::get_matrix(int InputCoords, int OutputCoords, MATRIX *M)
\brief Gets matrix from the LTA that will convert the given type of
input coordinate to the given type of output coordinate. If the
coordinates are of the same type, the it just runs LTAchangeType() and
returns that matrix.  Recomputes the matrix each time, ie, does not
use cached matrices. Does not change the LTA itself. The Coords must
be one of FS_COORDS_TKREG_RAS (1), FS_COORDS_SCANNER_RAS (2), and/or
FS_COORDS_VOXEL(3).
*/
MATRIX *LINEAR_TRANSFORM_ARRAY::get_matrix(int InputCoords, int OutputCoords, MATRIX *M)
{
  LTA *lta = LTAcopy(this,NULL); // so that the change types don't affect the LTA
  if(lta==NULL){
    printf("ERROR: LTA::GetMatrix: LTA is NULL\n");
    return(NULL);
  }
  if(lta->num_xforms > 1){
    printf("ERROR: LTA::GetMatrix: num_xforms > 1\n");
    return(NULL);
  }
  printf("LTA::GetMatrix: %d %d\n",InputCoords,OutputCoords);
  if(InputCoords == FS_COORDS_SCANNER_RAS && OutputCoords == FS_COORDS_SCANNER_RAS) {
    // scannerRAS to scannerRAS
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    M = MatrixCopy(lta->xforms[0].m_L,M);
  }
  else if(InputCoords == FS_COORDS_SCANNER_RAS && OutputCoords == FS_COORDS_VOXEL) {
    // scannerRAS to Vox
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    MATRIX *dst_ras2vox = lta->xforms[0].dst.get_RAS2Vox();
    M = MatrixMultiplyD(dst_ras2vox,lta->xforms[0].m_L,NULL);
    MatrixPrint(stdout,dst_ras2vox);
    MatrixPrint(stdout,lta->xforms[0].m_L);
    MatrixFree(&dst_ras2vox);
  }
  else if(InputCoords == FS_COORDS_SCANNER_RAS && OutputCoords == FS_COORDS_TKREG_RAS) {
    // scannerRAS to tkRAS
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    MATRIX *dst_ras2tkras = lta->xforms[0].dst.get_RAS2TkregRAS();
    M = MatrixMultiplyD(dst_ras2tkras,lta->xforms[0].m_L,NULL);
    MatrixFree(&dst_ras2tkras);
  }
  else if(InputCoords == FS_COORDS_VOXEL && OutputCoords == FS_COORDS_SCANNER_RAS) {
    // Vox to scannerRAS
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    MATRIX *src_vox2ras = lta->xforms[0].src.get_Vox2RAS();
    M = MatrixMultiplyD(lta->xforms[0].m_L,src_vox2ras,NULL);
    MatrixFree(&src_vox2ras);
  }
  else if(InputCoords == FS_COORDS_VOXEL && OutputCoords == FS_COORDS_VOXEL){
    // Vox to Vox
    LTAchangeType(lta,LINEAR_VOX_TO_VOX);
    M = MatrixCopy(lta->xforms[0].m_L,M);
  }
  else if(InputCoords == FS_COORDS_VOXEL && OutputCoords == FS_COORDS_TKREG_RAS){
    // Vox to tkRAS
    LTAchangeType(lta,LINEAR_VOX_TO_VOX);
    MATRIX *dst_vox2tkras = lta->xforms[0].dst.get_Vox2TkregRAS();
    M = MatrixMultiplyD(dst_vox2tkras,lta->xforms[0].m_L,NULL);
    MatrixFree(&dst_vox2tkras);
  }
  else if(InputCoords == FS_COORDS_TKREG_RAS && OutputCoords == FS_COORDS_SCANNER_RAS){
    // tkRAS to scannerRAS
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
    MATRIX *src_tkras2ras = lta->xforms[0].src.get_TkregRAS2RAS();
    M = MatrixMultiplyD(lta->xforms[0].m_L,src_tkras2ras,NULL);
    MatrixFree(&src_tkras2ras);
  }
  else if(InputCoords == FS_COORDS_TKREG_RAS && OutputCoords == FS_COORDS_VOXEL){
    // tkRAS to Vox
    LTAchangeType(lta,LINEAR_VOX_TO_VOX);
    MATRIX *src_tkras2vox = lta->xforms[0].src.get_TkregRAS2Vox();
    M = MatrixMultiplyD(lta->xforms[0].m_L,src_tkras2vox,NULL);
    MatrixFree(&src_tkras2vox);
  }
  else if(InputCoords == FS_COORDS_TKREG_RAS && OutputCoords == FS_COORDS_TKREG_RAS) {
    // tkRAS to tkRAS
    LTAchangeType(lta,REGISTER_DAT);
    M = MatrixInverse(lta->xforms[0].m_L,M); // have to invert
  }
  else {
    printf("ERROR: LTA::GetMatrix: intype=%d, outtype=%d not recognized\n",InputCoords,OutputCoords);
    LTAfree(&lta);
    return(NULL);
  }
  LTAfree(&lta);
  return(M);
}

/*
  \fn LTA *LTAcopy(LTA *lta, LTA *ltacp);
  \brief Copys lta to ltacp. If ltacp is NULL, allocs a new lta.
   If ltacp is not NULL, num_xforms must match that of lta.
 */
LTA *LTAcopy(const LTA *lta, LTA *ltacp)
{
  int i;

  if (ltacp == NULL) ltacp = LTAalloc(lta->num_xforms, NULL);
  if (lta->num_xforms != ltacp->num_xforms) {
    printf("LTAcopy(): ERROR: number of xforms does not match (%d,%d)\n", lta->num_xforms, ltacp->num_xforms);
    return (NULL);
  }

  ltacp->type = lta->type;
  ltacp->fscale = lta->fscale;
  strcpy(ltacp->subject, lta->subject);

  for (i = 0; i < lta->num_xforms; i++) {
    LTcopy(&lta->xforms[i], &ltacp->xforms[i]);
    LTcopy(&lta->inv_xforms[i], &ltacp->inv_xforms[i]);
  }
  return (ltacp);
}

/*
  \fn LINEAR_TRANSFORM *LTcopy(LT *lt, LT *ltcp)
  \brief Copys lt to ltcp. If ltcp cannot be NULL.
 */
LINEAR_TRANSFORM *LTcopy(const LT *lt, LT *ltcp)
{
  if (ltcp == NULL) {
    printf("ERROR: LTcopy() destination LT cannot be NULL\n");
    return (NULL);
  }
  ltcp->x0 = lt->x0;
  ltcp->y0 = lt->y0;
  ltcp->z0 = lt->z0;
  ltcp->sigma = lt->sigma;
  MatrixCopy(lt->m_L, ltcp->m_L);
  MatrixCopy(lt->m_dL, ltcp->m_dL);
  MatrixCopy(lt->m_last_dL, ltcp->m_last_dL);
  memcpy(&ltcp->src, &lt->src, sizeof(VOL_GEOM));
  memcpy(&ltcp->dst, &lt->dst, sizeof(VOL_GEOM));
  return (ltcp);
}

/*
  \fn int LTAdiff(LTA *lta1, LTA *lta2, double thresh)
  \brief Checks whether two LTAs are different by checking:
   num_xforms, type, src volume geometry, dst volume geometry, and
   transform matrix (m_L).  An element in the transform matrix must be
   different by more than thresh for the difference to be
   registered. Returns 1 if different, 0 if the same. Does not check
   the inverse transform, subject, or other parameters.
 */
int LTAdiff(LTA *lta1, LTA *lta2, double thresh)
{
  int i, ret, c, r, CheckInverse = 0;
  double d;

  ret = 0;
  if (lta1->num_xforms != lta2->num_xforms) {
    printf("LTAdiff(): number of xforms does not match (%d,%d)\n", lta1->num_xforms, lta2->num_xforms);
    ret = 1;
  }
  if (lta1->type != lta2->type) {
    printf("LTAdiff(): types do not match (%d,%d)\n", lta1->type, lta2->type);
    ret = 1;
  }
  if (ret) return (ret);

  for (i = 0; i < lta1->num_xforms; i++) {
    if (!vg_isEqual(&lta1->xforms[i].src, &lta2->xforms[i].src)) {
      printf("LTAdiff(): i=%d src does not match\n", i);
      ret = 1;
    }
    if (!vg_isEqual(&lta1->xforms[i].dst, &lta2->xforms[i].dst)) {
      printf("LTAdiff(): i=%d dst does not match\n", i);
      ret = 1;
    }
    for (r = 1; r <= 4; r++) {
      for (c = 1; c <= 4; c++) {
        d = (lta1->xforms[i].m_L->rptr[r][c] - lta2->xforms[i].m_L->rptr[r][c]);
        if (fabs(d) > thresh) {
          printf("LTAdiff(): i=%d L(%d,%d) does not match (%g,%g) %g\n",
                 i,
                 r,
                 c,
                 lta1->xforms[i].m_L->rptr[r][c],
                 lta2->xforms[i].m_L->rptr[r][c],
                 d);
          ret = 1;
        }
      }
    }

    if (CheckInverse) {
      if (!vg_isEqual(&lta1->inv_xforms[i].src, &lta2->inv_xforms[i].src)) {
        printf("LTAdiff(): i=%d inv src does not match\n", i);
        ret = 1;
      }
      if (!vg_isEqual(&lta1->inv_xforms[i].dst, &lta2->inv_xforms[i].dst)) {
        printf("LTAdiff(): i=%d inv dst does not match\n", i);
        ret = 1;
      }
      for (r = 1; r <= 4; r++) {
        for (c = 1; c <= 4; c++) {
          d = (lta1->inv_xforms[i].m_L->rptr[r][c] - lta2->inv_xforms[i].m_L->rptr[r][c]);
          if (fabs(d) > thresh) {
            printf("LTAdiff(): i=%d inv L(%d,%d) does not match (%g,%g) %g\n",
                   i,
                   r,
                   c,
                   lta1->inv_xforms[i].m_L->rptr[r][c],
                   lta2->inv_xforms[i].m_L->rptr[r][c],
                   d);
            ret = 1;
          }
        }
      }
    }
  }
  return (ret);
}

// this function is obsolete. Use vg->vgprint() instead.
void vg_print(const VOL_GEOM *vg)
{
  if (vg->valid == 1) {
    fprintf(stderr, "volume geometry:\n");
    fprintf(stderr, "extent  : (%d, %d, %d)\n", vg->width, vg->height, vg->depth);
    fprintf(stderr, "voxel   : (%7.4f, %7.4f, %7.4f)\n", vg->xsize, vg->ysize, vg->zsize);
    fprintf(stderr, "x_(ras) : (%7.4f, %7.4f, %7.4f)\n", vg->x_r, vg->x_a, vg->x_s);
    fprintf(stderr, "y_(ras) : (%7.4f, %7.4f, %7.4f)\n", vg->y_r, vg->y_a, vg->y_s);
    fprintf(stderr, "z_(ras) : (%7.4f, %7.4f, %7.4f)\n", vg->z_r, vg->z_a, vg->z_s);
    fprintf(stderr, "c_(ras) : (%7.4f, %7.4f, %7.4f)\n", vg->c_r, vg->c_a, vg->c_s);
    fprintf(stderr, "file    : %s\n", vg->fname);
  }
  else
    fprintf(stderr, "volume geometry: info is either not contained or not valid.\n");

  fflush(stderr);
}

// what should be the initialized value?
// I guess make it the same as COR standard.
void initVolGeom(VOL_GEOM *vg)
{
  vg->valid = 0;
  vg->width = 256;
  vg->height = 256;
  vg->depth = 256;
  vg->xsize = 1;
  vg->ysize = 1;
  vg->zsize = 1;
  vg->x_r = -1.;
  vg->x_a = 0.;
  vg->x_s = 0.;
  vg->y_r = 0.;
  vg->y_a = 0.;
  vg->y_s = -1.;
  vg->z_r = 0.;
  vg->z_a = 1.;
  vg->z_s = 0.;
  vg->c_r = 0.;
  vg->c_a = 0.;
  vg->c_s = 0.;
  strcpy(vg->fname, "unknown");  // initialized to be "unknown"
}

/*!
\fn void getVolGeom(const MRI *src, VOL_GEOM *dst)
\brief Copy Volume Geometry from MRI
 */
void getVolGeom(const MRI *src, VOL_GEOM *dst)
{
  if (!src) ErrorExit(ERROR_BADPARM, "must have a valid MRI (src)");
  if (!dst) ErrorExit(ERROR_BADPARM, "must have a valid VOL_GEOM (dst)");

  dst->valid = 1;
  dst->width = src->width;
  dst->height = src->height;
  dst->depth = src->depth;
  dst->xsize = src->xsize;
  dst->ysize = src->ysize;
  dst->zsize = src->zsize;
  dst->x_r = src->x_r;
  dst->x_a = src->x_a;
  dst->x_s = src->x_s;
  dst->y_r = src->y_r;
  dst->y_a = src->y_a;
  dst->y_s = src->y_s;
  dst->z_r = src->z_r;
  dst->z_a = src->z_a;
  dst->z_s = src->z_s;
  dst->c_r = src->c_r;
  dst->c_a = src->c_a;
  dst->c_s = src->c_s;
  strcpy(dst->fname, src->fname);
}

/*
\fn MRI *MRIallocFromVolGeom(VOL_GEOM *vg, int type, int nframes, int HeaderOnly)
\brief Creates an MRI from a VOL_GEOM, copying the geometry info
*/
MRI *MRIallocFromVolGeom(VOL_GEOM *vg, int type, int nframes, int HeaderOnly)
{
  MRI *mri;
  if (HeaderOnly)
    mri = MRIallocHeader(vg->width, vg->height, vg->depth, type, nframes);
  else
    mri = MRIallocSequence(vg->width, vg->height, vg->depth, type, nframes);
  if (mri) useVolGeomToMRI(vg, mri);
  return (mri);
}

/*
\fn void useVolGeomToMRI(const VOL_GEOM *src, MRI *dst)
\brief Copy geometry info from VOL_GEOM to an MRI structure
*/
void useVolGeomToMRI(const VOL_GEOM *src, MRI *dst)
{
  if (!src) ErrorExit(ERROR_BADPARM, "must have a valid VOL_GEOM (src)");
  if (!dst) ErrorExit(ERROR_BADPARM, "must have a valid MRI (dst)");

  dst->ras_good_flag = 1;
  dst->width = src->width;
  dst->height = src->height;
  dst->depth = src->depth;
  dst->xsize = src->xsize;
  dst->ysize = src->ysize;
  dst->zsize = src->zsize;
  dst->x_r = src->x_r;
  dst->x_a = src->x_a;
  dst->x_s = src->x_s;
  dst->y_r = src->y_r;
  dst->y_a = src->y_a;
  dst->y_s = src->y_s;
  dst->z_r = src->z_r;
  dst->z_a = src->z_a;
  dst->z_s = src->z_s;
  dst->c_r = src->c_r;
  dst->c_a = src->c_a;
  dst->c_s = src->c_s;
  strcpy(dst->fname, src->fname);
  // now we cache transform and thus we have to do the following whenever
  // we change direction cosines
  MRIreInitCache(dst);
}

int TransformCopyVolGeomToMRI(TRANSFORM *transform, MRI *mri)
{
  LTA *lta;
  GCA_MORPH *gcam;

  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *)transform->xform;
  }
  else {
    lta = (LTA *)transform->xform;
    MRIcopyVolGeomToMRI(mri, &lta->xforms[0].dst);
  }
  return (NO_ERROR);
}

void copyVolGeom(const VOL_GEOM *src, VOL_GEOM *dst)
{
  dst->valid = src->valid;
  dst->width = src->width;
  dst->height = src->height;
  dst->depth = src->depth;
  dst->xsize = src->xsize;
  dst->ysize = src->ysize;
  dst->zsize = src->zsize;
  dst->x_r = src->x_r;
  dst->x_a = src->x_a;
  dst->x_s = src->x_s;
  dst->y_r = src->y_r;
  dst->y_a = src->y_a;
  dst->y_s = src->y_s;
  dst->z_r = src->z_r;
  dst->z_a = src->z_a;
  dst->z_s = src->z_s;
  dst->c_r = src->c_r;
  dst->c_a = src->c_a;
  dst->c_s = src->c_s;
  strcpy(dst->fname, src->fname);
}

void writeVolGeom(FILE *fp, const VOL_GEOM *vg)
{
  if (vg->valid == 0)
    fprintf(fp, "valid = %d  # volume info invalid\n", vg->valid);
  else
    fprintf(fp, "valid = %d  # volume info valid\n", vg->valid);
  fprintf(fp, "filename = %s\n", vg->fname);
  fprintf(fp, "volume = %d %d %d\n", vg->width, vg->height, vg->depth);
  fprintf(fp, "voxelsize = %.15e %.15e %.15e\n", vg->xsize, vg->ysize, vg->zsize);
  fprintf(fp, "xras   = %.15e %.15e %.15e\n", vg->x_r, vg->x_a, vg->x_s);
  fprintf(fp, "yras   = %.15e %.15e %.15e\n", vg->y_r, vg->y_a, vg->y_s);
  fprintf(fp, "zras   = %.15e %.15e %.15e\n", vg->z_r, vg->z_a, vg->z_s);
  fprintf(fp, "cras   = %.15e %.15e %.15e\n", vg->c_r, vg->c_a, vg->c_s);
}

void readVolGeom(FILE *fp, VOL_GEOM *vg)
{
  char line[STRLEN + 16];
  char param[64];
  char eq[2];
  char buf[STRLEN];
  int vgRead = 0;
  char *p = 0;
  int counter = 0;
  long pos = 0;
  int fail = 0;
  while ((p = fgets(line, sizeof(line), fp)) && counter < 8) {
    if (strlen(p) == 0) break;
    sscanf(line, "%s %s %*s", param, eq);
    if (!strcmp(param, "valid")) {
      sscanf(line, "%s %s %d \n", param, eq, &vg->valid);
      vgRead = 1;
      counter++;
    }
    else if (!strcmp(param, "filename")) {
      // 1023 = STRLEN - 1
      if (sscanf(line, "%s %s %1023s%*s\n", param, eq, buf) >= 3) {
        strcpy(vg->fname, buf);
      }
      counter++;
    }
    else if (!strcmp(param, "volume")) {
      // rescan again
      sscanf(line, "%s %s %d %d %d\n", param, eq, &vg->width, &vg->height, &vg->depth);
      counter++;
    }
    else if (!strcmp(param, "voxelsize")) {
      // rescan again
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->xsize, &vg->ysize, &vg->zsize);
      counter++;
    }
    else if (!strcmp(param, "xras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->x_r, &vg->x_a, &vg->x_s);
      counter++;
    }
    else if (!strcmp(param, "yras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->y_r, &vg->y_a, &vg->y_s);
      counter++;
    }
    else if (!strcmp(param, "zras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->z_r, &vg->z_a, &vg->z_s);
      counter++;
    }
    else if (!strcmp(param, "cras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->c_r, &vg->c_a, &vg->c_s);
      counter++;
    }
    // remember the current position
    pos = ftell(fp);  // if fail = 0, then ok
  };
  if (p)  // we read one more line
  {
    if (pos > 0)                        // if success in getting pos, then
      fail = fseek(fp, pos, SEEK_SET);  // restore the position
    // note that this won't allow compression using pipe
  }
  if (!vgRead) {
    fprintf(stderr, "INFO: volume info was not present.\n");
    initVolGeom(vg);
  }
}

// scanner space vox2ras from vol geom
// call to this function is equivalent to vg->get_Vox2RAS();
MATRIX *vg_i_to_r(const VOL_GEOM *vg)
{
  MATRIX *mat = extract_i_to_r(vg);
  return mat;
}
// call to this function is equivalent to vg->get_RAS2Vox();
MATRIX *vg_r_to_i(const VOL_GEOM *vg)
{
  MATRIX *mat = extract_r_to_i(vg);
  return mat;
}
// tkregister space vox2ras from vol geom
// call to this function is equivalent to vg->get_Vox2TkregRAS();
MATRIX *TkrVox2RASfromVolGeom(const VOL_GEOM *vg)
{
  MATRIX *mat = MRIxfmCRS2XYZtkreg(vg);
  return (mat);
}
// tkregister space ras2vox from vol geom
// call to this function is equivalent to vg->get_TkregRAS2Vox();
MATRIX *TkrRAS2VoxfromVolGeom(const VOL_GEOM *vg)
{
  MATRIX *mat = NULL;
  mat = TkrVox2RASfromVolGeom(vg);
  mat = MatrixInverse(mat, mat);
  return (mat);
}
/*
 \fn MATRIX *VGtkreg2RAS(VOL_GEOM *vg, MATRIX *tkreg2ras)
 \brief Returns a matrix that converts from tkreg coords to ras
*/
MATRIX *VGtkreg2RAS(VOL_GEOM *vg, MATRIX *tkreg2ras)
{
  // Create a dummy MRI
  MRI *mri = MRIallocFromVolGeom(vg, MRI_INT, 1, 1);
  // Extract the matrix from the MRI
  tkreg2ras = RASFromSurfaceRAS_(mri,tkreg2ras);
  MRIfree(&mri);
  return(tkreg2ras);
}
/*
 \fn MATRIX *VGras2tkreg(VOL_GEOM *vg, MATRIX *ras2tkreg)
 \brief Returns a matrix that converts from ras coords to tkreg coords
*/
MATRIX *VGras2tkreg(VOL_GEOM *vg, MATRIX *ras2tkreg)
{
  MATRIX *M = VGtkreg2RAS(vg, NULL);
  MATRIX *P = MatrixInverse(M,NULL);
  MatrixFree(&M);
  return(P);
}

int vg_isEqual(const VOL_GEOM *vg1, const VOL_GEOM *vg2)
{
  int rt;
  extern double vg_isEqual_Threshold;
  // rt = vg_isEqualThresh(vg1, vg2, FLT_EPSILON);
  rt = vg_isNotEqualThresh(vg1, vg2, vg_isEqual_Threshold);
  if (rt == 0)
    return (1);
  else
    return (0);
}

int vg_isNotEqualThresh(const VOL_GEOM *vg1, const VOL_GEOM *vg2, const double thresh)
{
  if (vg1->valid != vg2->valid) return (1);
  if (vg1->width != vg2->width) return (2);
  if (vg1->height != vg2->height) return (3);
  if (vg1->depth != vg2->depth) return (4);
  if (!FZEROTHR(vg1->xsize - vg2->xsize, thresh)) return (5);
  if (!FZEROTHR(vg1->ysize - vg2->ysize, thresh)) return (6);
  if (!FZEROTHR(vg1->zsize - vg2->zsize, thresh)) return (7);
  if (!FZEROTHR(vg1->x_r - vg2->x_r, thresh)) return (8);
  if (!FZEROTHR(vg1->x_a - vg2->x_a, thresh)) return (9);
  if (!FZEROTHR(vg1->x_s - vg2->x_s, thresh)) return (10);
  if (!FZEROTHR(vg1->y_r - vg2->y_r, thresh)) return (11);
  if (!FZEROTHR(vg1->y_a - vg2->y_a, thresh)) return (12);
  if (!FZEROTHR(vg1->y_s - vg2->y_s, thresh)) return (13);
  if (!FZEROTHR(vg1->z_r - vg2->z_r, thresh)) return (14);
  if (!FZEROTHR(vg1->z_a - vg2->z_a, thresh)) return (15);
  if (!FZEROTHR(vg1->z_s - vg2->z_s, thresh)) return (16);
  if (!FZEROTHR(vg1->c_r - vg2->c_r, thresh)) return (17);
  if (!FZEROTHR(vg1->c_a - vg2->c_a, thresh)) return (18);
  if (!FZEROTHR(vg1->c_s - vg2->c_s, thresh)) return (19);
  return (0);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
LINEAR_TRANSFORM_ARRAY *LTAalloc(int nxforms, MRI *mri)
{
  int i;
  LINEAR_TRANSFORM_ARRAY *lta;
  MRI_REGION bbox;
  float x0, y0, z0;

  if (mri) {
    MRIboundingBox(mri, 70, &bbox);
    x0 = bbox.x + bbox.dx / 2;
    y0 = bbox.y + bbox.dy / 2;
    z0 = bbox.z + bbox.dz / 2;
  }
  else
    x0 = y0 = z0 = 0;

  lta = (LINEAR_TRANSFORM_ARRAY *)calloc(1, sizeof(LTA));
  if (!lta) ErrorExit(ERROR_NOMEMORY, "LTAalloc(%d): could not allocate LTA", nxforms);
  lta->num_xforms = nxforms;
  lta->xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT));
  if (!lta->xforms) ErrorExit(ERROR_NOMEMORY, "LTAalloc(%d): could not allocate xforms", nxforms);
  lta->inv_xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT));
  if (!lta->inv_xforms) ErrorExit(ERROR_NOMEMORY, "LTAalloc(%d): could not allocate inverse xforms", nxforms);
  for (i = 0; i < nxforms; i++) {
    lta->xforms[i].x0 = x0;
    lta->xforms[i].y0 = y0;
    lta->xforms[i].z0 = z0;
    lta->xforms[i].sigma = 10000.0f;
    lta->xforms[i].m_L = MatrixIdentity(4, NULL);
    lta->xforms[i].m_dL = MatrixAlloc(4, 4, MATRIX_REAL);
    lta->xforms[i].m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL);
    initVolGeom(&lta->xforms[i].src);
    initVolGeom(&lta->xforms[i].dst);
    lta->xforms[i].type = UNKNOWN;
    ;

    lta->inv_xforms[i].x0 = x0;
    lta->inv_xforms[i].y0 = y0;
    lta->inv_xforms[i].z0 = z0;
    lta->inv_xforms[i].sigma = 10000.0f;
    lta->inv_xforms[i].m_L = MatrixIdentity(4, NULL);
    lta->inv_xforms[i].m_dL = MatrixAlloc(4, 4, MATRIX_REAL);
    lta->inv_xforms[i].m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL);
    initVolGeom(&lta->inv_xforms[i].src);
    initVolGeom(&lta->inv_xforms[i].dst);
    lta->inv_xforms[i].type = UNKNOWN;
    ;
  }
  return (lta);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int LTAwrite(LTA *lta, const char *fname)
{
  return (LTAwriteEx(lta, fname));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

LTA *LTAread(const char *fname)
{
  int type;
  LTA *lta;
  MATRIX *V, *W, *m_tmp;

  return (LTAreadEx(fname));  // no reason not to always use it
  type = TransformFileNameType(fname);
  switch (type) {
    case REGISTER_DAT:
      lta = ltaReadRegisterDat(fname, NULL, NULL);
      if (!lta) return (NULL);

      lta->type = LINEAR_CORONAL_RAS_TO_CORONAL_RAS;
      break;
    case MNI_TRANSFORM_TYPE:
      lta = ltaMNIread(fname);
      if (!lta) return (NULL);

      /* by default convert MNI files to voxel coords.
        Sorry, I know this shouldn't be done here, particularly since we
        don't know enough to convert to scanner RAS coords, but I don't want
        to risk breaking the Talairach code by mucking around with it (BRF).
      */
      /* convert to voxel coords */
      V = MatrixAlloc(4, 4, MATRIX_REAL); /* world to voxel transform */
      W = MatrixAlloc(4, 4, MATRIX_REAL); /* voxel to world transform */
      *MATRIX_RELT(V, 1, 1) = -1;
      *MATRIX_RELT(V, 1, 4) = 128;
      *MATRIX_RELT(V, 2, 3) = -1;
      *MATRIX_RELT(V, 2, 4) = 128;
      *MATRIX_RELT(V, 3, 2) = 1;
      *MATRIX_RELT(V, 3, 4) = 128;
      *MATRIX_RELT(V, 4, 4) = 1;

      *MATRIX_RELT(W, 1, 1) = -1;
      *MATRIX_RELT(W, 1, 4) = 128;
      *MATRIX_RELT(W, 2, 3) = 1;
      *MATRIX_RELT(W, 2, 4) = -128;
      *MATRIX_RELT(W, 3, 2) = -1;
      *MATRIX_RELT(W, 3, 4) = 128;
      *MATRIX_RELT(W, 4, 4) = 1;

      m_tmp = MatrixMultiply(lta->xforms[0].m_L, W, NULL);
      MatrixMultiply(V, m_tmp, lta->xforms[0].m_L);
      MatrixFree(&V);
      MatrixFree(&W);
      MatrixFree(&m_tmp);
      lta->type = LINEAR_VOX_TO_VOX;
      break;
    case LINEAR_VOX_TO_VOX:
    case LINEAR_RAS_TO_RAS:
    case TRANSFORM_ARRAY_TYPE:
    default:
      lta = ltaReadFile(fname);
      break;
  }
  return (lta);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static LTA *ltaReadFile(const char *fname)
{
  FILE *fp;
  LINEAR_TRANSFORM *lt;
  int i, nxforms, type, skip = 0;
  char line[STRLEN], *cp;
  LTA *lta;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't open file", fname));
  cp = fgetl(line, STRLEN - 1, fp);
  if (cp == NULL) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data", fname));
  }
  sscanf(cp, "type      = %d\n", &type);
  cp = fgetl(line, STRLEN - 1, fp);
  sscanf(cp, "nxforms   = %d\n", &nxforms);
  lta = LTAalloc(nxforms, NULL);
  lta->type = type;
  for (i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    if (skip == 0) cp = fgetl(line, STRLEN - 1, fp);
    sscanf(cp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0);
    cp = fgetl(line, STRLEN - 1, fp);
    sscanf(cp, "sigma     = %f\n", &lt->sigma);
    MatrixAsciiReadFrom(fp, lt->m_L);
    cp = fgetl(line, STRLEN - 1, fp);
    if (strncmp(cp, "label", 5) == 0)  // not all files have the label tag
    {
      sscanf(cp, "label     = %d\n", &lt->label);
      skip = 0;
    }
    else
      skip = 1;
  }
  fclose(fp);
  return (lta);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int LTAfree(LTA **plta)
{
  LTA *lta;
  int i;

  lta = *plta;
  *plta = NULL;
  for (i = 0; i < lta->num_xforms; i++) {
    MatrixFree(&lta->xforms[i].m_L);
    MatrixFree(&lta->xforms[i].m_dL);
    MatrixFree(&lta->xforms[i].m_last_dL);

    MatrixFree(&lta->inv_xforms[i].m_L);
    MatrixFree(&lta->inv_xforms[i].m_dL);
    MatrixFree(&lta->inv_xforms[i].m_last_dL);
  }
  free(lta->xforms);
  free(lta->inv_xforms);
  free(lta);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  split each transform into 8 and space them
  evenly.
  ------------------------------------------------------*/
int LTAdivide(LTA *lta, MRI *mri)
{
  MRI_REGION bbox;
  int oldi, i, nxforms, row_size, x, y, z;
  LT *new_xforms, *lt;
  float sigma, dx, dy, dz, len, x0, y0, z0;

  MRIboundingBox(mri, 130, &bbox);
  dx = bbox.dx;
  dy = bbox.dy;
  dz = bbox.dz;
  nxforms = lta->num_xforms * 8;
  len = pow((double)nxforms, 1.0 / 3.0); /* # along each dimension */
  row_size = nint(len);                  /* # of nodes in each dimension */
  dx = dx / (len + 1);                   /* x spacing between nodes */
  dy = dy / (len + 1);                   /* y spacing between nodes */
  dz = dz / (len + 1);                   /* z spacing between nodes */
  new_xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT));

  sigma = 1.0 * (dx + dy + dz) / 3.0f; /* average node spacing */

  /*
    distribute the new transforms uniformly throughout the volume.
    Next, for each new transform, calculate the old transform at
    that point, and use it as the initial value for the new transform.
  */
  for (i = x = 0; x < row_size; x++) {
    x0 = dx / 2 + x * dx;
    for (y = 0; y < row_size; y++) {
      y0 = dy / 2 + y * dy;
      for (z = 0; z < row_size; z++, i++) {
        z0 = dz / 2 + z * dz;
        lt = &new_xforms[i];
        lt->x0 = x0;
        lt->y0 = y0;
        lt->z0 = z0;
        lt->sigma = sigma;
        lt->m_L = LTAtransformAtPoint(lta, x0, y0, z0, NULL);
        lt->m_dL = MatrixAlloc(4, 4, MATRIX_REAL);
        lt->m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL);
      }
    }
  }

  /* free old ones */
  for (oldi = 0; oldi < lta->num_xforms; oldi++) {
    MatrixFree(&lta->xforms[oldi].m_L);
    MatrixFree(&lta->xforms[oldi].m_dL);
    MatrixFree(&lta->xforms[oldi].m_last_dL);
  }
  free(lta->xforms);

  /* update lta structure with new info */
  lta->xforms = new_xforms;
  lta->num_xforms = nxforms;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Apply a transform array to an image.
  ------------------------------------------------------*/
MRI *LTAtransform(MRI *mri_src, MRI *mri_dst, LTA *lta)
{
  return LTAtransformInterp(mri_src, mri_dst, lta, SAMPLE_TRILINEAR);
}

MRI *LTAtransformInterp(MRI *mri_src, MRI *mri_dst, LTA *lta, int interp)
{
  int y1, y2, y3, width, height, depth, xi, yi, zi, f;
  VECTOR *v_X, *v_Y; /* original and transformed coordinate systems */
  double x1, x2, x3;
  MATRIX *m_L, *m_L_inv;
  LT *tran = &lta->xforms[0];
  MATRIX *r2i = 0;
  MATRIX *i2r = 0;
  MATRIX *tmp = 0;
  // MATRIX *v2v = 0;
  MRI *resMRI = 0;

  if (lta->type == REGISTER_DAT) {
    printf("warning: changing input transform type from REGISTER_DAT to VOX_TO_VOX\n");
    LTAchangeType(lta, LINEAR_VOX_TO_VOX);
  }
  if (lta->num_xforms == 1) {
    /////////////////////////////////////////////////////////////////////////
    //  The transform was created using src and dst
    //           src vox ----> RAS
    //            |             |
    //            V             V
    //           dst vox ----> RAS
    //
    //  Note: RAS is the physical space.  vox is the way to
    // embed voxel in physical space
    //
    //  You can take arbitray embedding of the dst but keeping the RAS same.
    //
    //           src vox ----> RAS
    //            | v2v         | r2r
    //            V             V
    //           dst vox ----> RAS
    //            |             || identity
    //            V             ||
    //           dst' vox ---> RAS
    //
    //  where dst->dst' map is given by V2V' = r2i(dst')*i2r(dst)
    //  (here dst is from the lta, and dst' from mri_dst)
    //
    //  Note that in order to obtain src->dst' with r2r, you
    // "don't need any info from dst (only from dst')"
    //  since
    //          src->dst'  = r2i(dst')* r2r * i2r(src)
    //
    //  However, using v2v, you need the information from dst
    //  since
    //          src->dst'  = r2i(dst')* i2r(dst) * v2v
    //
    //  Similarly src geometry can be different, therefore we convert
    //  any passed v2v lta into r2r first and apply it on possible different
    //  image geometries src and dst.
    //
    ////////////////////////////////////////////////////////////////////////

    // when the mri_dst volume is not given
    if (!mri_dst) {
      if (tran->dst.valid == 1)  // transform dst is valid
      {
        // modify dst geometry using the transform dst value,
        // i.e. put the head in the same position in the dst'
        // volume as the transform dst was
        if (DIAG_VERBOSE_ON)
          fprintf(stderr,
                  "INFO: Modifying dst geometry, "
                  "using the transform dst\n");
        // allocate dst space (take type from src and geometry from transform):
        mri_dst = MRIallocSequence(tran->dst.width, tran->dst.height, tran->dst.depth, mri_src->type, mri_src->nframes);
        // cp rest of header information from src:
        MRIcopyHeader(mri_src, mri_dst);
        // make sure the geometry is taken from the transform, not from src:
        useVolGeomToMRI(&tran->dst, mri_dst);
      }
      else if (getenv("USE_AVERAGE305")) {
        fprintf(stderr,
                "INFO: Environmental variable "
                "USE_AVERAGE305 set\n");
        fprintf(stderr,
                "INFO: Modifying dst c_(r,a,s), "
                "using average_305 values\n");
        // use the same volume size as the src
        mri_dst = MRIclone(mri_src, NULL);
        // reset talairach transform file name:
        mri_dst->transform_fname[0] = '\0';
        mri_dst->c_r = -0.0950;
        mri_dst->c_a = -16.5100;
        mri_dst->c_s = 9.7500;
        mri_dst->ras_good_flag = 1;
        // maye one should set also the other geometry entries
        // from the average ???
        //
        // now we cache transform and thus we have to
        // do the following whenever
        // we change direction cosines
        MRIreInitCache(mri_dst);
      }
      else {
        fprintf(stderr,
                "INFO: Transform dst volume "
                "info is not used (valid flag = 0).\n");
        // use the same volume size as the src
        mri_dst = MRIclone(mri_src, NULL);
        // reset talairach transform file name:
        mri_dst->transform_fname[0] = '\0';
        // maybe also reset or concatenate the actual transform
        // if available in mri_dst->transform (not yet implemented) ...
      }
    }

    ////////////////////////////////////////////////////////////////////////
    if (lta->type == LINEAR_RAS_TO_RAS) {
      // don't need any info from dst(in lta) only from mri_dst:
      return (MRIapplyRASlinearTransformInterp(mri_src, mri_dst, lta->xforms[0].m_L, interp));
    }
    else if (lta->type == LINEAR_VOX_TO_VOX)  // vox-to-vox
    {
      // convert to ras_to_ras using xforms from lta if available
      // this strips geometry information and allows to use the
      // lta on volumes with different geometries
      // if lta geomery information is missing (should not happen)
      // we use the geometry from the passed source and target image
      if (lta->xforms[0].dst.valid) {
        i2r = vg_i_to_r(&lta->xforms[0].dst);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA dst geometry information missing!\n"
                "      We assume that the dst volume passed is the\n"
                "      same as the dst for the transform.\n");
        i2r = extract_i_to_r(mri_dst);
      }
      tmp = MatrixMultiply(i2r, lta->xforms[0].m_L, NULL);
      if (lta->xforms[0].src.valid) {
        r2i = vg_r_to_i(&lta->xforms[0].src);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA src geometry information missing!\n"
                "      We assume that the src volume passed is the\n"
                "      same as the src for the transform.\n");
        r2i = extract_r_to_i(mri_src);
      }
      tmp = MatrixMultiply(tmp, r2i, tmp);

      resMRI = MRIapplyRASlinearTransformInterp(mri_src, mri_dst, tmp, interp);
      MatrixFree(&i2r);
      MatrixFree(&r2i);
      MatrixFree(&tmp);
      return resMRI;

      // old code only treated different dst not different src geometries
      /*if (lta->xforms[0].dst.valid)
      {
        i2r = vg_i_to_r(&lta->xforms[0].dst); // allocated
        r2i = extract_r_to_i(mri_dst);
        tmp = MatrixMultiply(i2r, lta->xforms[0].m_L, NULL);
        v2v = MatrixMultiply(r2i, tmp, NULL);
        resMRI = MRIlinearTransformInterp(mri_src, mri_dst, v2v, interp);
        MatrixFree(&v2v);
        v2v = 0;
        MatrixFree(&i2r);
        i2r = 0;
        MatrixFree(&r2i);
        r2i = 0;
        MatrixFree(&tmp);
        tmp = 0;
        return resMRI;
      }
      else
      {
        fprintf(stderr, "INFO: assumes that the dst volume "
                "given is the same as the dst for the transform\n");
        return(MRIlinearTransformInterp(mri_src, mri_dst,
                                        lta->xforms[0].m_L, interp)) ;
      }*/
    }
    else if (lta->type == LINEAR_PHYSVOX_TO_PHYSVOX) {
      // must have both transform src and dst geometry information
      LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      return (MRIapplyRASlinearTransformInterp(mri_src, mri_dst, lta->xforms[0].m_L, interp));
    }
    else
      ErrorExit(ERROR_BADPARM, "LTAtransform: unknown linear transform\n");
  }
  fprintf(stderr, "applying octree transform to image...\n");
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL);
  v_Y->rptr[4][1] = 1.0f;

  if (lta->num_xforms == 1) {
    LTAtransformAtPoint(lta, 0, 0, 0, m_L);
    if (MatrixInverse(m_L, m_L_inv) == NULL) {
      MatrixFree(&m_L);
      MatrixFree(&m_L_inv);
      ErrorReturn(NULL, (ERROR_BADPARM, "LTAtransform: could not invert matrix"));
    }
  }
  for (y3 = 0; y3 < depth; y3++) {
    DiagHeartbeat((float)y3 / (float)(depth - 1));
    V3_Z(v_Y) = y3;

    for (y2 = 0; y2 < height; y2++) {
      V3_Y(v_Y) = y2;
      for (y1 = 0; y1 < width; y1++) {
        V3_X(v_Y) = y1;

        /*
          this is not quite right. Should use the weighting at
          the point that maps to (y1,y2,y3), not (y1,y2,y3) itself,
          but this is much easier....
        */
        if (lta->num_xforms > 1) {
          LTAtransformAtPoint(lta, y1, y2, y3, m_L);
          if (MatrixInverse(m_L, m_L_inv) == NULL) continue;
        }
        MatrixMultiply(m_L_inv, v_Y, v_X);
        x1 = V3_X(v_X);
        x2 = V3_Y(v_X);
        x3 = V3_Z(v_X);

        xi = mri_dst->xi[nint(x1)];
        yi = mri_dst->yi[nint(x2)];
        zi = mri_dst->zi[nint(x3)];
        for (f = 0; f < mri_dst->nframes; f++)
          MRIsetVoxVal(mri_dst, y1, y2, y3, f, MRIgetVoxVal(mri_src, xi, yi, zi, f));
      }
    }
  }

  MatrixFree(&v_X);
  MatrixFree(&v_Y);
  MatrixFree(&m_L);
  MatrixFree(&m_L_inv);

  return (mri_dst);
}
MRI *LTAinverseTransformInterp(MRI *mri_src, MRI *mri_dst, LTA *lta, int interp)
{
  int y1, y2, y3, width, height, depth, xi, yi, zi, f;
  VECTOR *v_X, *v_Y; /* original and transformed coordinate systems */
  double x1, x2, x3;
  MATRIX *m_L, *m_L_inv;
  LT *tran = &lta->inv_xforms[0];
  MATRIX *r2i = 0;
  MATRIX *i2r = 0;
  MATRIX *tmp = 0;
  // MATRIX *v2v = 0;
  MRI *resMRI = 0;

  if (lta->num_xforms == 1) {
    /////////////////////////////////////////////////////////////////////////
    //  The transform was created using src and dst
    //           src vox ----> RAS
    //            |             |
    //            V             V
    //           dst vox ----> RAS
    //
    //  Note: RAS is the physical space.  vox is the way to
    // embed voxel in physical space
    //
    //  You can take arbitray embedding of the dst but keeping the RAS same.
    //
    //           src vox ----> RAS
    //            | v2v         | r2r
    //            V             V
    //           dst vox ----> RAS
    //            |             || identity
    //            V             ||
    //           dst' vox ---> RAS
    //
    //  where dst->dst' map is given by V2V' = r2i(dst')*i2r(dst)
    //
    //  Note that in order to obtain src->dst' with r2r, you
    // "don't need any info from dst"
    //  since
    //          src->dst'  = r2i(dst')* r2r * i2r(src)
    //
    //  However, using v2v, you need the information from dst
    //  since
    //          src->dst'  = r2i(dst')* i2r(dst) * v2v
    //
    ////////////////////////////////////////////////////////////////////////
    // when the dst volume is not given
    if (!mri_dst) {
      // use the same volume size as the src
      mri_dst = MRIclone(mri_src, NULL);
      if (tran->dst.valid == 1)  // transform dst is valid
      {
        // modify dst c_(r,a,s) using the transform dst value
        // to make the better positioning (i.e. put the
        // head in the same position in
        // the volume as the dst was)
        if (DIAG_VERBOSE_ON)
          fprintf(stderr,
                  "INFO: Modifying dst c_(r,a,s), "
                  "using the transform dst\n");
        mri_dst->c_r = tran->dst.c_r;
        mri_dst->c_a = tran->dst.c_a;
        mri_dst->c_s = tran->dst.c_s;
        mri_dst->ras_good_flag = 1;
        // now we cache transform and thus we have to do
        // the following whenever
        // we change direction cosines
        MRIreInitCache(mri_dst);
      }
      else if (getenv("USE_AVERAGE305")) {
        fprintf(stderr,
                "INFO: Environmental variable "
                "USE_AVERAGE305 set\n");
        fprintf(stderr,
                "INFO: Modifying dst c_(r,a,s), "
                "using average_305 values\n");
        mri_dst->c_r = -0.0950;
        mri_dst->c_a = -16.5100;
        mri_dst->c_s = 9.7500;
        mri_dst->ras_good_flag = 1;
        // now we cache transform and thus we have to
        // do the following whenever
        // we change direction cosines
        MRIreInitCache(mri_dst);
      }
      else
        fprintf(stderr,
                "INFO: Transform dst volume "
                "info is not used (valid flag = 0).\n");
    }
    ////////////////////////////////////////////////////////////////////////
    if (lta->type == LINEAR_RAS_TO_RAS) {
      // don't need any info from dst
      return (MRIapplyRASlinearTransformInterp(mri_src, mri_dst, lta->inv_xforms[0].m_L, interp));
    }
    else if (lta->type == LINEAR_VOX_TO_VOX)  // vox-to-vox
    {
      // convert to ras_to_ras using xforms from lta if available
      // this strips geometry information and allows to use the
      // lta on volumes with different geometries
      // if lta geomery information is missing (should not happen)
      // we use the geometry from the passed source and target image
      // WARNING: adoped from above, not tested here (using inv)
      if (lta->inv_xforms[0].dst.valid) {
        i2r = vg_i_to_r(&lta->inv_xforms[0].dst);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA inv dst geometry information missing!\n"
                "      We assume that the dst volume passed is the\n"
                "      same as the dst for the transform.\n");
        i2r = extract_i_to_r(mri_dst);
      }
      tmp = MatrixMultiply(i2r, lta->inv_xforms[0].m_L, NULL);
      if (lta->inv_xforms[0].src.valid) {
        r2i = vg_r_to_i(&lta->inv_xforms[0].src);  // allocated
      }
      else {
        fprintf(stderr,
                "INFO: LTA inv src geometry information missing!\n"
                "      We assume that the src volume passed is the\n"
                "      same as the src for the transform.\n");
        r2i = extract_r_to_i(mri_src);
      }
      tmp = MatrixMultiply(tmp, r2i, tmp);

      resMRI = MRIapplyRASlinearTransformInterp(mri_src, mri_dst, tmp, interp);
      MatrixFree(&i2r);
      MatrixFree(&r2i);
      MatrixFree(&tmp);
      return resMRI;

      // old code only treated different dst not different src geometries
      /*if (lta->xforms[0].dst.valid)
      {
        i2r = vg_i_to_r(&lta->inv_xforms[0].dst); // allocated
        r2i = extract_r_to_i(mri_dst);
        tmp = MatrixMultiply(i2r, lta->inv_xforms[0].m_L, NULL);
        v2v = MatrixMultiply(r2i, tmp, NULL);
        resMRI = MRIlinearTransformInterp(mri_src, mri_dst, v2v, interp);
        MatrixFree(&v2v);
        v2v = 0;
        MatrixFree(&i2r);
        i2r = 0;
        MatrixFree(&r2i);
        r2i = 0;
        MatrixFree(&tmp);
        tmp = 0;
        return resMRI;
      }
      else
      {
        fprintf(stderr, "INFO: assumes that the dst volume "
                "given is the same as the dst for the transform\n");
        return(MRIlinearTransformInterp(mri_src, mri_dst,
                                        lta->inv_xforms[0].m_L, interp)) ;
      }*/
    }
    else if (lta->type == LINEAR_PHYSVOX_TO_PHYSVOX) {
      // must have both transform src and dst geometry information
      LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      return (MRIapplyRASlinearTransformInterp(mri_src, mri_dst, lta->inv_xforms[0].m_L, interp));
    }
    else
      ErrorExit(ERROR_BADPARM, "LTAtransform: unknown linear transform\n");
  }
  fprintf(stderr, "applying octree transform to image...\n");
  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  v_X = VectorAlloc(4, MATRIX_REAL); /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL); /* transformed (dst) coordinates */
  m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL);
  v_Y->rptr[4][1] = 1.0f;

  if (lta->num_xforms == 1) {
    LTAtransformAtPoint(lta, 0, 0, 0, m_L_inv);
    if (m_L_inv == NULL) {
      ErrorReturn(NULL, (ERROR_BADPARM, "LTAinverseTransform: could not invert matrix"));
    }
  }
  for (y3 = 0; y3 < depth; y3++) {
    DiagHeartbeat((float)y3 / (float)(depth - 1));
    V3_Z(v_Y) = y3;

    for (y2 = 0; y2 < height; y2++) {
      V3_Y(v_Y) = y2;
      for (y1 = 0; y1 < width; y1++) {
        V3_X(v_Y) = y1;

        /*
          this is not quite right. Should use the weighting at
          the point that maps to (y1,y2,y3), not (y1,y2,y3) itself,
          but this is much easier....
        */
        if (lta->num_xforms > 1) {
          LTAtransformAtPoint(lta, y1, y2, y3, m_L);
          if (MatrixInverse(m_L, m_L_inv) == NULL) continue;
        }
        MatrixMultiply(m_L_inv, v_Y, v_X);
        x1 = V3_X(v_X);
        x2 = V3_Y(v_X);
        x3 = V3_Z(v_X);

        xi = mri_dst->xi[nint(x1)];
        yi = mri_dst->yi[nint(x2)];
        zi = mri_dst->zi[nint(x3)];
        for (f = 0; f < mri_dst->nframes; f++)
          MRIsetVoxVal(mri_dst, y1, y2, y3, f, MRIgetVoxVal(mri_src, xi, yi, zi, f));
      }
    }
  }

  MatrixFree(&v_X);
  MatrixFree(&v_Y);
  MatrixFree(&m_L);
  MatrixFree(&m_L_inv);

  return (mri_dst);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MATRIX *LTAtransformAtPoint(LTA *lta, float x, float y, float z, MATRIX *m_L)
{
  LT *lt;
  int i;
  double w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin;
  static MATRIX *m_tmp = NULL;

  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  else
    MatrixClear(m_L);

  if (m_tmp == NULL) m_tmp = MatrixAlloc(4, 4, MATRIX_REAL);

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    dx = lta->xforms[i].x0 - x;
    dy = lta->xforms[i].y0 - y;
    dz = lta->xforms[i].z0 - z;
    dsq = (dx * dx + dy * dy + dz * dz);
    sigma = lt->sigma;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq / (2 * sigma * sigma));
    wtotal += w_p[i];
  }

  if (DZERO(wtotal)) /* no transforms in range??? */
    MatrixCopy(lta->xforms[0].m_L, m_L);
    else /* now calculate linear combination of transforms at this point */
    {
      wmin = 0.1 / (double)lta->num_xforms;
      for (i = 0; i < lta->num_xforms; i++) {
        lt = &lta->xforms[i];
        w_k_p = w_p[i] / wtotal;
        if (w_k_p < wmin) /* optimization - ignore this transform */
          continue;
        MatrixScalarMul(lt->m_L, w_k_p, m_tmp);
        MatrixAdd(m_L, m_tmp, m_L);
      }
    }
  return (m_L);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MATRIX *LTAinverseTransformAtPoint(LTA *lta, float x, float y, float z, MATRIX *m_L)
{
  LT *lt;
  int i;
  double w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin;
  static MATRIX *m_tmp = NULL;
  MATRIX *m_inv;

  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  else
    MatrixClear(m_L);

  if (m_tmp == NULL) m_tmp = MatrixAlloc(4, 4, MATRIX_REAL);

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    dx = lta->xforms[i].x0 - x;
    dy = lta->xforms[i].y0 - y;
    dz = lta->xforms[i].z0 - z;
    dsq = (dx * dx + dy * dy + dz * dz);
    sigma = lt->sigma;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq / (2 * sigma * sigma));
    wtotal += w_p[i];
  }

  if (lta->num_xforms == 1) {
    m_L = MatrixInverse(lta->xforms[0].m_L, NULL);
  }
  else {
    if (DZERO(wtotal)) /* no transforms in range??? */
      MatrixIdentity(4, m_L);
    else /* now calculate linear combination of transforms at this point */
    {
      wmin = 0.1 / (double)lta->num_xforms;
      for (i = 0; i < lta->num_xforms; i++) {
        lt = &lta->xforms[i];
        w_k_p = w_p[i] / wtotal;
        if (w_k_p < wmin) /* optimization - ignore this transform */
          continue;

        m_inv = MatrixInverse(lt->m_L, NULL);
        if (!m_inv) continue;
        MatrixScalarMul(m_inv, w_k_p, m_tmp);
        MatrixAdd(m_L, m_tmp, m_L);
        MatrixFree(&m_inv);
      }
    }
  }
  return (m_L);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
VECTOR *LTAtransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  static MATRIX *m_L = NULL;

  m_L = LTAtransformAtPoint(lta, V3_X(v_X), V3_Y(v_X), V3_Z(v_X), m_L);
  v_Y = MatrixMultiply(m_L, v_X, v_Y);
  return (v_Y);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
VECTOR *LTAinverseTransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  static MATRIX *m_L = NULL;

  m_L = LTAinverseTransformAtPoint(lta, V3_X(v_X), V3_Y(v_X), V3_Z(v_X), m_L);
  v_Y = MatrixMultiply(m_L, v_X, v_Y);
  return (v_Y);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Same as LTAtransformPoint but also return weight normalization
  factor.
  ------------------------------------------------------*/
double LTAtransformPointAndGetWtotal(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  LT *lt;
  int i;
  double w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin, x, y, z;
  static MATRIX *m_tmp = NULL, *m_L = NULL;

  x = V3_X(v_X);
  y = V3_Y(v_X);
  z = V3_Z(v_X);
  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL);
  else
    MatrixClear(m_L);

  if (m_tmp == NULL) m_tmp = MatrixAlloc(4, 4, MATRIX_REAL);

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    dx = lta->xforms[i].x0 - x;
    dy = lta->xforms[i].y0 - y;
    dz = lta->xforms[i].z0 - z;
    dsq = (dx * dx + dy * dy + dz * dz);
    sigma = lt->sigma;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq / (2 * sigma * sigma));
    wtotal += w_p[i];
  }

  if (DZERO(wtotal)) /* no transforms in range??? */
    MatrixIdentity(4, m_L);
  else /* now calculate linear combination of transforms at this point */
  {
    wmin = 0.1 / (double)lta->num_xforms;
    for (i = 0; i < lta->num_xforms; i++) {
      lt = &lta->xforms[i];
      w_k_p = w_p[i] / wtotal;
      if (w_k_p < wmin) /* optimization - ignore this transform */
        continue;
      MatrixScalarMul(lt->m_L, w_k_p, m_tmp);
      MatrixAdd(m_L, m_tmp, m_L);
    }
  }

  MatrixMultiply(m_L, v_X, v_Y);
  return (wtotal);
}
/*-----------------------------------------------------*/
int TransformFileNameType(const char *fname)
{
  int file_type = TRANSFORM_ARRAY_TYPE;
  char *dot, buf[500], *number;

  strcpy(buf, fname);
  dot = strrchr(buf, '@');
  number = strchr(buf, '#');
  if (number) *number = 0; /* don't consider : part of extension */

  if (!dot) dot = strrchr(buf, '.');

  if (dot) {
    dot++;
    StrUpper(buf);
    if (!strcmp(dot, "M3D"))
      return (MORPH_3D_TYPE);
    else if (!strcmp(dot, "M3Z"))
      return (MORPH_3D_TYPE);
    else if (!strcmp(dot, "OCT"))
      return (TRANSFORM_ARRAY_TYPE);
    else if (!strcmp(dot, "XFM"))
      return (MNI_TRANSFORM_TYPE);
    else if (!strcmp(dot, "FSLMAT"))
      return (FSLREG_TYPE);
    else if (!strcmp(dot, "LTA"))
      return (TRANSFORM_ARRAY_TYPE);
    else if (!strcmp(dot, "DAT"))
      return (REGISTER_DAT);
    else if (!strcmp(dot, "REG"))
      return (REGISTER_DAT);
    else
      return (REGISTER_DAT);
  }

  return (file_type);
}

static int ltaFSLwrite(const LTA *lta, const char *fname)
{
  FILE *fp;
  MATRIX *m_L;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "ltaFSLwrite: could not open file %s", fname));

  // create shallow copy of LTA
  LTA *ltatmp = LTAalloc(1, NULL);
  ltatmp->xforms[0].m_L = MatrixCopy(lta->xforms[0].m_L, NULL);
  copyVolGeom(&lta->xforms[0].src, &ltatmp->xforms[0].src);
  copyVolGeom(&lta->xforms[0].dst, &ltatmp->xforms[0].dst);
  ltatmp->type = lta->type;

  // convert to FSL matrix if necessary
  if (ltatmp->type != FSLREG_TYPE) LTAchangeType(ltatmp, FSLREG_TYPE);
  if (ltatmp->type != FSLREG_TYPE)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "ltaFSLwrite: could not convert to FSL format"));

  // write matrix to file:
  m_L = ltatmp->xforms[0].m_L;
  int i, j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) fprintf(fp, "%13.8f ", m_L->rptr[i + 1][j + 1]);
    fprintf(fp, "\n");
  }
  fclose(fp);

  LTAfree(&ltatmp);
  return (NO_ERROR);
}

#include "minc.h"

static int ltaMNIwrite(const LTA *lta, const char *fname)
{
  FILE *fp;
  int row;
  MATRIX *m_L;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "ltaMNIwrite: could not open file %s", fname));

  fprintf(fp, "MNI Transform File\n");
  // now saves src and dst in comment line
  fprintf(fp, "%%Generated by %s src %s dst %s\n", Progname, lta->xforms[0].src.fname, lta->xforms[0].dst.fname);
  fprintf(fp, "\n");
  fprintf(fp, "Transform_Type = Linear;\n");
  fprintf(fp, "Linear_Transform =\n");

  // MNI_TRANSFORM_TYPE is used to force MNI output of otherwise RAS2RAS lta
  if (lta->type == LINEAR_RAS_TO_RAS || lta->type == MNI_TRANSFORM_TYPE) {
    m_L = lta->xforms[0].m_L;
    for (row = 1; row <= 3; row++) {
      fprintf(fp,
              "%13.8f %13.8f %13.8f %13.8f ",
              *MATRIX_RELT(m_L, row, 1),
              *MATRIX_RELT(m_L, row, 2),
              *MATRIX_RELT(m_L, row, 3),
              *MATRIX_RELT(m_L, row, 4));
      if (row == 3) fprintf(fp, ";");
      fprintf(fp, "\n");
    }
  }
  else if (lta->type == LINEAR_VOX_TO_VOX) {
    // we use src and dst info to create RAS_TO_RAS xfm
    MATRIX *voxFromRAS = 0;
    MATRIX *rasFromVoxel = 0;
    MATRIX *tmp = 0;
    MATRIX *rasToRAS = 0;
    LT *lt = &lta->xforms[0];
    voxFromRAS = extract_r_to_i(&lt->src);
    tmp = MatrixMultiply(lta->xforms[0].m_L, voxFromRAS, NULL);
    rasFromVoxel = extract_i_to_r(&lt->dst);
    rasToRAS = MatrixMultiply(rasFromVoxel, tmp, NULL);
    for (row = 1; row <= 3; row++) {
      fprintf(fp,
              "%13.8f %13.8f %13.8f %13.8f ",
              *MATRIX_RELT(rasToRAS, row, 1),
              *MATRIX_RELT(rasToRAS, row, 2),
              *MATRIX_RELT(rasToRAS, row, 3),
              *MATRIX_RELT(rasToRAS, row, 4));
      if (row == 3) fprintf(fp, ";");
      fprintf(fp, "\n");
    }
    MatrixFree(&voxFromRAS);
    MatrixFree(&rasFromVoxel);
    MatrixFree(&tmp);
    MatrixFree(&rasToRAS);
  }
  fclose(fp);
  return (NO_ERROR);
}
static LTA *ltaMNIread(const char *fname)
{
  LTA *lta;
  LINEAR_TRANSFORM *lt;
  char *cp, line[1000];
  FILE *fp;
  int row;
  MATRIX *m_L;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "ltaMNIread: could not open file %s", fname));

  lta = LTAalloc(1, NULL);
  lt = &lta->xforms[0];
  lt->sigma = 1.0f;
  lt->x0 = lt->y0 = lt->z0 = 0;

  fgetl(line, 900, fp);                        /* MNI Transform File */
  fgetl(line, 900, fp);                        /* Transform_type = Linear; */
  while (line[0] == '%') fgetl(line, 900, fp); /* variable # of comments */
  fgetl(line, 900, fp);

  m_L = lt->m_L;
  for (row = 1; row <= 3; row++) {
    cp = fgetl(line, 900, fp);
    if (!cp) {
      LTAfree(&lta);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "ltaMNIread: could not "
                   "read row %d from %s (%s)",
                   row,
                   fname,
                   line));
    }
    sscanf(cp,
           "%f %f %f %f",
           MATRIX_RELT(m_L, row, 1),
           MATRIX_RELT(m_L, row, 2),
           MATRIX_RELT(m_L, row, 3),
           MATRIX_RELT(m_L, row, 4));
  }
  fclose(fp);
  return (lta);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int LTAworldToWorld(LTA *lta, float x, float y, float z, float *px, float *py, float *pz)
{
  static VECTOR *v_X, *v_Y = NULL;

  if (v_Y == NULL) {
    v_X = VectorAlloc(4, MATRIX_REAL);
    v_Y = VectorAlloc(4, MATRIX_REAL);
  }
  /* world to voxel */
  v_X->rptr[4][1] = 1.0f;
  V3_X(v_X) = x;
  V3_Y(v_X) = y;
  V3_Z(v_X) = z;

  LTAtransformPoint(lta, v_X, v_Y);

/* voxel to world */
  *px = V3_X(v_Y);
  *py = V3_Y(v_Y);
  *pz = V3_Z(v_Y);

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  This is the same as LTAworldToWorld but doesn't
  do the weird (and mostly incorrect) voxel conversion.
  ------------------------------------------------------*/
int LTAworldToWorldEx(LTA *lta, float x, float y, float z, float *px, float *py, float *pz)
{
  static VECTOR *v_X, *v_Y = NULL;

  if (v_Y == NULL) {
    v_X = VectorAlloc(4, MATRIX_REAL);
    v_Y = VectorAlloc(4, MATRIX_REAL);
  }

  v_X->rptr[4][1] = 1.0f;
  V3_X(v_X) = x;
  V3_Y(v_X) = y;
  V3_Z(v_X) = z;

  LTAtransformPoint(lta, v_X, v_Y);

  *px = V3_X(v_Y);
  *py = V3_Y(v_Y);
  *pz = V3_Z(v_Y);

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int LTAinverseWorldToWorld(LTA *lta, float x, float y, float z, float *px, float *py, float *pz)
{
  static VECTOR *v_X, *v_Y = NULL;

  if (v_Y == NULL) {
    v_X = VectorAlloc(4, MATRIX_REAL);
    v_Y = VectorAlloc(4, MATRIX_REAL);
  }
  /* world to voxel */
  v_X->rptr[4][1] = 1.0f;
  V3_X(v_X) = 128.0 - x;
  V3_Z(v_X) = (y + 128.0);
  V3_Y(v_X) = (-z + 128.0);

  LTAinverseTransformPoint(lta, v_X, v_Y);

  /* voxel to world */
  *px = 128.0 - V3_X(v_Y);
  *py = V3_Z(v_Y) - 128.0;
  *pz = -(V3_Y(v_Y) - 128.0);

  return (NO_ERROR);
}

int LTAinverseWorldToWorldEx(LTA *lta, float x, float y, float z, float *px, float *py, float *pz)
{
  static VECTOR *v_X, *v_Y = NULL;

  if (v_Y == NULL) {
    v_X = VectorAlloc(4, MATRIX_REAL);
    v_Y = VectorAlloc(4, MATRIX_REAL);
  }
  /* world to voxel */
  v_X->rptr[4][1] = 1.0f;
  V3_X(v_X) = x;
  V3_Y(v_X) = y;
  V3_Z(v_X) = z;

  LTAinverseTransformPoint(lta, v_X, v_Y);

  /* voxel to world */
  *px = V3_X(v_Y);
  *py = V3_Y(v_Y);
  *pz = V3_Z(v_Y);

  return (NO_ERROR);
}

// this assumes that lta was ras to ras
int LTAtoVoxelCoords(LTA *lta, MRI *mri)
{
  MATRIX *m_L;
  int i;

  for (i = 0; i < lta->num_xforms; i++) {
    m_L = MRIrasXformToVoxelXform(mri, mri, lta->xforms[i].m_L, NULL);
    MatrixFree(&lta->xforms[0].m_L);
    lta->xforms[0].m_L = m_L;
  }
  lta->type = LINEAR_VOX_TO_VOX;
  return (NO_ERROR);
}

int LTAvoxelToRasXform(LTA *lta, MRI *mri_src, MRI *mri_dst)
{
  MATRIX *m_L;
  int i;

  for (i = 0; i < lta->num_xforms; i++) {
    if (mri_src == NULL) {
      MATRIX *m_source_r2v, *m_dst_v2r, *m_tmp;
      // Before 10/2018, VGget*To*Xform() returned the inverse of the
      // transform one would expect from the function name. This is now
      // fixed. It seems the problem was unnoticed here, however. To keep the
      // output of LTAvoxelToRasXform() unchanged, we swapped the following
      // two function invocations:
      m_source_r2v = VGgetVoxelToRasXform(&lta->xforms[i].src, NULL, 0);
      m_dst_v2r = VGgetRasToVoxelXform(&lta->xforms[i].dst, NULL, 0);
      m_tmp = MatrixMultiply(lta->xforms[i].m_L, m_source_r2v, NULL);
      m_L = MatrixMultiply(m_dst_v2r, m_tmp, NULL);
      MatrixFree(&m_tmp);
      MatrixFree(&m_dst_v2r);
      MatrixFree(&m_source_r2v);
    }
    else {
      m_L = MRIvoxelXformToRasXform(mri_src, mri_dst, lta->xforms[i].m_L, NULL);
    }
    MatrixFree(&lta->xforms[i].m_L);
    lta->xforms[i].m_L = m_L;
  }

  lta->type = LINEAR_RAS_TO_RAS;
  return (NO_ERROR);
}

int LTArasToVoxelXform(LTA *lta, MRI *mri_src, MRI *mri_dst)
{
  MATRIX *m_L;
  int i;

  if (lta->type == LINEAR_VOX_TO_VOX) return (NO_ERROR);
  for (i = 0; i < lta->num_xforms; i++) {
    m_L = MRIrasXformToVoxelXform(mri_src, mri_dst, lta->xforms[i].m_L, NULL);
    MatrixFree(&lta->xforms[0].m_L);
    lta->xforms[0].m_L = m_L;
  }

  LTAsetVolGeom(lta, mri_src, mri_dst);
  lta->type = LINEAR_VOX_TO_VOX;
  return (NO_ERROR);
}
int LTAvoxelTransformToCoronalRasTransform(LTA *lta)
{
  MATRIX *V, *W, *m_tmp;

  V = MatrixAlloc(4, 4, MATRIX_REAL); /* world to voxel transform */
  W = MatrixAlloc(4, 4, MATRIX_REAL); /* voxel to world transform */
  *MATRIX_RELT(V, 1, 1) = -1;
  *MATRIX_RELT(V, 1, 4) = 128;
  *MATRIX_RELT(V, 2, 3) = -1;
  *MATRIX_RELT(V, 2, 4) = 128;
  *MATRIX_RELT(V, 3, 2) = 1;
  *MATRIX_RELT(V, 3, 4) = 128;
  *MATRIX_RELT(V, 4, 4) = 1;

  *MATRIX_RELT(W, 1, 1) = -1;
  *MATRIX_RELT(W, 1, 4) = 128;
  *MATRIX_RELT(W, 2, 3) = 1;
  *MATRIX_RELT(W, 2, 4) = -128;
  *MATRIX_RELT(W, 3, 2) = -1;
  *MATRIX_RELT(W, 3, 4) = 128;
  *MATRIX_RELT(W, 4, 4) = 1;

  m_tmp = MatrixMultiply(lta->xforms[0].m_L, V, NULL);
  MatrixMultiply(W, m_tmp, lta->xforms[0].m_L);
  MatrixFree(&V);
  MatrixFree(&W);
  MatrixFree(&m_tmp);
  lta->type = LINEAR_RAS_TO_RAS;

  return (NO_ERROR);
}
/*-----------------------------------------------------------
  FixMNITal() - function to compute the "real" talairach coordinates
  from the MNI talaiarch coordinates. Can be done in-place. This is
  based on Matthew Brett's 10/8/98 transform. See
  http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
  -----------------------------------------------------------*/
int FixMNITal(float xmni, float ymni, float zmni, float *xtal, float *ytal, float *ztal)
{
  MATRIX *T, *xyzMNI, *xyzTal;

  T = MatrixAlloc(4, 4, MATRIX_REAL);
  if (zmni >= 0.0) {
    stuff_four_by_four(
        T, .9900, .0000, .0000, 0, .0000, .9688, .0460, 0, .0000, -.0485, .9189, 0, .0000, .0000, .0000, 1);
  }
  else {
    stuff_four_by_four(
        T, .9900, .0000, .0000, 0, .0000, .9688, .0420, 0, .0000, -.0485, .8390, 0, .0000, .0000, .0000, 1);
  }

  xyzMNI = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzMNI->rptr[1][1] = xmni;
  xyzMNI->rptr[2][1] = ymni;
  xyzMNI->rptr[3][1] = zmni;
  xyzMNI->rptr[4][1] = 1.0;

  xyzTal = MatrixMultiply(T, xyzMNI, NULL);

  *xtal = xyzTal->rptr[1][1];
  *ytal = xyzTal->rptr[2][1];
  *ztal = xyzTal->rptr[3][1];

  MatrixFree(&T);
  MatrixFree(&xyzMNI);
  MatrixFree(&xyzTal);

  return (0);
}
/*-----------------------------------------------------------------
  DevolveXFM() - change a tranformation matrix to work with tkreg RAS
  coordinates (ie, c_ras = 0). For example, when used on a surface
  xyz, will produce talairach coords. Assumes that the transformation
  matrix was computed with the orig-volume native/scanner vox2ras
  transform.  The XFM maps from native orig RAS to a native target
  RAS. Note: this has no effect when the orig volume has c_ras = 0. If
  XFM is empty, then xfmname is read from the subject's transforms
  directory.  If xfmname is empty, then talairach.xfm is used.

  Note: uses LTAvoxelTransformToCoronalRasTransform().
  -----------------------------------------------------------------*/
MATRIX *DevolveXFM(const char *subjid, MATRIX *XFM, const char *xfmname)
{
  return (DevolveXFMWithSubjectsDir(subjid, XFM, xfmname, NULL));
}

MATRIX *DevolveXFMWithSubjectsDir(const char *subjid, MATRIX *XFM, const char *xfmname, const char *sdir)
{
  MATRIX *Torig_tkreg, *invTorig_tkreg, *Torig_native, *Mfix;
  char dirname[2000], xfmpath[2000];
  const char *sd;
  MRI *mriorig;
  FILE *fp;
  LTA *lta;

  if (sdir)
    sd = sdir;
  else {
    sd = getenv("SUBJECTS_DIR");
    if (sd == NULL) {
      printf("ERROR: SUBJECTS_DIR not defined\n");
      return (NULL);
    }
  }

  /* Check that the subject exists */
  sprintf(dirname, "%s/%s", sd, subjid);
  if (!fio_IsDirectory(dirname)) {
    printf("ERROR: cannot find subject %s in %s\n", subjid, sd);
    return (NULL);
  }

  /* Load the orig header for the subject */
  sprintf(dirname, "%s/%s/mri/orig.mgz", sd, subjid);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Trying %s\n", dirname);
  if (fio_FileExistsReadable(dirname))
    mriorig = MRIreadHeader(dirname, MRI_MGH_FILE);
  else
    mriorig = NULL;

  if (mriorig == NULL) {
    sprintf(dirname, "%s/%s/mri/orig.mgh", sd, subjid);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Trying %s\n", dirname);
    if (fio_FileExistsReadable(dirname)) {
      mriorig = MRIreadHeader(dirname, MRI_MGH_FILE);
    }
    else
      mriorig = NULL;
    if (mriorig == NULL) {
      sprintf(dirname, "%s/%s/mri/orig", sd, subjid);
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("Trying %s\n", dirname);
      if (fio_IsDirectory(dirname)) {
        mriorig = MRIreadHeader(dirname, MRI_CORONAL_SLICE_DIRECTORY);
      }
      else
        mriorig = NULL;
      if (mriorig == NULL) {
        printf("ERROR: could not read header for %s\n", dirname);
        return (NULL);
      }
    }
  }

  if (XFM == NULL) {
    /* Read in the talairach.xfm matrix */
    if (xfmname == NULL) xfmname = "talairach.xfm";
    sprintf(xfmpath, "%s/%s/mri/transforms/%s", sd, subjid, xfmname);
    fp = fopen(xfmpath, "r");
    if (fp == NULL) {
      printf("ERROR: could not open %s for reading \n", xfmpath);
      return (NULL);
    }
    lta = LTAreadEx(xfmpath);
    if (lta == NULL) {
      printf("ERROR: reading %s\n", xfmpath);
      return (NULL);
    }
    if (lta->type == LINEAR_VOX_TO_VOX) LTAvoxelTransformToCoronalRasTransform(lta);
    XFM = MatrixCopy(lta->xforms[0].m_L, NULL);
    LTAfree(&lta);
    fclose(fp);
  }

  /* Mfix = Torig_native*inv(Torig_tkreg) */
  /* X2 = X*Mfix */
  Torig_tkreg = MRIxfmCRS2XYZtkreg(mriorig);
  Torig_native = MRIxfmCRS2XYZ(mriorig, 0);
  invTorig_tkreg = MatrixInverse(Torig_tkreg, NULL);
  Mfix = MatrixMultiply(Torig_native, invTorig_tkreg, NULL);
  MatrixMultiply(XFM, Mfix, XFM);

  MatrixFree(&Mfix);
  MatrixFree(&Torig_tkreg);
  MatrixFree(&invTorig_tkreg);
  MatrixFree(&Torig_native);
  MRIfree(&mriorig);

  return (XFM);
}
/*----------------------------------------------------------------*/
TRANSFORM *TransformRead(const char *fname)
{
  TRANSFORM *trans;
  GCA_MORPH *gcam = NULL;
  char fname_no_path[STRLEN];

  trans = (TRANSFORM *)calloc(1, sizeof(TRANSFORM));
  memset(trans, 0, sizeof(TRANSFORM));

  // firstly, check for filename 'identify.nofile', which does not exist
  // as a file, but instead is used to force creation of an identity
  // matrix of type linear vox2vox
  FileNameOnly(fname, fname_no_path);
  if (0 == strcmp(fname_no_path, "identity.nofile")) {
    trans->type = LINEAR_RAS_TO_RAS;
    LTA *lta = LTAalloc(1, NULL);
    lta->xforms[0].m_L = MatrixIdentity(4, NULL);
    lta->xforms[0].type = trans->type;
    lta->type = trans->type;
    trans->xform = (void *)lta;
    return trans;
  }
  // continue normal processing...

  trans->type = TransformFileNameType(fname);
  switch (trans->type) {
    case MNI_TRANSFORM_TYPE:
    case LINEAR_VOX_TO_VOX:
    case FSLREG_TYPE:
    case LINEAR_RAS_TO_RAS:
    case TRANSFORM_ARRAY_TYPE:
    default:
      trans->xform = (void *)LTAreadEx(fname);
      if (!trans->xform) {
        free(trans);
        return (NULL);
      }
      trans->type = ((LTA *)trans->xform)->type;
      break;
    case MORPH_3D_TYPE:
      gcam = GCAMread(fname);
      if (!gcam) {
        free(trans);
        return (NULL);
      }
      trans->xform = (void *)gcam;
      break;
  }
  return (trans);
}

int TransformFree(TRANSFORM **ptrans)
{
  int errCode = NO_ERROR;
  TRANSFORM *trans;

  trans = *ptrans;
  *ptrans = NULL;

  switch (trans->type) {
    default: {
      void *pvoid = (void *)&trans->xform;
      errCode = LTAfree((LTA **)pvoid);
      break;
    }
    case MORPH_3D_TYPE: {
      void *pvoid = (void *)&trans->xform;
      errCode = GCAMfree((GCA_MORPH **)pvoid);
      break;
    }
  }
  free(trans);

  return errCode;
}

/*
  take a voxel in MRI space and find the voxel in the
  gca/gcamorph space to which it maps. Note that the caller
  will scale this down depending on the spacing of
  the gca/gcamorph.
*/
// no range check is done here.   user must validate the range
int TransformSample(TRANSFORM *transform, float xv, float yv, float zv, float *px, float *py, float *pz)
{
  static VECTOR *v_input, *v_canon = NULL;
  float xt, yt, zt;
  LTA *lta;
  GCA_MORPH *gcam;
  int errCode = NO_ERROR, xi, yi, zi;
  double xd, yd, zd;

  *px = *py = *pz = 0;
  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *)transform->xform;

    if (!gcam->mri_xind)
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "TransformSample: gcam has not been inverted!"));

    // the following should not happen /////////////////
    if (xv < 0) xv = 0;
    if (xv >= gcam->mri_xind->width) xv = gcam->mri_xind->width - 1;
    if (yv < 0) yv = 0;
    if (yv >= gcam->mri_yind->height) yv = gcam->mri_yind->height - 1;
    if (zv < 0) zv = 0;
    if (zv >= gcam->mri_zind->depth) zv = gcam->mri_zind->depth - 1;

    xi = nint(xv);
    yi = nint(yv);
    zi = nint(zv);
    if (xi < 0)
      xi = 0;
    else if (xi >= gcam->mri_xind->width)
      xi = gcam->mri_xind->width - 1;
    if (yi < 0)
      yi = 0;
    else if (yi >= gcam->mri_yind->height)
      yi = gcam->mri_yind->height - 1;
    if (zi < 0)
      zi = 0;
    else if (zi >= gcam->mri_zind->depth)
      zi = gcam->mri_zind->depth - 1;

    xt = nint(MRIgetVoxVal(gcam->mri_xind, xi, yi, zi, 0) * gcam->spacing);
    yt = nint(MRIgetVoxVal(gcam->mri_yind, xi, yi, zi, 0) * gcam->spacing);
    zt = nint(MRIgetVoxVal(gcam->mri_zind, xi, yi, zi, 0) * gcam->spacing);
    MRIsampleVolume(gcam->mri_xind, xv, yv, zv, &xd);
    MRIsampleVolume(gcam->mri_yind, xv, yv, zv, &yd);
    MRIsampleVolume(gcam->mri_zind, xv, yv, zv, &zd);
    xt = (float)xd * gcam->spacing;
    yt = (float)yd * gcam->spacing;
    zt = (float)zd * gcam->spacing;
  }
  else {
    lta = (LTA *)transform->xform;
    if (lta->type != LINEAR_VOXEL_TO_VOXEL) {
      int i;
      printf("Converting to LTA type LINEAR_VOXEL_TO_VOXEL...\n");
      lta = LTAchangeType(lta, LINEAR_VOXEL_TO_VOXEL);
      printf("After conversion:\n");
      for (i = 0; i < lta->num_xforms; i++) {
        LINEAR_TRANSFORM *lt = &lta->xforms[i];
        MatrixAsciiWriteInto(stdout, lt->m_L);
      }
    }
    if (!v_canon) {
      v_input = VectorAlloc(4, MATRIX_REAL);
      v_canon = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_input, 4, 1) = 1.0;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0;
    }
    V3_X(v_input) = xv;
    V3_Y(v_input) = yv;
    V3_Z(v_input) = zv;
    MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon);
    xt = V3_X(v_canon);
    yt = V3_Y(v_canon);
    zt = V3_Z(v_canon);

  }
  *px = xt;
  *py = yt;
  *pz = zt;

  return errCode;
}
int TransformSampleReal(TRANSFORM *transform, float xv, float yv, float zv, float *px, float *py, float *pz)
{
  static VECTOR *v_input, *v_canon = NULL;
  float xt, yt, zt;
  LTA *lta;
  GCA_MORPH *gcam;
  int errCode = NO_ERROR, xi, yi, zi;

  *px = *py = *pz = 0;
  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *)transform->xform;
    if (!gcam->mri_xind)
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "TransformSample: gcam has not been inverted!"));

    // the following should not happen /////////////////
    if (xv < 0) xv = 0;
    if (xv >= gcam->mri_xind->width) xv = gcam->mri_xind->width - 1;
    if (yv < 0) yv = 0;
    if (yv >= gcam->mri_yind->height) yv = gcam->mri_yind->height - 1;
    if (zv < 0) zv = 0;
    if (zv >= gcam->mri_zind->depth) zv = gcam->mri_zind->depth - 1;

    xi = nint(xv);
    yi = nint(yv);
    zi = nint(zv);

    xt = MRIgetVoxVal(gcam->mri_xind, xi, yi, zi, 0) * gcam->spacing;
    yt = MRIgetVoxVal(gcam->mri_yind, xi, yi, zi, 0) * gcam->spacing;
    zt = MRIgetVoxVal(gcam->mri_zind, xi, yi, zi, 0) * gcam->spacing;
  }
  else {
    lta = (LTA *)transform->xform;
    if (lta->type != LINEAR_VOXEL_TO_VOXEL) {
      int i;
      printf("Converting to LTA type LINEAR_VOXEL_TO_VOXEL...\n");
      lta = LTAchangeType(lta, LINEAR_VOXEL_TO_VOXEL);
      printf("After conversion:\n");
      for (i = 0; i < lta->num_xforms; i++) {
        LINEAR_TRANSFORM *lt = &lta->xforms[i];
        MatrixAsciiWriteInto(stdout, lt->m_L);
      }
    }
    if (!v_canon) {
      v_input = VectorAlloc(4, MATRIX_REAL);
      v_canon = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_input, 4, 1) = 1.0;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0;
    }
    V3_X(v_input) = xv;
    V3_Y(v_input) = yv;
    V3_Z(v_input) = zv;
    MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon);
    xt = V3_X(v_canon);
    yt = V3_Y(v_canon);
    zt = V3_Z(v_canon);

    if (xt < 0) xt = 0;
    if (yt < 0) yt = 0;
    if (zt < 0) zt = 0;

    if (!v_canon) {
      VectorFree(&v_input);
      VectorFree(&v_canon);
    }
  }
  *px = xt;
  *py = yt;
  *pz = zt;

  return errCode;
}

// with interpolation
int TransformSampleReal2(TRANSFORM *transform, float xv, float yv, float zv, float *px, float *py, float *pz)
{
  static VECTOR *v_input, *v_canon = NULL;
  // float           xt, yt, zt ;
  double xt, yt, zt;
  LTA *lta;
  GCA_MORPH *gcam;
  int errCode = NO_ERROR;  //, xi, yi, zi;

  *px = *py = *pz = 0;
  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *)transform->xform;
    if (!gcam->mri_xind)
      ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "TransformSample: gcam has not been inverted!"));

    // the following should not happen /////////////////
    if (xv < 0) xv = 0;
    if (xv >= gcam->mri_xind->width) xv = gcam->mri_xind->width - 1;
    if (yv < 0) yv = 0;
    if (yv >= gcam->mri_yind->height) yv = gcam->mri_yind->height - 1;
    if (zv < 0) zv = 0;
    if (zv >= gcam->mri_zind->depth) zv = gcam->mri_zind->depth - 1;

    /*xi = nint(xv) ;
    yi = nint(yv) ;
    zi = nint(zv) ;

    xt = MRIgetVoxVal(gcam->mri_xind, xi, yi, zi, 0)*gcam->spacing ;
    yt = MRIgetVoxVal(gcam->mri_yind, xi, yi, zi, 0)*gcam->spacing ;
    zt = MRIgetVoxVal(gcam->mri_zind, xi, yi, zi, 0)*gcam->spacing ;*/

    MRIsampleVolume(gcam->mri_xind, xv, yv, zv, &xt);
    MRIsampleVolume(gcam->mri_yind, xv, yv, zv, &yt);
    MRIsampleVolume(gcam->mri_zind, xv, yv, zv, &zt);
    xt *= gcam->spacing;
    yt *= gcam->spacing;
    zt *= gcam->spacing;
  }
  else {
    lta = (LTA *)transform->xform;
    if (lta->type != LINEAR_VOXEL_TO_VOXEL) {
      int i;
      printf("Converting to LTA type LINEAR_VOXEL_TO_VOXEL...\n");
      lta = LTAchangeType(lta, LINEAR_VOXEL_TO_VOXEL);
      printf("After conversion:\n");
      for (i = 0; i < lta->num_xforms; i++) {
        LINEAR_TRANSFORM *lt = &lta->xforms[i];
        MatrixAsciiWriteInto(stdout, lt->m_L);
      }
    }
    if (!v_canon) {
      v_input = VectorAlloc(4, MATRIX_REAL);
      v_canon = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_input, 4, 1) = 1.0;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0;
    }
    V3_X(v_input) = xv;
    V3_Y(v_input) = yv;
    V3_Z(v_input) = zv;
    MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon);
    xt = V3_X(v_canon);
    yt = V3_Y(v_canon);
    zt = V3_Z(v_canon);

    if (xt < 0) xt = 0;
    if (yt < 0) yt = 0;
    if (zt < 0) zt = 0;

    if (!v_canon) {
      VectorFree(&v_input);
      VectorFree(&v_canon);
    }
  }
  *px = (float)xt;
  *py = (float)yt;
  *pz = (float)zt;

  return errCode;
}

/*
  take a voxel in gca/gcamorph space and find the MRI voxel to which
  it maps
*/
int TransformSampleInverse(TRANSFORM *transform, int xv, int yv, int zv, float *px, float *py, float *pz)
{
  static VECTOR *v_input, *v_canon = NULL;
  static MATRIX *m_L_inv;
  float xt, yt, zt;
  int xn, yn, zn;
  LTA *lta;
  GCA_MORPH *gcam;
  GCA_MORPH_NODE *gcamn;
  int errCode = NO_ERROR;

  if (transform->type == MORPH_3D_TYPE) {
    gcam = (GCA_MORPH *)transform->xform;
    xn = nint(xv / gcam->spacing);
    yn = nint(yv / gcam->spacing);
    zn = nint(zv / gcam->spacing);

    if (xn >= gcam->width) xn = gcam->width - 1;
    if (yn >= gcam->height) yn = gcam->height - 1;
    if (zn >= gcam->depth) zn = gcam->depth - 1;
    if (xn < 0) xn = 0;
    if (yn < 0) yn = 0;
    if (zn < 0) zn = 0;

    gcamn = &gcam->nodes[xn][yn][zn];
    xt = gcamn->x;
    yt = gcamn->y;
    zt = gcamn->z;
    // if marked invalid, then return error
    if (gcamn->invalid) errCode = ERROR_BADPARM;
  }
  else {
    lta = (LTA *)transform->xform;
    if (lta->type != LINEAR_VOXEL_TO_VOXEL) {
      int i;
      printf("Converting to LTA type LINEAR_VOXEL_TO_VOXEL...\n");
      lta = LTAchangeType(lta, LINEAR_VOXEL_TO_VOXEL);
      printf("After conversion:\n");
      for (i = 0; i < lta->num_xforms; i++) {
        LINEAR_TRANSFORM *lt = &lta->xforms[i];
        MatrixAsciiWriteInto(stdout, lt->m_L);
      }
    }
    if (!v_canon) {
      v_input = VectorAlloc(4, MATRIX_REAL);
      v_canon = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_input, 4, 1) = 1.0;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0;
      m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL);
    }

    V3_X(v_canon) = (float)xv;
    V3_Y(v_canon) = (float)yv;
    V3_Z(v_canon) = (float)zv;
    MatrixMultiply(lta->inv_xforms[0].m_L, v_canon, v_input);
    xt = V3_X(v_input);
    yt = V3_Y(v_input);
    zt = V3_Z(v_input);
    // here I cannot get access to width, height, depth values
    // thus I cannot judge the point is good or bad
    // errCode remains to be valid
  }
  *px = xt;
  *py = yt;
  *pz = zt;

  return errCode;
}

int TransformSampleInverseFloat(const TRANSFORM *transform, float xv, float yv,
                                float zv, float *px, float *py, float *pz)
{
  static VECTOR *v_input, *v_canon = NULL;
  static MATRIX *m_L_inv;
  LTA *lta;
  int errCode = NO_ERROR;

  if (transform->type == MORPH_3D_TYPE) {
    // Return error if out of bounds instead of closest valid coordinates: when
    // interpolating, this will prevent things like repeating voxels at the tip
    // of the nose until reaching the edge of the image.
    return GCAMsampleMorph((GCAM *)transform->xform, xv, yv, zv, px, py, pz);
  }
  else {
    lta = (LTA *)transform->xform;
    if (lta->type != LINEAR_VOXEL_TO_VOXEL) {
      int i;
      printf("Converting to LTA type LINEAR_VOXEL_TO_VOXEL...\n");
      lta = LTAchangeType(lta, LINEAR_VOXEL_TO_VOXEL);
      printf("After conversion:\n");
      for (i = 0; i < lta->num_xforms; i++) {
        LINEAR_TRANSFORM *lt = &lta->xforms[i];
        MatrixAsciiWriteInto(stdout, lt->m_L);
      }
    }
    if (!v_canon) {
      v_input = VectorAlloc(4, MATRIX_REAL);
      v_canon = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_input, 4, 1) = 1.0;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0;
      m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL);
    }

    V3_X(v_canon) = (float)xv;
    V3_Y(v_canon) = (float)yv;
    V3_Z(v_canon) = (float)zv;
    MatrixMultiply(lta->inv_xforms[0].m_L, v_canon, v_input);
    *px = V3_X(v_input);
    *py = V3_Y(v_input);
    *pz = V3_Z(v_input);
    // here I cannot get access to width, height, depth values
    // thus I cannot judge the point is good or bad
    // errCode remains to be valid
  }

  return errCode;
}
int TransformSampleInverseVoxel(
    TRANSFORM *transform, int width, int height, int depth, int xv, int yv, int zv, int *px, int *py, int *pz)
{
  float xf, yf, zf;
  int errCode = NO_ERROR;

  errCode = TransformSampleInverse(transform, xv, yv, zv, &xf, &yf, &zf);

  xv = nint(xf);
  yv = nint(yf);
  zv = nint(zf);
  if (xv < 0) {
    errCode = ERROR_BADPARM;
    xv = 0;
  }
  if (xv >= width) {
    errCode = ERROR_BADPARM;
    xv = width - 1;
  }
  if (yv < 0) {
    errCode = ERROR_BADPARM;
    yv = 0;
  }
  if (yv >= height) {
    errCode = ERROR_BADPARM;
    yv = height - 1;
  }
  if (zv < 0) {
    errCode = ERROR_BADPARM;
    zv = 0;
  }
  if (zv >= depth) {
    errCode = ERROR_BADPARM;
    zv = depth - 1;
  }
  *px = xv;
  *py = yv;
  *pz = zv;

  return errCode;
}

TRANSFORM *TransformAlloc(int type, MRI *mri)
{
  TRANSFORM *transform;

  switch (type) {
    default:
      transform = (TRANSFORM *)calloc(1, sizeof(TRANSFORM));
      transform->type = type;
      transform->xform = (void *)LTAalloc(1, mri);
      break;
    case MORPH_3D_TYPE:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "TransformAlloc(MORPH_3D_TYPE): unsupported"));
      break;
  }
  return (transform);
}

int TransformSwapInverse(TRANSFORM *transform)
{
  LT *lt;
  LTA *lta;

  if (transform->type == MORPH_3D_TYPE) {
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED, "TransformSwapInverse: MORPH_3D_TYPE not supported"));
  }
  else {
    lta = (LTA *)(transform->xform);
    lt = lta->xforms;
    lta->xforms = lta->inv_xforms;
    lta->inv_xforms = lt;
  }
  return (NO_ERROR);
}

int TransformInvert(TRANSFORM *transform, MRI *mri)
{
  LTA *lta;
  GCA_MORPH *gcam;

  switch (transform->type) {
    default:
      lta = (LTA *)transform->xform;
      LTAfillInverse(lta);
      break;
    case MORPH_3D_TYPE:
      if (!mri) return (NO_ERROR);
      gcam = (GCA_MORPH *)transform->xform;
      GCAMinvert(gcam, mri);
      break;
  }
  return (NO_ERROR);
}
MRI *TransformApplyType(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst, int interp_type)
{
  LTA *lta = 0;
  switch (transform->type) {
    case MORPH_3D_TYPE:
      mri_dst = GCAMmorphToAtlas(mri_src, (GCA_MORPH *)transform->xform, mri_dst, -1, interp_type);
      break;
    default:
      // this does not work for RAS-to-RAS
      // mri_dst = MRIlinearTransformInterp(mri_src, NULL,
      //          ((LTA *)transform->xform)->xforms[0].m_L, interp_type);
      lta = (LTA *)transform->xform;
      mri_dst = LTAtransformInterp(mri_src, mri_dst, lta, interp_type);
      break;
  }
  return (mri_dst);
}

MRI *TransformApply(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  LTA *lta = 0;
  switch (transform->type) {
    case MORPH_3D_TYPE:
      // does take care of the dst c_ras position using atlas information
      mri_dst = GCAMmorphToAtlas(mri_src, (GCA_MORPH *)transform->xform, NULL, -1, SAMPLE_TRILINEAR);
      break;
    default:
      // now assumes that this is the LTA type
      // the following does the wrong thing for RAS-to-RAS
      // mri_dst = MRIlinearTransform(mri_src, NULL,
      //       ((LTA *)transform->xform)->xforms[0].m_L);
      // the following take care of dst c_ras poisiton
      lta = (LTA *)transform->xform;
      mri_dst = LTAtransform(mri_src, NULL, lta);
      break;
  }
  return (mri_dst);
}

MRI *TransformApplyInverse(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  LTA *lta = 0;
  switch (transform->type) {
    case MORPH_3D_TYPE:
      mri_dst = GCAMmorphFromAtlas(mri_src, (GCA_MORPH *)transform->xform, NULL, SAMPLE_NEAREST);
      break;
    default:
      // the following does not work for ras-to-ras
      // mri_dst = MRIinverseLinearTransform(mri_src, NULL,
      //      ((LTA *)transform->xform)->xforms[0].m_L);
      lta = (LTA *)transform->xform;
      LTAfillInverse(lta);
      mri_dst = LTAinverseTransformInterp(mri_src, mri_dst, lta, SAMPLE_TRILINEAR);
      break;
  }
  return (mri_dst);
}
MRI *TransformApplyInverseType(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst, int interp_type)
{
  LTA *lta = 0;
  switch (transform->type) {
    case MORPH_3D_TYPE:
      mri_dst = GCAMmorphFromAtlas(mri_src, (GCA_MORPH *)transform->xform, NULL, interp_type);
      break;
    default:
      // the following does not work for ras-to-ras
      // mri_dst = MRIinverseLinearTransform(mri_src, NULL,
      //      ((LTA *)transform->xform)->xforms[0].m_L);
      lta = (LTA *)transform->xform;
      LTAfillInverse(lta);
      mri_dst = LTAinverseTransformInterp(mri_src, mri_dst, lta, interp_type);
      break;
  }
  return (mri_dst);
}

#include "stats.h"
LTA *ltaReadRegisterDat(const char *fname, const char *mov, const char *ref)
{
  LTA *lta;
  char *tmpstr;
  float ipr, bpr;
  int err, float2int;
  MATRIX *R;

  lta = LTAalloc(1, NULL);
  lta->xforms[0].sigma = 1.0f;
  // (mr) does not seem to make a difference, so I removed it:
  // lta->xforms[0].x0 = lta->xforms[0].y0 = lta->xforms[0].z0 = 0 ;
  err = regio_read_register((char *)fname, &tmpstr, &ipr, &bpr, &lta->fscale, &R, &float2int);
  if (err > 0) ErrorReturn(NULL, (ERROR_BADPARM, "ltaReadRegisterDat: could not read %s", fname));
  strcpy(lta->subject, tmpstr);
  free(tmpstr);
  MatrixCopy(R, lta->xforms[0].m_L);
  MatrixFree(&R);

  /* The definitions of mov=src and ref=dst are consistent with
     tkregister2, LTAchangeType() and ltaReadRegisterDat(). This is an
     unfortunate definition because the registration matrix actually
     does from ref to mov. But this was an error introduced a long
     time ago and the rest of the code base has built up around it. */
  if (mov != NULL) {
    MRI *mritmp;
    mritmp = MRIreadHeader(mov, MRI_VOLUME_TYPE_UNKNOWN);
    if (mritmp == NULL) return (NULL);
    getVolGeom(mritmp, &lta->xforms[0].src);
    MRIfree(&mritmp);
  }
  else {
    // at least set what is known
    lta->xforms[0].src.valid = 0;
    lta->xforms[0].src.xsize = ipr;
    lta->xforms[0].src.ysize = ipr;
    lta->xforms[0].src.zsize = bpr;
  }
  if (ref != NULL) {
    MRI *mritmp;
    mritmp = MRIreadHeader(ref, MRI_VOLUME_TYPE_UNKNOWN);
    if (mritmp == NULL) return (NULL);
    getVolGeom(mritmp, &lta->xforms[0].dst);
    MRIfree(&mritmp);
  }
  else {
    lta->xforms[0].dst.valid = 0;
    lta->xforms[0].dst.xsize = 1;
    lta->xforms[0].dst.ysize = 1;
    lta->xforms[0].dst.zsize = 1;
  }

  // (mr) I dont think CORONAL_RAS_TO_CORONAL_RAS is the same here,
  // it did not work for me so I changed it to the REGISTER_DAT
  // if you change it back, make sure lta_convert does not break
  // lta->type = LINEAR_CORONAL_RAS_TO_CORONAL_RAS ;
  lta->type = REGISTER_DAT;
  return (lta);
}

// find volumes which created the transform.
// if buffer == NULL, then it will allocate memory
// parsing with strtok.  not thread safe.
int mincFindVolume(const char *line, const char *line2, char **srcVol, char **dstVol)
{
  static int count = 0;
  char buf[1024];
  char *pch;

  int mncorig = 0;
  // if not MGH way, try MNC way
  if (!strstr(line, "%Generated by")) {
    /*
    // There are two ways of getting .xfm files.
    // The most common one is by mritotal which has the
    // second line of xfm file looks like
    // %Wed Jan  7 13:23:45 2004>>> /usr/pubsw/packages/mni/1.1/
    //         bin/minctracc \
    // -clobber /tmp/mritotal_7539/orig_8_dxyz.mnc \
    // /usr/pubsw/packages/mni/1.1/share/mni_autoreg/average_305_8_dxyz.mnc \
    // transforms/talairach.xfm -transformation
    //      /tmp/mritotal_7539/orig_8tmp2c.xfm \
    // -lsq9 -xcorr -model_mask /usr/pubsw/packages
    //          /mni/1.1/share/mni_autoreg/average_305_8_mask.mnc \
    // -center 3.857440 25.165291 -28.701599 \
    //     -step 4 4 4 -tol 0.004 -simplex 2
    //
    // Another type is produced by GUI registration tool which
    //  has the following two lines
    // trick is that the first one is the target and the
    //     second one is the src
    // %Volume: /autofs/space/sound_003/users/ex_vivo_recons
    //      /I007/mri/T1/T1.mnc
    // %Volume: /space/solo/12/users/tmp/flash25_down.mnc
    //
    */
    // if minctracc way, then
    if (strstr(line, "minctracc")) {
      pch = strtok((char *)line, " ");
      while (pch != NULL) {
        strcpy(buf, pch);
        if (strstr(buf, ".mnc"))  // first src mnc volume
        {
          mncorig = 1;
          if (count == 0)  // src
          {
            count++;
            *srcVol = (char *)malloc(strlen(pch) + 1);
            strcpy(*srcVol, pch);
            if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Src volume %s\n", *srcVol);
          }
          else if (count == 1)
          // this is the second one must be dest volume
          {
            *dstVol = (char *)malloc(strlen(pch) + 1);
            strcpy(*dstVol, pch);
            if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Target volume %s\n", *dstVol);
            count = 0;  // restore for another case
            return 1;
          }
        }
        pch = strtok(NULL, " ");
      }
    }
    else  // let us assume MINC GUI transform
    {
      pch = strtok((char *)line, " ");  // points to %Volume
      pch = strtok(NULL, " ");          // now points to filename
      if (pch != NULL) {
        strcpy(buf, pch);
        if (strstr(buf, ".mnc"))  // if it is a minc file.  it is dst
        {
          *dstVol = (char *)malloc(strlen(pch) + 1);
          strcpy(*dstVol, pch);
        }
        pch = strtok((char *)line2, " ");  // points to %Volume
        pch = strtok(NULL, " ");           // now points to filename
        if (pch != NULL) {
          strcpy(buf, pch);
          if (strstr(buf, ".mnc"))  // if it is a minc file   it is src
          {
            *srcVol = (char *)malloc(strlen(pch) + 1);
            strcpy(*srcVol, pch);
          }
        }
      }
      if (*srcVol) {
        if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Src volume %s\n", *srcVol);
      }
      if (*dstVol) {
        if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Target volume %s\n", *dstVol);
      }
    }
  }
  else {
    // now MGH way  line has %Generated by ... src ... dst ...
    pch = strtok((char *)line, " ");
    while (pch != NULL) {
      strcpy(buf, pch);
      if (strstr(buf, "src"))  // first src mnc volume
      {
        // get next token
        pch = strtok(NULL, " ");
        *srcVol = (char *)malloc(strlen(pch) + 1);
        strcpy(*srcVol, pch);
        if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Src volume %s\n", *srcVol);
      }
      else if (strstr(buf, "dst")) {
        // get next token
        pch = strtok(NULL, " ");
        *dstVol = (char *)malloc(strlen(pch) + 1);
        strcpy(*dstVol, pch);
        if (DIAG_VERBOSE_ON) fprintf(stdout, "INFO: Target volume %s\n", *dstVol);
        return 1;
      }
      pch = strtok(NULL, " ");
    }
  }

  // neither MNC nor MGH
  return 0;
}

// find the volume and get the information
void mincGetVolumeInfo(const char *srcVol, VOL_GEOM *vgSrc)
{
  MRI *mri = 0;
  struct stat stat_buf;
  int ret;

  if (srcVol != 0) {
    // check the existence of a file
    ret = stat(srcVol, &stat_buf);
    if (ret != 0) {
      // annoying useless message commented-out:
      // fprintf(stderr, "INFO: Volume %s cannot be found.\n", srcVol);
      // now check whether it is average_305
      if (strstr(srcVol, "average_305")) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("INFO: The transform was made with average_305.mnc.\n");
        // average_305 value
        vgSrc->width = 172;
        vgSrc->height = 220;
        vgSrc->depth = 156;
        vgSrc->xsize = 1;
        vgSrc->ysize = 1;
        vgSrc->zsize = 1;
        vgSrc->x_r = 1;
        vgSrc->x_a = 0;
        vgSrc->x_s = 0;
        vgSrc->y_r = 0;
        vgSrc->y_a = 1;
        vgSrc->y_s = 0;
        vgSrc->z_r = 0;
        vgSrc->z_a = 0;
        vgSrc->z_s = 1;
        vgSrc->c_r = -0.0950;
        vgSrc->c_a = -16.5100;
        vgSrc->c_s = 9.7500;
        vgSrc->valid = 1;
      }
      else {
        // printf("INFO: Set Volume %s to the
        //      standard COR type.\n", srcVol);
        initVolGeom(vgSrc);  // valid = 0; so no need to give info
      }
    }
    else  // file exists
    {
      // note that both mri volume but also gca can be read
      mri = MRIreadHeader((char *)srcVol, MRI_VOLUME_TYPE_UNKNOWN);
      if (mri)  // find the MRI volume
        getVolGeom(mri, vgSrc);
      else  // cound not find the volume
        initVolGeom(vgSrc);
    }
  }
  else {
    initVolGeom(vgSrc);
  }
  // copy filename even if file does not exist to keep track
  if (srcVol)
    strcpy(vgSrc->fname, srcVol);
  else
    strcpy(vgSrc->fname, "unknown");

  // free memory
  if (mri) MRIfree(&mri);

  return;
}

void mincGetVolInfo(const char *infoline, const char *infoline2, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
{
  char *psrcVol = 0;
  char *pdstVol = 0;
  int retVal;

  retVal = mincFindVolume(infoline, infoline2, &psrcVol, &pdstVol);
  mincGetVolumeInfo(psrcVol, vgSrc);  // src may not be found
  mincGetVolumeInfo(pdstVol, vgDst);  // dst can be found
  if (vgDst->valid == 0) {
    if (getenv("USE_AVERAGE305")) {
      // average_305 value
      fprintf(stderr, "INFO: using average_305 info, since \n");
      fprintf(stderr, "INFO: environment var USE_AVERAGE305 set\n");
      vgDst->width = 172;
      vgDst->height = 220;
      vgDst->depth = 156;
      vgDst->xsize = 1;
      vgDst->ysize = 1;
      vgDst->zsize = 1;
      vgDst->x_r = 1;
      vgDst->x_a = 0;
      vgDst->x_s = 0;
      vgDst->y_r = 0;
      vgDst->y_a = 1;
      vgDst->y_s = 0;
      vgDst->z_r = 0;
      vgDst->z_a = 0;
      vgDst->z_s = 1;
      vgDst->c_r = -0.0950;
      vgDst->c_a = -16.5100;
      vgDst->c_s = 9.7500;
      vgDst->valid = 1;
    }
  }

  // check if they need to be freed
  if (psrcVol)
    free(psrcVol);
  if (pdstVol)
    free(pdstVol);
}

LTA *ltaMNIreadEx(const char *fname)
{
  LTA *lta = 0;
  LINEAR_TRANSFORM *lt;
  char *cp, line[2048], infoline[2048], infoline2[2048];
  FILE *fp;
  int row;
  MATRIX *m_L;
  int no_volinfo = 0;

  line[0] = 0;
  infoline[0] = 0;
  infoline2[0] = 0;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "ltMNIreadEx: could not open file %s", fname));

  lta = LTAalloc(1, NULL);
  lt = &lta->xforms[0];
  lt->sigma = 1.0f;
  lt->x0 = lt->y0 = lt->z0 = 0;

  fgetl(line, 900, fp); /* MNI Transform File */
  if (strncmp("MNI Transform File", line, 18)) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_NOFILE, "ltMNIreadEx:%s does not start as 'MNI Transform File'", fname));
  }

  fgetl(line, 900, fp); /* fileinfo line */
  if (line[0] == '%')
    strcpy(infoline, line);
  else {
    no_volinfo = 1;
    if (!strncmp("Transform_Type", line, 14)) {
      fgetl(line, 900, fp);
      goto get_transform;
    }
  }
  // second line in %
  fgetl(line, 900, fp);
  if (line[0] == '%') {
    strcpy(infoline2, line);
    while (line[0] == '%') fgetl(line, 900, fp); /* variable # of comments */
    fgetl(line, 900, fp);
    if (!strncmp("Transform_Type", line, 14)) {
      fgetl(line, 900, fp);
      goto get_transform;
    }
  }
  else {
    if (!strncmp("Transform_Type", line, 14)) {
      fgetl(line, 900, fp);
      goto get_transform;
    }
    while (line[0] == '%') fgetl(line, 900, fp); /* variable # of comments */
  }

get_transform:
  m_L = lt->m_L;
  for (row = 1; row <= 3; row++) {
    cp = fgetl(line, 900, fp);
    if (!cp) {
      fclose(fp);
      LTAfree(&lta);
      ErrorReturn(NULL, (ERROR_BADFILE, "ltMNIreadEx: could not read row %d from %s (%s)", row, fname, line));
    }
    sscanf(cp,
           "%f %f %f %f",
           MATRIX_RELT(m_L, row, 1),
           MATRIX_RELT(m_L, row, 2),
           MATRIX_RELT(m_L, row, 3),
           MATRIX_RELT(m_L, row, 4));
  }
  if (!lta) {
    fclose(fp);
    return NULL;
  }
  fclose(fp);

  // add original src and dst information
  if (no_volinfo == 0) mincGetVolInfo(infoline, infoline2, &lta->xforms[0].src, &lta->xforms[0].dst);
  lta->type = LINEAR_RAS_TO_RAS;
  return lta;
}

LTA *ltaReadFileEx(const char *fname)
{
  FILE *fp;
  LINEAR_TRANSFORM *lt;
  int i, nxforms, type;
  char line[STRLEN], *cp;
  LTA *lta;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't open file", fname));
  cp = fgetl(line, STRLEN - 1, fp);
  if (cp == NULL) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data", fname));
  }
  sscanf(cp, "type      = %d\n", &type);
  cp = fgetl(line, STRLEN - 1, fp);
  sscanf(cp, "nxforms   = %d\n", &nxforms);
  lta = LTAalloc(nxforms, NULL);
  lta->type = type;
  for (i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    fscanf(fp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0);
    fscanf(fp, "sigma     = %f\n", &lt->sigma);
    MatrixAsciiReadFrom(fp, lt->m_L);
  }
  // oh, well this is the added part
  for (i = 0; i < lta->num_xforms; i++) {
    if (fgets(line, STRLEN - 1, fp)) {
      if (strncmp(line, "src volume info", 15) == 0) {
        char *p;
        if (DIAG_VERBOSE_ON) fprintf(stderr, "INFO: src volume info present\n");
        readVolGeom(fp, &lta->xforms[i].src);
        p = fgets(line, STRLEN - 1, fp);
        if (strncmp(line, "dst volume info", 15) == 0) {
          if (DIAG_VERBOSE_ON) fprintf(stderr, "INFO: dst volume info present\n");

          readVolGeom(fp, &lta->xforms[i].dst);
        }
      }
    }
  }

  // these are extras for tkregister2 that may or may not be in the file
  lta->subject[0] = 0;
  lta->fscale = .15;
  while (fgets(line, STRLEN - 1, fp)) {
    // printf("reading extra input line %s", line) ;
    if (strncmp(line, "subject", 7) == 0)
      sscanf(line, "%*s %s", lta->subject);
    else if (strncmp(line, "fscale", 6) == 0)
      sscanf(line, "%*s %f", &lta->fscale);
  }

  fclose(fp);
  return (lta);
}

LTA *LTAreadExType(const char *fname, int type)
{
  char fname_no_path[STRLEN];
  LTA *lta = NULL;

  FileNameOnly(fname, fname_no_path);
  // firstly, check for filename 'identify.nofile', which does not exist
  // as a file, but instead is used to force creation of an identity
  // matrix of type linear vox2vox
  if (0 == strcmp(fname_no_path, "identity.nofile")) {
    LTA *lta = LTAalloc(1, NULL);
    lta->type = LINEAR_RAS_TO_RAS;
    lta->xforms[0].m_L = MatrixIdentity(4, NULL);
    lta->xforms[0].type = lta->type;
    return lta;
  }
  // continue normal processing...

  switch (type) {
    case FSLREG_TYPE:
      lta = ltaFSLread(fname);
      break;
    case REGISTER_DAT:
      printf(
          "INFO: This REGISTER_DAT transform "
          "is valid only for volumes between "
          " COR types with c_(r,a,s) = 0.\n");
      lta = ltaReadRegisterDat((char *)fname, NULL, NULL);
      if (!lta) return (NULL);

// (mr) I dont think CORONAL_RAS_TO_CORONAL_RAS is the same here,
// it did not work for me so I changed it to the REGISTER_DAT
// if you change it back, make sure lta_convert does not break
// lta->type = LINEAR_CORONAL_RAS_TO_CORONAL_RAS ;
      break;

    case MNI_TRANSFORM_TYPE:
      // this routine collect info on original src and dst volume
      // so that you can use the information to modify c_(r,a,s)
      // we no longer convert the transform to vox-to-vox
      // the transform is ras-to-ras
      lta = ltaMNIreadEx(fname);
      break;

    case LINEAR_VOX_TO_VOX:
    case LINEAR_RAS_TO_RAS:
    case TRANSFORM_ARRAY_TYPE:
    default:
      // get the original src and dst information
      lta = ltaReadFileEx(fname);
      break;
  }
  return (lta);
}

LTA *LTAreadEx(const char *fname)
{
  int type;
  type = TransformFileNameType((char *)fname);
  return LTAreadExType(fname, type);
}

// write out the information on src and dst volume
// in addition to the transform
int LTAprint(FILE *fp, const LTA *lta)
{
  int i, c, r;
  LT *lt;

  // fprintf(fp, "type      = %d\n", lta->type) ;
  fprintf(fp, "type      = %d ", lta->type);
  if (lta->type == LINEAR_VOX_TO_VOX) fprintf(fp, "# LINEAR_VOX_TO_VOX");
  else if (lta->type == LINEAR_RAS_TO_RAS) fprintf(fp, "# LINEAR_RAS_TO_RAS");
  else if (lta->type == REGISTER_DAT) fprintf(fp, "# REGISTER_DAT");
  fprintf(fp, "\n");
  fprintf(fp, "nxforms   = %d\n", lta->num_xforms);
  for (i = 0; i < lta->num_xforms; i++) {
    lt = &lta->xforms[i];
    fprintf(fp, "mean      = %6.4f %6.4f %6.4f\n", lt->x0, lt->y0, lt->z0);
    fprintf(fp, "sigma     = %6.4f\n", lt->sigma);
    // MatrixAsciiWriteInto(fp, lt->m_L) ;
    fprintf(fp, "1 4 4\n");  // Matrix size
    for (r = 1; r <= 4; r++) {
      for (c = 1; c <= 4; c++) fprintf(fp, "%18.15le ", (double)lt->m_L->rptr[r][c]);
      fprintf(fp, "\n");
    }
  }
  // write out src and dst volume info if there is one
  // note that this info may or may not be valid depending on vg->valid value
  // oh, well this is the addition
  for (i = 0; i < lta->num_xforms; ++i) {
    fprintf(fp, "src volume info\n");
    writeVolGeom(fp, &lta->xforms[i].src);
    fprintf(fp, "dst volume info\n");
    writeVolGeom(fp, &lta->xforms[i].dst);
  }

  // tkregister2 stuff
  {
    if (strlen(lta->subject) > 0) fprintf(fp, "subject %s\n", lta->subject);
    if (lta->fscale > 0) fprintf(fp, "fscale %f\n", lta->fscale);
  }
  return (NO_ERROR);
}

int  // (NOT TRUE ANYMORE: OK means 1, BAD means 0 )  use standard ERROR returns
LTAwriteEx(const LTA *lta, const char *fname)
{
  FILE *fp;
  const char *user;
  char ext[STRLEN];

  if (!stricmp(FileNameExtension((char *)fname, ext), "XFM") || lta->type == MNI_TRANSFORM_TYPE) {
    // someone defined NO_ERROR to be 0 and thus I have to change it
    return (ltaMNIwrite((LTA *)lta, (char *)fname));
  }
  else if (!stricmp(FileNameExtension((char *)fname, ext), "DAT") ||
           !stricmp(FileNameExtension((char *)fname, ext), "REG") || 
	   (lta->type == REGISTER_DAT && stricmp(FileNameExtension((char *)fname, ext), "LTA"))) {
    // If extension is "lta", then don't save it as a .dat
    int err;
    err = regio_write_register((char *)fname,
                               (char *)lta->subject,
                               lta->xforms[0].src.xsize,
                               lta->xforms[0].src.zsize,
                               lta->fscale,
                               lta->xforms[0].m_L,
                               FLT2INT_ROUND);
    if (err == 0)
      return (NO_ERROR);
    else
      return (ERROR_NOFILE);
  }
  else if (!stricmp(FileNameExtension((char *)fname, ext), "FSLMAT") || lta->type == FSLREG_TYPE) {
    return ltaFSLwrite(lta, fname);
  }

  fp = fopen(fname, "w");
  if (fp == NULL) ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "LTAwrite(%s): can't create file", fname));
  user = getenv("USER");
  if (!user) user = getenv("LOGNAME");
  if (!user) user = "UNKNOWN";
  fprintf(fp, "# transform file %s\n# created by %s on %s\n\n", fname, user, currentDateTime().c_str());
  LTAprint(fp, lta);
  fclose(fp);

  return (NO_ERROR);
}

// add src and voxel information
int LTAvoxelXformToRASXform(const MRI *src, const MRI *dst, LT *voxTran, LT *rasTran)
{
  MRIvoxelXformToRasXform((MRI *)src, (MRI *)dst, voxTran->m_L, rasTran->m_L);
  getVolGeom(src, &voxTran->src);
  getVolGeom(dst, &voxTran->dst);
  return 1;
}

static LTA *ltaFSLread(const char *fname)
{
  LTA *lta;
  LINEAR_TRANSFORM *lt;
  char *cp, line[1000];
  FILE *fp;
  int row;
  MATRIX *m_L;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NOFILE, "ltFSLread: could not open file %s", fname));

  lta = LTAalloc(1, NULL);
  lt = &lta->xforms[0];
  lt->sigma = 1.0f;
  lt->x0 = lt->y0 = lt->z0 = 0;

  m_L = lt->m_L;
  for (row = 1; row <= 3; row++) {
    cp = fgetl(line, 900, fp);
    if (!cp) {
      LTAfree(&lta);
      ErrorReturn(NULL, (ERROR_BADFILE, "ltFSLread: could not read row %d from %s", row, fname));
    }
    sscanf(cp,
           "%f %f %f %f",
           MATRIX_RELT(m_L, row, 1),
           MATRIX_RELT(m_L, row, 2),
           MATRIX_RELT(m_L, row, 3),
           MATRIX_RELT(m_L, row, 4));
  }
  fclose(fp);
  // (mr) I dont think physvox_to_physvox is the same here,
  // it did not work for me so I changed it to the FSLREG_TYPE
  // if you change it back, make sure lta_convert does not break
  // lta->type = LINEAR_PHYSVOX_TO_PHYSVOX;
  lta->type = FSLREG_TYPE;

  return (lta);
}

/*
  \fn LTA *LTAconcat2(LTA *lta1, LTA *lta2, int Reduce)
  \brief Same as LTAconcat() but for exactly two LTAs (which I
  seem to be doing a lot of).
 */
LTA *LTAconcat2(LTA *lta1, LTA *lta2, int Reduce)
{
  LTA *ltaArray[2], *ltaout;
  ltaArray[0] = lta1;
  ltaArray[1] = lta2;
  ltaout = LTAconcat(ltaArray, 2, Reduce);
  return (ltaout);
}

/*
  \fn LTA *LTAconcat(LTA **ltaArray, int nLTAs, int Reduce)
  \brief Concatenates all the transforms in all the LTAs into a single
   LTA. The LTAs can be of different types.  The LTAs do not need to
   be in the proper direction (direction is figured out from the
   volume geometry). If Reduce=0, all the transforms appear in the
   output LTA. If Reduce=1, all the transforms are combined into a
   single matrix using LTAreduce().
 */
LTA *LTAconcat(LTA **ltaArray, int nLTAs, int Reduce)
{
  int ntot, n, m, k, nx, DoInv = 0;
  LTA *lta, *lta1, *ltatmp;

  ntot = 0;
  for (n = 0; n < nLTAs; n++) {
    ntot += ltaArray[n]->num_xforms;
    LTAfillInverse(ltaArray[n]);  // just in case
  }
  lta = LTAalloc(ntot, NULL);

  // Copy the xforms from the first LTA into the new LTA
  k = 0;
  lta1 = ltaArray[0];
  for (m = 0; m < lta1->num_xforms; m++) {
    LTcopy(&lta1->xforms[m], &lta->xforms[k]);
    k++;
  }
  // Type, subject, and fscale get that of the first LTA
  lta->type = lta1->type;
  lta->fscale = lta1->fscale;
  strcpy(lta->subject, lta1->subject);

  // Go through the rest of the array
  for (n = 1; n < nLTAs; n++) {
    lta1 = ltaArray[n];
    if (lta1->type != lta->type) LTAchangeType(lta1, lta->type);
    nx = lta1->num_xforms;
    if (vg_isEqual(&lta->xforms[k - 1].dst, &lta1->xforms[0].src))
      DoInv = 0;
    else if (vg_isEqual(&lta->xforms[k - 1].dst, &lta1->xforms[nx - 1].dst)) {
      DoInv = 1;
      printf("WARNING: LTAconcat(): inverting LTA %d to match geometry\n", n);
    }
    else {
      printf("ERROR: LTAconcat(): LTAs %d and %d do not match\n", n - 1, n);
      printf("LTA %d -------------------\n", n - 1);
      LTAprint(stdout, ltaArray[n - 1]);
      printf("LTA %d -------------------\n", n);
      LTAprint(stdout, ltaArray[n]);
      printf("--------------------------\n");
      return (NULL);
    }
    if (DoInv == 0) {
      for (m = 0; m < nx; m++) {
        LTcopy(&lta1->xforms[m], &lta->xforms[k]);
        k++;
      }
    }
    else {
      for (m = nx - 1; m >= 0; m--) {
        LTcopy(&lta1->inv_xforms[m], &lta->xforms[k]);
        k++;
      }
    }
  }
  if (Reduce) {
    ltatmp = LTAreduce(lta);
    LTAfree(&lta);
    lta = ltatmp;
  }

  return (lta);
}

/*
  \fn LTA *LTAreduce(LTA *lta0)
  \brief Reduces a LTA with multiple transforms into a LTA with a
   single transform. The component LTs do not need to be in the
   correct direction (direction is figured out from the volume
   geometry).
 */
LTA *LTAreduce(const LTA *lta0)
{
  LTA *ltar, *lta;
  LT *lt0, *lt, *ltinv, *ltprev;
  int n, nx, DoInv = 0;
  MATRIX *M;

  lta = LTAcopy(lta0, NULL);
  LTAfillInverse(lta);

  ltar = LTAalloc(1, NULL);
  ltar->type = lta->type;
  ltar->fscale = lta->fscale;
  strcpy(ltar->subject, lta->subject);

  // Copy the first LT into the LT of the new LTA
  lt0 = &ltar->xforms[0];
  LTcopy(&lta->xforms[0], lt0);

  nx = lta->num_xforms;
  lt = lt0;  // in case only 1 xform
  for (n = 1; n < lta->num_xforms; n++) {
    ltprev = &lta->xforms[n - 1];
    lt = &lta->xforms[n];
    ltinv = &lta->inv_xforms[n];
    if (vg_isEqual(&ltprev->dst, &lt->src))
      DoInv = 0;
    else if (vg_isEqual(&ltprev->dst, &lt->dst)) {
      DoInv = 1;
      printf("WARNING: LTAreduce(): inverting LT %d to match geometry\n", n);
    }
    else {
      printf("ERROR: LTAreduce(): LTs %d and %d do not match\n", n - 1, n);
      return (NULL);
    }
    // printf("LTAreduce: n=%d, DoInv = %d\n",n,DoInv);
    if (DoInv == 0)
      M = lt->m_L;
    else
      M = ltinv->m_L;
    MatrixMultiply(M, lt0->m_L, lt0->m_L);  // new = M*old, M must be on left
  }
  if (DoInv == 0)
    memcpy(&lt0->dst, &lt->dst, sizeof(VOL_GEOM));
  else
    memcpy(&lt0->dst, &lt->src, sizeof(VOL_GEOM));

  LTAfree(&lta);
  return (ltar);
}

/*
  \fn LTA *LTAinvert(LTA *lta, LTA *ltainv)
  \breif Inverts the LTA by copying it and filling in the inverse
  using LTAfillInverse() then swaps xforms and inv_xforms. Note: there
  was a function with this name that would simply fill the inverse
  without performing the inverse. That function is now (Feb 2014)
  called LTAfillInverse().  The passed lta is not changed.
 */
LTA *LTAinvert(LTA *lta, LTA *ltainv)
{
  int i, j;
  LINEAR_TRANSFORM *lt;
  if (lta != ltainv) ltainv = LTAcopy(lta, ltainv);
  if (ltainv == NULL) return (NULL);
  LTAfillInverse(ltainv);
  
  // Because inv(A*B) = inv(B)*inv(A), we have to flip first. Flip xforms, too,
  // such that inv_xforms[n] stays the inverse of xforms[n]. Deep copy
  // unnecessary here.
  lt = (LT *)calloc(1, sizeof(LT));
  if (!lt) ErrorExit(ERROR_NOMEMORY, "ERROR LTAinvert(): no memory");
  for (i=0, j=lta->num_xforms-1; i<(lta->num_xforms+1)/2; i++, j--) {
    *lt = ltainv->inv_xforms[i];
    ltainv->inv_xforms[i] = ltainv->inv_xforms[j];
    ltainv->inv_xforms[j] = *lt;
    *lt = ltainv->xforms[i];
    ltainv->xforms[i] = ltainv->xforms[j];
    ltainv->xforms[j] = *lt;
  }
  free(lt);
  
  lt = ltainv->inv_xforms;
  ltainv->inv_xforms = ltainv->xforms;
  ltainv->xforms = lt;
  return (ltainv);
}
/*
  \fn LTA *LTAfillInverse(LTA *lta)
  \brief Computes and fills the inverse transforms (inv_xforms) and
  swaps vol_geom for the LTA. Note: this does not change the value of
  the forward m_L transform to be the inverse! It merely fills in the
  inverse xforms. This function used to be called LTAinvert() but the
  name was changed to better reflect what it actually does.  The
  LTAinvert() function now (Feb 2014) actually inverts the LTA.  I
  think the passed LTA is changed.
*/
LTA *LTAfillInverse(LTA *lta)
{
  int i;

  for (i = 0; i < lta->num_xforms; ++i) {
    if (MatrixInverse(lta->xforms[i].m_L, lta->inv_xforms[i].m_L) == NULL)
      ErrorExit(ERROR_BADPARM, "TransformInvert: xform noninvertible");
    memmove(&lta->inv_xforms[i].src, &lta->xforms[i].dst, sizeof(VOL_GEOM));
    memmove(&lta->inv_xforms[i].dst, &lta->xforms[i].src, sizeof(VOL_GEOM));
  }
  return lta;
}

// verify lta stored src and dst are the same as src and dst
int LTAmodifySrcDstGeom(LTA *lta, MRI *src, MRI *dst)
{
  LINEAR_TRANSFORM *lt = 0;  // work pointer
  int i;
  int resSrc = 0;
  int resDst = 0;
  int countValidSrc = 0;
  int countValidDst = 0;
  int res;
  VOL_GEOM svg, dvg;

  for (i = 0; i < lta->num_xforms; ++i) {
    lt = &lta->xforms[i];
    if (src)  // src is not NULL
    {
      if (lt->src.valid) {
        countValidSrc++;
        getVolGeom(src, &svg);
        res = vg_isEqual(&lt->src, &svg);
        resSrc += res;
        if (res == 0) {
          fprintf(stderr,
                  "INFO: src volume info "
                  "differs from the one stored in lta. "
                  "gets modified now.\n");
          svg.vgprint();
          lt->src.vgprint();
          getVolGeom(src, &lt->src);
        }
      }
      else if (lt->src.valid == 0)  // if not valid, just copy
      {
        getVolGeom(src, &lt->src);
      }
    }
    if (dst)  // dst is not NULL
    {
      if (lt->dst.valid) {
        countValidDst++;
        getVolGeom(dst, &dvg);
        res = vg_isEqual(&lt->dst, &dvg);
        resDst += res;
        if (res == 0) {
          fprintf(stderr,
                  "INFO: dst volume info differs "
                  "from the one stored in lta.  gets modified now.\n");
          dvg.vgprint();
          lt->dst.vgprint();
          getVolGeom(dst, &lt->dst);
        }
      }
      else if (lt->dst.valid == 0)  // if not valid, just copy
      {
        getVolGeom(dst, &lt->dst);
      }
    }
  }
  // currently just do nothing for error checking
  return 0;
}

// this one does not need to allocate sI2R and dR2I
// compared with MRIgetRasToVoxelXform().
static void LTAgetV2V(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
{
  //           sI2R
  //     src -------> RAS
  //      |?           | mod (input)
  //      V            V
  //     dst <------- RAS
  //           dR2I
  MATRIX *sI2R = vg_i_to_r(vgSrc);
  MATRIX *dR2I = vg_r_to_i(vgDst);
  MATRIX *tmp = 0;
  if (sI2R == 0 || dR2I == 0)
    ErrorExit(ERROR_BADPARM,
              "LTAgetV2V: passed volumes did "
              "not have the info on i_to_r or r_to_i.");
  tmp = MatrixMultiply(mod, sI2R, NULL);
  MatrixMultiply(dR2I, tmp, mod);  // m_L gets modified -> lta gets modified
  MatrixFree(&tmp);
  MatrixFree(&sI2R);
  MatrixFree(&dR2I);
}

// this one does not need to allocate sR2I and dI2R
// compared with MRIgetVoxToRasXform().
static void LTAgetR2R(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
{
  //           sR2I
  //     src <------- RAS
  //      | mod        | ?
  //      V            V
  //     dst -------> RAS
  //           dI2R
  MATRIX *sR2I = vg_r_to_i(vgSrc);
  MATRIX *dI2R = vg_i_to_r(vgDst);
  MATRIX *tmp = 0;
  if (sR2I == 0 || dI2R == 0)
    ErrorExit(ERROR_BADPARM,
              "LTAgetR2R: passed volumes did "
              "not have the info on r_to_i or i_to_r");
  tmp = MatrixMultiply(mod, sR2I, NULL);
  MatrixMultiply(dI2R, tmp, mod);  // m_L gets modified -> lta gets modified
  MatrixFree(&tmp);
  MatrixFree(&sR2I);
  MatrixFree(&dI2R);
}

/*
  \fn LTA *LTAchangeType(LTA *lta, int ltatype)
  \brief Changes the transform type. The LTA itself is changed.
 */
LTA *LTAchangeType(LTA *lta, int ltatype)
{
  LINEAR_TRANSFORM *lt;  // work pointer
  MATRIX *m_L;           // work pointer
  MATRIX *sISize = 0;
  MATRIX *sSize = 0;
  MATRIX *dISize = 0;
  MATRIX *dSize = 0;
  MATRIX *tmp = 0;
  int i;
  MRI *movmri = 0;
  MRI *refmri = 0;
  MATRIX *mreg = 0;
  MATRIX *mfsl = 0;
  // if it is the same, don't do anything
  if (lta->type == ltatype) return lta;

  // verify both src and dst have valid geometry
  for (i = 0; i < lta->num_xforms; ++i) {
    lt = &lta->xforms[i];
    if (lt->src.valid == 0) ErrorExit(ERROR_BADPARM, "LTAchangeType: src geometry must be valid\n");
    if (lt->dst.valid == 0) ErrorExit(ERROR_BADPARM, "LTAchangeType: dst geometry must be valid\n");
  }

  // ras-to-ras
  if (lta->type == LINEAR_RAS_TO_RAS) {
    switch (ltatype) {
      case LINEAR_VOX_TO_VOX:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          LTAgetV2V(m_L, &lt->src, &lt->dst);  // m_L gets modified
        }
        lta->type = LINEAR_VOX_TO_VOX;
        break;
      case LINEAR_PHYSVOX_TO_PHYSVOX:
        //       src' (regarded as 1 mm unit voxel)
        //        |sISize
        //        V
        //       src  -------> RAS
        //        |             |
        //        V             V
        //       dst  -------> RAS
        //        |dSize
        //        V
        //       dst' (regarded as 1 mm unix voxel)
        //  we need src' to dst'
        sISize = MatrixIdentity(4, 0);
        dSize = MatrixIdentity(4, 0);
        tmp = MatrixIdentity(4, NULL);
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          *MATRIX_RELT(sISize, 1, 1) = 1. / lt->src.xsize;
          *MATRIX_RELT(sISize, 2, 2) = 1. / lt->src.ysize;
          *MATRIX_RELT(sISize, 3, 3) = 1. / lt->src.zsize;
          *MATRIX_RELT(dSize, 1, 1) = lt->dst.xsize;
          *MATRIX_RELT(dSize, 2, 2) = lt->dst.ysize;
          *MATRIX_RELT(dSize, 3, 3) = lt->dst.zsize;
          LTAgetV2V(m_L, &lt->src, &lt->dst);  // m_L gets modified to be V2V
          tmp = MatrixMultiply(m_L, sISize, NULL);
          MatrixMultiply(dSize, tmp, m_L);  // modified to physvox to physvox
        }
        MatrixFree(&sISize);
        sISize = 0;
        MatrixFree(&dSize);
        dSize = 0;
        MatrixFree(&tmp);
        tmp = 0;
        lta->type = LINEAR_PHYSVOX_TO_PHYSVOX;
        break;
      case REGISTER_DAT:
        // from LINEAR_RAS_TO_RAS to REGISTER_DAT:
        for (i = 0; i < lta->num_xforms; ++i) {
          /* The definitions of mov=src and ref=dst are consistent
             with tkregister2, LTAchangeType() and
             ltaReadRegisterDat(). This is an unfortunate definition
             because the TKR registration matrix actually goes from
             ref to mov. All applications should recognize this and
             invert the matrix as needed. This was an error introduced
             a long, long time ago and the rest of the code base has built
             up around it. */
          lt = &lta->xforms[i];  // movsrc->refdst
          m_L = lt->m_L;
          mreg = MRItkRegMtx(&lt->dst, &lt->src, m_L);  // refdst->movsrc
          MatrixCopy(mreg, m_L);
          MatrixFree(&mreg);
        }
        lta->type = REGISTER_DAT;
        break;
      case FSLREG_TYPE:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          mreg = MRItkRegMtx(&lt->dst, &lt->src, m_L);
          mfsl = MRItkreg2FSL(&lt->dst, &lt->src, mreg);
          MatrixCopy(mfsl, m_L);
          MatrixFree(&mreg);
          MatrixFree(&mfsl);
        }
        lta->type = FSLREG_TYPE;
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType: you are "
                  "requesting ras-to-ras to %d ",
                  ltatype);
        break;
    }
  }
  else if (lta->type == LINEAR_VOX_TO_VOX) {
    switch (ltatype) {
      case LINEAR_RAS_TO_RAS:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          LTAgetR2R(m_L, &lt->src, &lt->dst);  // m_L gets modified to be R2R
        }
        lta->type = LINEAR_RAS_TO_RAS;
        break;
      case LINEAR_PHYSVOX_TO_PHYSVOX:
        //       src' (regarded as 1 mm unit voxel)
        //        |sISize
        //        V
        //       src
        //        |
        //        V
        //       dst
        //        |dSize
        //        V
        //       dst' (regarded as 1 mm unix voxel)
        //  we need src' to dst'
        sISize = MatrixIdentity(4, 0);
        dSize = MatrixIdentity(4, 0);
        tmp = MatrixIdentity(4, NULL);
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          *MATRIX_RELT(sISize, 1, 1) = 1. / lt->src.xsize;
          *MATRIX_RELT(sISize, 2, 2) = 1. / lt->src.ysize;
          *MATRIX_RELT(sISize, 3, 3) = 1. / lt->src.zsize;
          *MATRIX_RELT(dSize, 1, 1) = lt->dst.xsize;
          *MATRIX_RELT(dSize, 2, 2) = lt->dst.ysize;
          *MATRIX_RELT(dSize, 3, 3) = lt->dst.zsize;
          tmp = MatrixMultiply(m_L, sISize, NULL);
          MatrixMultiply(dSize, tmp, m_L);  // modified to physvox to physvox
        }
        MatrixFree(&sISize);
        sISize = 0;
        MatrixFree(&dSize);
        dSize = 0;
        MatrixFree(&tmp);
        tmp = 0;
        lta->type = LINEAR_PHYSVOX_TO_PHYSVOX;
        break;
      case REGISTER_DAT:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, REGISTER_DAT);
        break;
      case FSLREG_TYPE:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, FSLREG_TYPE);
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType: you are "
                  "requesting vox-to-vox to %d ",
                  ltatype);
        break;
    }
  }
  else if (lta->type == LINEAR_PHYSVOX_TO_PHYSVOX) {
    switch (ltatype) {
      case LINEAR_VOX_TO_VOX:
        //       src' (regarded as 1 mm unit voxel)
        //        |sISize
        //        V
        //       src
        //        | ?
        //        V
        //       dst
        //        |dSize
        //        V
        //       dst' (regarded as 1 mm unix voxel)
        //  we need src to dst.  thus
        //      X = dISize*M*sSize
        sSize = MatrixIdentity(4, 0);
        dISize = MatrixIdentity(4, 0);
        tmp = MatrixIdentity(4, NULL);
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          *MATRIX_RELT(sSize, 1, 1) = lt->src.xsize;
          *MATRIX_RELT(sSize, 2, 2) = lt->src.ysize;
          *MATRIX_RELT(sSize, 3, 3) = lt->src.zsize;
          *MATRIX_RELT(dISize, 1, 1) = 1. / lt->dst.xsize;
          *MATRIX_RELT(dISize, 2, 2) = 1. / lt->dst.ysize;
          *MATRIX_RELT(dISize, 3, 3) = 1. / lt->dst.zsize;
          tmp = MatrixMultiply(m_L, sSize, NULL);
          MatrixMultiply(dISize, tmp, m_L);  // modified to physvox to physvox
        }
        MatrixFree(&sSize);
        sSize = 0;
        MatrixFree(&dISize);
        dISize = 0;
        MatrixFree(&tmp);
        tmp = 0;
        lta->type = LINEAR_VOX_TO_VOX;
        break;
      case LINEAR_RAS_TO_RAS:
        sSize = MatrixIdentity(4, 0);
        dISize = MatrixIdentity(4, 0);
        tmp = MatrixIdentity(4, NULL);
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          *MATRIX_RELT(sSize, 1, 1) = lt->src.xsize;
          *MATRIX_RELT(sSize, 2, 2) = lt->src.ysize;
          *MATRIX_RELT(sSize, 3, 3) = lt->src.zsize;
          *MATRIX_RELT(dISize, 1, 1) = 1. / lt->dst.xsize;
          *MATRIX_RELT(dISize, 2, 2) = 1. / lt->dst.ysize;
          *MATRIX_RELT(dISize, 3, 3) = 1. / lt->dst.zsize;
          LTAgetV2V(m_L, &lt->src, &lt->dst);
          tmp = MatrixMultiply(m_L, sSize, NULL);
          MatrixMultiply(dISize, tmp, m_L);  // modified to physvox to physvox
        }
        MatrixFree(&sSize);
        sSize = 0;
        MatrixFree(&dISize);
        dISize = 0;
        MatrixFree(&tmp);
        tmp = 0;
        lta->type = LINEAR_RAS_TO_RAS;
        break;
      case REGISTER_DAT:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, REGISTER_DAT);
        break;
      case FSLREG_TYPE:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, FSLREG_TYPE);
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType: you are"
                  " requesting physvox-to-physvox to %d ",
                  ltatype);
        break;
    }
  }
  else if (lta->type == REGISTER_DAT) {
    switch (ltatype) {
      case LINEAR_RAS_TO_RAS:
        // from REGISTER_DAT to LINEAR_RAS_TO_RAS
        for (i = 0; i < lta->num_xforms; ++i) {
          /* The definitions of mov=src and ref=dst are consistent with
             tkregister2, LTAchangeType() and ltaReadRegisterDat(). This is an
             unfortunate definition because the registration matrix actually
             does from ref to mov. But this was an error introduced a long
             time ago and the rest of the code base has built up around it. */
          lt = &lta->xforms[0];  // movsrc->refdst
          m_L = lt->m_L;
          movmri = MRIallocHeader(lt->src.width, lt->src.height, lt->src.depth, MRI_UCHAR, 1);
          refmri = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR, 1);
          MRIcopyVolGeomToMRI(movmri, &lt->src);
          MRIcopyVolGeomToMRI(refmri, &lt->dst);
          mreg = MRItkReg2Native(refmri, movmri, m_L);  // movsrc->refdst
          MatrixCopy(mreg, m_L);
          MatrixFree(&mreg);
          MRIfree(&movmri);
          MRIfree(&refmri);
        }
        lta->type = LINEAR_RAS_TO_RAS;
        break;
      case LINEAR_VOX_TO_VOX:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);
        break;
      case FSLREG_TYPE:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          mfsl = MRItkreg2FSL(&lt->dst, &lt->src, m_L);
          MatrixCopy(mfsl, m_L);
          MatrixFree(&mfsl);
        }
        lta->type = FSLREG_TYPE;
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType unsupported: you are"
                  " requesting REGISTER_DAT to %d ",
                  ltatype);
        break;
    }
    if (Gdiag_no > 0) {
      printf("transformed matrix:\n");
      MatrixPrint(Gstdout, lta->xforms[0].m_L);
    }
  }
  else if (lta->type == LINEAR_CORONAL_RAS_TO_CORONAL_RAS) {
    MATRIX *m_sras2ras;
    MRI *mri_tmp;
    switch (ltatype) {
      case LINEAR_RAS_TO_RAS:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          mri_tmp = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR, 1);
          MRIcopyVolGeomToMRI(mri_tmp, &lt->dst);
          m_sras2ras = RASFromSurfaceRAS_(mri_tmp,NULL);
          m_L = MatrixMultiply(m_sras2ras, m_L, m_L);
          MatrixFree(&m_sras2ras);
          MRIfree(&mri_tmp);
        }
        lta->type = LINEAR_RAS_TO_RAS;
        break;
      case REGISTER_DAT:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, REGISTER_DAT);
        break;
      case FSLREG_TYPE:
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        lta = LTAchangeType(lta, FSLREG_TYPE);
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType unsupported: you are"
                  " requesting COR_RAS_TO_COR_RAS to %d ",
                  ltatype);
        break;
    }
    if (Gdiag_no > 0) {
      printf("transformed matrix:\n");
      MatrixPrint(Gstdout, lta->xforms[0].m_L);
    }
  }
  else if (lta->type == FSLREG_TYPE) {
    switch (ltatype) {
      case LINEAR_RAS_TO_RAS:
        lta = LTAchangeType(lta, REGISTER_DAT);
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
        break;
      case LINEAR_VOX_TO_VOX:
        lta = LTAchangeType(lta, REGISTER_DAT);
        lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);
        break;
      case REGISTER_DAT:
        for (i = 0; i < lta->num_xforms; ++i) {
          lt = &lta->xforms[i];
          m_L = lt->m_L;
          mreg = MRIfsl2TkReg(&lt->dst, &lt->src, m_L);
          MatrixCopy(mreg, m_L);
          MatrixFree(&mreg);
        }
        lta->type = REGISTER_DAT;
        break;
      default:
        ErrorExit(ERROR_BADPARM,
                  "LTAchangeType unsupported: you are"
                  " requesting FSLREG_TYPE to %d ",
                  ltatype);
        break;
    }
    if (Gdiag_no > 0) {
      printf("transformed matrix:\n");
      MatrixPrint(Gstdout, lta->xforms[0].m_L);
    }
  }

  // fill inverse part
  LTAfillInverse(lta);
  return lta;
}

// lta is the transform from src to dst for ras2ras or vox2vox
MATRIX *surfaceRASFromSurfaceRAS_(MRI *dst, MRI *src, LTA *lta)
{
  MATRIX *res = 0;
  MATRIX *tmp = 0;
  MATRIX *surf2src = 0;
  MATRIX *dst2surf = 0;
  LT *lt = 0;
  int ltaabsent = 0;
  // this is the combined operation
  //          surfaceRAS(src)
  //               |
  //               V
  //              src ---> RAS
  //               |        |
  //               |        | lta
  //               |        |
  //               V        V
  //              dst----> RAS
  //               |
  //               V
  //           surfaceRAS(dst)
  if (lta == 0) {
    ltaabsent = 1;
    fprintf(stderr, "INFO: assumes the identity RAS2RAS transform\n");
    lta = LTAalloc(1, NULL);
    lta->type = LINEAR_RAS_TO_RAS;
  }
  lt = &lta->xforms[0];
  if (lta->type == LINEAR_PHYSVOX_TO_PHYSVOX) {
    LTAchangeType(lta, LINEAR_RAS_TO_RAS);
  }
  if (lta->type == LINEAR_RAS_TO_RAS) {
    surf2src = RASFromSurfaceRAS_(src,NULL);
    dst2surf = surfaceRASFromRAS_(dst);
  }
  else if (lta->type == LINEAR_VOX_TO_VOX) {
    surf2src = voxelFromSurfaceRAS_(src);
    dst2surf = surfaceRASFromVoxel_(dst);
  }
  tmp = MatrixMultiply(lt->m_L, surf2src, NULL);
  res = MatrixMultiply(dst2surf, tmp, NULL);
  // memory management
  MatrixFree(&tmp);
  MatrixFree(&surf2src);
  MatrixFree(&dst2surf);
  //
  if (ltaabsent == 1) LTAfree(&lta);

  return res;
}

MRI *TransformCreateDensityMap(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  if (transform->type == MORPH_3D_TYPE) {
    GCA_MORPH *gcam = (GCA_MORPH *)(transform->xform);
    mri_dst = GCAMmorphToAtlasWithDensityCorrection(mri_src, gcam, mri_dst, 0);
  }
  else /* compute determinant of jacobian and apply it everywhere */
  {
    double det;

    mri_dst = TransformApply(transform, mri_src, mri_dst);
    det = MatrixDeterminant(((LTA *)(transform->xform))->xforms[0].m_L);
    printf("scaling volume by %2.3f...\n", det);
    MRIscalarMul(mri_dst, mri_dst, det);
  }
  return (mri_dst);
}

int TransformWrite(TRANSFORM *transform, const char *fname)
{
  switch (transform->type) {
    case MORPH_3D_TYPE:
      return (GCAMwrite((GCA_MORPH *)(transform->xform), fname));
      break;
    default:
      return (LTAwriteEx((LTA *)(transform->xform), fname));
      break;
  }
  return (NO_ERROR); /* will never get here */
}

// Allocates memory for new transform if destination NULL.
TRANSFORM *TransformCopy(const TRANSFORM *tsrc, TRANSFORM *tdst)
{
  GCAM *gcam_src, *gcam_dst;
  LTA *lta_src, *lta_dst;
  
  if (!tdst) tdst = (TRANSFORM *)calloc(1, sizeof(TRANSFORM));
  tdst->type = tsrc->type;
  tdst->mri_xn = tsrc->mri_xn ? MRIcopy(tsrc->mri_xn, tdst->mri_xn) : NULL;
  tdst->mri_yn = tsrc->mri_yn ? MRIcopy(tsrc->mri_yn, tdst->mri_yn) : NULL;
  tdst->mri_zn = tsrc->mri_zn ? MRIcopy(tsrc->mri_zn, tdst->mri_zn) : NULL;
  switch (tsrc->type) {
    case MORPH_3D_TYPE:
      gcam_src = (GCAM *)tsrc->xform;
      gcam_dst = (GCAM *)tdst->xform;
      tdst->xform = (void *)GCAMcopy(gcam_src, gcam_dst);
      break;
    default:
      lta_src = (LTA *)(tsrc->xform);
      lta_dst = (LTA *)(tdst->xform);
      tdst->xform = (void *)LTAcopy(lta_src, lta_dst);
      break;
  }
  return (tdst);
}

/*-------------------------------------------------------------
  LTAtransformTypeName() - just returns the name of the transform type
  -------------------------------------------------------------*/
const char *LTAtransformTypeName(int ltatype)
{
  switch (ltatype) {
    case LINEAR_VOX_TO_VOX:
      return ("linear_vox_to_vox");
      break;
    case LINEAR_RAS_TO_RAS:
      return ("linear_ras_to_ras");
      break;
    case LINEAR_PHYSVOX_TO_PHYSVOX:
      return ("linear_physvox_to_physvox");
      break;
    case TRANSFORM_ARRAY_TYPE:
      return ("transform_array");
      break;
    case MORPH_3D_TYPE:
      return ("morph_3d");
      break;
    case MNI_TRANSFORM_TYPE:
      return ("mni_transform");
      break;
    case MATLAB_ASCII_TYPE:
      return ("matlab_ascii");
      break;
    case LINEAR_COR_TO_COR:
      return ("linear_cor_to_cor");
      break;
    case REGISTER_DAT:
      return ("register.dat");
      break;
    case FSLREG_TYPE:
      return ("FSL");
      break;
  }
  return ("unknown");
}
/*-------------------------------------------------------------------
  LTAdumpVolGeom() - prints volume geometry to a stream
  -------------------------------------------------------------------*/
int LTAdumpVolGeom(FILE *fp, VG *vg)
{
  fprintf(fp, "valid  %d\n", vg->valid);
  fprintf(fp, "width  %d\n", vg->width);
  fprintf(fp, "height %d\n", vg->height);
  fprintf(fp, "depth  %d\n", vg->depth);
  fprintf(fp, "xsize  %f\n", vg->xsize);
  fprintf(fp, "ysize  %f\n", vg->ysize);
  fprintf(fp, "zsize  %f\n", vg->zsize);
  fprintf(fp, "xdc    %f %f %f\n", vg->x_r, vg->x_a, vg->x_s);
  fprintf(fp, "ydc    %f %f %f\n", vg->y_r, vg->y_a, vg->y_s);
  fprintf(fp, "zdc    %f %f %f\n", vg->z_r, vg->z_a, vg->z_s);
  fprintf(fp, "fname  %s\n", vg->fname);
  return (0);
}
/*-------------------------------------------------------------------
  LTAdumpLinearTransform() - prints LT to stream
  -------------------------------------------------------------------*/
int LTAdumpLinearTransform(FILE *fp, LT *lt)
{
  fprintf(fp, "x0 %f\n", lt->x0);
  fprintf(fp, "y0 %f\n", lt->y0);
  fprintf(fp, "z0 %f\n", lt->z0);
  fprintf(fp, "sigma %f\n", lt->sigma);
  fprintf(fp, "type %d\n", lt->type);
  fprintf(fp, "typename %s\n", LTAtransformTypeName(lt->type));
  fprintf(fp, "label %d\n", lt->label);

  if (lt->m_L) {
    fprintf(fp, "m_L Matrix ----------\n");
    MatrixPrint(fp, lt->m_L);
    fprintf(fp, "---------------------\n");
  }
  fprintf(fp, "Source Geometry ----------\n");
  LTAdumpVolGeom(fp, &lt->src);
  fprintf(fp, "Destination Geometry ----------\n");
  LTAdumpVolGeom(fp, &lt->dst);

  return (0);
}

/*-------------------------------------------------------------------
  LTAdump() - prints LTA to stream
  -------------------------------------------------------------------*/
int LTAdump(FILE *fp, LTA *lta)
{
  int nthxform;

  fprintf(fp, "num_xforms %d\n", lta->num_xforms);
  fprintf(fp, "type %d\n", lta->type);
  fprintf(fp, "typename %s\n", LTAtransformTypeName(lta->type));

  for (nthxform = 0; nthxform < lta->num_xforms; nthxform++) {
    fprintf(fp, "nthxform %d =====================================\n", nthxform);
    LTAdumpLinearTransform(fp, &lta->xforms[nthxform]);
  }
  return (0);
}

int LTAsetVolGeom(LTA *lta, MRI *mri_src, MRI *mri_dst)
{
  int i;

  for (i = 0; i < lta->num_xforms; i++) {
    MRIcopyVolGeomFromMRI(mri_src, &lta->xforms[i].src);
    MRIcopyVolGeomFromMRI(mri_dst, &lta->xforms[i].dst);
  }
  return (NO_ERROR);
}

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

  Note: to compute the matrix with respect to the first voxel being
  at CRS 1,1,1 instead of 0,0,0, then set base = 1. This is
  necessary with SPM matrices.

  See also: MRIxfmCRS2XYZtkreg, MRItkReg2Native, extract_i_to_r().
  surfaceRASFromVoxel_(MRI *mri), voxelFromSurfaceRAS_().

  Note: MRIgetVoxelToRasXform is #defined to be extract_i_to_r().
  ----------------------------------------------------------------*/
// this function is same as MATRIX *MRIxfmCRS2XYZ(const VOL_GEOM *mri, int base)
// call to this function is equivalent to vg->get_Vox2RAS(base);
MATRIX *VGgetVoxelToRasXform(VOL_GEOM *vg, MATRIX *m, int base)
{
  MATRIX *PxyzOffset, *Pcrs;

  if (m == NULL) m = MatrixAlloc(4, 4, MATRIX_REAL);

  /* direction cosine between columns scaled by
     distance between colums */
  *MATRIX_RELT(m, 1, 1) = vg->x_r * vg->xsize;
  *MATRIX_RELT(m, 2, 1) = vg->x_a * vg->xsize;
  *MATRIX_RELT(m, 3, 1) = vg->x_s * vg->xsize;

  /* direction cosine between rows scaled by
     distance between rows */
  *MATRIX_RELT(m, 1, 2) = vg->y_r * vg->ysize;
  *MATRIX_RELT(m, 2, 2) = vg->y_a * vg->ysize;
  *MATRIX_RELT(m, 3, 2) = vg->y_s * vg->ysize;

  /* direction cosine between slices scaled by
     distance between slices */
  *MATRIX_RELT(m, 1, 3) = vg->z_r * vg->zsize;
  *MATRIX_RELT(m, 2, 3) = vg->z_a * vg->zsize;
  *MATRIX_RELT(m, 3, 3) = vg->z_s * vg->zsize;

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
  *MATRIX_RELT(Pcrs, 1, 1) = vg->width / 2.0 + base;
  *MATRIX_RELT(Pcrs, 2, 1) = vg->height / 2.0 + base;
  *MATRIX_RELT(Pcrs, 3, 1) = vg->depth / 2.0 + base;
  *MATRIX_RELT(Pcrs, 4, 1) = 1.0;

  /* XYZ offset the first Col, Row, and Slice from Center */
  /* PxyzOffset = Mdc*D*PcrsCenter */
  PxyzOffset = MatrixMultiply(m, Pcrs, NULL);

  /* XYZ at the Center of the Volume is vg->c_r, c_a, c_s  */

  /* The location of the center of the voxel at CRS = (0,0,0)*/
  *MATRIX_RELT(m, 1, 4) = vg->c_r - PxyzOffset->rptr[1][1];
  *MATRIX_RELT(m, 2, 4) = vg->c_a - PxyzOffset->rptr[2][1];
  *MATRIX_RELT(m, 3, 4) = vg->c_s - PxyzOffset->rptr[3][1];

  MatrixFree(&Pcrs);
  MatrixFree(&PxyzOffset);

  return (m);
}

// call to this function is equivalent to vg->get_RAS2Vox(base)
MATRIX *VGgetRasToVoxelXform(VOL_GEOM *vg, MATRIX *m, int base)
{
  MATRIX *m_inv;

  m_inv = VGgetVoxelToRasXform(vg, NULL, base);
  m = MatrixInverse(m_inv, m);
  MatrixFree(&m_inv);
  return (m);
}
/*!
  \fn LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R)
  \brief Converts a tkregister-style registration matrix to LTA.
   The LTA will be LINEAR_VOX_TO_VOX that maps target vox to
   mov vox. If R=NULL, then computes the LTA based on header geometry.
*/
LTA *TransformRegDat2LTA(MRI *targ, MRI *mov, MATRIX *R)
{
  LTA *lta;
  MATRIX *vox2vox;  // Targ->Mov
  MATRIX *Ttarg, *Tmov, *invTmov;
  int freeR = 0;

  if (R == NULL) {
    Ttarg = MRIxfmCRS2XYZ(targ, 0);
    Tmov = MRIxfmCRS2XYZ(mov, 0);
    R = MatrixIdentity(4, NULL);
    freeR = 1;
  }
  else {
    Ttarg = MRIxfmCRS2XYZtkreg(targ);
    Tmov = MRIxfmCRS2XYZtkreg(mov);
  }
  invTmov = MatrixInverse(Tmov, NULL);

  // vox2vox = invTmov * R * Ttarg
  vox2vox = MatrixMultiply(invTmov, R, NULL);
  MatrixMultiply(vox2vox, Ttarg, vox2vox);

  lta = LTAalloc(1, NULL);
  lta->type = LINEAR_VOX_TO_VOX;
  lta->xforms[0].type = LINEAR_VOX_TO_VOX;
  getVolGeom(targ, &lta->xforms[0].src);
  getVolGeom(mov, &lta->xforms[0].dst);
  lta->xforms[0].m_L = MatrixCopy(vox2vox, NULL);

  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&invTmov);
  MatrixFree(&vox2vox);
  if (freeR) MatrixFree(&R);

  return (lta);
}
/*!
  \fn LTA *TransformLTARegDat(LTA *lta)
  \brief Converts the first xform from LTA to a tkregister-style
         registration matrix. Assumes LTA is vox2vox.
*/
MATRIX *TransformLTA2RegDat(LTA *lta)
{
  MATRIX *Vox2Vox = NULL;           // Mov->Targ
  MATRIX *invVox2Vox = NULL;        // Targ->Mov
  MATRIX *Ttarg, *Tmov, *invTtarg;  // tkreg space vox2ras
  MATRIX *R;                        // TkRegMat

  if (lta->type != LINEAR_RAS_TO_RAS && lta->type != LINEAR_VOX_TO_VOX && lta->type != REGISTER_DAT) {
    printf("ERROR: TransformLTA2RegDat(): type = %d, must be %d or %d or %d\n",
           lta->type,
           LINEAR_RAS_TO_RAS,
           LINEAR_VOX_TO_VOX,
           REGISTER_DAT);

    return (NULL);
  }

  if (lta->type == REGISTER_DAT) {
    R = MatrixCopy(lta->xforms[0].m_L, NULL);
    return (R);
  }

  // Get MovVox-to-TargVox (Vox2Vox)
  if (lta->type == LINEAR_VOX_TO_VOX) Vox2Vox = MatrixCopy(lta->xforms[0].m_L, NULL);
  if (lta->type == LINEAR_RAS_TO_RAS) {
    MATRIX *M;                        // Scanner Space MovRAS-to-TargRAS
    MATRIX *Starg, *Smov, *invStarg;  // Scanner Space vox2ras
    M = MatrixCopy(lta->xforms[0].m_L, NULL);
    Starg = vg_i_to_r(&lta->xforms[0].dst);
    Smov = vg_i_to_r(&(lta->xforms[0].src));
    invStarg = MatrixInverse(Starg, NULL);
    // MovVox2TargVox = inv(Starg)*M*Smov
    Vox2Vox = MatrixMultiply(invStarg, M, NULL);
    Vox2Vox = MatrixMultiply(Vox2Vox, Smov, Vox2Vox);
    MatrixFree(&M);
    MatrixFree(&Smov);
    MatrixFree(&Starg);
    MatrixFree(&invStarg);
  }
  invVox2Vox = MatrixInverse(Vox2Vox, NULL);

  // TkReg Space vox2ras
  Ttarg = TkrVox2RASfromVolGeom(&lta->xforms[0].dst);
  Tmov = TkrVox2RASfromVolGeom(&(lta->xforms[0].src));
  invTtarg = MatrixInverse(Ttarg, NULL);

  // R = Tmov * invVox2Vox * invTtarg
  R = MatrixMultiply(Tmov, invVox2Vox, NULL);
  R = MatrixMultiply(R, invTtarg, R);

  if (Gdiag_no > 0) {
    printf("TransformLTA2RegDat() -----------");
    printf("src/targ Vol Geom");
    lta->xforms[0].src.vgprint();
    printf("dst/mov  Vol Geom");
    lta->xforms[0].dst.vgprint();
    printf("Vox2Vox---------------------------\n");
    MatrixPrint(stdout, Vox2Vox);
    printf("Tmov ---------------------------\n");
    MatrixPrint(stdout, Tmov);
    printf("invTtarg ---------------------------\n");
    MatrixPrint(stdout, invTtarg);
    printf("---------------------------\n");
    MatrixPrint(stdout, R);
    printf("---------------------------\n");
  }

  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&invTtarg);
  MatrixFree(&Vox2Vox);
  MatrixFree(&invVox2Vox);

  return (R);
}
int TransformRas2Vox(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  if (transform->type == MORPH_3D_TYPE)
    return (GCAMrasToVox((GCA_MORPH *)(transform->xform), mri_src));
  else {
    transform->type = LINEAR_VOX_TO_VOX;
    return (LTArasToVoxelXform((LTA *)(transform->xform), mri_src, mri_dst));
  }
}
int TransformVox2Ras(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  if (transform->type == MORPH_3D_TYPE)
    return (GCAMvoxToRas((GCA_MORPH *)(transform->xform)));
  else {
    transform->type = LINEAR_RAS_TO_RAS;
    return (LTAvoxelToRasXform((LTA *)(transform->xform), mri_src, mri_dst));
  }
}
int TransformSampleDirection(TRANSFORM *transform,
                             float x0,
                             float y0,
                             float z0,
                             float nx,
                             float ny,
                             float nz,
                             float *pnx,
                             float *pny,
                             float *pnz)
{
  float xa0, ya0, za0, xa1, ya1, za1, mag;

  TransformSampleReal(transform, x0, y0, z0, &xa0, &ya0, &za0);
  TransformSampleReal(transform, x0 + 2 * nx, y0 + 2 * ny, z0 + 2 * nz, &xa1, &ya1, &za1);
  nx = xa1 - xa0;
  ny = ya1 - ya0;
  nz = za1 - za0;
  mag = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= mag;
  ny /= mag;
  nz /= mag;

  *pnx = nx;
  *pny = ny;
  *pnz = nz;
  return (NO_ERROR);
}

/*!
  \fn MATRIX *MRIangles2RotMat(double *angles)
  \brief Convert 3 euler angles into a 4x4 rotation matrix
  \param angles is a 3x1 vector in radians.
    angles[0] - pitch - rotation about x or LR axis (gamma)
    angles[1] - yaw   - rotation about y or AP axis (beta)
    angles[2] - roll  - rotation about z or SI axis (alpha)
  Ref: Craig, Intro to Robotics
*/
MATRIX *MRIangles2RotMat(double *angles)
{
  double gamma, beta, alpha;
  int r, c;
  MATRIX *R, *R3, *Rx, *Ry, *Rz;

  gamma = angles[0];
  beta = angles[1];
  alpha = angles[2];

  // printf("angles %g %g %g\n",angles[0],angles[1],angles[2]);

  Rx = MatrixZero(3, 3, NULL);
  Rx->rptr[1][1] = +1;
  Rx->rptr[2][2] = +cos(gamma);
  Rx->rptr[2][3] = -sin(gamma);
  Rx->rptr[3][2] = +sin(gamma);
  Rx->rptr[3][3] = +cos(gamma);
  // printf("Rx ----------------\n");
  // MatrixPrint(stdout,Rx);

  Ry = MatrixZero(3, 3, NULL);
  Ry->rptr[1][1] = +cos(beta);
  Ry->rptr[1][3] = +sin(beta);
  Ry->rptr[2][2] = 1;
  Ry->rptr[3][1] = -sin(beta);
  Ry->rptr[3][3] = +cos(beta);
  // printf("Ry ----------------\n");
  // MatrixPrint(stdout,Ry);

  Rz = MatrixZero(3, 3, NULL);
  Rz->rptr[1][1] = +cos(alpha);
  Rz->rptr[1][2] = -sin(alpha);
  Rz->rptr[2][1] = +sin(alpha);
  Rz->rptr[2][2] = +cos(alpha);
  Rz->rptr[3][3] = +1;
  // printf("Rz ----------------\n");
  // MatrixPrint(stdout,Rz);

  // This will be a 3x3 matrix
  R3 = MatrixMultiply(Rz, Ry, NULL);
  R3 = MatrixMultiply(R3, Rx, R3);

  // Stuff 3x3 into a 4x4 matrix, with (4,4) = 1
  R = MatrixZero(4, 4, NULL);
  for (c = 1; c <= 3; c++) {
    for (r = 1; r <= 3; r++) {
      R->rptr[r][c] = R3->rptr[r][c];
    }
  }
  R->rptr[4][4] = 1;

  MatrixFree(&Rx);
  MatrixFree(&Ry);
  MatrixFree(&Rz);
  MatrixFree(&R3);

  // printf("R ----------------\n");
  // MatrixPrint(stdout,R);

  return (R);
}

/*!
  \fn double *SegRegCost(MRI *regseg, MRI *f, double *costs)
  \brief Compute cost function for segmentation-based registration.
  \param regseg - segmentation used to identify WM and Ctx voxels.
    Must be uchar, int, or float
  \param f - volume of intensities used to compute the cost function.
    Must be float (usually is because of resampling).
  \param costs - 8 element array with cost measures. If NULL, it will
    be allocated.
*/
double *SegRegCost(MRI *regseg, MRI *f, double *costs)
{
  double wmsum, wmsum2, wmmean, wmstd;
  double ctxsum, ctxsum2, ctxmean, ctxstd;
  double t, cost;
  int nwmhits, nctxhits;

  if ((regseg->type != MRI_INT) && (regseg->type != MRI_UCHAR) && (regseg->type != MRI_FLOAT)) {
    printf("ERROR: SegRegCost(): regseg type must be int, uchar, or float\n");
    return (NULL);
  }

  if (f->type != MRI_FLOAT) {
    printf("ERROR: SegRegCost(): f type must be float\n");
    return (NULL);
  }

  // Should check that f and regseg have consistent dims

  if (costs == NULL) costs = (double *)calloc(sizeof(double), 8);

  nwmhits = 0;
  nctxhits = 0;
  wmsum = 0;
  wmsum2 = 0;
  ctxsum = 0;
  ctxsum2 = 0;

  for (unsigned int c = 0; c < f->width; c++) {
    for (unsigned int r = 0; r < f->height; r++) {
      for (unsigned int s = 0; s < f->depth; s++) {
        double vf = MRIgetVoxVal(f, c, r, s, 0);
        if (vf == 0) {
          continue;
        }
        double vseg = MRIgetVoxVal(regseg, c, r, s, 0);
        if (vseg == 2 || vseg == 41) {
          // white matter
          wmsum += vf;
          wmsum2 += (vf * vf);
          nwmhits++;
        }
        if (vseg == 3 || vseg == 42) {
          // cortex
          ctxsum += vf;
          ctxsum2 += (vf * vf);
          nctxhits++;
        }
      }
    }
  }

  // printf("wmsum2 = %lf ctxsum2 = %lf\n",wmsum2,ctxsum2);

  wmmean = wmsum / nwmhits;
  wmstd = sum2stddev(wmsum, wmsum2, nwmhits);
  // wmstd = sqrt( (wmsum2 - 2*wmmean*wmsum + nwmhits*wmmean*wmmean)/(nwmhits-1) );

  ctxmean = ctxsum / nctxhits;
  ctxstd = sum2stddev(ctxsum, ctxsum2, nctxhits);
  // ctxstd = sqrt( (ctxsum2 - 2*ctxmean*ctxsum + nctxhits*ctxmean*ctxmean)/nctxhits );

  t = fabs(ctxmean - wmmean) / sqrt(ctxstd * ctxstd + wmstd * wmstd);
  cost = 1 / t;

  // printf("WM: %6d %6.1f %6.1f   CTX: %6d %6.1f %6.1f  Cost: %g\n",
  // nwmhits,wmmean,wmstd, nctxhits,ctxmean,ctxstd, cost);

  costs[0] = nwmhits;
  costs[1] = wmmean;
  costs[2] = wmstd;
  costs[3] = nctxhits;
  costs[4] = ctxmean;
  costs[5] = ctxstd;
  costs[6] = t;
  costs[7] = cost;

  return (0);
}

/*!
  \fn MRI *MRIaffineDisplacment(MRI *mri, MATRIX *R)
  \brief Computes the displacment magnitude at each voxel due to
  applying R (which is a tkreg matrix). The output is a single-frame
  MRI with the unsigned displacment.
*/
MRI *MRIaffineDisplacment(MRI *mri, MATRIX *R)
{
  MRI *disp;
  int c, r, s;
  MATRIX *Pcrs, *Pras, *Pras2, *Vox2RAS, *Vox2RAS2;
  double dx, dy, dz, d;

  disp = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 1);
  MRIcopyHeader(mri, disp);

  Pcrs = MatrixAlloc(4, 1, MATRIX_REAL);
  Pcrs->rptr[4][1] = 1;
  Pras = MatrixAlloc(4, 1, MATRIX_REAL);
  Pras2 = MatrixAlloc(4, 1, MATRIX_REAL);
  Vox2RAS = MRIxfmCRS2XYZtkreg(disp);
  Vox2RAS2 = MatrixMultiply(R, Vox2RAS, NULL);

  for (c = 0; c < disp->width; c++) {
    for (r = 0; r < disp->height; r++) {
      for (s = 0; s < disp->depth; s++) {
        Pcrs->rptr[1][1] = c;
        Pcrs->rptr[2][1] = r;
        Pcrs->rptr[3][1] = s;
        Pras = MatrixMultiply(Vox2RAS, Pcrs, Pras);
        Pras2 = MatrixMultiply(Vox2RAS2, Pcrs, Pras2);
        dx = Pras->rptr[1][1] - Pras2->rptr[1][1];
        dy = Pras->rptr[2][1] - Pras2->rptr[2][1];
        dz = Pras->rptr[3][1] - Pras2->rptr[3][1];
        d = sqrt(dx * dx + dy * dy + dz * dz);
        MRIsetVoxVal(disp, c, r, s, 0, d);
      }
    }
  }

  MatrixFree(&Pcrs);
  MatrixFree(&Pras);
  MatrixFree(&Pras2);
  MatrixFree(&Vox2RAS);
  MatrixFree(&Vox2RAS2);

  return (disp);
}

int TransformGetSrcVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)
{
  GCAM *gcam;
  LTA *lta;

  switch (transform->type) {
    case MORPH_3D_TYPE:
      gcam = (GCA_MORPH *)transform->xform;
      *vg = *(&gcam->image);
      break;
    default:  // linear tranforms
      lta = (LTA *)transform->xform;
      *vg = *(&lta->xforms[0].src);
      break;
  }
  return (NO_ERROR);
}
int TransformGetDstVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)
{
  GCAM *gcam;
  LTA *lta;

  switch (transform->type) {
    case MORPH_3D_TYPE:
      gcam = (GCA_MORPH *)transform->xform;
      *vg = *(&gcam->atlas);
      break;
    default:  // linear tranforms
      lta = (LTA *)transform->xform;
      *vg = *(&lta->xforms[0].dst);
      break;
  }
  return (NO_ERROR);
}
int TransformSetMRIVolGeomToSrc(const TRANSFORM *transform, MRI *mri)
{
  VOL_GEOM vg;

  TransformGetSrcVolGeom(transform, &vg);
  MRIcopyVolGeomToMRI(mri, &vg);
  return (NO_ERROR);
}
int TransformSetMRIVolGeomToDst(const TRANSFORM *transform, MRI *mri)
{
  VOL_GEOM vg;

  TransformGetDstVolGeom(transform, &vg);
  MRIcopyVolGeomToMRI(mri, &vg);
  return (NO_ERROR);
}

LTA *LTAcompose(LTA *lta_src, MATRIX *m_left, MATRIX *m_right, LTA *lta_dst)
{
  if (lta_dst == NULL) {
    lta_dst = LTAalloc(1, NULL);
    copyVolGeom(&lta_src->xforms[0].src, &lta_dst->xforms[0].src);
    copyVolGeom(&lta_src->xforms[0].dst, &lta_dst->xforms[0].dst);
  }
  if (m_left) MatrixMultiply(m_left, lta_src->xforms[0].m_L, lta_dst->xforms[0].m_L);
  if (m_right) MatrixMultiply(lta_src->xforms[0].m_L, m_right, lta_dst->xforms[0].m_L);
  return (lta_dst);
}

TRANSFORM *TransformCompose(TRANSFORM *t_src, MATRIX *m_left, MATRIX *m_right, TRANSFORM *t_dst)
{
  LTA *lta_src, *lta_dst;

  if (t_dst == NULL) t_dst = TransformAlloc(t_src->type, NULL);

  switch (t_src->type) {
    default:
      lta_src = (LTA *)(t_src->xform);
      lta_dst = (LTA *)(t_dst->xform);
      LTAcompose(lta_src, m_left, m_right, lta_dst);
      break;
    case MORPH_3D_TYPE:
      ErrorReturn(NULL, (ERROR_UNSUPPORTED, "TransformCompose: MORPH_3D_TYPE unsupported"));
      break;
  }
  return (t_dst);
}

int TransformSourceVoxelToAtlas(
    TRANSFORM *transform, MRI *mri, int xv, int yv, int zv, double *px, double *py, double *pz)
{
  float xt, yt, zt;
  LTA *lta;

  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, px, py, pz);
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "RFAsourceVoxelToNode: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
    *px = (double)xt;
    *py = (double)yt;
    *pz = (double)zt;
  }
  return (NO_ERROR);
}

/*
  \fn int LTAmriIsSource(const LTA *lta, const MRI *mri)
  \brief Returns 1 if the geometry in mri matches that in
  lta->xforms[0].src, otherwise returns 0.
 */
int LTAmriIsSource(const LTA *lta, const MRI *mri)
{
  VOL_GEOM mrivg;
  int IsSource;
  getVolGeom(mri, &mrivg);
  IsSource = vg_isEqual(&lta->xforms[0].src, &mrivg);
  return (IsSource);
}
/*!
  \fn int LTAmriIsTarget(const LTA *lta, const MRI *mri)
  \brief Returns 1 if the geometry in mri matches that in
  lta->xforms[num_xforms-1].dst, otherwise returns 0.
 */
int LTAmriIsTarget(const LTA *lta, const MRI *mri)
{
  VOL_GEOM mrivg;
  int IsTarget, nxforms;
  getVolGeom(mri, &mrivg);
  nxforms = lta->num_xforms;
  IsTarget = vg_isEqual(&lta->xforms[nxforms - 1].dst, &mrivg);
  return (IsTarget);
}
/*!
  \fn LTA *LTAcreate(MRI *src, MRI *dst, MATRIX *T, int type)
  \brief Create an LTA of the given type with the given matrix.
  Note: when the matrix is of type REGISTER_DAT, the matrix
  should point from destination to source. 
 */
LTA *LTAcreate(MRI *src, MRI *dst, MATRIX *T, int type)
{
  LTA *lta;
  lta = LTAalloc(1, NULL);
  LTAsetVolGeom(lta, src, dst);
  lta->xforms[0].m_L = MatrixCopy(T, NULL);
  lta->xforms[0].type = type;
  lta->type = type;
  if(type == REGISTER_DAT){
    // To help keep me from going crazy
    printf("  ... using LTAcreate() with REGISTER_DAT, make sure matrix points in the right direction\n");
  }
  return (lta);
}

/*!
  \fn LTA *LTAmat2RotMat(LTA *lta)
  \brief Convert LTA into a pure rotation matrix
 */
int LTAmat2RotMat(LTA *lta)
{
  int ltatype = lta->type;
  if(lta->type != LINEAR_RAS_TO_RAS)
    LTAchangeType(lta,LINEAR_RAS_TO_RAS);
  LTAinvert(lta,lta);
  double par[12];
  TranformExtractAffineParams(lta->xforms[0].m_L,par);
  // Only keep the rotational components
  int k;
  for(k=0; k <  3; k++) par[k] = 0; // trans
  for(k=6; k <  9; k++) par[k] = 1; // scale must stay 1
  for(k=9; k < 12; k++) par[k] = 0; // shear
  MATRIX *T = TranformAffineParams2Matrix(par, NULL);
  MatrixCopy(T, lta->xforms[0].m_L);
  MatrixFree(&T);
  LTAinvert(lta,lta);
  if(ltatype != LINEAR_RAS_TO_RAS) LTAchangeType(lta,ltatype);
  return(0);
}

/*!
  \fn double RMSregDiffMJ(MATRIX *T1, MATRIX *T2, double radius)
  \brief Computes the average RMS differences between the two registration
  matrices within a sphere of the given radius. Suggest radius=70.
  Based on Jenkinson, et al, NI 2002 and tech report tr99mj1.pdf.
  Typical usage:
    LTAchangeType(lta1,REGISTER_DAT);
    LTAchangeType(lta2,REGISTER_DAT);
    rms = RMSregDiffMJ(lta1->xforms[0].m_L, lta2->xforms[0].m_L, 70);
 */
double RMSregDiffMJ(MATRIX *T1, MATRIX *T2, double radius)
{
  MATRIX *Q, *M, *t, *Mt, *MtM, *tt, *ttt;
  int c, r;
  double rms, MtMtrace;

  Q = MatrixSubtract(T1, T2, NULL);
  M = MatrixAlloc(3, 3, MATRIX_REAL);
  t = MatrixAlloc(3, 1, MATRIX_REAL);
  for (r = 1; r <= 3; r++) {
    for (c = 1; c <= 3; c++) M->rptr[r][c] = Q->rptr[r][c];
    t->rptr[r][1] = Q->rptr[r][4];
  }
  Mt = MatrixTranspose(M, NULL);
  MtM = MatrixMultiplyD(Mt, M, NULL);
  MtMtrace = MatrixTrace(MtM);

  tt = MatrixTranspose(t, NULL);
  ttt = MatrixMultiplyD(tt, t, NULL);

  rms = sqrt(radius * radius * MtMtrace / 5.0 + ttt->rptr[1][1]);

  MatrixFree(&Q);
  MatrixFree(&M);
  MatrixFree(&Mt);
  MatrixFree(&MtM);
  MatrixFree(&t);
  MatrixFree(&tt);
  MatrixFree(&ttt);

  return (rms);
}

// Concatenate any combination of LTAs and GCAMs. New nemory is allocated.
// Inverse of GCAMs will be freed. LTAs will be inverted if geometries don't
// match, or error. The first transform would be applied to images first
// (to coordinates last).
TRANSFORM *TransformConcat(TRANSFORM** trxArray, unsigned numTrx)
{
  GCAM *gcam;
  LTA *lta;
  TRANSFORM *out;
  TRANSFORM *next;
  if (numTrx == 0) {
    ErrorExit(ERROR_BADPARM, "TransformConcat(): no transform passed\n");
  }

  next = trxArray[--numTrx];
  out = TransformCopy(next, NULL);
  
  while (numTrx > 0) {
    next = trxArray[--numTrx];
    if (out->type == MORPH_3D_TYPE) {
      gcam = (GCAM *)out->xform;
      if (next->type == MORPH_3D_TYPE) {
        GCAMconcat2(/*gcam1*/(GCAM *)next->xform, /*gcam2*/gcam, /*out*/gcam);
      }
      else {
        lta = (LTA *)next->xform;
        out->xform = (void *)GCAMconcat3(lta, gcam, /*lta2*/NULL, /*out*/NULL);
        GCAMfree(&gcam);
      }
      continue;
    }
    
    lta = (LTA *)out->xform;
    if (next->type == MORPH_3D_TYPE) {
      gcam = (GCAM *)next->xform;
      out->xform = (void *)GCAMconcat3(/*lta1*/NULL, gcam, lta, /*out*/NULL);
      out->type = MORPH_3D_TYPE;
    }
    else {
      out->xform = (void *)LTAconcat2((LTA *)next->xform, lta, /*Reduce*/0);
      if (!out->xform) {
        ErrorExit(ERROR_BADPARM, "ERROR: TransformConcat(): LTAs do not match");
      }
    }
    LTAfree(&lta);
  }
  return (out);
}

// Inverts transform in-place, i.e. the original transform is replaced. MRI not
// needed for LTAs.
void TransformInvertReplace(TRANSFORM *transform, const MRI *mri)
{
  LTA *lta;
  GCAM *gcam;
  switch (transform->type) {
    case MORPH_3D_TYPE:
      gcam = (GCAM *)transform->xform;
      if (mri) getVolGeom(mri, /*to*/&gcam->image); // For GCAMfillInverse.
      transform->xform = (void *)GCAMfillInverse(gcam);
      GCAMfree(&gcam);
      break;
    default:
      lta = (LTA *)transform->xform;
      transform->xform = (void *)LTAinvert(lta, NULL);
      LTAfree(&lta);
      break;
  }
}

/*!
  \fn int LTAinversionNeeded(const MRI *src, const MRI *dst, const LTA *lta)
  \brief Determines whether the LTA needs to be inverted or not given a source
  and destination MRI. This is useful when loading in two volumes and a
  transform to make sure the transform goes in the right direction.
  \return
    0 = no inversion needed
    1 = inversion needed
    2 = error: geometries are equal, cannot determine whether inversion needed
    3 = error: source MRI does not match geometry of either LTA source or destination
    4 = error: destination MRI does not match geometry of either LTA source or destination
  Note: it is highly likely that there are going to be some very small differences 
  in the geometry, so it is suggested to vg_isEqual_Threshold = 10e-4;
  See also int LTAinvertIfNeeded()
 */
int LTAinversionNeeded(const MRI *src, const MRI *dst, const LTA *lta)
{
  int InversionNeeded=0;
  VOL_GEOM vgsrc, vgdst;

  getVolGeom(src, &vgsrc);
  getVolGeom(dst, &vgdst);

  if(Gdiag_no > 0){
    printf("LTAinversionNeeded(): vg_isEqual_Threshold=%g\n",vg_isEqual_Threshold);
    printf("mri src ==================================\n");
    writeVolGeom(stdout, &vgsrc);
    printf("mri dst ==================================\n");
    writeVolGeom(stdout, &vgdst);
    printf("lta src ==================================\n");
    writeVolGeom(stdout, &(lta->xforms[0].src));
    printf("lta dst ==================================\n");
    writeVolGeom(stdout, &(lta->xforms[0].dst));
    printf("==================================\n");
  }

  if(vg_isEqual(&vgsrc, &vgdst)){
    printf("WARNING: LTAinversionNeeded(): src and dst MRIs have same vol geom.\n");
    printf("         Cannot determine if an inversion is needed ...\n");
    return(2);
  }
  if(!vg_isEqual(&vgsrc, &(lta->xforms[0].src))){
    if(!vg_isEqual(&vgsrc, &(lta->xforms[0].dst))){
      printf("ERROR: LTAinversionNeeded(): src MRI matches neither src nor dst in LTA (%g)\n",
	     vg_isEqual_Threshold);
      return(3);
    }
    InversionNeeded = 1;
  }
  if(!vg_isEqual(&vgdst, &(lta->xforms[0].dst))){
    if(!vg_isEqual(&vgdst, &(lta->xforms[0].src))){
      printf("ERROR: LTAinversionNeeded(): dst MRI matches neither src nor dst in LTA (%g)\n",
	     vg_isEqual_Threshold);
      return(4);
    }
    InversionNeeded = 1; // this should be redundant
  }
  
  return(InversionNeeded);
}

/*!
  \fn int LTAinvertIfNeeded(const MRI *src, const MRI *dst, LTA *lta)
  \brief Determines whether the LTA needs to be inverted or not given a source
  and destination MRI and performs the inversion. This is useful when loading in two 
  volumes and a transform to make sure the transform goes in the right direction.
  \return
    0 = no inversion needed
    1 = inversion needed
    2 = error: geometries are equal, cannot determine whether inversion needed
    3 = error: source MRI does not match geometry of either LTA source or destination
    4 = error: destination MRI does not match geometry of either LTA source or destination
  Note: it is highly likely that there are going to be some very small differences 
  in the geometry, so it is suggested to vg_isEqual_Threshold = 10e-4;
  See also int LTAinversionNeeded()
 */
int LTAinvertIfNeeded(const MRI *src, const MRI *dst, LTA *lta)
{
  int InversionNeeded = LTAinversionNeeded(src, dst, lta);
  if(InversionNeeded != 1) return(InversionNeeded);
  printf("LTAinvertIfNeeded(): Inverting LTA\n");
  LTAinvert(lta,lta);
  return(InversionNeeded);  
}


/*!
  \fn int TransformCRS2MNI305(const MRI *mri, 
      const double col, const double row, const double slice, 
      const LTA *talxfm, double *R, double *A, double *S)
  \brief Convert a (col, row, slice) from the given volume into an RAS
         in MNI305. The LTA must be the talairach.xfm (as read by
         LTAreadEx(). The volume must have the same geometry as that
         used to generate the talairach.xfm (eg, orig.mgz). Note that
         the RAS is MNI305, not "talairach". The output agrees with
         tkmedit MNI305 coords (but not currently (4/2020) FV, I think
         there is a bug in FV).
 */
int TransformCRS2MNI305(const MRI *mri, const double col, const double row, const double slice, 
			const LTA *talxfm, double *R, double *A, double *S)
{
  MATRIX *Norig = MRIxfmCRS2XYZ(mri,0);

  // M = XFM*Norig
  //MatrixPrint(stdout,Norig);
  //MatrixPrint(stdout,talxfm->xforms[0].m_L);
  MATRIX *M = MatrixMultiplyD(talxfm->xforms[0].m_L,Norig,NULL);
  MATRIX *crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[1][1] = col;
  crs->rptr[2][1] = row;
  crs->rptr[3][1] = slice;
  crs->rptr[4][1] = 1;
  MATRIX *mni305ras = MatrixMultiplyD(M,crs,NULL);
  *R = mni305ras->rptr[1][1];
  *A = mni305ras->rptr[2][1];
  *S = mni305ras->rptr[3][1];

  MatrixFree(&Norig);
  MatrixFree(&crs);
  MatrixFree(&mni305ras);
  MatrixFree(&M);

  //printf("%g %g %g    %g %g %g\n",col,row,slice,*R,*A,*S);

  return(0);
}


/*!
  \fn MATRIX *TranformAffineParams2Matrix(double *p, MATRIX *M)
  \brief Computes a RAS-to-RAS transformation matrix given
    the 12 parameters. p[0-2] translation, p[3-5] rotation
    in degrees, p[6-8] scale, p[9-11] shear. 
    M = T*R1*R2*R3*SCALE*SHEAR
    Consistent with TranformExtractAffineParams()
*/
MATRIX *TranformAffineParams2Matrix(double *p, MATRIX *M)
{
  MATRIX *T, *R1, *R2, *R3, *R, *ZZ, *S;

  // translations
  T = MatrixIdentity(4,NULL);
  T->rptr[1][4] = p[0];
  T->rptr[2][4] = p[1];
  T->rptr[3][4] = p[2];

  // rotations
  R1 = MatrixIdentity(4,NULL);
  R1->rptr[2][2] = cos(p[3]*M_PI/180);
  R1->rptr[2][3] = sin(p[3]*M_PI/180);
  R1->rptr[3][2] = -sin(p[3]*M_PI/180);
  R1->rptr[3][3] = cos(p[3]*M_PI/180);

  R2 = MatrixIdentity(4,NULL);
  R2->rptr[1][1] = cos(p[4]*M_PI/180);
  R2->rptr[1][3] = sin(p[4]*M_PI/180);
  R2->rptr[3][1] = -sin(p[4]*M_PI/180);
  R2->rptr[3][3] = cos(p[4]*M_PI/180);

  R3 = MatrixIdentity(4,NULL);
  R3->rptr[1][1] = cos(p[5]*M_PI/180);
  R3->rptr[1][2] = sin(p[5]*M_PI/180);
  R3->rptr[2][1] = -sin(p[5]*M_PI/180);
  R3->rptr[2][2] = cos(p[5]*M_PI/180);

  // R = R1*R2*R3
  R = MatrixMultiplyD(R1,R2,NULL);
  MatrixMultiplyD(R,R3,R);

  // scale, use ZZ because some idiot #defined Z
  ZZ = MatrixIdentity(4,NULL);
  ZZ->rptr[1][1] = p[6];
  ZZ->rptr[2][2] = p[7];
  ZZ->rptr[3][3] = p[8];
  ZZ->rptr[4][4] = 1;

  // shear
  S = MatrixIdentity(4,NULL);
  S->rptr[1][2] = p[9];
  S->rptr[1][3] = p[10];
  S->rptr[2][3] = p[11];

  // M = T*R*ZZ*S
  M = MatrixMultiplyD(T,R,M);
  MatrixMultiplyD(M,ZZ,M);
  MatrixMultiplyD(M,S,M);
  //MatrixPrint(stdout,M);

  MatrixFree(&T);
  MatrixFree(&R1);
  MatrixFree(&R2);
  MatrixFree(&R3);
  MatrixFree(&R);
  MatrixFree(&ZZ);
  MatrixFree(&S);

  return(M);
}
/*!
  \fn double *TranformExtractAffineParams(MATRIX *M, double *p)
  \brief Extracts parameters from a 12 dof (ie, affine) transformation
  matrix.  This is consistent with TranfromAffineParams2Matrix().
  Angles are in degrees. Handles matrices with a negative determinant
  by making the first scale parameter negative; this is arbitrary
  as any of the scale parameters could have been negated.
 */
double *TranformExtractAffineParams(MATRIX *M, double *p)
{
  int c,r;
  if(M==NULL){
    printf("ERROR: TranformExtractAffineParams() input matrix is NULL\n");
    return(NULL);
  }

  if(p==NULL) p = (double*) calloc(12,sizeof(double));

  // Translation
  p[0] = M->rptr[1][4];
  p[1] = M->rptr[2][4];
  p[2] = M->rptr[3][4];

  // Decompose
  MATRIX *Q = MatrixAlloc(4,4,MATRIX_REAL);
  MATRIX *R = MatrixAlloc(4,4,MATRIX_REAL);
  MatrixQRdecomposition(M,Q,R);

  // Compute a new decomposition such that the diagonal of the
  // R matrix is always positive (these are the scale parameters)
  MATRIX *P = MatrixAlloc(4,4,MATRIX_REAL);
  for(c=1; c <=4; c++){
    for(r=1; r <=4; r++){
      if(r != c) continue;
      P->rptr[r][c] = 1;
      if(R->rptr[r][c] < 0) P->rptr[r][c] = -1;
    }
  }

  // Have to do something when the det<0 because the scales will be >1
  // and so the scale matrix will have det>0 and the rotation matrix
  // will have det=1, so just change the sign of the first scale to
  // preserve the neg determinant. Most real transforms will have det>0
  // anyway; this handles some corner case.
  double det = MatrixDeterminant(M);
  if(det < 0) P->rptr[1][1] *= -1;

  // M = (Q*P)*(P*R) where P*P=I
  //New Q = Q*P
  MatrixMultiply(Q,P,Q);
  // New R = P*R
  MatrixMultiply(P,R,R);
  MatrixFree(&P);

  // Rotation angles (in degrees)
  p[3] = atan2(Q->rptr[2][3],Q->rptr[3][3])*180/M_PI;
  p[4] = atan2(Q->rptr[1][3],sqrt(pow(Q->rptr[2][3],2.0) + pow(Q->rptr[3][3],2.0)))*180/M_PI;
  p[5] = atan2(Q->rptr[1][2],Q->rptr[1][1])*180/M_PI;

  // Scale is the diagonal of the R matrix. Already always positive
  p[6] = R->rptr[1][1];
  p[7] = R->rptr[2][2];
  p[8] = R->rptr[3][3];

  // Shear (this can fail if one of the scale parameters is 0)
  p[ 9] = R->rptr[1][2]/R->rptr[1][1];
  p[10] = R->rptr[1][3]/R->rptr[1][1];
  p[11] = R->rptr[2][3]/R->rptr[2][2];

  MatrixFree(&Q);
  MatrixFree(&R);
  return(p);
}

/*!
  \fn double TransformAffineParamTest(int niters, double thresh)
  \brief This is a test of TranformAffineParams2Matrix() and
  TranformExtractAffineParams() to assure consistency. A random set of
  affine parameters is generated and used to generate a matrix.  The
  affine parameters are then extracted from the matrix and compared
  with the original. The maximum abs diff across all iterations is
  returned. If the max diff for a given iteration exceeds thresh, then
  info is printed out. The consistency is not perfect due to lack of 
  numerical precision (but it is very close).
 */
double TransformAffineParamTest(int niters, double thresh)
{
  MATRIX *M  = MatrixAlloc(4,4,MATRIX_REAL);
  MATRIX *M2 = MatrixAlloc(4,4,MATRIX_REAL);
  double p[12],p2[12],emax,en;
  int k,n,r,c;

  emax = 0;
  for(n=0; n < niters; n++){
    for(k=0; k <  3; k++) p[k] = 10*(drand48()-0.5); // translations (-5 to +5)
    for(k=3; k <  6; k++) p[k] = 2*(drand48()-0.5)*180; // rotations (-180 to +180)
    // Those with prod(scale)<0 will have neg determinant
    for(k=6; k <  9; k++) p[k] = 4*(drand48()-0.5); // scale (-2 to +2)
    for(k=9; k < 12; k++) p[k] = 2*(drand48()-0.5); // shear (-1 to + 1)
    // dont let scale get close to 0
    for(k=6; k <=8; k++) if(fabs(p[k]) < .01) p[k] = 1;
    TranformAffineParams2Matrix(p, M);
    TranformExtractAffineParams(M, p2);
    TranformAffineParams2Matrix(p2, M2);
    en = 0;
    //for(k=0; k < 12; k++) if(en < fabs(p[k]-p2[k])) en = fabs(p[k]-p2[k]);
    for(r=1; r<=3; r++){
      for(c=1; c<=4; c++){
	double en0 = fabs(M->rptr[r][c]-M2->rptr[r][c]);
	if(en < en0) en = en0;
      }
    }
    if(emax < en) emax = en;
    if(en > thresh){
      // This is just a random threshold
      printf("n = %d, en = %g, det = %g ============\n",n,en,MatrixDeterminant(M));
      printf("M = [\n");
      MatrixPrint(stdout,M);
      printf("]\n");
      printf("M2 = [\n");
      MatrixPrint(stdout,M2);
      printf("]\n");
      for(k=0; k < 12; k++) printf("%g ",p[k]);
      printf("\n");
      for(k=0; k < 12; k++) printf("%g ",p2[k]);
      printf("\n");
    }
  }
  printf("n = %d, emax = %g\n",n,emax);
  MatrixFree(&M);
  return(emax);
}


int RegLandmarks::ReadCoords(void)
{
  // stvp is created by mris_apply_reg --stvp via MRISapplyReg() 
  if(txyzfile.compare("stvp") == 0){
    int err = ReadSTVPairFile(sxyzfile);
    return(err);
  }
  if(sxyzfile.find("json") > 0 && txyzfile.find("json") > 0){
    int err = ReadJsonPair(sxyzfile,txyzfile);
    return(err);
  }

  printf("ERROR: reglandmarks only accepts stvp files right now\n");
  return(1);
}
LTA *RegLandmarks::ComputeLTA(void)
{

  if(mrisrc==NULL){
    mrisrc = MRIreadHeader(mrisrcfile.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
    if(mrisrc==NULL) return(NULL);
  }
  if(mritrg==NULL){
    mritrg = MRIreadHeader(mritrgfile.c_str(),MRI_VOLUME_TYPE_UNKNOWN);
    if(mritrg==NULL) return(NULL);
  }
  int coordtype = -100;
  if(coordtypename.compare("RAS")==0) coordtype = LINEAR_RAS_TO_RAS;
  if(coordtypename.compare("VOX")==0) coordtype = LINEAR_VOX_TO_VOX;
  if(coordtypename.compare("TKR")==0) coordtype = REGISTER_DAT;
  if(coordtype == -100){
    printf("ERROR: unrecognized coord type name %s\n",coordtypename.c_str());
    printf("   Expecting RAS, VOX, or TKR\n");
    return(NULL);
  }

  int err = ReadCoords();
  if(err) return(NULL);
  MATRIX *R = ComputeReg(NULL);
  if(coordtype == REGISTER_DAT){
    R = MatrixInverse(R,R);
    if(R==NULL) return(NULL);
  }

  LTA *lta = LTAcreate(mrisrc, mritrg, R, coordtype);
  LTAchangeType(lta,LINEAR_RAS_TO_RAS);

  return(lta);
}
MATRIX *RegLandmarks::ComputeReg(MATRIX *R)
{
  int npoints = sxyz.size();
  MATRIX *X = MatrixAlloc(npoints,4,MATRIX_REAL);
  MATRIX *y = MatrixAlloc(npoints,3,MATRIX_REAL);
  int r,c;
  for(r=0; r < npoints; r++){
    X->rptr[r+1][4] = 1;
    for(c=0; c < 3; c++){
      X->rptr[r+1][c+1] = sxyz[r][c];
      y->rptr[r+1][c+1] = txyz[r][c];
    }
  }
  MATRIX *Xnorm = MatrixNormalizeCol(X,NULL,NULL);
  double Xcond = MatrixNSConditionNumber(Xnorm);
  MatrixFree(&Xnorm);
  printf("normalized condition = %g\n",Xcond);
  MATRIX *Xt = MatrixTranspose(X,NULL);
  MATRIX *XtX = MatrixMultiply(Xt,X,NULL);
  MATRIX *iXtX = MatrixInverse(XtX,NULL);
  if(iXtX==NULL) return(NULL);
  MATRIX *iXtXXt = MatrixMultiply(iXtX,Xt,NULL);
  MATRIX *beta = MatrixMultiply(iXtXXt,y,NULL);
  if(R==NULL) R = MatrixAlloc(4,4,MATRIX_REAL);
  R->rptr[4][4] = 1.0;
  for(r=0; r < 4; r++){
    for(c=0; c < 4; c++){
      if(r < 3) R->rptr[r+1][c+1] = beta->rptr[c+1][r+1];
    }
  }
  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&iXtXXt);
  MatrixFree(&beta);
  return(R);
}
int RegLandmarks::PrintXYZ(FILE *fp, std::vector<std::vector<double>> xyz)
{
  for ( const auto &row : xyz ) {
    for ( const auto &s : row ) {
      fprintf(fp,"%lf ",s);
    }
    fprintf(fp,"\n");
  }
  return(0);
}
int RegLandmarks::ReadSTVPairFile(std::string stvpairfile)
{
  // stvp is created by mris_apply_reg --stvp via MRISapplyReg() 
  std::ifstream ifs;
  ifs.open(stvpairfile.c_str());
  if(ifs.fail()){
    printf("ReadSTVPairFile(): %s %s\n",strerror(errno),stvpairfile.c_str());
    return(1);
  }
  std::string line;
  int sno, tno;
  double sx, sy, sz, tx, ty, tz;
  while(! ifs.eof()){
    std::getline(ifs,line);
    if(line.length()==0) break;
    //printf("%s\n",line.c_str());
    sscanf(line.c_str(),"%d %lf %lf %lf %d %lf %lf %lf",&sno,&sx,&sy,&sz,&tno,&tx,&ty,&tz);
    svtxno.push_back(sno);
    tvtxno.push_back(tno);
    std::vector<double> vrow;
    vrow.push_back(sx);
    vrow.push_back(sy);
    vrow.push_back(sz);
    sxyz.push_back(vrow);
    vrow.clear();
    vrow.push_back(tx);
    vrow.push_back(ty);
    vrow.push_back(tz);
    txyz.push_back(vrow);
  }
  return(0);
}
int RegLandmarks::ReadJsonPair(std::string jsonfile1,std::string jsonfile2)
{
  fsPointSet ps1 = loadfsPointSet(jsonfile1);
  fsPointSet ps2 = loadfsPointSet(jsonfile2);

  if(ps1.npoints() != ps2.npoints()){
    printf("ERROR: RegLandmarks::ReadJsonPair(): dim mismatch %d %d\n",ps1.npoints(),ps2.npoints());
    return(1);
  }

  for(int n=0; n < ps1.npoints(); n++){
    fsPointSet::Point p1 = ps1.get(n);
    fsPointSet::Point p2 = ps2.get(n);
    svtxno.push_back(n);
    tvtxno.push_back(n);
    std::vector<double> vrow;
    vrow.push_back(p1.x);
    vrow.push_back(p1.y);
    vrow.push_back(p1.z);
    sxyz.push_back(vrow);
    vrow.clear();
    vrow.push_back(p2.x);
    vrow.push_back(p2.y);
    vrow.push_back(p2.z);
    txyz.push_back(vrow);
  }
  return(0);
}

int RegLandmarks::WriteCoords(const char *fname, const LTA *lta)
{
  FILE *fp = fopen(fname,"w");
  if(fp==NULL) return(1);
  for(int n=0; n < sxyz.size(); n++){
    fprintf(fp,"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf",
	    sxyz[n][0],sxyz[n][1],sxyz[n][2],txyz[n][0],txyz[n][1],txyz[n][2]);
    if(lta == NULL){
      fprintf(fp,"\n");
      continue;
    }
    MATRIX *ras = MatrixAlloc(4,1,MATRIX_REAL); ras->rptr[4][1] = 1;
    for(int k=0; k < 3; k++) ras->rptr[k+1][1] = sxyz[n][k];
    MATRIX *newras = MatrixMultiplyD(lta->xforms[0].m_L,ras,NULL);
    fprintf(fp," %9.4lf %9.4lf %9.4lf\n",newras->rptr[1][1],newras->rptr[2][1],newras->rptr[3][1]);
    MatrixFree(&ras);
    MatrixFree(&newras);
  }
  return(0);
}
fsPointSet RegLandmarks::Coords2PointSet(const LTA *lta)
{
  fsPointSet ps;
  for(int n=0; n < sxyz.size(); n++){
    MATRIX *ras = MatrixAlloc(4,1,MATRIX_REAL); ras->rptr[4][1] = 1;
    for(int k=0; k < 3; k++) ras->rptr[k+1][1] = sxyz[n][k];
    MATRIX *newras = MatrixMultiplyD(lta->xforms[0].m_L,ras,NULL);
    fsPointSet::Point p;
    p.x = newras->rptr[1][1];
    p.y = newras->rptr[2][1];
    p.z = newras->rptr[3][1];
    ps.add(p);
    MatrixFree(&ras);
    MatrixFree(&newras);
  }
  return(ps);
}

