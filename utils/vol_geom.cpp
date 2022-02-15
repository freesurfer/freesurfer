#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "log.h"
#include "vol_geom.h"

// utils/transform.cpp:void initVolGeom(VOL_GEOM *vg)
VOL_GEOM::VOL_GEOM(const char *srcVol)
{
  vox2ras = NULL;
  ras2vox = NULL;
  tkregvox2ras = NULL;
  tkregras2vox = NULL;

  if (srcVol == NULL)
  {
    init();
  }
  else
  {
    struct stat stat_buf;

    // check the existence of a file
    int ret = stat(srcVol, &stat_buf);
    if (ret != 0)
    {
      // srcVol doesn't exist
      if (strstr(srcVol, "average_305"))
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) printf("INFO: The transform was made with average_305.mnc.\n");

        init_average_305();
      }
      else
        init();
    }
    else  // srcVol exists
    {
      // note that both mri volume but also gca can be read
      MRI *mri = MRIreadHeader((char *)srcVol, MRI_VOLUME_TYPE_UNKNOWN);
      if (mri)  // find the MRI volume
        copyFromMRI(mri);
      else  // cound not find the volume
        init();
    }

    // copy filename even if file does not exist to keep track
    strcpy(fname, srcVol);
  }
}


VOL_GEOM::~VOL_GEOM()
{
  if (vox2ras != NULL)
    MatrixFree(&vox2ras);

  if (ras2vox != NULL)
    MatrixFree(&ras2vox);

  if (tkregvox2ras != NULL)
    MatrixFree(&tkregvox2ras);

  if (tkregras2vox != NULL)
    MatrixFree(&tkregras2vox);
}

// utils/transform.cpp:void initVolGeom(VOL_GEOM *vg)
void VOL_GEOM::init()
{
  valid = 0;
  width = 256;
  height = 256;
  depth = 256;
  xsize = 1;
  ysize = 1;
  zsize = 1;
  x_r = -1.;
  x_a = 0.;
  x_s = 0.;
  y_r = 0.;
  y_a = 0.;
  y_s = -1.;
  z_r = 0.;
  z_a = 1.;
  z_s = 0.;
  c_r = 0.;
  c_a = 0.;
  c_s = 0.;
  strcpy(fname, "unknown");  // initialized to be "unknown"
}

void VOL_GEOM::init_average_305()
{
  // average_305 value
  width = 172;
  height = 220;
  depth = 156;
  xsize = 1;
  ysize = 1;
  zsize = 1;
  x_r = 1;
  x_a = 0;
  x_s = 0;
  y_r = 0;
  y_a = 1;
  y_s = 0;
  z_r = 0;
  z_a = 0;
  z_s = 1;
  c_r = -0.0950;
  c_a = -16.5100;
  c_s = 9.7500;
  valid = 1;
}

// utils/transform.cpp:void vg_print(const VOL_GEOM *vg)
void VOL_GEOM::print()
{
  if (valid == 1) {
    fprintf(stderr, "volume geometry:\n");
    fprintf(stderr, "extent  : (%d, %d, %d)\n", width, height, depth);
    fprintf(stderr, "voxel   : (%7.4f, %7.4f, %7.4f)\n", xsize, ysize, zsize);
    fprintf(stderr, "x_(ras) : (%7.4f, %7.4f, %7.4f)\n", x_r, x_a, x_s);
    fprintf(stderr, "y_(ras) : (%7.4f, %7.4f, %7.4f)\n", y_r, y_a, y_s);
    fprintf(stderr, "z_(ras) : (%7.4f, %7.4f, %7.4f)\n", z_r, z_a, z_s);
    fprintf(stderr, "c_(ras) : (%7.4f, %7.4f, %7.4f)\n", c_r, c_a, c_s);
    fprintf(stderr, "file    : %s\n", fname);
  }
  else
    fprintf(stderr,
            "volume geometry info is either not contained "
            "or not valid.\n");
  fflush(stderr);
}

// utils/transform.cpp:void getVolGeom(const MRI *src, VOL_GEOM *dst)
void VOL_GEOM::copyFromMRI(const MRI *src)
{
  if (!src) ErrorExit(ERROR_BADPARM, "must have a valid MRI (src)");

  valid = 1;
  width = src->width;
  height = src->height;
  depth = src->depth;
  xsize = src->xsize;
  ysize = src->ysize;
  zsize = src->zsize;
  x_r = src->x_r;
  x_a = src->x_a;
  x_s = src->x_s;
  y_r = src->y_r;
  y_a = src->y_a;
  y_s = src->y_s;
  z_r = src->z_r;
  z_a = src->z_a;
  z_s = src->z_s;
  c_r = src->c_r;
  c_a = src->c_a;
  c_s = src->c_s;
  strcpy(fname, src->fname);

  vox2ras = NULL; ras2vox = NULL;
  getVox2RAS();
  getRAS2Vox();
  
}

// utils/transform.cpp:void useVolGeomToMRI(const VOL_GEOM *src, MRI *dst)
void VOL_GEOM::copyToMRI(MRI *dst)
{
  if (!dst) ErrorExit(ERROR_BADPARM, "must have a valid MRI (dst)");

  dst->ras_good_flag = 1;
  dst->width = width;
  dst->height = height;
  dst->depth = depth;
  dst->xsize = xsize;
  dst->ysize = ysize;
  dst->zsize = zsize;
  dst->x_r = x_r;
  dst->x_a = x_a;
  dst->x_s = x_s;
  dst->y_r = y_r;
  dst->y_a = y_a;
  dst->y_s = y_s;
  dst->z_r = z_r;
  dst->z_a = z_a;
  dst->z_s = z_s;
  dst->c_r = c_r;
  dst->c_a = c_a;
  dst->c_s = c_s;
  strcpy(dst->fname, fname);
  // now we cache transform and thus we have to do the following whenever
  // we change direction cosines
  MRIreInitCache(dst);
}

// utils/transform.cpp:void writeVolGeom(FILE *fp, const VOL_GEOM *vg)
void VOL_GEOM::write(FILE *fp)
{
  if (valid == 0)
    fprintf(fp, "valid = %d  # volume info invalid\n", valid);
  else
    fprintf(fp, "valid = %d  # volume info valid\n", valid);
  fprintf(fp, "filename = %s\n", fname);
  fprintf(fp, "volume = %d %d %d\n", width, height, depth);
  fprintf(fp, "voxelsize = %.15e %.15e %.15e\n", xsize, ysize, zsize);
  fprintf(fp, "xras   = %.15e %.15e %.15e\n", x_r, x_a, x_s);
  fprintf(fp, "yras   = %.15e %.15e %.15e\n", y_r, y_a, y_s);
  fprintf(fp, "zras   = %.15e %.15e %.15e\n", z_r, z_a, z_s);
  fprintf(fp, "cras   = %.15e %.15e %.15e\n", c_r, c_a, c_s);
}

// utils/transform.cpp:void readVolGeom(FILE *fp, VOL_GEOM *vg)
void VOL_GEOM::read(FILE *fp)
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
      sscanf(line, "%s %s %d \n", param, eq, &valid);
      vgRead = 1;
      counter++;
    }
    else if (!strcmp(param, "filename")) {
      // 1023 = STRLEN - 1
      if (sscanf(line, "%s %s %1023s%*s\n", param, eq, buf) >= 3) {
        strcpy(fname, buf);
      }
      counter++;
    }
    else if (!strcmp(param, "volume")) {
      // rescan again
      sscanf(line, "%s %s %d %d %d\n", param, eq, &width, &height, &depth);
      counter++;
    }
    else if (!strcmp(param, "voxelsize")) {
      // rescan again
      sscanf(line, "%s %s %f %f %f\n", param, eq, &xsize, &ysize, &zsize);
      counter++;
    }
    else if (!strcmp(param, "xras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &x_r, &x_a, &x_s);
      counter++;
    }
    else if (!strcmp(param, "yras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &y_r, &y_a, &y_s);
      counter++;
    }
    else if (!strcmp(param, "zras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &z_r, &z_a, &z_s);
      counter++;
    }
    else if (!strcmp(param, "cras")) {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &c_r, &c_a, &c_s);
      counter++;
    }
    // remember the current position
    pos = ftell(fp);  // if fail = 0, then ok
  }
  if (p)  // we read one more line
  {
    if (pos > 0)                        // if success in getting pos, then
      fail = fseek(fp, pos, SEEK_SET);  // restore the position
    // note that this won't allow compression using pipe
  }
  if (!vgRead) {
    fprintf(stderr, "INFO: volume info was not present.\n");
    init();
  }
}

// utils/transform.cpp:MATRIX *vg_i_to_r(const VOL_GEOM *vg)
// #define vg_getVoxelToRasXform vg_i_to_r
// scanner space vox2ras from vol geom
#if 0
MATRIX* VOL_GEOM::getVox2RAS()
{ 
  MATRIX *mat = 0;
  MRI *tmp = 0;
  tmp = MRIallocHeader(width, height, depth, MRI_UCHAR, 1);
  copyToMRI(tmp);
  mat = extract_i_to_r(tmp);
  MRIfree(&tmp);
  return mat;
}

// utils/transform.cpp:MATRIX *vg_r_to_i(const VOL_GEOM *vg)
// #define vg_getRasToVoxelXform vg_r_to_i
MATRIX* VOL_GEOM::getRAS2Vox()
{
  MATRIX *mat = 0;
  MRI *tmp = 0;
  tmp = MRIallocHeader(width, height, depth, MRI_UCHAR, 1);
  copyToMRI(tmp);
  mat = extract_r_to_i(tmp);
  MRIfree(&tmp);
  return mat;
}

// utils/transform.cpp:MATRIX *TkrVox2RASfromVolGeom(const VOL_GEOM *vg)
// tkregister space vox2ras from vol geom
MATRIX* VOL_GEOM::getTkregVox2RAS()
{
  MATRIX *mat = NULL;
  MRI *tmp = 0;
  tmp = MRIallocHeader(width, height, depth, MRI_UCHAR, 1);
  copyToMRI(tmp);
  mat = MRIxfmCRS2XYZtkreg(tmp);
  MRIfree(&tmp);
  return (mat);
}

// utils/transform.cpp:MATRIX *TkrRAS2VoxfromVolGeom(const VOL_GEOM *vg)
// tkregister space ras2vox from vol geom
MATRIX* VOL_GEOM::getTkregRAS2Vox()
{ 
  MATRIX *mat = NULL;
  mat = getTkrVox2RAS();
  mat = MatrixInverse(mat, mat);
  return (mat);                                                 
}
#endif

// utils/transform.cpp:int vg_isEqual(const VOL_GEOM *vg1, const VOL_GEOM *vg2)
int VOL_GEOM::isEqual(const VOL_GEOM* vg2)
{
  int rt;
  extern double vg_isEqual_Threshold;
  // rt = vg_isEqualThresh(vg1, vg2, FLT_EPSILON);
  rt = isNotEqualThresh(vg2, vg_isEqual_Threshold);
  if (rt == 0)
    return (1);
  else
    return (0);
}

// utils/transform.cpp:int vg_isNotEqualThresh(const VOL_GEOM *vg1, const VOL_GEOM *vg2, const double thresh)
int VOL_GEOM::isNotEqualThresh(const VOL_GEOM *vg2, const double thresh)
{
  if (valid != vg2->valid) return (1);
  if (width != vg2->width) return (2);
  if (height != vg2->height) return (3);
  if (depth != vg2->depth) return (4);
  if (!FZEROTHR(xsize - vg2->xsize, thresh)) return (5);
  if (!FZEROTHR(ysize - vg2->ysize, thresh)) return (6);
  if (!FZEROTHR(zsize - vg2->zsize, thresh)) return (7);
  if (!FZEROTHR(x_r - vg2->x_r, thresh)) return (8);
  if (!FZEROTHR(x_a - vg2->x_a, thresh)) return (9);
  if (!FZEROTHR(x_s - vg2->x_s, thresh)) return (10);
  if (!FZEROTHR(y_r - vg2->y_r, thresh)) return (11);
  if (!FZEROTHR(y_a - vg2->y_a, thresh)) return (12);
  if (!FZEROTHR(y_s - vg2->y_s, thresh)) return (13);
  if (!FZEROTHR(z_r - vg2->z_r, thresh)) return (14);
  if (!FZEROTHR(z_a - vg2->z_a, thresh)) return (15);
  if (!FZEROTHR(z_s - vg2->z_s, thresh)) return (16);
  if (!FZEROTHR(c_r - vg2->c_r, thresh)) return (17);
  if (!FZEROTHR(c_a - vg2->c_a, thresh)) return (18);
  if (!FZEROTHR(c_s - vg2->c_s, thresh)) return (19);
  return (0);
}

// utils/transform.cpp:MATRIX *VGgetVoxelToRasXform(VOL_GEOM *vg, MATRIX *m, int base)
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
//MATRIX* VOL_GEOM::getVoxelToRasXform(int base)
MATRIX* VOL_GEOM::getVox2RAS(int base)
{
  if (vox2ras != NULL)
    return vox2ras;

  MATRIX *PxyzOffset, *Pcrs;

  vox2ras = MatrixAlloc(4, 4, MATRIX_REAL);

  /* direction cosine between columns scaled by
     distance between colums */
  *MATRIX_RELT(vox2ras, 1, 1) = (double)x_r * xsize;
  *MATRIX_RELT(vox2ras, 2, 1) = (double)x_a * xsize;
  *MATRIX_RELT(vox2ras, 3, 1) = (double)x_s * xsize;

  /* direction cosine between rows scaled by
     distance between rows */
  *MATRIX_RELT(vox2ras, 1, 2) = (double)y_r * ysize;
  *MATRIX_RELT(vox2ras, 2, 2) = (double)y_a * ysize;
  *MATRIX_RELT(vox2ras, 3, 2) = (double)y_s * ysize;

  /* direction cosine between slices scaled by
     distance between slices */
  *MATRIX_RELT(vox2ras, 1, 3) = (double)z_r * zsize;
  *MATRIX_RELT(vox2ras, 2, 3) = (double)z_a * zsize;
  *MATRIX_RELT(vox2ras, 3, 3) = (double)z_s * zsize;

  /* Preset the offsets to 0 */
  *MATRIX_RELT(vox2ras, 1, 4) = 0.0;
  *MATRIX_RELT(vox2ras, 2, 4) = 0.0;
  *MATRIX_RELT(vox2ras, 3, 4) = 0.0;

  /* Last row of matrix */
  *MATRIX_RELT(vox2ras, 4, 1) = 0.0;
  *MATRIX_RELT(vox2ras, 4, 2) = 0.0;
  *MATRIX_RELT(vox2ras, 4, 3) = 0.0;
  *MATRIX_RELT(vox2ras, 4, 4) = 1.0;

  /* At this point, m = Mdc * D */

  /* Col, Row, Slice at the Center of the Volume */
  Pcrs = MatrixAlloc(4, 1, MATRIX_REAL);
  *MATRIX_RELT(Pcrs, 1, 1) = (double)width / 2.0 + base;
  *MATRIX_RELT(Pcrs, 2, 1) = (double)height / 2.0 + base;
  *MATRIX_RELT(Pcrs, 3, 1) = (double)depth / 2.0 + base;
  *MATRIX_RELT(Pcrs, 4, 1) = 1.0;

  /* XYZ offset the first Col, Row, and Slice from Center */
  /* PxyzOffset = Mdc*D*PcrsCenter */
  PxyzOffset = MatrixMultiply(vox2ras, Pcrs, NULL);

  /* XYZ at the Center of the Volume is c_r, c_a, c_s  */

  /* The location of the center of the voxel at CRS = (0,0,0)*/
  *MATRIX_RELT(vox2ras, 1, 4) = (double)c_r - PxyzOffset->rptr[1][1];
  *MATRIX_RELT(vox2ras, 2, 4) = (double)c_a - PxyzOffset->rptr[2][1];
  *MATRIX_RELT(vox2ras, 3, 4) = (double)c_s - PxyzOffset->rptr[3][1];

  MatrixFree(&Pcrs);
  MatrixFree(&PxyzOffset);

  return (vox2ras);
}

//utils/transform.cpp:MATRIX *VGgetRasToVoxelXform(VOL_GEOM *vg, MATRIX *m, int base)
//MATRIX* VOL_GEOM::getRasToVoxelXform(int base)
MATRIX* VOL_GEOM::getRAS2Vox(int base)
{
  if (ras2vox != NULL)
    return ras2vox;

  if (vox2ras == NULL)
    vox2ras = getVox2RAS(base);

  ras2vox = MatrixInverse(vox2ras, ras2vox);

  return (ras2vox);
}


// tkregister space vox2ras from vol geom
// MATRIX *TkrVox2RASfromVolGeom(const VOL_GEOM *vg)
// MATRIX *MRIxfmCRS2XYZtkreg(const MRI *mri)
MATRIX* VOL_GEOM::getTkregVox2RAS(int base)
{
  if (tkregvox2ras != NULL)
    return tkregvox2ras;
  
  tkregvox2ras = MatrixAlloc(4, 4, MATRIX_REAL);

  /* Set tkregister defaults */
  /* column         row           slice          center      */
  float x_r_tkreg = -1;
  float y_r_tkreg = 0;
  float z_r_tkreg = 0;
  float c_r_tkreg = 0.0;
  float x_a_tkreg = 0;
  float y_a_tkreg = 0;
  float z_a_tkreg = 1;
  float c_a_tkreg = 0.0;
  float x_s_tkreg = 0;
  float y_s_tkreg = -1;
  float z_s_tkreg = 0;
  float c_s_tkreg = 0.0;

  /* direction cosine between columns scaled by
     distance between colums */
  *MATRIX_RELT(tkregvox2ras, 1, 1) = (double)x_r_tkreg * xsize;
  *MATRIX_RELT(tkregvox2ras, 2, 1) = (double)x_a_tkreg * xsize;
  *MATRIX_RELT(tkregvox2ras, 3, 1) = (double)x_s_tkreg * xsize;

  /* direction cosine between rows scaled by
     distance between rows */
  *MATRIX_RELT(tkregvox2ras, 1, 2) = (double)y_r_tkreg * ysize;
  *MATRIX_RELT(tkregvox2ras, 2, 2) = (double)y_a_tkreg * ysize;
  *MATRIX_RELT(tkregvox2ras, 3, 2) = (double)y_s_tkreg * ysize;

  /* direction cosine between slices scaled by
     distance between slices */
  *MATRIX_RELT(tkregvox2ras, 1, 3) = (double)z_r_tkreg * zsize;
  *MATRIX_RELT(tkregvox2ras, 2, 3) = (double)z_a_tkreg * zsize;
  *MATRIX_RELT(tkregvox2ras, 3, 3) = (double)z_s_tkreg * zsize;

  /* Preset the offsets to 0 */
  *MATRIX_RELT(tkregvox2ras, 1, 4) = 0.0;
  *MATRIX_RELT(tkregvox2ras, 2, 4) = 0.0;
  *MATRIX_RELT(tkregvox2ras, 3, 4) = 0.0;

  /* Last row of matrix */
  *MATRIX_RELT(tkregvox2ras, 4, 1) = 0.0;
  *MATRIX_RELT(tkregvox2ras, 4, 2) = 0.0;
  *MATRIX_RELT(tkregvox2ras, 4, 3) = 0.0;
  *MATRIX_RELT(tkregvox2ras, 4, 4) = 1.0;

  /* At this point, m = Mdc * D */

  /* Col, Row, Slice at the Center of the Volume */
  MATRIX *Pcrs = MatrixAlloc(4, 1, MATRIX_REAL);
  *MATRIX_RELT(Pcrs, 1, 1) = (double)width / 2.0 + base;
  *MATRIX_RELT(Pcrs, 2, 1) = (double)height / 2.0 + base;
  *MATRIX_RELT(Pcrs, 3, 1) = (double)depth / 2.0 + base;
  *MATRIX_RELT(Pcrs, 4, 1) = 1.0;

  /* XYZ offset the first Col, Row, and Slice from Center */
  /* PxyzOffset = Mdc*D*PcrsCenter */
  MATRIX *PxyzOffset = MatrixMultiply(tkregvox2ras, Pcrs, NULL);

  /* XYZ at the Center of the Volume is c_r, c_a, c_s  */

  /* The location of the center of the voxel at CRS = (0,0,0)*/
  *MATRIX_RELT(tkregvox2ras, 1, 4) = (double)c_r_tkreg - PxyzOffset->rptr[1][1];
  *MATRIX_RELT(tkregvox2ras, 2, 4) = (double)c_a_tkreg - PxyzOffset->rptr[2][1];
  *MATRIX_RELT(tkregvox2ras, 3, 4) = (double)c_s_tkreg - PxyzOffset->rptr[3][1];

  MatrixFree(&Pcrs);
  MatrixFree(&PxyzOffset);

  return (tkregvox2ras);
}

//utils/transform.cpp:MATRIX *VGgetRasToVoxelXform(VOL_GEOM *vg, MATRIX *m, int base)
//MATRIX* VOL_GEOM::getRasToVoxelXform(int base)
MATRIX* VOL_GEOM::getTkregRAS2Vox(int base)
{
  if (tkregras2vox != NULL)
    return tkregras2vox;

  if (tkregvox2ras == NULL)
    tkregvox2ras = getTkregVox2RAS(base);

  tkregras2vox = MatrixInverse(tkregvox2ras, tkregras2vox);

  return tkregras2vox;
}

// utils/transform.cpp:MRI *MRIallocFromVolGeom(VOL_GEOM *vg, int type, int nframes, int HeaderOnly)
// Creates an MRI from this VOL_GEOM, copying the geometry info
MRI* VOL_GEOM::allocMRI(int type, int nframes, int HeaderOnly)
{
  MRI *mri;
  if (HeaderOnly)
    mri = MRIallocHeader(width, height, depth, type, nframes);
  else
    mri = MRIallocSequence(width, height, depth, type, nframes);
  if (mri) copyToMRI(mri);
  return (mri);
}

// utils/transform.cpp:void copyVolGeom(const VOL_GEOM *src, VOL_GEOM *dst)
void  VOL_GEOM::copy(VOL_GEOM *dst)
{
  dst->valid = valid;
  dst->width = width;
  dst->height = height;
  dst->depth = depth;
  dst->xsize = xsize;
  dst->ysize = ysize;
  dst->zsize = zsize;
  dst->x_r = x_r;
  dst->x_a = x_a;
  dst->x_s = x_s;
  dst->y_r = y_r;
  dst->y_a = y_a;
  dst->y_s = y_s;
  dst->z_r = z_r;
  dst->z_a = z_a;
  dst->z_s = z_s;
  dst->c_r = c_r;
  dst->c_a = c_a;
  dst->c_s = c_s;
  strcpy(dst->fname, fname);

  // copy transform matrix too???
  // ...
}

// VOL_GEOM created with srcVol
// utils/transform.cpp:int mincFindVolume(const char *line, const char *line2, char **srcVol, char **dstVol)
// utils/transform.cpp:void mincGetVolumeInfo(const char *srcVol, VOL_GEOM *vgSrc)
// utils/transform.cpp:void mincGetVolInfo(const char *infoline, const char *infoline2, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)

// utils/transform.cpp:static void LTAgetV2V(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
// utils/transform.cpp:static void LTAgetR2R(MATRIX *mod, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)



// utils/transform.cpp:int TransformGetSrcVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)
// utils/transform.cpp:int TransformGetDstVolGeom(const TRANSFORM *transform, VOL_GEOM *vg)

// utils/mri.cpp:int MRIcopyVolGeomToMRI(MRI *mri, const VOL_GEOM *vg)
// utils/mri.cpp:int MRIcopyVolGeomFromMRI(const MRI *mri, VOL_GEOM *vg)


// mri_modify/mri_modify.cpp:int get_option(int argc, char *argv[], VOL_GEOM &vg)

// utils/gca.cpp:void GCAsetVolGeom(GCA *gca, VOL_GEOM *vg)
