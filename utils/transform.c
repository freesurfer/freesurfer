/*
 *       FILE NAME:   transform.c
 *
 *       DESCRIPTION: utilities for linear transforms
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        11/98
 *
*/

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "gcamorph.h"
#include "mri.h"
#include "mrinorm.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "fio.h"
#include "proto.h"
#include "matrix.h"
#include "transform.h"

#define MAX_TRANSFORMS (1024*4)

static LTA  *ltaReadRegisterDat(char *fname) ;
static LTA  *ltaMNIread(char *fname) ;
static int  ltaMNIwrite(LTA *lta, char *fname) ;
static LTA  *ltaReadFile(char *fname) ;

static LTA *ltaMNIreadEx(const char *fname);
static LTA *ltaReadFileEx(const char *fname);

void initVolGeom(VOL_GEOM *vg);
void writeVolGeom(FILE *fp, const VOL_GEOM *vg);
void readVolGeom(FILE *fp, VOL_GEOM *vg);

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
  vg->x_r = -1.; vg->x_a = 0.; vg->x_s =  0.;
  vg->y_r =  0.; vg->y_a = 0.; vg->y_s = -1.;
  vg->z_r =  0.; vg->z_a = 1.; vg->z_s =  0.;
  vg->c_r =  0.; vg->c_a = 0.; vg->c_s =  0.;
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
  dst->x_r = src->x_r; dst->x_a = src->x_a; dst->x_s = src->x_s;
  dst->y_r = src->y_r; dst->y_a = src->y_a; dst->y_s = src->y_s;
  dst->z_r = src->z_r; dst->z_a = src->z_a; dst->z_s = src->z_s;
  dst->c_r = src->c_r; dst->c_a = src->c_a; dst->c_s = src->c_s;
}
 
void writeVolGeom(FILE *fp, const VOL_GEOM *vg)
{
  if (vg->valid==0)
    fprintf(fp, "valid = %d  # volume info invalid\n", vg->valid);
  else
    fprintf(fp, "valid = %d  # volume info valid\n", vg->valid);
  fprintf(fp, "volume = %d %d %d\n", vg->width, vg->height, vg->depth);
  fprintf(fp, "voxelsize = %.4f %.4f %.4f\n", vg->xsize, vg->ysize, vg->zsize);
  fprintf(fp, "xras   = %.4f %.4f %.4f\n", vg->x_r, vg->x_a, vg->x_s);
  fprintf(fp, "yras   = %.4f %.4f %.4f\n", vg->y_r, vg->y_a, vg->y_s);
  fprintf(fp, "zras   = %.4f %.4f %.4f\n", vg->z_r, vg->z_a, vg->z_s);
  fprintf(fp, "cras   = %.4f %.4f %.4f\n", vg->c_r, vg->c_a, vg->c_s);
}

void readVolGeom(FILE *fp, VOL_GEOM *vg)
{
  char line[256];
  char param[64];
  char eq[2];

  
  if (fgets(line, sizeof(line), fp))
  {
    sscanf(line, "%s %s %d \n", param, eq, &vg->valid);
    if (strcmp(param, "valid")!=0)
    {
      fprintf(stderr, "INFO: parameter was %s. should be 'valid'\n", param);
      initVolGeom(vg);
      return;
    }
    if (fgets(line, sizeof(line), fp))
    {
      sscanf(line, "%s %s %d %d %d\n", param, eq, &vg->width, &vg->height, &vg->depth);
      if (strcmp(param, "volume")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'volume'\n", param);
	initVolGeom(vg);
	return;
      }
    }
    if (fgets(line, sizeof(line), fp))  
    {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->xsize, &vg->ysize, &vg->zsize);
      if (strcmp(param, "voxelsize")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'voxelsize'\n", param);
	initVolGeom(vg);
	return;
      }
    }
    if (fgets(line, sizeof(line), fp))  
    {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->x_r, &vg->x_a, &vg->x_s);
      if (strcmp(param, "xras")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'xras'\n", param);
	initVolGeom(vg);
	return;
      }
    }
    if (fgets(line, sizeof(line), fp))  
    {
      sscanf(line, "%s %s %f %f %f\n", param, eq,  &vg->y_r, &vg->y_a, &vg->y_s);
      if (strcmp(param, "yras")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'yras'\n", param);
	initVolGeom(vg);
	return;
      }
    }
    if (fgets(line, sizeof(line), fp))  
    {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->z_r, &vg->z_a, &vg->z_s);
      if (strcmp(param, "zras")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'zras\n", param);
	initVolGeom(vg);
	return;
      }
    }
    if (fgets(line, sizeof(line), fp))  
    {
      sscanf(line, "%s %s %f %f %f\n", param, eq, &vg->c_r, &vg->c_a, &vg->c_s);
      if (strcmp(param, "cras")!=0)
      {
	fprintf(stderr, "INFO: parameter was %s. should be 'cras'\n", param);
	initVolGeom(vg);
	return;
      }
    }
  }
  else
    fprintf(stderr, "INFO: no volume info stored\n");
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
LINEAR_TRANSFORM_ARRAY *
LTAalloc(int nxforms, MRI *mri)
{
  int                    i ;
  LINEAR_TRANSFORM_ARRAY *lta ;
  MRI_REGION             bbox ;
  float                  x0, y0, z0 ;

  if (mri)
  {
    MRIboundingBox(mri, 70, &bbox) ;
    x0 = bbox.x + bbox.dx/2 ; 
    y0 = bbox.y + bbox.dy/2 ; 
    z0 = bbox.z + bbox.dz/2 ; 
  }
  else
    x0 = y0 = z0 = 0 ;

  lta = (LINEAR_TRANSFORM_ARRAY *)calloc(1, sizeof(LTA)) ;
  if (!lta)
    ErrorExit(ERROR_NOMEMORY, "LTAalloc(%d): could not allocate LTA",nxforms);
  lta->num_xforms = nxforms ;
  lta->xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT)) ;
  if (!lta->xforms)
    ErrorExit(ERROR_NOMEMORY, "LTAalloc(%d): could not allocate xforms",
              nxforms);
  lta->inv_xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT)) ;
  if (!lta->inv_xforms)
    ErrorExit(ERROR_NOMEMORY, 
              "LTAalloc(%d): could not allocate inverse xforms",
              nxforms);
  for (i = 0 ; i < nxforms ; i++)
  {
    lta->xforms[i].x0 = x0 ; 
    lta->xforms[i].y0 = y0 ; 
    lta->xforms[i].z0 = z0 ;
    lta->xforms[i].sigma = 10000.0f ;
    lta->xforms[i].m_L = MatrixIdentity(4, NULL) ;
    lta->xforms[i].m_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
    lta->xforms[i].m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
    initVolGeom(&lta->xforms[i].src);
    initVolGeom(&lta->xforms[i].dst);
    lta->xforms[i].type = UNKNOWN;;

    lta->inv_xforms[i].x0 = x0 ; 
    lta->inv_xforms[i].y0 = y0 ; 
    lta->inv_xforms[i].z0 = z0 ;
    lta->inv_xforms[i].sigma = 10000.0f ;
    lta->inv_xforms[i].m_L = MatrixIdentity(4, NULL) ;
    lta->inv_xforms[i].m_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
    lta->inv_xforms[i].m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
    initVolGeom(&lta->inv_xforms[i].src);
    initVolGeom(&lta->inv_xforms[i].dst);
    lta->inv_xforms[i].type = UNKNOWN;;
  }
  return(lta) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LTAwrite(LTA *lta, char *fname)
{
  FILE             *fp;
  time_t           tt ;
  char             *user, *time_str ;
  LINEAR_TRANSFORM *lt ;
  int              i ;
  char             ext[STRLEN] ;

  if (!stricmp(FileNameExtension(fname, ext), "XFM"))
    return(ltaMNIwrite(lta, fname)) ;

  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "LTAwrite(%s): can't create file",fname));
  user = getenv("USER") ;
  if (!user)
    user = getenv("LOGNAME") ;
  if (!user)
    user = "UNKNOWN" ;
  tt = time(&tt) ;
  time_str = ctime(&tt) ;
  fprintf(fp, "# transform file %s\n# created by %s on %s\n",
          fname, user, time_str) ;
  fprintf(fp, "type      = %d\n", lta->type) ;
  fprintf(fp, "nxforms   = %d\n", lta->num_xforms) ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fprintf(fp, "mean      = %2.3f %2.3f %2.3f\n", lt->x0, lt->y0, lt->z0) ;
    fprintf(fp, "sigma     = %2.3f\n", lt->sigma) ;
    MatrixAsciiWriteInto(fp, lt->m_L) ;
  }
  fclose(fp) ;
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/

LTA *
LTAread(char *fname)
{
  int       type ;
  LTA       *lta ;
  MATRIX    *V, *W, *m_tmp ;

  type = TransformFileNameType(fname);
  switch (type)
  {
  case REGISTER_DAT:
    lta = ltaReadRegisterDat(fname) ;
    if (!lta) return(NULL) ;

    V = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* world to voxel transform */
    W = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* voxel to world transform */
    *MATRIX_RELT(V, 1, 1) = -1 ; *MATRIX_RELT(V, 1, 4) = 128 ;
    *MATRIX_RELT(V, 2, 3) = -1 ; *MATRIX_RELT(V, 2, 4) = 128 ;
    *MATRIX_RELT(V, 3, 2) = 1 ;  *MATRIX_RELT(V, 3, 4) = 128 ;
    *MATRIX_RELT(V, 4, 4) = 1 ;
    
    *MATRIX_RELT(W, 1, 1) = -1 ; *MATRIX_RELT(W, 1, 4) = 128 ;
    *MATRIX_RELT(W, 2, 3) = 1 ; *MATRIX_RELT(W, 2, 4) = -128 ;
    *MATRIX_RELT(W, 3, 2) = -1 ;  *MATRIX_RELT(W, 3, 4) = 128 ;
    *MATRIX_RELT(W, 4, 4) = 1 ;
    
    m_tmp = MatrixMultiply(lta->xforms[0].m_L, W, NULL) ;
    MatrixMultiply(V, m_tmp, lta->xforms[0].m_L) ;
    MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
    lta->type = LINEAR_VOX_TO_VOX ;
		break ;
  case MNI_TRANSFORM_TYPE:
    lta = ltaMNIread(fname) ;
    if (!lta)
      return(NULL) ;
    
    /* by default convert MNI files to voxel coords.
       Sorry, I know this shouldn't be done here, particularly since we
       don't know enough to convert to scanner RAS coords, but I don't want 
       to risk breaking the Talairach code by mucking around with it (BRF).
    */
    /* convert to voxel coords */
    V = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* world to voxel transform */
    W = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* voxel to world transform */
    *MATRIX_RELT(V, 1, 1) = -1 ; *MATRIX_RELT(V, 1, 4) = 128 ;
    *MATRIX_RELT(V, 2, 3) = -1 ; *MATRIX_RELT(V, 2, 4) = 128 ;
    *MATRIX_RELT(V, 3, 2) = 1 ;  *MATRIX_RELT(V, 3, 4) = 128 ;
    *MATRIX_RELT(V, 4, 4) = 1 ;
    
    *MATRIX_RELT(W, 1, 1) = -1 ; *MATRIX_RELT(W, 1, 4) = 128 ;
    *MATRIX_RELT(W, 2, 3) = 1 ; *MATRIX_RELT(W, 2, 4) = -128 ;
    *MATRIX_RELT(W, 3, 2) = -1 ;  *MATRIX_RELT(W, 3, 4) = 128 ;
    *MATRIX_RELT(W, 4, 4) = 1 ;
    
    m_tmp = MatrixMultiply(lta->xforms[0].m_L, W, NULL) ;
    MatrixMultiply(V, m_tmp, lta->xforms[0].m_L) ;
    MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
    lta->type = LINEAR_VOX_TO_VOX ;
    break ;
  case LINEAR_VOX_TO_VOX:
  case LINEAR_RAS_TO_RAS:
  case TRANSFORM_ARRAY_TYPE:
  default:
    lta = ltaReadFile(fname) ;
    break ;
  }
  return(lta) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static LTA *
ltaReadFile(char *fname)
{
  FILE             *fp;
  LINEAR_TRANSFORM *lt ;
  int              i, nxforms, type ;
  char             line[STRLEN], *cp ;
  LTA              *lta ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(NULL,
                (ERROR_BADFILE, "ltaReadFile(%s): can't open file",fname));
  cp = fgetl(line, 199, fp) ; 
  if (cp == NULL)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data",fname));
  }
  sscanf(cp, "type      = %d\n", &type) ;
  cp = fgetl(line, 199, fp) ; sscanf(cp, "nxforms   = %d\n", &nxforms) ;
  lta = LTAalloc(nxforms, NULL) ;
  lta->type = type ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fscanf(fp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0) ;
    fscanf(fp, "sigma     = %f\n", &lt->sigma) ;
    MatrixAsciiReadFrom(fp, lt->m_L) ;
  }
  fclose(fp) ;
  return(lta) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LTAfree(LTA **plta)
{
  LTA *lta ;
  int i ;

  lta = *plta ;
  *plta = NULL ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    MatrixFree(&lta->xforms[i].m_L) ;
    MatrixFree(&lta->xforms[i].m_dL) ;
    MatrixFree(&lta->xforms[i].m_last_dL) ;

    MatrixFree(&lta->inv_xforms[i].m_L) ;
    MatrixFree(&lta->inv_xforms[i].m_dL) ;
    MatrixFree(&lta->inv_xforms[i].m_last_dL) ;
  }
  free(lta->xforms) ;
  free(lta) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
            split each transform into 8 and space them
            evenly.
------------------------------------------------------*/
int
LTAdivide(LTA *lta, MRI *mri)
{
  MRI_REGION  bbox ;
  int         oldi, i, nxforms, row_size, x, y, z ;
  LT          *new_xforms, *lt ;
  float       sigma, dx, dy, dz, len,x0, y0, z0 ;

  MRIboundingBox(mri, 130, &bbox) ;
  dx = bbox.dx ; dy = bbox.dy ; dz = bbox.dz ;
  nxforms = lta->num_xforms*8 ;
  len = pow((double)nxforms, 1.0/3.0) ;  /* # along each dimension */
  row_size = nint(len) ;  /* # of nodes in each dimension */
  dx = dx / (len+1) ;  /* x spacing between nodes */
  dy = dy / (len+1) ;  /* y spacing between nodes */
  dz = dz / (len+1) ;  /* z spacing between nodes */
  new_xforms = (LINEAR_TRANSFORM *)calloc(nxforms, sizeof(LT)) ;

  sigma = 1.0 * (dx + dy + dz) / 3.0f ;   /* average node spacing */

  /*
    distribute the new transforms uniformly throughout the volume.
    Next, for each new transform, calculate the old transform at
    that point, and use it as the initial value for the new transform.
  */
  for (i = x = 0 ; x < row_size ; x++)
  {
    x0 = dx/2 + x * dx ;
    for (y = 0 ; y < row_size ; y++)
    {
      y0 = dy/2 + y * dy ;
      for (z = 0 ; z < row_size ; z++, i++)
      {
        z0 = dz/2 + z * dz ;
        lt = &new_xforms[i] ;
        lt->x0 = x0 ; lt->y0 = y0 ; lt->z0 = z0 ; lt->sigma = sigma ;
        lt->m_L = LTAtransformAtPoint(lta, x0, y0, z0, NULL) ;
        lt->m_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
        lt->m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
      }
    }
  }

  /* free old ones */
  for (oldi = 0 ; oldi < lta->num_xforms ; oldi++)
  {
    MatrixFree(&lta->xforms[oldi].m_L) ;
    MatrixFree(&lta->xforms[oldi].m_dL) ;
    MatrixFree(&lta->xforms[oldi].m_last_dL) ;
  }
  free(lta->xforms) ;

  /* update lta structure with new info */
  lta->xforms = new_xforms ; lta->num_xforms = nxforms ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Apply a transform array to an image.
------------------------------------------------------*/
MRI *
LTAtransform(MRI *mri_src, MRI *mri_dst, LTA *lta)
{
  int         y1, y2, y3, width, height, depth, xi, yi, zi ;
  VECTOR      *v_X, *v_Y ;/* original and transformed coordinate systems */
  Real        x1, x2, x3 ;
  MATRIX      *m_L, *m_L_inv ;
  LT *tran = &lta->xforms[0];

  if (lta->num_xforms == 1)
  {
    if (lta->type == LINEAR_RAS_TO_RAS)
    {
      if (!mri_dst) 
	mri_dst = MRIclone(mri_src, NULL) ;
      if (tran->dst.valid == 1) // transform dst is valid
      {
	// modify dst c_(r,a,s) using the transform dst value
	// to make the better positioning
	printf("INFO: Modifying dst c_(r,a,s), using the transform dst\n");
	mri_dst->c_r = tran->dst.c_r;
	mri_dst->c_a = tran->dst.c_a;
	mri_dst->c_s = tran->dst.c_s;
      }
      else
	printf("INFO: Transform dst volume info is not used (valid flag = 0).\n");

      return(MRIapplyRASlinearTransform(mri_src, mri_dst, lta->xforms[0].m_L)) ;
    }
    else
      return(MRIlinearTransform(mri_src, mri_dst, lta->xforms[0].m_L)) ;
  }

  fprintf(stderr, "applying octree transform to image...\n") ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  v_X     = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y     = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */
  m_L     = MatrixAlloc(4, 4, MATRIX_REAL) ;
  m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL) ;
  v_Y->rptr[4][1] = 1.0f ;

  if (lta->num_xforms == 1)
  {
    LTAtransformAtPoint(lta, 0, 0, 0, m_L) ;
    if (MatrixInverse(m_L, m_L_inv) == NULL)
    {
      MatrixFree(&m_L) ; MatrixFree(&m_L_inv) ; 
      ErrorReturn(NULL,
                  (ERROR_BADPARM, "LTAtransform: could not invert matrix")) ;
    }
  }
  for (y3 = 0 ; y3 < depth ; y3++)
  {
    DiagHeartbeat((float)y3 / (float)(depth-1)) ;
    V3_Z(v_Y) = y3 ;

    for (y2 = 0 ; y2 < height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < width ; y1++)
      {
        V3_X(v_Y) = y1 ;

        /*
          this is not quite right. Should use the weighting at
          the point that maps to (y1,y2,y3), not (y1,y2,y3) itself,
          but this is much easier....
        */
        if (lta->num_xforms > 1)
        {
          LTAtransformAtPoint(lta, y1, y2, y3, m_L) ;
#if 0
          if (MatrixSVDInverse(m_L, m_L_inv) == NULL)
            continue ;
#else
          if (MatrixInverse(m_L, m_L_inv) == NULL)
            continue ;
#endif
        }
        MatrixMultiply(m_L_inv, v_Y, v_X) ;
        x1 = V3_X(v_X) ;
        x2 = V3_Y(v_X) ;
        x3 = V3_Z(v_X) ;

        xi = mri_dst->xi[nint(x1)] ; 
        yi = mri_dst->yi[nint(x2)] ; 
        zi = mri_dst->zi[nint(x3)] ; 
        MRIvox(mri_dst, y1, y2, y3) = MRIvox(mri_src, xi, yi, zi) ;
      }
    }
#if 0
    if (y3 > 10)
      exit(0) ;
#endif
  }

  MatrixFree(&v_X) ; MatrixFree(&v_Y) ;
  MatrixFree(&m_L) ; MatrixFree(&m_L_inv) ;

  return(mri_dst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MATRIX *
LTAtransformAtPoint(LTA *lta, float x, float y, float z, MATRIX *m_L)
{
  LT       *lt ;
  int      i ;
  double   w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin ;
  static MATRIX  *m_tmp = NULL ;


  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL) ;
  else
    MatrixClear(m_L) ;

  if (m_tmp == NULL)
    m_tmp = MatrixAlloc(4, 4, MATRIX_REAL) ;

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    dx = lta->xforms[i].x0 - x ;
    dy = lta->xforms[i].y0 - y ;
    dz = lta->xforms[i].z0 - z ;
    dsq = (dx*dx+dy*dy+dz*dz) ;
    sigma = lt->sigma ;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq/(2*sigma*sigma));
    wtotal += w_p[i] ;
  }

  if (DZERO(wtotal))   /* no transforms in range??? */
#if 0
    MatrixIdentity(4, m_L) ;
#else
  MatrixCopy(lta->xforms[0].m_L, m_L) ;
#endif
  else /* now calculate linear combination of transforms at this point */
  {
    wmin = 0.1 / (double)lta->num_xforms ;
    for (i = 0 ; i < lta->num_xforms ; i++)
    {
      lt = &lta->xforms[i] ;
      w_k_p = w_p[i] / wtotal ;
      if (w_k_p < wmin)   /* optimization - ignore this transform */
        continue ;
      MatrixScalarMul(lt->m_L, w_k_p, m_tmp) ;
      MatrixAdd(m_L, m_tmp, m_L) ;
    }
  }
  return(m_L) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MATRIX *
LTAinverseTransformAtPoint(LTA *lta, float x, float y, float z, MATRIX *m_L)
{
  LT       *lt ;
  int      i ;
  double   w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin ;
  static MATRIX  *m_tmp = NULL ;
  MATRIX   *m_inv ;


  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL) ;
  else
    MatrixClear(m_L) ;

  if (m_tmp == NULL)
    m_tmp = MatrixAlloc(4, 4, MATRIX_REAL) ;

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    dx = lta->xforms[i].x0 - x ;
    dy = lta->xforms[i].y0 - y ;
    dz = lta->xforms[i].z0 - z ;
    dsq = (dx*dx+dy*dy+dz*dz) ;
    sigma = lt->sigma ;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq/(2*sigma*sigma));
    wtotal += w_p[i] ;
  }


  if (lta->num_xforms == 1)
  {
    m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  }
  else
  {
    if (DZERO(wtotal))   /* no transforms in range??? */
      MatrixIdentity(4, m_L) ;
    else /* now calculate linear combination of transforms at this point */
    {
      wmin = 0.1 / (double)lta->num_xforms ;
      for (i = 0 ; i < lta->num_xforms ; i++)
      {
        lt = &lta->xforms[i] ;
        w_k_p = w_p[i] / wtotal ;
        if (w_k_p < wmin)   /* optimization - ignore this transform */
          continue ;
        
        m_inv = MatrixInverse(lt->m_L, NULL) ;
        if (!m_inv)
          continue ;
        MatrixScalarMul(m_inv, w_k_p, m_tmp) ;
        MatrixAdd(m_L, m_tmp, m_L) ;
        MatrixFree(&m_inv) ;
      }
    }
  }
  return(m_L) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
VECTOR *
LTAtransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  static MATRIX *m_L = NULL ;
  
  m_L = LTAtransformAtPoint(lta, V3_X(v_X), V3_Y(v_X), V3_Z(v_X), m_L) ;
  v_Y = MatrixMultiply(m_L, v_X, v_Y) ;
  return(v_Y) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
VECTOR *
LTAinverseTransformPoint(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  static MATRIX *m_L = NULL ;
  
  m_L = LTAinverseTransformAtPoint(lta, V3_X(v_X), V3_Y(v_X), V3_Z(v_X), m_L);
  v_Y = MatrixMultiply(m_L, v_X, v_Y) ;
  return(v_Y) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
         Same as LTAtransformPoint but also return weight normalization
         factor.
------------------------------------------------------*/
double
LTAtransformPointAndGetWtotal(LTA *lta, VECTOR *v_X, VECTOR *v_Y)
{
  LT       *lt ;
  int      i ;
  double   w_p[MAX_TRANSFORMS], wtotal, dsq, sigma, dx, dy, dz, w_k_p, wmin,
           x, y, z ;
  static MATRIX  *m_tmp = NULL, *m_L = NULL ;

  x = V3_X(v_X) ; y = V3_Y(v_X) ; z = V3_Z(v_X) ;
  if (m_L == NULL)
    m_L = MatrixAlloc(4, 4, MATRIX_REAL) ;
  else
    MatrixClear(m_L) ;

  if (m_tmp == NULL)
    m_tmp = MatrixAlloc(4, 4, MATRIX_REAL) ;

  /* first compute normalized weights */
  for (wtotal = 0.0, i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    dx = lta->xforms[i].x0 - x ;
    dy = lta->xforms[i].y0 - y ;
    dz = lta->xforms[i].z0 - z ;
    dsq = (dx*dx+dy*dy+dz*dz) ;
    sigma = lt->sigma ;
    w_p[i] = /*(1 / (sigma*sqrt(2.0*M_PI))) * */ exp(-dsq/(2*sigma*sigma));
    wtotal += w_p[i] ;
  }

  if (DZERO(wtotal))   /* no transforms in range??? */
    MatrixIdentity(4, m_L) ;
  else /* now calculate linear combination of transforms at this point */
  {
    wmin = 0.1 / (double)lta->num_xforms ;
    for (i = 0 ; i < lta->num_xforms ; i++)
    {
      lt = &lta->xforms[i] ;
      w_k_p = w_p[i] / wtotal ;
      if (w_k_p < wmin)   /* optimization - ignore this transform */
        continue ;
      MatrixScalarMul(lt->m_L, w_k_p, m_tmp) ;
      MatrixAdd(m_L, m_tmp, m_L) ;
    }
  }

  MatrixMultiply(m_L, v_X, v_Y) ;
  return(wtotal) ;
}
/*-----------------------------------------------------*/
int TransformFileNameType(char *fname)
{
  int file_type = TRANSFORM_ARRAY_TYPE ;
  char *dot, buf[500], *number ;

  strcpy(buf, fname) ;
  dot = strrchr(buf, '@') ;
  number = strchr(buf, '#') ;
  if(number)  *number = 0 ;  /* don't consider : part of extension */

  if(!dot)   dot = strrchr(buf, '.') ;

  if(dot)
  {
    dot++ ;
    StrUpper(buf) ;
    if (!strcmp(dot, "M3D"))
      return(MORPH_3D_TYPE) ;
    else if (!strcmp(dot, "LTA"))
      return(TRANSFORM_ARRAY_TYPE) ;
    else if (!strcmp(dot, "DAT"))
      return(REGISTER_DAT) ;
    else if (!strcmp(dot, "OCT"))
      return(TRANSFORM_ARRAY_TYPE) ;
    else if (!strcmp(dot, "XFM"))
      return(MNI_TRANSFORM_TYPE) ;
  }

  return(file_type) ;
}

#include "volume_io.h"

static int  
ltaMNIwrite(LTA *lta, char *fname)
{
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE, 
                (ERROR_NOFILE, "ltMNIwrite: could not open file %s",fname));

  fprintf(fp, "MNI Transform File\n") ;
  fprintf(fp, "%%Generated by %s\n", Progname) ;
  fprintf(fp, "\n") ;
  fprintf(fp, "Transform_Type = Linear;\n") ;
  fprintf(fp, "Linear_Transform =\n") ;

  m_L = lta->xforms[0].m_L ;
  for (row = 1 ; row <= 3 ; row++)
  {
    fprintf(fp, "      %f       %f       %f       %f",
           *MATRIX_RELT(m_L,row,1), *MATRIX_RELT(m_L,row,2), 
           *MATRIX_RELT(m_L,row,3), *MATRIX_RELT(m_L,row,4)) ;
    if (row == 3)
      fprintf(fp, ";") ;
    fprintf(fp, "\n") ;
  }
#if 0
  m_tmp = MatrixMultiply(lt->m_L, W, NULL) ;
  MatrixMultiply(V, m_tmp, lt->m_L) ;
  /*  MatrixAsciiWriteInto(stderr, lt->m_L) ;*/
  MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
#endif
  fclose(fp) ;
  return(NO_ERROR) ;
}
static LTA  *
ltaMNIread(char *fname)
{
  LTA              *lta ;
  LINEAR_TRANSFORM *lt ;
  char             *cp, line[1000] ;
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NOFILE, "ltMNIread: could not open file %s",fname));

  lta = LTAalloc(1, NULL) ;
  lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;

  fgetl(line, 900, fp) ;   /* MNI Transform File */
  fgetl(line, 900, fp) ;   /* Transform_type = Linear; */
  while (line[0] == '%')
    fgetl(line, 900, fp) ; /* variable # of comments */
  fgetl(line, 900, fp) ;

  m_L = lt->m_L ;
  for (row = 1 ; row <= 3 ; row++)
  {
    cp = fgetl(line, 900, fp) ;
    if (!cp)
    {
      LTAfree(&lta) ;
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "ltMNIread: could not read row %d from %s",
                   row, fname)) ;
    }
    sscanf(cp, "%f %f %f %f",
           MATRIX_RELT(m_L,row,1), MATRIX_RELT(m_L,row,2), 
           MATRIX_RELT(m_L,row,3), MATRIX_RELT(m_L,row,4)) ;
  }
  fclose(fp) ;
  return(lta) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LTAworldToWorld(LTA *lta, float x, float y, float z, float *px, float *py, 
                float *pz)
{
  static VECTOR *v_X, *v_Y = NULL ;

  if (v_Y == NULL)
  {
    v_X = VectorAlloc(4, MATRIX_REAL) ;
    v_Y = VectorAlloc(4, MATRIX_REAL) ;
  }
  /* world to voxel */
  v_X->rptr[4][1] = 1.0f ;
  V3_X(v_X) = 128.0 - x ;
  V3_Z(v_X) = (y + 128.0) ;
  V3_Y(v_X) = (-z + 128.0) ;

  LTAtransformPoint(lta, v_X, v_Y) ;

  /* voxel to world */
  *px = 128.0 - V3_X(v_Y) ;
  *py = V3_Z(v_Y)  - 128.0 ;
  *pz = -(V3_Y(v_Y) - 128.0) ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
LTAinverseWorldToWorld(LTA *lta, float x, float y, float z, float *px, 
                       float *py, float *pz)
{
  static VECTOR *v_X, *v_Y = NULL ;

  if (v_Y == NULL)
  {
    v_X = VectorAlloc(4, MATRIX_REAL) ;
    v_Y = VectorAlloc(4, MATRIX_REAL) ;
  }
  /* world to voxel */
  v_X->rptr[4][1] = 1.0f ;
  V3_X(v_X) = 128.0 - x ;
  V3_Z(v_X) = (y + 128.0) ;
  V3_Y(v_X) = (-z + 128.0) ;

  LTAinverseTransformPoint(lta, v_X, v_Y) ;

  /* voxel to world */
  *px = 128.0 - V3_X(v_Y) ;
  *py = V3_Z(v_Y)  - 128.0 ;
  *pz = -(V3_Y(v_Y) - 128.0) ;

  return(NO_ERROR) ;
}

int
LTAtoVoxelCoords(LTA *lta, MRI *mri)
{
  MATRIX *m_L ;
  int    i ;

  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    m_L = MRIrasXformToVoxelXform(mri, mri, lta->xforms[i].m_L, NULL) ;
    MatrixFree(&lta->xforms[0].m_L) ;
    lta->xforms[0].m_L = m_L ;
  }
  lta->type = LINEAR_VOX_TO_VOX ;
  return(NO_ERROR) ;
}


int
LTAvoxelToRasXform(LTA *lta, MRI *mri_src, MRI *mri_dst)
{

  MATRIX *m_L ;
  int    i ;

  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    m_L = MRIvoxelXformToRasXform(mri_src, mri_dst, lta->xforms[i].m_L, NULL) ;
    MatrixFree(&lta->xforms[0].m_L) ;
    lta->xforms[0].m_L = m_L ;
  }  

  lta->type = LINEAR_RAS_TO_RAS ;
  return(NO_ERROR) ;
}

int
LTArasToVoxelXform(LTA *lta, MRI *mri_src, MRI *mri_dst)
{

  MATRIX *m_L ;
  int    i ;

  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    m_L = MRIrasXformToVoxelXform(mri_src, mri_dst, lta->xforms[i].m_L, NULL) ;
    MatrixFree(&lta->xforms[0].m_L) ;
    lta->xforms[0].m_L = m_L ;
  }  

  lta->type = LINEAR_VOX_TO_VOX ;
  return(NO_ERROR) ;
}
int
LTAvoxelTransformToCoronalRasTransform(LTA *lta)
{
  MATRIX   *V, *W, *m_tmp ;

  V = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* world to voxel transform */
  W = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* voxel to world transform */
  *MATRIX_RELT(V, 1, 1) = -1 ; *MATRIX_RELT(V, 1, 4) = 128 ;
  *MATRIX_RELT(V, 2, 3) = -1 ; *MATRIX_RELT(V, 2, 4) = 128 ;
  *MATRIX_RELT(V, 3, 2) = 1 ;  *MATRIX_RELT(V, 3, 4) = 128 ;
  *MATRIX_RELT(V, 4, 4) = 1 ;
  
  *MATRIX_RELT(W, 1, 1) = -1 ; *MATRIX_RELT(W, 1, 4) = 128 ;
  *MATRIX_RELT(W, 2, 3) = 1 ; *MATRIX_RELT(W, 2, 4) = -128 ;
  *MATRIX_RELT(W, 3, 2) = -1 ;  *MATRIX_RELT(W, 3, 4) = 128 ;
  *MATRIX_RELT(W, 4, 4) = 1 ;

  m_tmp = MatrixMultiply(lta->xforms[0].m_L, V, NULL) ;
  MatrixMultiply(W, m_tmp, lta->xforms[0].m_L) ;
  MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
  lta->type = LINEAR_RAS_TO_RAS ;

  return(NO_ERROR) ;
}
/*-----------------------------------------------------------
  FixMNITal() - function to compute the "real" talairach coordinates
  from the MNI talaiarch coordinates. Can be done in-place. This is
  based on Matthew Brett's 10/8/98 transform. See
  http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
  -----------------------------------------------------------*/
int FixMNITal(float  xmni, float  ymni, float  zmni,
        float *xtal, float *ytal, float *ztal)
{
  MATRIX *T, *xyzMNI, *xyzTal;


  T = MatrixAlloc(4, 4, MATRIX_REAL);
  if(zmni >= 0.0){
    stuff_four_by_four(T, 
           .9900,  .0000, .0000, 0,
           .0000,  .9688, .0460, 0,
           .0000, -.0485, .9189, 0,
           .0000,  .0000, .0000, 1);
  }
  else {
    stuff_four_by_four(T, 
           .9900,  .0000, .0000, 0,
           .0000,  .9688, .0420, 0,
           .0000, -.0485, .8390, 0,
           .0000,  .0000, .0000, 1);
  }

  xyzMNI = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzMNI->rptr[1][1] = xmni;
  xyzMNI->rptr[2][1] = ymni;
  xyzMNI->rptr[3][1] = zmni;
  xyzMNI->rptr[4][1] = 1.0;

  xyzTal = MatrixMultiply(T,xyzMNI,NULL);

  *xtal = xyzTal->rptr[1][1];
  *ytal = xyzTal->rptr[2][1];
  *ztal = xyzTal->rptr[3][1];

  MatrixFree(&T);
  MatrixFree(&xyzMNI);
  MatrixFree(&xyzTal);

  return(0);
}
/*-----------------------------------------------------------------
  DevolveXFM() - change a tranformation matrix to work with tkreg RAS
  coordinates (ie, c_ras = 0). Assumes that the transformation matrix
  was computed with the orig-volume native vox2ras transform (or any
  volume with the same geometry). The XFM maps from native orig RAS to
  a native target RAS. Note: this has no effect when the orig volume
  has c_ras = 0. If XFM is empty, then xfmname is read from the
  subject's transforms directory.  If xfmname is empty, then
  talairach.xfm is used.  
  Note: uses LTAvoxelTransformToCoronalRasTransform(). 
  -----------------------------------------------------------------*/
MATRIX *DevolveXFM(char *subjid, MATRIX *XFM, char *xfmname)
{
  MATRIX *Torig_tkreg, *invTorig_tkreg, *Torig_native, *Mfix;
  char dirname[2000], xfmpath[2000], *sd;
  MRI *mriorig;
  FILE *fp;
  LTA    *lta;

  sd = getenv("SUBJECTS_DIR") ;
  if(sd==NULL){
    printf("ERROR: SUBJECTS_DIR not defined\n");
    return(NULL);
  }

  /* Check that the subject exists */
  sprintf(dirname,"%s/%s",sd,subjid);
  if(!fio_IsDirectory(dirname)){
    printf("ERROR: cannot find subject %s in %s\n",subjid,sd);
    return(NULL);
  }

  /* Load the orig header for the subject */
  sprintf(dirname,"%s/%s/mri/orig",sd,subjid);
  mriorig = MRIreadHeader(dirname,MRI_CORONAL_SLICE_DIRECTORY);
  if(mriorig == NULL){
    printf("ERROR: could not read header for %s\n",dirname);
    return(NULL);
  }

  if(XFM==NULL){
    /* Read in the talairach.xfm matrix */
    if(xfmname == NULL) xfmname = "talairach.xfm";
    sprintf(xfmpath,"%s/%s/mri/transforms/%s",sd,subjid,xfmname);
    fp = fopen(xfmpath,"r");
    if(fp == NULL){
      printf("ERROR: could not open %s for reading \n",xfmpath);
      return(NULL);
    }
    lta = LTAread(xfmpath);
    if(lta == NULL){
      printf("ERROR: reading %s\n",xfmpath);
      return(NULL);
    }
    if(lta->type == LINEAR_VOX_TO_VOX)
      LTAvoxelTransformToCoronalRasTransform(lta);
    XFM = MatrixCopy(lta->xforms[0].m_L,NULL);
    LTAfree(&lta);
  }

  /* Mfix = inv(Torig_tkreg)*Torig_native */
  /* X2 = Mfix*X */
  Torig_tkreg  = MRIxfmCRS2XYZtkreg(mriorig);
  Torig_native = MRIxfmCRS2XYZ(mriorig,0);
  invTorig_tkreg = MatrixInverse(Torig_tkreg,NULL);
  Mfix = MatrixMultiply(Torig_native,invTorig_tkreg,NULL);
  MatrixMultiply(XFM,Mfix,XFM);

  MatrixFree(&Mfix);
  MatrixFree(&Torig_tkreg);
  MatrixFree(&invTorig_tkreg);
  MatrixFree(&Torig_native);
  MRIfree(&mriorig);

  return(XFM);
}
/*----------------------------------------------------------------*/
TRANSFORM *
TransformRead(char *fname)
{
  TRANSFORM *trans ;
  GCA_MORPH *gcam ;

  trans = (TRANSFORM *)calloc(1, sizeof(TRANSFORM)) ;

  trans->type = TransformFileNameType(fname) ;
  switch (trans->type)
  {
  case MNI_TRANSFORM_TYPE:
  case LINEAR_VOX_TO_VOX:
  case LINEAR_RAS_TO_RAS:
  case TRANSFORM_ARRAY_TYPE:
  default:
    trans->xform = (void *)LTAread(fname) ;
    if (!trans->xform)
    {
      free(trans) ;
      return(NULL) ;
    }
    trans->type = ((LTA *)trans->xform)->type ;
    break ;
  case MORPH_3D_TYPE:
    gcam = GCAMread(fname) ;
    if (!gcam)
    {
      free(trans) ;
      return(NULL) ;
    }
    trans->xform = (void *)gcam ;
    break ;
  }
  return(trans) ;
}


int
TransformFree(TRANSFORM **ptrans)
{
  TRANSFORM *trans ;
  
  trans = *ptrans ; *ptrans = NULL ;

  switch (trans->type)
  {
  default:
    LTAfree((LTA **)(&trans->xform)) ;
    break ;
  case MORPH_3D_TYPE:
    GCAMfree((GCA_MORPH **)(&trans->xform)) ;
    break ;
  }
  free(trans) ;
    
  return(NO_ERROR) ;
}

/*
  take a voxel in MRI space and find the voxel in the gca/gcamorph space to which
  it maps. Note that the caller will scale this down depending on the spacing of
  the gca/gcamorph.
*/
int
TransformSample(TRANSFORM *transform, int xv, int yv, int zv, float *px, float *py, float *pz)
{
  static VECTOR   *v_input, *v_canon = NULL ;
  float           xt, yt, zt ;
  LTA             *lta ;
  GCA_MORPH       *gcam ;

  *px = *py = *pz = 0 ;
  if (transform->type == MORPH_3D_TYPE)
  {
    gcam = (GCA_MORPH *)transform->xform ;
    if (!gcam->mri_xind)
      ErrorReturn(ERROR_UNSUPPORTED, 
                  (ERROR_UNSUPPORTED, "TransformSample: gcam has not been inverted!")) ;
    if (xv < 0)
      xv = 0 ;
    if (xv >= gcam->mri_xind->width)
      xv = gcam->mri_xind->width-1 ;
    if (yv < 0)
      yv = 0 ;
    if (yv >= gcam->mri_yind->height)
      yv = gcam->mri_yind->height-1 ;
    if (zv < 0)
      zv = 0 ;
    if (zv >= gcam->mri_zind->depth)
      zv = gcam->mri_zind->depth-1 ;
    xt = (int)MRISvox(gcam->mri_xind, xv, yv, zv)*gcam->spacing ;
    yt = (int)MRISvox(gcam->mri_yind, xv, yv, zv)*gcam->spacing ;
    zt = (int)MRISvox(gcam->mri_zind, xv, yv, zv)*gcam->spacing ;
  }
  else
  {
    lta = (LTA *)transform->xform ;
    if (!v_canon)
    {
      v_input = VectorAlloc(4, MATRIX_REAL) ;
      v_canon = VectorAlloc(4, MATRIX_REAL) ;
      *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
    }
    
    V3_X(v_input) = (float)xv;
    V3_Y(v_input) = (float)yv;
    V3_Z(v_input) = (float)zv;
    MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon) ;
    xt = V3_X(v_canon) ; yt = V3_Y(v_canon) ; zt = V3_Z(v_canon) ; 
  }
  *px = xt ; *py = yt ; *pz = zt ;
  return(NO_ERROR) ;
}
/*
  take a voxel in gca/gcamorph space and find the MRI voxel to which
  it maps
*/
int
TransformSampleInverse(TRANSFORM *transform, int xv, int yv, int zv, float *px, float *py, float *pz)
{
  static VECTOR  *v_input, *v_canon = NULL ;
  static MATRIX  *m_L_inv ;
  float          xt, yt, zt ;
  int            xn, yn, zn ;
  LTA            *lta ;
  GCA_MORPH      *gcam ;
  GCA_MORPH_NODE *gcamn ;

  if (transform->type == MORPH_3D_TYPE)
  {
    gcam = (GCA_MORPH *)transform->xform ;
    xn = nint(xv/gcam->spacing) ; yn = nint(yv/gcam->spacing) ; zn = nint(zv/gcam->spacing) ;
    if (xn >= gcam->width)
      xn = gcam->width-1 ;
    if (yn >= gcam->height)
      yn = gcam->height-1 ;
    if (zn >= gcam->depth)
      yn = gcam->depth-1 ;
    gcamn = &gcam->nodes[xn][yn][zn] ;
    xt = gcamn->x ; yt = gcamn->y ; zt = gcamn->z ; 
  }
  else
  {
    lta = (LTA *)transform->xform ;
    if (!v_canon)
    {
      v_input = VectorAlloc(4, MATRIX_REAL) ;
      v_canon = VectorAlloc(4, MATRIX_REAL) ;
      *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
      *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
      m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL) ;
    }
    
    V3_X(v_canon) = (float)xv; V3_Y(v_canon) = (float)yv; V3_Z(v_canon) = (float)zv;
#if 0
    MatrixInverse(lta->xforms[0].m_L, m_L_inv) ;
    MatrixMultiply(m_L_inv, v_canon, v_input) ;
#else
    MatrixMultiply(lta->inv_xforms[0].m_L, v_canon, v_input) ;
#endif
    xt = V3_X(v_input) ; yt = V3_Y(v_input) ; zt = V3_Z(v_input) ; 
  }
  *px = xt ; *py = yt ; *pz = zt ;
  return(NO_ERROR) ;
}
int
TransformSampleInverseVoxel(TRANSFORM *transform, int width, int height, int depth,
                            int xv, int yv, int zv, 
                            int *px, int *py, int *pz)
{
  float   xf, yf, zf ;

  TransformSampleInverse(transform, xv, yv, zv, &xf, &yf, &zf) ;
  xv = nint(xf) ; yv = nint(yf); zv = nint(zf);
  if (xv < 0) xv = 0 ;
  if (xv >= width) xv = width-1 ;
  if (yv < 0) yv = 0 ;
  if (yv >= height) yv = height-1 ;
  if (zv < 0) zv = 0 ;
  if (zv >= depth) zv = depth-1 ;
  *px = xv ; *py = yv ; *pz = zv ;
  return(NO_ERROR) ;
}


TRANSFORM *
TransformAlloc(int type, MRI *mri)
{
  TRANSFORM *transform ;

  
  switch (type)
  {
  default:
    transform = (TRANSFORM *)calloc(1, sizeof(TRANSFORM)) ;
    transform->type = type ;
    transform->xform = (void *)LTAalloc(1, mri) ;
    break ;
  case MORPH_3D_TYPE:
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "TransformAlloc(MORPH_3D_TYPE): unsupported")) ;
    break ;
  }
  return(transform) ;
}

int
TransformInvert(TRANSFORM *transform, MRI *mri)
{
  LTA       *lta ;
  GCA_MORPH *gcam ;

  switch (transform->type)
  {
  default:
    lta = (LTA *)transform->xform ;
    if (MatrixInverse(lta->xforms[0].m_L, lta->inv_xforms[0].m_L) == NULL)
      ErrorExit(ERROR_BADPARM, "TransformInvert: xform noninvertible") ;
    break ;
  case MORPH_3D_TYPE:
    if (!mri)
      return(NO_ERROR) ;
    gcam = (GCA_MORPH *)transform->xform ;
    GCAMinvert(gcam, mri) ;
    break ;
  }
  return(NO_ERROR) ;
}
MRI *
TransformApply(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  switch (transform->type)
  {
  case MORPH_3D_TYPE:
    mri_dst = GCAMmorphToAtlas(mri_src, (GCA_MORPH*)transform->xform, NULL) ;
    break ;
  default:
    mri_dst = MRIlinearTransform(mri_src, NULL, ((LTA *)transform->xform)->xforms[0].m_L);
    break ;
  }
  return(mri_dst) ;
}

MRI *
TransformApplyInverse(TRANSFORM *transform, MRI *mri_src, MRI *mri_dst)
{
  switch (transform->type)
  {
  case MORPH_3D_TYPE:
    mri_dst = GCAMmorphFromAtlas(mri_src,(GCA_MORPH*)transform->xform,NULL);
    break ;
  default:
    mri_dst = MRIinverseLinearTransform(mri_src, NULL, ((LTA *)transform->xform)->xforms[0].m_L);
    break ;
  }
  return(mri_dst) ;
}

#include "stats.h"
static LTA *
ltaReadRegisterDat(char *fname)
{
  LTA        *lta ;
  fMRI_REG   *reg ;

  lta = LTAalloc(1, NULL) ;
  lta->xforms[0].sigma = 1.0f ;
  lta->xforms[0].x0 = lta->xforms[0].y0 = lta->xforms[0].z0 = 0 ;
  reg = StatReadRegistration(fname) ;
  MatrixCopy(reg->mri2fmri, lta->xforms[0].m_L) ;
  StatFreeRegistration(&reg) ;

  lta->type = LINEAR_VOXEL_TO_VOXEL ;
  return(lta) ;
}


// find volumes which created the transform.
// if buffer == NULL, then it will allocate memory
int mincFindVolume(const char *line, char **srcVol, char **dstVol)
{
  struct stat stat_buf;
  int ret;
  static int count = 0;
  char buf[1000];
  char *pch = strtok((char *) line, " ");
  while (pch != NULL)
  {
    strcpy(buf, pch);
    if (strstr(buf, ".mnc")) // first src mnc volume
    {
      if (count ==0) // src
      {
	count++;
	*srcVol = (char *) malloc(strlen(pch)+1);
	strcpy(*srcVol, pch);
	// check the existence of the file
	ret = stat(*srcVol, &stat_buf);
	if (ret != 0)
	{
	  fprintf(stderr, "INFO: Src volume %s cannot be found. No problem.\n",*srcVol);
	  free(*srcVol);
	  *srcVol = 0;
	}
	else
	  fprintf(stdout, "INFO: Src volume %s\n", *srcVol);
      }	
      else if (count == 1) // this is the second one must be dest volume
      {
	*dstVol = (char *) malloc(strlen(pch)+1);
	strcpy(*dstVol, pch);
	// check the existence of the file
	ret = stat(*dstVol, &stat_buf);
	if (ret != 0)
	{
	  fprintf(stderr, "INFO: Target volume %s cannot be found.\n", *dstVol);
	  free(*dstVol);
	  *dstVol = 0;
	  return 0;
	}
	else
	  fprintf(stdout, "INFO: Target volume %s\n", *dstVol);
	return 1;
      }
    }
    pch = strtok(NULL, " ");
  }
  return 0;
}

// find the volume and get the information
void mincGetVolumeInfo(const char *srcVol, VOL_GEOM *vgSrc)
{
  Volume vol;
  Status status;
  char *dim_names[4];
  int dim_sizes[4];
  int ndims;
  volume_input_struct input_info;
  Real voxel[4];
  Real separations[4];
  Real worldr, worlda, worlds;

  dim_names[0] = MIxspace;
  dim_names[1] = MIyspace;
  dim_names[2] = MIzspace;
  dim_names[3] = MItime;

  if (srcVol != 0)
  {
    status = start_volume_input((char *)srcVol, 0, dim_names, NC_UNSPECIFIED, 0, 0, 0, 
				TRUE, &vol, NULL, &input_info);
    if(status != OK)
    {
      // now check whether it is average_305
      if (strstr(srcVol, "average_305"))
      {
	printf("INFO: The transform was made with average_305.mnc.\n");
	// average_305 value
	vgSrc->width = 172; vgSrc->height = 220; vgSrc->depth = 156;
	vgSrc->xsize = 1; vgSrc->ysize = 1; vgSrc->zsize = 1;
	vgSrc->x_r = 1; vgSrc->x_a = 0; vgSrc->x_s = 0;
	vgSrc->y_r = 0; vgSrc->y_a = 1; vgSrc->y_s = 0;
	vgSrc->z_r = 0; vgSrc->z_a = 0; vgSrc->z_s = 1;      
	vgSrc->c_r = -0.0950;
	vgSrc->c_a = -16.5100;
	vgSrc->c_s = 9.7500;
	vgSrc->valid = 1;
      }
      else
      {
	printf("INFO: Could not find the original src volume %s\n", srcVol);
	printf("INFO: Set the src volume to the standard cor file.\n");
	initVolGeom(vgSrc);
      }
      return;
    }
    ndims = get_volume_n_dimensions(vol);
    if(ndims != 3 && ndims != 4)
      ErrorExit(ERROR_BADFILE, "ERROR:Dimension expecting 3 or 4 but got %d", ndims);
    
    /* ----- get the dimension sizes ----- */
    get_volume_sizes(vol, dim_sizes);
    
    /* --- one time point if there are only three dimensions in the file --- */
    if(ndims == 3) 
      dim_sizes[3] = 1;
    
    vgSrc->width = dim_sizes[0];
    vgSrc->height = dim_sizes[1];
    vgSrc->depth = dim_sizes[2];
    
    get_volume_separations(vol, separations);
    vgSrc->xsize = fabs(separations[0]);
    vgSrc->ysize = fabs(separations[1]);
    vgSrc->zsize = fabs(separations[2]);

    vgSrc->x_r = vol->direction_cosines[0][0];  
    vgSrc->x_a = vol->direction_cosines[0][1];  
    vgSrc->x_s = vol->direction_cosines[0][2];
    
    vgSrc->y_r = vol->direction_cosines[1][0];  
    vgSrc->y_a = vol->direction_cosines[1][1];  
    vgSrc->y_s = vol->direction_cosines[1][2];
    
    vgSrc->z_r = vol->direction_cosines[2][0];  
    vgSrc->z_a = vol->direction_cosines[2][1];  
    vgSrc->z_s = vol->direction_cosines[2][2];

    if(separations[0] < 0)
    {
      vgSrc->x_r = -vgSrc->x_r;  vgSrc->x_a = -vgSrc->x_a;  vgSrc->x_s = -vgSrc->x_s;
    }
    if(separations[1] < 0)
    {
      vgSrc->y_r = -vgSrc->y_r;  vgSrc->y_a = -vgSrc->y_a;  vgSrc->y_s = -vgSrc->y_s;
    }
    if(separations[2] < 0)
    {
      vgSrc->z_r = -vgSrc->z_r;  vgSrc->z_a = -vgSrc->z_a;  vgSrc->z_s = -vgSrc->z_s;
    }
    // get our definition of c_(r,a,s)
    voxel[0] = dim_sizes[0] / 2.0;
    voxel[1] = dim_sizes[1] / 2.0;
    voxel[2] = dim_sizes[2] / 2.0;
    voxel[3] = 0.0;
    convert_voxel_to_world(vol, voxel, &worldr, &worlda, &worlds);
    vgSrc->c_r = worldr; vgSrc->c_a = worlda; vgSrc->c_s = worlds;

    vgSrc->valid = 1;
    printf("INFO: Got the dst volume info.\n");
  }
  else
  {
    initVolGeom(vgSrc);
  }
  return;
}

void mincGetVolInfo(const char *xfmfile, VOL_GEOM *vgSrc, VOL_GEOM *vgDst)
{
  char *psrcVol=0;
  char *pdstVol=0;
  char line[1000];
  int retVal;

  FILE *fp = fopen(xfmfile, "r");
  if (!fp)
  {
    ErrorExit(ERROR_BADPARM, "Could not open xform") ;    
  }
  // skip one line at the header
  fgetl(line, 900, fp);
  fgetl(line, 900, fp);
  retVal = mincFindVolume(line, &psrcVol, &pdstVol);
  mincGetVolumeInfo(psrcVol, vgSrc); // src may not be found
  mincGetVolumeInfo(pdstVol, vgDst); // dst can be found
  free(psrcVol);
  free(pdstVol);
}

LTA *ltaMNIreadEx(const char *fname)
{
  LTA *lta = 0;
  lta = ltaMNIread((char *) fname);
  if (!lta)
    return NULL;
  // add original src and dst information
  mincGetVolInfo(fname, &lta->xforms[0].src, &lta->xforms[0].dst);   
  lta->type = LINEAR_RAS_TO_RAS;
  return lta;
}

LTA *ltaReadFileEx(const char *fname)
{
  FILE             *fp;
  LINEAR_TRANSFORM *lt ;
  int              i, nxforms, type ;
  char             line[STRLEN], *cp ;
  LTA              *lta ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(NULL,
                (ERROR_BADFILE, "ltaReadFile(%s): can't open file",fname));
  cp = fgetl(line, 199, fp) ; 
  if (cp == NULL)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data",fname));
  }
  sscanf(cp, "type      = %d\n", &type) ;
  cp = fgetl(line, 199, fp) ; sscanf(cp, "nxforms   = %d\n", &nxforms) ;
  lta = LTAalloc(nxforms, NULL) ;
  lta->type = type ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fscanf(fp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0) ;
    fscanf(fp, "sigma     = %f\n", &lt->sigma) ;
    MatrixAsciiReadFrom(fp, lt->m_L) ;
  }
  // oh, well this is the added part
  for (i=0; i < lta->num_xforms; i++)
  {
    if (fgets(line, 199, fp))
    {
      if (strncmp(line, "src volume info", 15)==0)
	readVolGeom(fp, &lta->xforms[i].src);
      fgets(line, 199, fp);
      if (strncmp(line, "dst volume info", 15)==0)
	readVolGeom(fp, &lta->xforms[i].dst);
    }
  }
  fclose(fp) ;
  return(lta) ;
}


LTA *
LTAreadEx(const char *fname)
{
  int       type ;
  LTA       *lta=0 ;
  MATRIX *V, *W, *m_tmp;

  type = TransformFileNameType((char *) fname) ;
  switch (type)
  {
  case REGISTER_DAT:
    printf("INFO: This REGISTER_DAT transform is valid only for volumes between "
	   " COR types with c_(r,a,s) = 0.");
    lta = ltaReadRegisterDat((char *) fname);
    if (!lta)
      return(NULL) ;

    V = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* world to voxel transform */
    W = MatrixAlloc(4, 4, MATRIX_REAL) ;  /* voxel to world transform */
    *MATRIX_RELT(V, 1, 1) = -1 ; *MATRIX_RELT(V, 1, 4) = 128 ;
    *MATRIX_RELT(V, 2, 3) = -1 ; *MATRIX_RELT(V, 2, 4) = 128 ;
    *MATRIX_RELT(V, 3, 2) = 1 ;  *MATRIX_RELT(V, 3, 4) = 128 ;
    *MATRIX_RELT(V, 4, 4) = 1 ;
    
    *MATRIX_RELT(W, 1, 1) = -1 ; *MATRIX_RELT(W, 1, 4) = 128 ;
    *MATRIX_RELT(W, 2, 3) = 1 ; *MATRIX_RELT(W, 2, 4) = -128 ;
    *MATRIX_RELT(W, 3, 2) = -1 ;  *MATRIX_RELT(W, 3, 4) = 128 ;
    *MATRIX_RELT(W, 4, 4) = 1 ;
    
    m_tmp = MatrixMultiply(lta->xforms[0].m_L, W, NULL) ;
    MatrixMultiply(V, m_tmp, lta->xforms[0].m_L) ;
    MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
    lta->type = LINEAR_VOX_TO_VOX ;
    break ;

  case MNI_TRANSFORM_TYPE:
    // this routine collect info on original src and dst volume
    // so that you can use the information to modify c_(r,a,s)
    // we no longer convert the transform to vox-to-vox
    // the transform is ras-to-ras
    lta = ltaMNIreadEx(fname) ;
    break ;

  case LINEAR_VOX_TO_VOX:
  case LINEAR_RAS_TO_RAS:
  case TRANSFORM_ARRAY_TYPE:
  default:
    // get the original src and dst information
    lta = ltaReadFileEx(fname);
    break ;
  }
  return(lta) ;
}

// write out the information on src and dst volume 
// in addition to the transform
int LTAprint(FILE *fp, const LTA *lta)
{
  int i;
  LT *lt;

  fprintf(fp, "type      = %d\n", lta->type) ;
  fprintf(fp, "nxforms   = %d\n", lta->num_xforms) ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fprintf(fp, "mean      = %2.3f %2.3f %2.3f\n", lt->x0, lt->y0, lt->z0) ;
    fprintf(fp, "sigma     = %2.3f\n", lt->sigma) ;
    MatrixAsciiWriteInto(fp, lt->m_L) ;
  }
  // write out src and dst volume info if there is one
  // note that this info may or may not be valid depending on vg->valid value
  // oh, well this is the addition
  for (i = 0; i < lta->num_xforms; ++i)
  {
    fprintf(fp, "src volume info\n");
    writeVolGeom(fp, &lta->xforms[i].src);
    fprintf(fp, "dst volume info\n");
    writeVolGeom(fp, &lta->xforms[i].dst);
  }
  return 1;
}

int
LTAwriteEx(const LTA *lta, const char *fname)
{
  FILE             *fp;
  time_t           tt ;
  char             *user, *time_str ;
  char             ext[STRLEN] ;

  if (!stricmp(FileNameExtension((char *) fname, ext), "XFM"))
  {
    fprintf(stderr, "ERROR: XFM format no longer supported as output.  Please use lta extension");
    return 0;
  }

  fp = fopen(fname,"w");
  if (fp==NULL) 
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "LTAwrite(%s): can't create file",fname));
  user = getenv("USER") ;
  if (!user)
    user = getenv("LOGNAME") ;
  if (!user)
    user = "UNKNOWN" ;
  tt = time(&tt) ;
  time_str = ctime(&tt) ;
  fprintf(fp, "# transform file %s\n# created by %s on %s\n",
          fname, user, time_str) ;
  LTAprint(fp, lta);
  fclose(fp) ;

  return 1;
}
