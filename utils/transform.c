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

#include "mrimorph.h"
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

static LTA  *ltaMNIread(char *fname) ;

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
  for (i = 0 ; i < nxforms ; i++)
  {
    lta->xforms[i].x0 = x0 ; 
    lta->xforms[i].y0 = y0 ; 
    lta->xforms[i].z0 = z0 ;
    lta->xforms[i].sigma = 10000.0f ;
    lta->xforms[i].m_L = MatrixIdentity(4, NULL) ;
    lta->xforms[i].m_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
    lta->xforms[i].m_last_dL = MatrixAlloc(4, 4, MATRIX_REAL) ;
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
  fprintf(fp, "type      = %d\n", LTA_TYPE) ;
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
  FILE             *fp;
  LINEAR_TRANSFORM *lt ;
  int              i, type, nxforms ;
  char             line[200], *cp ;
  LTA              *lta ;

  type = TransformFileNameType(fname) ;
  if (type == MNI_TRANSFORM_TYPE)
    return(ltaMNIread(fname)) ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(NULL,
                (ERROR_BADFILE, "LTAread(%s): can't open file",fname));
  cp = fgetl(line, 199, fp) ; sscanf(cp, "type      = %d\n", &type) ;
  cp = fgetl(line, 199, fp) ; sscanf(cp, "nxforms   = %d\n", &nxforms) ;
  lta = LTAalloc(nxforms, NULL) ;
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

  if (lta->num_xforms == 1)
    return(MRIlinearTransform(mri_src, mri_dst, lta->xforms[0].m_L)) ;

  fprintf(stderr, "applying octree transform to image...\n") ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  v_X     = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y     = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */
  m_L     = MatrixAlloc(4, 4, MATRIX_REAL) ;
  m_L_inv = MatrixAlloc(4, 4, MATRIX_REAL) ;
  v_Y->rptr[4][1] = 1.0f / mri_src->thick ;

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
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
         Same as LTAtransformPoint but also return weight normalization
         factor.
------------------------------------------------------*/
int
TransformFileNameType(char *fname)
{
  int file_type = TRANSFORM_ARRAY_TYPE ;
  char *dot, buf[500], *number ;

  strcpy(buf, fname) ;
  dot = strrchr(buf, '@') ;
  number = strchr(buf, '#') ;
  if (number)
    *number = 0 ;  /* don't consider : part of extension */

  if (!dot)
    dot = strrchr(buf, '.') ;

  if (dot)
  {
    dot++ ;
    StrUpper(buf) ;
    if (!strcmp(dot, "M3D"))
      return(MORPH_3D_TYPE) ;
    else if (!strcmp(dot, "LTA"))
      return(TRANSFORM_ARRAY_TYPE) ;
    else if (!strcmp(dot, "OCT"))
      return(TRANSFORM_ARRAY_TYPE) ;
    else if (!strcmp(dot, "XFM"))
      return(MNI_TRANSFORM_TYPE) ;
  }

  return(file_type) ;
}

#include "volume_io.h"

static LTA  *
ltaMNIread(char *fname)
{
  LTA              *lta ;
  LINEAR_TRANSFORM *lt ;
  char             *cp, line[1000] ;
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L, *V, *W, *m_tmp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL, 
                (ERROR_NOFILE, "ltMNIread: could not open file %s",fname));
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

  lta = LTAalloc(1, NULL) ;
  lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;

  fgetl(line, 900, fp) ;
  fgetl(line, 900, fp) ;
  fgetl(line, 900, fp) ;
  fgetl(line, 900, fp) ;
  /*  fgetl(line, 900, fp) ;*/

  m_L = lt->m_L ;
  for (row = 1 ; row <= 3 ; row++)
  {
    cp = fgetl(line, 900, fp) ;
    sscanf(cp, "%f %f %f %f",
           MATRIX_RELT(m_L,row,1), MATRIX_RELT(m_L,row,2), 
           MATRIX_RELT(m_L,row,3), MATRIX_RELT(m_L,row,4)) ;
  }
  m_tmp = MatrixMultiply(lt->m_L, W, NULL) ;
  MatrixMultiply(V, m_tmp, lt->m_L) ;
  /*  MatrixAsciiWriteInto(stderr, lt->m_L) ;*/
  MatrixFree(&V) ; MatrixFree(&W) ; MatrixFree(&m_tmp) ;
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
  V3_Z(v_X) = (y - 128.0) ;
  V3_Y(v_X) = (-z + 128.0) ;

  LTAinverseTransformPoint(lta, v_X, v_Y) ;

  /* voxel to world */
  *px = 128.0 - V3_X(v_Y) ;
  *py = V3_Z(v_Y)  - 128.0 ;
  *pz = -(V3_Y(v_Y) - 128.0) ;

  return(NO_ERROR) ;
}

