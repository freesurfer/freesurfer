#define COMPILING_MRISURF_METRIC_PROPERTIES
#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_metricProperties.h"

#include "mrisurf_MRIS.h"
#include "mrisurf_MRIS_MP.h"
#include "mrisurf_MRISPV.h"

#include "mrisurf_SurfaceFromMRIS_generated.h"
#include "mrisurf_SurfaceFromMRISPV_generated.h"
#include "mrisurf_SurfaceFromMRIS_MP_generated.h"

#include "face_barycentric_coords.h"
#include "mrisutils.h"
#include "surfgrad.h"
#include "surfcluster.h"

#include "mrisurf_base.h"


static int int_compare(const void* lhs_ptr, const void* rhs_ptr) {
   int lhs = *(int*)lhs_ptr;
   int rhs = *(int*)rhs_ptr;
   return lhs - rhs;
}



//==================================================================================================================
// Support code that accelerates finding the vertices and faces needed during defect correction.
//
int mrisurf_activeRealmTreesSize;
int mrisurf_orig_clock;
  // To do this, it must be able to tell when vertex orig[xyz] are changed.
  // Such changes, when relevant, need to be reported via noteVnoMovedInActiveRealmTrees.
  // To test it is correct, the code can scan all vertices of an mris and verify their origxyz are what was expected.


//==================================================================================================================
// Simple properties
//
// xyz are set in over 100 places
//
static void MRISsetXYZwkr1(MRIS *mris, const char * file, int line, bool* laterTime) {
  if (mris->dist_alloced_flags & 1) {
    copeWithLogicProblem2(laterTime, NULL,"dist not freed when setXYZ called.  Add call to MRISfreeDistsButNotOrig(MRIS*)",file,line,"<unknown>");
  }
}

static void MRISsetXYZwkr2(MRIS *mris, int vno, float x, float y, float z) {
  cheapAssertValidVno(mris,vno);
  VERTEX * v = &mris->vertices[vno];

  const float * pcx = &v->x;  float * px = (float*)pcx; *px = x;
  const float * pcy = &v->y;  float * py = (float*)pcy; *py = y;
  const float * pcz = &v->z;  float * pz = (float*)pcz; *pz = z;
}

void MRISsetXYZwkr(MRIS *mris, int vno, float x, float y, float z, const char * file, int line, bool* laterTime)  {
  MRISsetXYZwkr1(mris, file, line, laterTime);
  MRISsetXYZwkr2(mris, vno, x, y, z);
}

void MRISmemalignNFloats(size_t n, float** ppx, float** ppy, float** ppz) {
  // cache aligned to improve the performance of loops that use the vectors
  if (ppx) cheapAssert(0 == posix_memalign((void**)ppx, 64, n*sizeof(float)));
  if (ppy) cheapAssert(0 == posix_memalign((void**)ppy, 64, n*sizeof(float)));
  if (ppz) cheapAssert(0 == posix_memalign((void**)ppz, 64, n*sizeof(float)));
}

void MRISexportXYZ(MRIS *mris, float** ppx, float** ppy, float** ppz) {
  int const nvertices = mris->nvertices;

  float *px, *py, *pz;
  MRISmemalignNFloats(nvertices, &px, &py, &pz);
  
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];
    px[vno] = v->x; py[vno] = v->y; pz[vno] = v->z;
  }
  
  *ppx = px;
  *ppy = py;
  *ppz = pz;
}

#define SURFACE_DIMENSION_CALC_INIT \
  float xlo = 10000, ylo = xlo, zlo = ylo, xhi = -xlo, yhi = -ylo, zhi = -zlo;

#define SURFACE_DIMENSION_CALC_ITER \
    if (x > xhi) xhi = x; \
    if (x < xlo) xlo = x; \
    if (y > yhi) yhi = y; \
    if (y < ylo) ylo = y; \
    if (z > zhi) zhi = z; \
    if (z < zlo) zlo = z; \
    // end of macro

#define SURFACE_DIMENSION_CALC_FINI \
  mris->xlo = xlo; \
  mris->xhi = xhi; \
  mris->ylo = ylo; \
  mris->yhi = yhi; \
  mris->zlo = zlo; \
  mris->zhi = zhi; \
  \
  mris->xctr = 0.5f * (float)((double)xlo + (double)xhi); \
  mris->yctr = 0.5f * (float)((double)ylo + (double)yhi); \
  mris->zctr = 0.5f * (float)((double)zlo + (double)zhi); \
  // end of macro


void MRISimportXYZ(MRIS *mris, const float* const    px, const float* const    py, const float* const    pz) 
{
  static bool laterTime;
  MRISsetXYZwkr1(mris, __FILE__, __LINE__, &laterTime);

  int const nvertices = mris->nvertices;

  SURFACE_DIMENSION_CALC_INIT
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    float x = px[vno], y = py[vno], z = pz[vno];

    MRISsetXYZwkr2(mris, vno, x,y,z);

    SURFACE_DIMENSION_CALC_ITER
  }
  
  SURFACE_DIMENSION_CALC_FINI
}


void MRIScopyXYZ(MRIS *mris, MRIS* mris_from) {
  static bool laterTime;
  MRISsetXYZwkr1(mris, __FILE__, __LINE__, &laterTime);
  
  int const nvertices = mris->nvertices;
  cheapAssert(nvertices == mris_from->nvertices);
  
  SURFACE_DIMENSION_CALC_INIT
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX * v = &mris_from->vertices[vno];
    float x = v->x, y = v->y, z = v->z;
    
    MRISsetXYZwkr2(mris, vno, x, y, z);
    
    SURFACE_DIMENSION_CALC_ITER
  }
  
  SURFACE_DIMENSION_CALC_FINI
}


//  orig[xyz] are set during the creation of a surface
//  and during deformations that create an improved 'original' surface
//
void MRISsetOriginalXYZwkr(MRIS *mris, int vno, float origx, float origy, float origz, const char * file, int line, bool* laterTime) 
{
  if (mris->dist_alloced_flags & 2) {
    copeWithLogicProblem2(laterTime, NULL,"dist_orig not freed when setOriginalXYZ called.  Add call to MRISfreeDistOrigs(MRIS*)",file,line,"<unknown>");
  }

  cheapAssertValidVno(mris,vno);
  VERTEX * v = &mris->vertices[vno];
    
  const float * pcx = &v->origx;  float * px = (float*)pcx; *px = origx;
  const float * pcy = &v->origy;  float * py = (float*)pcy; *py = origy;
  const float * pcz = &v->origz;  float * pz = (float*)pcz; *pz = origz;

  if (hasActiveRealmTrees()) {
    noteVnoMovedInActiveRealmTrees(mris, vno);
  }
}


void MRISsetOriginalXYZfromXYZ(MRIS *mris) {
  
  MRISfreeDistOrigs(mris);  
    // Old values no longer valid
  
  mris->origxyz_status = mris->status;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];
    const float * pcx = &v->origx;  float * px = (float*)pcx; *px = v->x;
    const float * pcy = &v->origy;  float * py = (float*)pcy; *py = v->y;
    const float * pcz = &v->origz;  float * pz = (float*)pcz; *pz = v->z;
  }
  if (hasActiveRealmTrees()) {
    noteVnoMovedInActiveRealmTrees(mris, vno);
  }
}


/*-------------------------------------------------------------
  MRISavgVetexRadius() - computes the average and stddev of
  the distance from the origin to each vertex. If StdDev is NULL,
  it is ignored.
  -------------------------------------------------------------*/
double MRISavgVetexRadius(MRIS *mris, double *stdDev)
{
  double sum = 0, sum2 = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    double d = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    sum  += d;
    sum2 += (d * d);
  }

  int const N = mris->nvertices;
  double const avg = sum / N;
  if (stdDev != NULL) {
    *stdDev = sqrt(N * (sum2 / N - avg * avg) / (N - 1));
  }

  // printf("\n\nN = %ld, sum = %g, sum2 = %g, avg=%g, std = %g\n\n",
  // N,sum,sum2,avg,*stdDev);

  return avg;
}


int MRISimportVertexCoords(MRIS * const mris, float * locations[3], int const which)
{
  int const nvertices = mris->nvertices;

  if (which == CURRENT_VERTICES) MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...
   
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    switch (which) {
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISimportVertexCoords: which %d not supported", which);
        break;
      case CURRENT_VERTICES:
        v->x = locations[0][vno];
        v->y = locations[1][vno];
        v->z = locations[2][vno];
        break;
    }
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------*/
/*!
  \fn int MRISreverseCoords(MRIS *mris, int which_reverse, int reverse_face_order, int which_coords)
  \brief Reverse sign of one of the dimensions of the surface coords.
  If reversing X, the order of the verticies is also reversed.
*/
int MRISreverseCoords(MRIS *mris, int which_direction, int reverse_face_order, int which_coords)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    float x = 0, y = 0, z = 0;

    switch (which_coords) {
      case CURRENT_VERTICES:
        x = v->x;
        y = v->y;
        z = v->z;
        break;
      case CANONICAL_VERTICES:
        x = v->cx;
        y = v->cy;
        z = v->cz;
        break;
      case ORIGINAL_VERTICES:
        x = v->origx;
        y = v->origy;
        z = v->origz;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISreverseCoords: unsupported which_vertices %d", which_coords);
    }

    switch (which_direction) {
      default:
      case REVERSE_X:
        x = -x;
        break;
      case REVERSE_Y:
        y = -y;
        break;
      case REVERSE_Z:
        z = -z;
        break;
    }

    switch (which_coords) {
      case CURRENT_VERTICES:
        v->x = x;
        v->y = y;
        v->z = z;
        break;
      case CANONICAL_VERTICES:
        v->cx = x;
        v->cy = y;
        v->cz = z;
        break;
      case ORIGINAL_VERTICES:
        MRISsetOriginalXYZ(mris, vno, x,y,z);
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISreverseCoords: unsupported which_vertices %d", which_coords);
    }
  }
  if (which_direction == REVERSE_X && reverse_face_order)  // swap order of faces
  {
    MRISreverseFaceOrder(mris);
  }

  return (NO_ERROR);
}



/*-----------------------------------------------------*/
/*!
  \fn int MRISreverse(MRIS *mris, int which)
  \brief Reverse sign of one of the dimensions of the surface coords.
  If reversing X, the order of the verticies is also reversed.
*/
int MRISreverse(MRIS *mris, int which, int reverse_face_order)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    float x = v->x;
    float y = v->y;
    float z = v->z;
    switch (which) {
      default:
      case REVERSE_X:
        x = -x;
        break;
      case REVERSE_Y:
        y = -y;
        break;
      case REVERSE_Z:
        z = -z;
        break;
    }
    v->x = x;
    v->y = y;
    v->z = z;
  }

  if (which == REVERSE_X && reverse_face_order)  // swap order of faces
  {
    MRISreverseFaceOrder(mris);
  }

  return (NO_ERROR);
}


int MRISscale(MRIS *mris, double scale)
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...
   
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    v->x *= scale;
    v->y *= scale;
    v->z *= scale;
  }

  // Emulate mrisComputeSurfaceDimensions(mris)
  //
  mris->xlo *= scale;
  mris->xhi *= scale;
  mris->ylo *= scale;
  mris->yhi *= scale;
  mris->zlo *= scale;
  mris->zhi *= scale;
  
  if (scale < 0) {
    float xlo = mris->xlo; mris->xlo = mris->xhi; mris->xlo = xlo;
    float ylo = mris->ylo; mris->ylo = mris->yhi; mris->ylo = ylo;
    float zlo = mris->zlo; mris->zlo = mris->zhi; mris->zlo = zlo;
  }
  
  mris->xctr *= scale;
  mris->yctr *= scale;
  mris->zctr *= scale;
  
  return (NO_ERROR);
}


int mrisFlipPatch(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    v->x = -v->x;
  }

  float xlo = mris->xlo; mris->xlo = -mris->xhi; mris->xhi = -xlo;
  mris->xctr = -mris->xctr;
  
  return (NO_ERROR);
}



int MRIStranslate(MRIS *mris, float dx, float dy, float dz)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    v->x += dx;
    v->y += dy;
    v->z += dz;
  }

  // Emulate mrisComputeSurfaceDimensions(mris)
  //
  mris->xlo += dx;
  mris->xhi += dx;
  mris->ylo += dy;
  mris->yhi += dy;
  mris->zlo += dz;
  mris->zhi += dz;
  
  mris->xctr += dx;
  mris->yctr += dy;
  mris->zctr += dz;

  return (NO_ERROR);
}


void MRISscaleThenTranslate (MRIS *mris, double sx, double sy, double sz, double dx, double dy, double dz) {
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...
   
  //
  // This uses double because mri_brain_volume was using double,
  // and because the combined scaling and adding could be much less accurate in float.
  //
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    v->x = v->x*sx + dx;
    v->y = v->y*sy + dy;
    v->z = v->z*sz + dz;
  }

  // Emulate mrisComputeSurfaceDimensions(mris)
  //
  double xlo = mris->xlo;
  double xhi = mris->xhi;
  double ylo = mris->ylo;
  double yhi = mris->yhi;
  double zlo = mris->zlo;
  double zhi = mris->zhi;
  double xctr = mris->xctr;
  double yctr = mris->yctr;
  double zctr = mris->zctr;
  
  xlo *= sx;
  xhi *= sx;
  ylo *= sy;
  yhi *= sy;
  zlo *= sz;
  zhi *= sz;
  
  if (sx < 0) { double t = xlo; xlo = xhi; xlo = t; }
  if (sy < 0) { double t = ylo; ylo = yhi; ylo = t; }
  if (sz < 0) { double t = zlo; zlo = zhi; zlo = t; }
  
  xctr *= sx;
  yctr *= sy;
  zctr *= sz;
  
  xlo += dx;
  xhi += dx;
  ylo += dy;
  yhi += dy;
  zlo += dz;
  zhi += dz;
  
  xctr += dx;
  yctr += dy;
  zctr += dz;

  mris->xlo = (float)xlo;
  mris->xhi = (float)xhi;
  mris->ylo = (float)ylo;
  mris->yhi = (float)yhi;
  mris->zlo = (float)zlo;
  mris->zhi = (float)zhi;
  mris->xctr = (float)xctr;
  mris->yctr = (float)yctr;
  mris->zctr = (float)zctr;
}


int MRISanisotropicScale(MRIS *mris, float sx, float sy, float sz)
{
  mrisComputeSurfaceDimensions(mris);
 
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...
   
  if (sx == 1.0 && sy == 1.0 && sz == 1.0) return (NO_ERROR);
    // Not just for speed.   (v - d) + d  !=  v   in floating point arithmetic
    // and this has been seen to cause slight differences in outputs between code 
    // that did and did not scale by 1x 1x 1x before doing something else

  /* scale around the center */
  float const x0 = mris->xctr;
  float const y0 = mris->yctr;
  float const z0 = mris->zctr;

  int k;
  for (k = 0; k < mris->nvertices; k++) {
    VERTEX *v = &mris->vertices[k];
    v->x = (v->x - x0) * sx + x0;
    v->y = (v->y - y0) * sy + y0;
    v->z = (v->z - z0) * sz + z0;
  }

  // emulate mrisComputeSurfaceDimensions(mris)
  //
  float xlo = (mris->xlo - x0) * sx + x0;
  float xhi = (mris->xhi - x0) * sx + x0;
  float ylo = (mris->ylo - y0) * sy + y0;
  float yhi = (mris->yhi - y0) * sy + y0;
  float zlo = (mris->zlo - z0) * sz + z0;
  float zhi = (mris->zhi - z0) * sz + z0;

  if (xlo < xhi) { mris->xlo = xlo; mris->xhi = xhi; } else { mris->xlo = xhi; mris->xhi = xlo; }
  if (ylo < yhi) { mris->ylo = ylo; mris->yhi = yhi; } else { mris->ylo = yhi; mris->yhi = ylo; }
  if (zlo < zhi) { mris->zlo = zlo; mris->zhi = zhi; } else { mris->zlo = zhi; mris->zhi = zlo; }
  
  mris->xctr = 0.5f * (float)((double)xlo + (double)xhi);
  mris->yctr = 0.5f * (float)((double)ylo + (double)yhi);
  mris->zctr = 0.5f * (float)((double)zlo + (double)zhi);
  
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  The following function is broken, since it is applying a
  transform for surfaceRAS space.  If transform has c_(ras)
  values, the result would be different.

  conformed -----> surfaceRAS
  |                |        [ 1 Csrc]
  V                V        [ 0  1  ]
  src    ----->    RAS
  |                |
  |                |  Xfm
  V                V
  talvol   ----->  talRAS
  |                |        [ 1 -Ctal]
  |                |        [ 0  1   ]
  V                V
  conformed -----> surfaceRAS

  Thus

  surfRASToTalSurfRAS = [ 1 -Ctal ]*[ R  T ]*[ 1 Csrc ]=
  [ R  T + R*Csrc - Ctal ]
  [ 0   1   ] [ 0  1 ] [ 0   1  ]  [ 0         1          ]

  We need to know the Csrc and Ctal values
  ------------------------------------------------------*/
MRIS* MRIStalairachTransform(MRIS* mris_src, MRIS* mris)
{
  if (!mris) {
    mris = MRISclone(mris_src);
  }

  if (!mris_src->lta) {
    return mris;
  }

  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  SURFACE_DIMENSION_CALC_INIT

  int vno;
  for (vno = 0; vno < mris_src->nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];

    double x,y,z;
    TransformWithMatrix(mris_src->SRASToTalSRAS_, v->x, v->y, v->z, &x, &y, &z);
    v->x = x;
    v->y = y;
    v->z = z;

    SURFACE_DIMENSION_CALC_ITER
  }

  SURFACE_DIMENSION_CALC_FINI

  return mris;
}

template <class Surface>
void MRISrotate(Surface surface, float alpha, float beta, float gamma)
{
  surface.freeDistsButNotOrig();

  float const sa = sin(alpha);
  float const sb = sin(beta);
  float const sg = sin(gamma);
  float const ca = cos(alpha);
  float const cb = cos(beta);
  float const cg = cos(gamma);

  float const cacb = ca * cb;
  float const cacgsb = ca * cg * sb;
  float const sasg = sa * sg;
  float const cgsa = cg * sa;
  float const casbsg = ca * sb * sg;
  float const cbsa = cb * sa;
  float const cgsasb = cg * sa * sb;
  float const casg = ca * sg;
  float const cacg = ca * cg;
  float const sasbsg = sa * sb * sg;
  float const cbcg = cb * cg;
  float const cbsg = cb * sg;

  auto const nvertices = surface.nvertices();

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) schedule(guided)
#endif
  for (int vno = 0; vno < nvertices; vno++) {
    ROMP_PFLB_begin

    if (vno == Gdiag_no) {
      DiagBreak();
    }

    auto vertex = surface.vertices(vno);
    float x = vertex.x();
    float y = vertex.y();
    float z = vertex.z();
    float xp =  x * cacb + z * (-cacgsb - sasg) + y * (cgsa - casbsg);
    float yp = -x * cbsa + z * ( cgsasb - casg) + y * (cacg + sasbsg);
    float zp =  z * cbcg + x * sb + y * cbsg;
    vertex.set_x(xp);
    vertex.set_y(yp);
    vertex.set_z(zp);

    ROMP_PFLB_end
  }
  ROMP_PF_end
}

void MRISrotate(MRIS* mris, float alpha, float beta, float gamma)
{
    SurfaceFromMRIS::XYZPositionM::Surface surface(mris);
    MRISrotate(surface, alpha, beta, gamma);
}

void MRISrotate(MRIS_MP* mris, float alpha, float beta, float gamma)
{
    SurfaceFromMRIS_MP::XYZPositionM::Surface surface(mris);
    MRISrotate(surface, alpha, beta, gamma);
}

MRIS* MRISrotate(MRIS* mris_src, MRIS* mris_dst, float alpha, float beta, float gamma)
{
  if (!mris_dst) {
    mris_dst = MRISclone(mris_src);
  }
  MRISrotate(mris_dst, alpha, beta, gamma);
  return (mris_dst);
}


/*------------------------------------------------------------------------
  MRISmatrixMultiply() - simply multiplies matrix M by the vertex xyz.
  See also MRIStransform() and MRISltaMultiply().
  ------------------------------------------------------------------------*/
int MRISmatrixMultiply(MRIS *mris, MATRIX *M)
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  int vno;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX *v;
    MATRIX *xyz, *Mxyz;
    xyz = MatrixAlloc(4, 1, MATRIX_REAL);
    xyz->rptr[4][1] = 1.0;
    Mxyz = MatrixAlloc(4, 1, MATRIX_REAL);
    v = &mris->vertices[vno];
    xyz->rptr[1][1] = v->x;
    xyz->rptr[2][1] = v->y;
    xyz->rptr[3][1] = v->z;
    MatrixMultiply(M, xyz, Mxyz);

    v->x = Mxyz->rptr[1][1];
    v->y = Mxyz->rptr[2][1];
    v->z = Mxyz->rptr[3][1];
    MatrixFree(&xyz);
    MatrixFree(&Mxyz);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (0);
}


/*-------------------------------------------------------*/
/* Replaces x,y,z with theta,phi,radius. Assumes
    that the surface xyz are already on the sphere.
    Note: this is not related to the surface-based
    spherical coords.
 */
int MRISsphericalCoords(MRIS* mris)
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  int k;
  double x, y, z, d2, d, r, theta, phi;

  for (k = 0; k < mris->nvertices; k++) {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    d2 = x * x + y * y;
    d = sqrt(d2);
    r = sqrt(d2 + z * z);
    theta = atan2(y, x);
    phi = atan2(z, d);
    mris->vertices[k].x = theta;
    mris->vertices[k].y = phi;
    mris->vertices[k].z = r;
  }
  
  return (0);
}



void mrisDisturbVertices(MRIS* mris, double amount)
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  int vno ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX *v = &mris->vertices[vno] ;
    v->x += randomNumber(-amount, amount) ;
    v->y += randomNumber(-amount, amount) ;
  }

  MRIScomputeMetricProperties(mris) ;
}


void MRISmoveOrigin(MRIS* mris, float x0, float y0, float z0)
{
  mrisComputeSurfaceDimensions(mris);
  MRIStranslate(mris, x0 - mris->xctr, y0 - mris->yctr, z0 - mris->zctr);
  mris->xctr = x0;  // adjust for tiny numeric errors in MRIStranslate
  mris->yctr = y0;
  mris->zctr = z0;
}


MRIS* MRIScenter(MRIS* mris_src, MRIS* mris_dst)
{
  cheapAssert(mris_src == mris_dst);
  MRISmoveOrigin(mris_dst, 0.0f, 0.0f, 0.0f);
  return mris_dst;
}

void MRIScalculateCenterCOG2(MRIS *mris, double *xCOG, double *yCOG, double *zCOG)
{
  double x = 0.0, y = 0.0, z = 0.0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    x += v->x;
    y += v->y;
    z += v->z;
  }
  (*xCOG) = x / mris->nvertices;
  (*yCOG) = y / mris->nvertices;
  (*zCOG) = z / mris->nvertices;
}

// translate the COG of a surface to (0,0,0)
void MRIScenterCOG2(MRIS *mris, double *xCOG, double *yCOG, double *zCOG)
{
  double x, y, z;
  MRIScalculateCenterCOG2(mris, &x,&y,&z);
  MRIStranslate(mris, -x, -y, -z);

  if (xCOG || yCOG || zCOG) {
    (*xCOG) = x;
    (*yCOG) = y;
    (*zCOG) = z;
  }
  /*       fprintf(stderr,"\nCOG Centered at x=%f y=%f z=%f",
     (float)x,(float)y,(float)z);*/
}

void MRIScenterCOG(MRIS *mris)
{ 
  MRIScenterCOG2(mris, NULL, NULL, NULL); 
}



MRIS* MRISprojectOntoTranslatedSphere(
    MRIS *mris_src, MRIS * mris_dst, 
    double r,
    double x0, double y0, double z0) {
    
  if (FZERO(r))
    r = DEFAULT_RADIUS ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  if ((mris_dst->status != MRIS_SPHERE) &&
      (mris_dst->status != MRIS_PARAMETERIZED_SPHERE))
    MRIScenter(mris_dst, mris_dst) ;

  mris_dst->radius = r ;

  MRISfreeDistsButNotOrig(mris_dst);  // it is either this or adjust them...

  int vno ;
  for (vno = 0 ; vno < mris_dst->nvertices ; vno++) {
    VERTEX* const v = &mris_dst->vertices[vno];

    if (v->ripflag)  /* shouldn't happen */
      continue ;
    
    double x = v->x, x2 = x*x;
    double y = v->y, y2 = y*y;
    double z = v->z, z2 = z*z;

    double dist = sqrt(x2+y2+z2) ;
    
    double d = FZERO(dist) ? 0.0 : (1 - r / dist);

    v->x = x - d*x;
    v->y = y - d*y;
    v->z = z - d*z;

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z))
      DiagBreak() ;
  }

  MRIStranslate(mris_dst, x0,y0,z0);

  MRISupdateEllipsoidSurface(mris_dst) ;

  mris_dst->status = mris_src->status == MRIS_PARAMETERIZED_SPHERE ?
                     MRIS_PARAMETERIZED_SPHERE : MRIS_SPHERE ;

  return(mris_dst) ;
}


void MRISblendXYZandTXYZ(MRIS* mris, float xyzScale, float txyzScale) {

  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    v->x = xyzScale*v->x + txyzScale*v->tx;
    v->y = xyzScale*v->y + txyzScale*v->ty;
    v->z = xyzScale*v->z + txyzScale*v->tz;
  }

  // current only user did not have this, but did immediately call MRIScomputeMetricProperties(mris)
  //
  // mrisComputeSurfaceDimensions(mris);
}


void MRISblendXYZandNXYZ(MRIS* mris, float nxyzScale) 
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    v->x = v->x + nxyzScale*v->nx;
    v->y = v->y + nxyzScale*v->ny;
    v->z = v->z + nxyzScale*v->nz;
  }

  // old form did not have this
  //
  // mrisComputeSurfaceDimensions(mris);
}


template <class SurfaceOut, class SurfaceInp>
void MRIStranslate_along_vertex_dxdydz(SurfaceOut surfaceOut, SurfaceInp surfaceInp, double dt)
{
  surfaceOut.freeDistsButNotOrig();

  auto const nvertices = surfaceInp.nvertices();
  cheapAssert(nvertices == surfaceOut.nvertices());

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) schedule(guided)
#endif
  for (int vno = 0; vno < nvertices; vno++) {
    ROMP_PFLB_begin
      
    auto vertexInp = surfaceInp.vertices(vno);
    if (vertexInp.ripflag()) continue;
    
    float x = vertexInp.x();
    float y = vertexInp.y();
    float z = vertexInp.z();

    x += dt * vertexInp.dx();
    y += dt * vertexInp.dy();
    z += dt * vertexInp.dz();


    auto vertexOut = surfaceOut.vertices(vno);
    vertexOut.set_x(x);
    vertexOut.set_y(y);
    vertexOut.set_z(z);

    ROMP_PFLB_end
  }
  ROMP_PF_end
}

void MRIStranslate_along_vertex_dxdydz(MRIS*    dst, MRIS*    src, double dt)
{
    SurfaceFromMRIS::Distort::Surface      surfaceInp(src);
    SurfaceFromMRIS::XYZPositionM::Surface surfaceOut(dst);
    MRIStranslate_along_vertex_dxdydz(surfaceOut, surfaceInp, dt);
}

void MRIStranslate_along_vertex_dxdydz(MRIS_MP* dst, MRIS_MP* src, double dt)
{
    SurfaceFromMRIS_MP::Distort::Surface      surfaceInp(src);
    SurfaceFromMRIS_MP::XYZPositionM::Surface surfaceOut(dst);
    MRIStranslate_along_vertex_dxdydz(surfaceOut, surfaceInp, dt);
}


static void MRISaverageVertexPositionsWkr_part1(MRIS* mris, int navgs, float scale, bool includeCenter, int start, int step, float **ppx, float **ppy, float **ppz)
{
  MRISfreeDistsButNotOrig(mris);  // it is either this or adjust them...

  MRISexportXYZ(mris, ppx, ppy, ppz);

  int const nvertices = mris->nvertices;

  // first inputs
  //
  float *px = *ppx, *py = *ppy, *pz = *ppz;
  
  // first outputs
  //
  float *qx, *qy, *qz;
  MRISmemalignNFloats(nvertices, &qx, &qy, &qz);
  
  // iterate
  //
  while (navgs-- > 0) {
    int vno;
    for (vno = start; vno < nvertices; vno += step) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      float x, y, z;
      int count;
      if (includeCenter) {
        x = px[vno], y = py[vno], z = pz[vno];
        count = 1;
      } else {
        x = y = z = 0.0f;
        count = 0;    
      }
      if (!v->ripflag) { 
        int n;
        for (n = 0; n < vt->vnum; n++) {
          int vno2 = vt->v[n];
          VERTEX const * const vn = &mris->vertices[vno2];
          if (vn->ripflag) continue;
          count++;
          x += px[vno2]; y += py[vno2]; z += pz[vno2];
        }
      }
      x /= count;
      y /= count;
      z /= count;
      if (scale == 1.0f) {
        qx[vno] = x;
        qy[vno] = y;
        qz[vno] = z;
      } else {
        qx[vno] = px[vno] + scale * (x - px[vno]);
        qy[vno] = py[vno] + scale * (y - py[vno]);
        qz[vno] = pz[vno] + scale * (z - pz[vno]);
      }
    }
    
    // swap for next round
    float* t;
    t = px; px = qx; qx = t;
    t = py; py = qy; qy = t;
    t = pz; pz = qz; qz = t;
  }
  
  free(qx); free(qy), free(qz);
  
  *ppx = px; *ppy = py; *ppz = pz;
}

static void MRISaverageVertexPositionsWkr_part2(MRIS* mris, float *px, float *py, float *pz, bool checking)
{
  if (!checking) MRISimportXYZ(mris, px, py, pz);
  else {
    int const nvertices = mris->nvertices;
    int vno;
    for (vno = 0; vno < nvertices; vno++) {
      VERTEX const * const v  = &mris->vertices[vno];
      cheapAssert(px[vno] == v->x);
      cheapAssert(py[vno] == v->y);
      cheapAssert(pz[vno] == v->z);
    }

    // current only user did not have this
    // and using MRISimportXYZ may cause a problem because it does...
    //
    // mrisComputeSurfaceDimensions(mris);
  }
  
  free(px); free(py), free(pz);
  

}

int MRISaverageVertexPositions(MRIS *mris, int navgs)
{
  // This passes mris_inflate test
  //
  bool checking = false;
  if (checking) {
    switch (copeWithLogicProblem("FREESURFER_fix_MRISaverageVertexPositions","MRISaverageVertexPositions was inefficient")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      checking = false;
    }
  }
  
  float *px = NULL, *py = NULL, *pz = NULL;
  
  MRISaverageVertexPositionsWkr_part1(mris, navgs, 1.0f, true, 0, 1, &px, &py, &pz);
  
  if (checking) {
    int i, vno, vnb, vnum;
    float x, y, z, num;

    for (i = 0; i < navgs; i++) {
      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX                * const v  = &mris->vertices         [vno];
        if (v->ripflag) {
          continue;
        }
        x = v->x;
        y = v->y;
        z = v->z;
        int const * pnb  = vt->v;
        vnum = vt->vnum;
        for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
          VERTEX * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
          if (vn->ripflag) {
            continue;
          }
          num++;
          x += vn->x;
          y += vn->y;
          z += vn->z;
        }
        num++; /* account for central vertex */
        v->tdx = x / num;
        v->tdy = y / num;
        v->tdz = z / num;
      }
      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX * const v = &mris->vertices[vno];
        if (v->ripflag) {
          continue;
        }
        v->x = v->tdx;
        v->y = v->tdy;
        v->z = v->tdz;
      }
    }
  }
    
  MRISaverageVertexPositionsWkr_part2(mris, px, py, pz, checking);

  return (NO_ERROR);
}


// smooth a surface 'niter' times with a step (should be around 0.5)
//
void MRISsmoothSurface(MRIS *mris, int niter, float step)
{
  bool useOldBehaviour = false;
  switch (copeWithLogicProblem("FREESURFER_fix_MRISsmoothSurface1","MRISsmoothSurface was calling MRIScomputeMetricProperties too soon")) {
  case LogicProblemResponse_old: 
    break;
  case LogicProblemResponse_fix:
    useOldBehaviour = false;
  }

  bool checking = false;
  if (checking) {
    switch (copeWithLogicProblem("FREESURFER_fix_MRISsmoothSurface2","MRISsmoothSurface was not using common efficient code")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      checking = false;
    }
  }
  
  if (step > 1) step = 1.0f;

  float *px = NULL, *py = NULL, *pz = NULL;
  
  MRISaverageVertexPositionsWkr_part1(mris, niter, step, false, 0, 1, &px, &py, &pz);
  
  if (checking) {
    int iter;
    for (iter = 0; iter < niter; iter++) {
      if (useOldBehaviour) MRIScomputeMetricProperties(mris);

      int k;
      for (k = 0; k < mris->nvertices; k++) {
        VERTEX * v = &mris->vertices[k];
        v->tx = v->x;
        v->ty = v->y;
        v->tz = v->z;
      }

      for (k = 0; k < mris->nvertices; k++) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
        VERTEX                * const v  = &mris->vertices         [k];
        int n = 0;
        float x = 0, y = 0, z = 0;
        int m;
        for (m = 0; m < vt->vnum; m++) {
          x += mris->vertices[vt->v[m]].tx;
          y += mris->vertices[vt->v[m]].ty;
          z += mris->vertices[vt->v[m]].tz;
          n++;
        }
        x /= n;
        y /= n;
        z /= n;

        v->x = v->x + step * (x - v->x);
        v->y = v->y + step * (y - v->y);
        v->z = v->z + step * (z - v->z);
      }
    }
  }
  
  MRISaverageVertexPositionsWkr_part2(mris, px, py, pz, checking);
  
  if (!useOldBehaviour) MRIScomputeMetricProperties(mris);
}


/* Center the surface mris at location (cx,cy,cz) with a radius r
   such the energy sum((x-cx)^2+(y-cy)^2+(z-cz)^2-r^2)^2 is minimized 
*/

#define DEBUG_CENTER_SURFACE 0

/* this set of functions centers the surface mris,
   assuming that its vertices
   (using ORIGINAL_VERTICES) are sampled onto a
   sphere of radius R to be determined */

static double estimateNRG          (MRIS *mris, double cx, double cy, double cz, double R2);
static double estimateSquaredRadius(MRIS *mris, double cx, double cy, double cz);
static void   computeGradient      (MRIS *mris, double cx, double cy, double cz, double R2, double *gx, double *gy, double *gz);

/* compute the NRG associated with the center
   (cx,cy,cz) and the squared radius R2
   NRG = sum((x-cx)^2+(y-cy)^2+(z-cz)^2-R2)^2 */
static double estimateNRG(MRIS *mris, double cx, double cy, double cz, double R2)
{
  double NRG = 0.0;
  int n;
  for (n = 0; n < mris->nvertices; n++) {
    NRG += SQR(SQR(mris->vertices[n].x - cx) + SQR(mris->vertices[n].y - cy) + SQR(mris->vertices[n].z - cz) - R2);
  }
  return NRG;
}

/* estimate the squared radius that
   minimizes the energy defined above for a given
   center (cx,cy,cz) */
static double estimateSquaredRadius(MRIS *mris, double cx, double cy, double cz)
{
  double R2 = 0.0;

  int n;
  for (n = 0; n < mris->nvertices; n++) {
    R2 += SQR(mris->vertices[n].x - cx) + SQR(mris->vertices[n].y - cy) + SQR(mris->vertices[n].z - cz);
    if (!std::isfinite(R2)) {
      DiagBreak();
    }
  }
  return (R2 / (double)mris->nvertices);
}

/* compute the gradient of the energy
   defined above at location (cx,cy,cz,R2) */
static void computeGradient(MRIS *mris, double cx, double cy, double cz, double R2, double *gx, double *gy, double *gz)
{
  double tx = 0.0, ty = 0.0, tz = 0.0;
  int n;
  for (n = 0; n < mris->nvertices; n++) {
    tx += (SQR(mris->vertices[n].x - cx) + SQR(mris->vertices[n].y - cy) + SQR(mris->vertices[n].z - cz) - R2) *
          (mris->vertices[n].x - cx);
    ty += (SQR(mris->vertices[n].x - cx) + SQR(mris->vertices[n].y - cy) + SQR(mris->vertices[n].z - cz) - R2) *
          (mris->vertices[n].y - cy);
    tz += (SQR(mris->vertices[n].x - cx) + SQR(mris->vertices[n].y - cy) + SQR(mris->vertices[n].z - cz) - R2) *
          (mris->vertices[n].z - cz);
  }
  (*gx) = tx / (double)mris->nvertices;
  (*gy) = ty / (double)mris->nvertices;
  (*gz) = tz / (double)mris->nvertices;
}

static void MRIScenterSphere_old(MRIS* mris);
static void MRIScenterSphere_new(MRIS* mris);

void MRIScenterSphere(MRIS* mris)
{
  ROMP_SCOPE_begin
  MRIScenterSphere_old(mris);
  ROMP_SCOPE_end
}

static void MRIScenterSphere_new(MRIS* mris)
{
  cheapAssert(!"MRIScenterSphere_new TBD");
}


static void MRIScenterSphere_old(MRIS* mris)
{
  VERTEX *vertex;
  int n, niters;
  /* sphere parameters */
  double x, y, z, R2, xhi, xlo, yhi, ylo, zhi, zlo, cx, cy, cz, radius, scale;
  /* NRG parameters */
  double NRG, last_NRG;
  /* gradient parameters */
  double dx, dy, dz, d, epsilon;

  fprintf(WHICH_OUTPUT, "Finding true center and radius of Spherical Surface...");

  /* compute an initial estimate of the the center (x,y,z) */
  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (n = 0; n < mris->nvertices; n++) {
    vertex = &mris->vertices[n];
    x = (double)vertex->x;
    y = (double)vertex->y;
    z = (double)vertex->z;

    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
  }
  x = (xlo + xhi) / 2.0f;
  y = (ylo + yhi) / 2.0f;
  z = (zlo + zhi) / 2.0f;

  /* estimate the corresponding squared radius */
  R2 = estimateSquaredRadius(mris, x, y, z);

  /* compute the initial NRG */
  NRG = estimateNRG(mris, x, y, z, R2);

  /* verify if the center (0,0,0) with radius 100.0 is a better candidate */
  if (estimateNRG(mris, 0.0, 0.0, 0.0, 10000.0) < NRG) {
    x = 0.0;
    y = 0.0;
    z = 0.0;
    R2 = 10000.0;
    NRG = estimateNRG(mris, x, y, z, R2);
  }

  if (estimateNRG(mris, 0, 0, 0, estimateSquaredRadius(mris, 0, 0, 0)) < NRG) {
    x = y = z = 0.0;
    R2 = estimateSquaredRadius(mris, x, y, z);
    NRG = estimateNRG(mris, x, y, z, R2);
  }

#if DEBUG_CENTER_SURFACE
  fprintf(
      WHICH_OUTPUT, "\nInitial Configuration: NRG=%lf, R=%f and ( %f , %f , %f )", 100000.0 * NRG, sqrt(R2), x, y, z);
#endif

  /* iteratively minize the NRG */
  last_NRG = NRG + 1.0;
  niters = 0;
  while (NRG < last_NRG) {
    niters++;
    if (niters > 100) {
      break;
    }
    last_NRG = NRG;

    /* compute the gradient */
    computeGradient(mris, x, y, z, R2, &dx, &dy, &dz);

    /* bound gradient by displacement of 1.0 mm */
    d = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
    if (d > 1.0) {
      dx /= d;
      dy /= d;
      dz /= d;
    }

    /*   fprintf(WHICH_OUTPUT,"\n gradient:(%f,%f,%f)",dx,dy,dz); */

    epsilon = 2.0;
    while (NRG >= last_NRG) {
      epsilon /= 2.0;
      NRG = estimateNRG(mris,
                        x + epsilon * dx,
                        y + epsilon * dy,
                        z + epsilon * dz,
                        estimateSquaredRadius(mris, x + epsilon * dx, y + epsilon * dy, z + epsilon * dz));
      d = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
      if (epsilon * d < 0.00000000001)  // FZERO(epsilon*(SQR(dx)+SQR(dy)+SQR(dz))))
      {
        break;
      }
    }

    if (NRG < last_NRG) {
      x = x + epsilon * dx;
      y = y + epsilon * dy;
      z = z + epsilon * dz;

      R2 = estimateSquaredRadius(mris, x, y, z);

#if DEBUG_CENTER_SURFACE
      fprintf(WHICH_OUTPUT,
              "\nNew Minimum found: NRG=%lf, R=%f and ( %f , %f , %f )",
              100000.0 * NRG,
              sqrt(R2),
              10000. * x,
              10000. * y,
              10000. * z);
#endif
    }
    else {
      NRG = estimateNRG(mris, x, y, z, R2);
    }
  }

  /* now centering the surface at position (0,0,0) with radius 100.0 */
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;
  radius = 100.0;
  scale = sqrt(SQR(radius) / R2);
  for (n = 0; n < mris->nvertices; n++) {
    vertex = &mris->vertices[n];
    vertex->x = cx + scale * (vertex->x - x);
    vertex->y = cy + scale * (vertex->y - y);
    vertex->z = cz + scale * (vertex->z - z);
  }

  // scaling onto sphere with the exact right radius DEFAULT_RADIUS
  for (n = 0; n < mris->nvertices; n++) {
    vertex = &mris->vertices[n];
    R2 = SQR(vertex->x) + SQR(vertex->y) + SQR(vertex->z);
    scale = DEFAULT_RADIUS / sqrt(R2);
    vertex->x *= scale;
    vertex->y *= scale;
    vertex->z *= scale;
    R2 = SQR(vertex->x) + SQR(vertex->y) + SQR(vertex->z);
    scale = DEFAULT_RADIUS / sqrt(R2);
    vertex->x *= scale;
    vertex->y *= scale;
    vertex->z *= scale;
  }

#if DEBUG_CENTER_SURFACE
  fprintf(WHICH_OUTPUT,
          "\nFinal Minimum found: NRG=%lf, R=%f and ( %f , %f , %f )\n",
          100000.0 * estimateNRG(mris, cx, cy, cz, 10000.0),
          sqrt(estimateSquaredRadius(mris, cx, cy, cz)),
          10000. * cx,
          10000. * cy,
          10000. * cz);
#endif

  fprintf(WHICH_OUTPUT, "done\nSurface centered at (0,0,0) with radius 100.0 in %d iterations\n", niters);
}



void MRISrecenter(MRIS *mris, int which_move, int which_target) 
{
  // MRISrecenter(mris, CURRENT_VERTICES, CANONICAL_VERTICES) ;
  //
  double  xt, yt, zt, tx, ty, tz, xm, ym, zm, radius, r ;
  int     n, vno ;
  VERTEX  *v;

  MRIScomputeMetricProperties(mris) ;
  radius = sqrt(mris->total_area / (4.0*M_PI)) ; // preserve total surface area

  MRISsaveVertexPositions(mris, which_move) ;
  MRISrestoreVertexPositions(mris, which_target) ;
  MRIScenterSphere(mris) ;
  printf("rescaling canonical coordinates to use r=%2.1f\n", radius) ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISvertexCoord2XYZ_double(v, which_target, &xt, &yt, &zt) ;
    r = sqrt(xt*xt + yt*yt + zt*zt) ;
    v->x = radius*(v->x/r) ; v->y = radius*(v->y/r) ; v->z = radius*(v->z/r) ;
  }
  tx = ty = tz = 0 ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISvertexCoord2XYZ_double(v, which_move, &xm, &ym, &zm) ;
    MRISvertexCoord2XYZ_double(v, which_target, &xt, &yt, &zt) ;
    n++ ;
    tx += xt-xm ; ty += yt-ym ; tz += zt-zm ;
  }

  tx /= n ;  ty /= n ; tz /= n ;
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    switch (which_target)
    {
    case CANONICAL_VERTICES:
      v->cx += tx ; v->cy += ty ; v->cz += tz ;
      break ;
    default:
      ErrorExit(ERROR_UNSUPPORTED, "MRISrecenter: movable vertices %d not supported", which_move);
      break ;
    }
  }
  MRISrestoreVertexPositions(mris, which_move) ;
  MRIScomputeMetricProperties(mris) ;
}


// Distances
//
//      Note:   dist_orig is the distance between the origxyz of the vertices, 
//              which may not be the same as the distance between xyz when they were copied to origxyz, because 
//                  origxyz may have changed
//                  dist    might not have been calculated then
//                  nsize   might have changed since then
//
int MRISclearOrigDistances(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    int n;
    for (n = 0; n < vt->vtotal; n++) {
      v->dist_orig[n] = 0;
    }
  }
  return (NO_ERROR);
}


int MRISscaleDistances(MRIS *mris, float scale)
{
  // This operation seems wrong, since this is a derived metric
  // but there are two users...
  //
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if(v->ripflag) continue;
    int vsize = mrisVertexVSize(mris, vno);
    int n;
    for (n = 0; n < vsize; n++) {
      v->dist[n] *= scale;
    }
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------------------
  Calculate distances between each vertex and each of its neighbors.
  ----------------------------------------------------------------*/
static bool mrisComputeVertexDistancesWkr        (MRIS *mris, int new_dist_nsize, bool check);
static bool mrisComputeOriginalVertexDistancesWkr(MRIS *mris, int new_dist_nsize, bool check);

bool mrisCheckDist(MRIS const * mris) {
  if (!mris->dist_nsize) return true;
  return mrisComputeVertexDistancesWkr((MRIS*)mris, mris->dist_nsize, true);
}

bool mrisCheckDistOrig(MRIS const * mris) {
  if (!mris->dist_orig_nsize) return true;
  return mrisComputeOriginalVertexDistancesWkr((MRIS*)mris, mris->dist_orig_nsize, true);
}

int mrisComputeVertexDistances(MRIS *mris) {
  mrisComputeVertexDistancesWkr(mris, mris->nsize, false);
  mris->dist_nsize = mris->nsize;
  if (false) mrisCheckDist(mris);
  return NO_ERROR;
}

int mrisComputeOriginalVertexDistances(MRIS *mris) {
  mrisComputeOriginalVertexDistancesWkr(mris, mris->nsize, false);
  mris->dist_orig_nsize = mris->nsize;
  if (false) mrisCheckDistOrig(mris);
  return NO_ERROR;
}

void mrisComputeOriginalVertexDistancesIfNecessaryWkr(MRIS *mris, bool* laterTime, const char* file, int line)
{
  if (mris->dist_alloced_flags&2) return;
  
  bool useOldBehaviour = false;
  switch (copeWithLogicProblem("FREESURFER_fix_missing_dist_orig","dist_orig not already computed")) {
  case LogicProblemResponse_old: 
    useOldBehaviour = true;
    break;
  case LogicProblemResponse_fix:
    useOldBehaviour = false;
  }

  // The old code did not compute this distance, but instead had zero's loaded into the already allocated dist_orig
  // Computing it here may change the result, so the default is to zero the values after computing them
  //    thereby allocating the correct size, checking the calculation, but reverting to the old behaviour
  //
  mrisComputeOriginalVertexDistances(mris);

  if (!useOldBehaviour) return;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * v = &mris->vertices[vno];
    float* dist_orig = v->dist_orig;
    bzero(dist_orig, v->dist_orig_capacity*sizeof(float));
  }
}

// The only difference between these should be
// whether they use     xyz and write dist
//                  origxyz and write dist_orig
// but they diverged when these two did not share code.
//
#define FUNCTION_NAME mrisComputeVertexDistancesWkr
#define INPUT_STATUS status
#define INPUT_X x
#define INPUT_Y y
#define INPUT_Z z
#define OUTPUT_DIST dist
#define OUTPUT_MAKER MRISmakeDist
#include "mrisComputeVertexDistancesWkr_extracted.h"

#define FUNCTION_NAME mrisComputeOriginalVertexDistancesWkr
#define INPUT_STATUS origxyz_status
#define INPUT_X origx
#define INPUT_Y origy
#define INPUT_Z origz
#define OUTPUT_DIST dist_orig
#define OUTPUT_MAKER MRISmakeDistOrig
#include "mrisComputeVertexDistancesWkr_extracted.h"


double MRISpercentDistanceError(MRIS *mris)
{
  double dist_scale;
  if (mris->patch) {
    dist_scale = 1.0;
  } else {
    cheapAssert(mris->total_area > 0.0);
    cheapAssert(mris->orig_area  > 0.0);
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }

  double mean_dist  = 0.0;
  double mean_odist = 0.0;
  double mean_error = 0.0;
  double pct        = 0.0;
  int    nnbrs      = 0;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    int n;
    for (n = 0; n < vt->vtotal; n++) {
      double dist  = !v->dist      ? 0.0 : v->dist     [n] * dist_scale;
      double odist = !v->dist_orig ? 0.0 : v->dist_orig[n];
      
      nnbrs++;
      mean_dist  +=  dist;
      mean_odist += odist;
      mean_error += fabs(dist - odist);
      if (!FZERO(odist)) pct += fabs(dist - odist) / odist;
    }
  }

  if (mean_dist == 0.0f)  fs::debug() << "MRISpercentDistanceError called when all dist zero";
  if (mean_odist == 0.0f) fs::debug() << "MRISpercentDistanceError called when all dist_orig zero";
  
  if (nnbrs == 0) nnbrs = 1;
  
  mean_odist /= (double)nnbrs;
  mean_error /= (double)nnbrs;
  if (!FZERO(mean_odist)) {
    pct = mean_error / mean_odist;
  } else {
    pct = 1000.0; /* should never happen */
  }

  return (100.0 * pct);
}


// Area
//
int MRISclearOrigArea(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->origarea = 0;
  }
  return (NO_ERROR);
}

void MRISclearOrigAreaAndVal2(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->origarea = 0;
    v->val2 = 0;
  }
}


float mrisComputeArea(MRIS *mris, int fac, int n)
{
  int n0, n1;
  float v0[3], v1[3], d1, d2, d3, dot, area;

  n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;

  FACE const * const f = &mris->faces[fac];

  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;

  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;

  d1 = -v1[1] * v0[2] + v0[1] * v1[2];
  d2 =  v1[0] * v0[2] - v0[0] * v1[2];
  d3 = -v1[0] * v0[1] + v0[0] * v1[1];

  dot = mris->vertices[f->v[n]].x * d1 + mris->vertices[f->v[n]].y * d2 + mris->vertices[f->v[n]].z * d3;

  if (dot < 0.0f) /* not in same direction, area < 0 and reverse n */
  {
    area = -sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  }
  else {
    area =  sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  }

  return area;
}


float MRIScomputeOrigArea(MRIS* mris)
{
  float orig_area = 0.0f;
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE const * const f = &mris->faces[fno];
    if (f->ripflag) continue;
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
    orig_area += fNorm->orig_area;
  }
  return orig_area;
}


void MRISsetOrigArea(MRIS* mris) 
{
  mris->orig_area = MRIScomputeOrigArea(mris);
}


int MRISzeroNegativeAreas(MRIS *mris)
{
  int vno, fno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->area < 0) {
      v->area = 0;
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    if (face->area < 0.0f) {
      face->area = 0.0f;
    }
  }
  return (NO_ERROR);
}


static int count_MRIScomputeMetricProperties_calls = 0;

#include "mrisurf_metricProperties_slow.h"
#include "mrisurf_metricProperties_fast.h"

int MRIScomputeMetricProperties(MRIS *mris)
{
  if (!!getenv("FS_FASTER_MP")) {
    // this method has been shown to be 3-8x faster than
    // the "mrisurf_metricProperties_fast" implementation
    MRIScomputeMetricPropertiesFaster(mris);
    return NO_ERROR;
  }

  static int debug_count = 0;
  count_MRIScomputeMetricProperties_calls++;
  if (count_MRIScomputeMetricProperties_calls == debug_count) {
    fprintf(stdout, "b %s:%d   needed with debug_count:%d\n", __FILE__, __LINE__, count_MRIScomputeMetricProperties_calls);
  }

  bool const useOldBehaviour = !!getenv("FREESURFER_OLD_MRIScomputeMetricProperties");
  bool const useNewBehaviour = !!getenv("FREESURFER_NEW_MRIScomputeMetricProperties") || !useOldBehaviour;

  MRIS_MP mp = MRIS_MP();
  MRISMP_ctr(&mp);

  if (useNewBehaviour) {
    MRISfreeDistsButNotOrig(mris);     // So they can be stolen to avoid unnecessary mallocs and frees
    MRISMP_load(&mp, mris);            // Copy the input data before MRIScomputeMetricPropertiesWkr changes it
    MRIScomputeMetricProperties(&mp);  // It should not matter the order these are done in
  }
  
  if (useOldBehaviour) {
    MRIScomputeMetricPropertiesWkr(mris);
  }
  
  // Verify the two approaches got the same answers, or use the new answers 
  if (useNewBehaviour) MRISMP_unload(mris, &mp, useOldBehaviour);

  MRISMP_dtr(&mp);

  return NO_ERROR;
}

// Convenience functions
//
int load_orig_triangle_vertices(MRIS *mris, int fno, double U0[3], double U1[3], double U2[3])
{
  FACE const * const face = &mris->faces[fno];
  
  VERTEX *v;
  
  v = &mris->vertices[face->v[0]];
  U0[0] = v->origx;
  U0[1] = v->origy;
  U0[2] = v->origz;
  v = &mris->vertices[face->v[1]];
  U1[0] = v->origx;
  U1[1] = v->origy;
  U1[2] = v->origz;
  v = &mris->vertices[face->v[2]];
  U2[0] = v->origx;
  U2[1] = v->origy;
  U2[2] = v->origz;
  return (NO_ERROR);
}

int load_triangle_vertices(MRIS *mris, int fno, double U0[3], double U1[3], double U2[3], int which)
{
  VERTEX *v;
  FACE *face;

  face = &mris->faces[fno];
  switch (which) {
    default:
    case CURRENT_VERTICES:
      v = &mris->vertices[face->v[0]];
      U0[0] = v->x;
      U0[1] = v->y;
      U0[2] = v->z;
      v = &mris->vertices[face->v[1]];
      U1[0] = v->x;
      U1[1] = v->y;
      U1[2] = v->z;
      v = &mris->vertices[face->v[2]];
      U2[0] = v->x;
      U2[1] = v->y;
      U2[2] = v->z;
      break;
    case WHITE_VERTICES:
      v = &mris->vertices[face->v[0]];
      U0[0] = v->whitex;
      U0[1] = v->whitey;
      U0[2] = v->whitez;
      v = &mris->vertices[face->v[1]];
      U1[0] = v->whitex;
      U1[1] = v->whitey;
      U1[2] = v->whitez;
      v = &mris->vertices[face->v[2]];
      U2[0] = v->whitex;
      U2[1] = v->whitey;
      U2[2] = v->whitez;
      break;
    case PIAL_VERTICES:
      v = &mris->vertices[face->v[0]];
      U0[0] = v->pialx;
      U0[1] = v->pialy;
      U0[2] = v->pialz;
      v = &mris->vertices[face->v[1]];
      U1[0] = v->pialx;
      U1[1] = v->pialy;
      U1[2] = v->pialz;
      v = &mris->vertices[face->v[2]];
      U2[0] = v->pialx;
      U2[1] = v->pialy;
      U2[2] = v->pialz;
      break;
  }
  return (NO_ERROR);
}


int MRISextractVertexCoords(MRIS *mris, float *locations[3], int which)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    switch (which) {
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRISextractVertexCoords: which %d not supported", which);
        break;
      case CURRENT_VERTICES:
        locations[0][vno] = v->x;
        locations[1][vno] = v->y;
        locations[2][vno] = v->z;
        break;
      case TARGET_VERTICES:
        locations[0][vno] = v->targx;
        locations[1][vno] = v->targy;
        locations[2][vno] = v->targz;
        break;
      case WHITE_VERTICES:
        locations[0][vno] = v->whitex;
        locations[1][vno] = v->whitey;
        locations[2][vno] = v->whitez;
        break;
      case LAYERIV_VERTICES:
        locations[0][vno] = v->l4x;
        locations[1][vno] = v->l4y;
        locations[2][vno] = v->l4z;
        break;
      case PIAL_VERTICES:
        locations[0][vno] = v->pialx;
        locations[1][vno] = v->pialy;
        locations[2][vno] = v->pialz;
        break;
      case INFLATED_VERTICES:
        locations[0][vno] = v->infx;
        locations[1][vno] = v->infy;
        locations[2][vno] = v->infz;
        break;
      case FLATTENED_VERTICES:
        locations[0][vno] = v->fx;
        locations[1][vno] = v->fy;
        locations[2][vno] = v->fz;
        break;
      case CANONICAL_VERTICES:
        locations[0][vno] = v->cx;
        locations[1][vno] = v->cy;
        locations[2][vno] = v->cz;
        break;
      case ORIGINAL_VERTICES:
        locations[0][vno] = v->origx;
        locations[1][vno] = v->origy;
        locations[2][vno] = v->origz;
        break;
      case TMP2_VERTICES:
        locations[0][vno] = v->t2x;
        locations[1][vno] = v->t2y;
        locations[2][vno] = v->t2z;
        break;
      case TMP_VERTICES:
        locations[0][vno] = v->tx;
        locations[1][vno] = v->ty;
        locations[2][vno] = v->tz;
        break;
    }
  }
  return (NO_ERROR);
}


int MRISscaleBrainArea(MRIS *mris)
{
  float const scale = sqrt(mris->orig_area / (mris->total_area + mris->neg_area));
  
  MRISscaleBrain(mris, mris, scale);
  
  MRIScomputeMetricProperties(mris);
  
  return (NO_ERROR);
}


/*!
  \fn int MRISltaMultiply(MRIS *surf, LTA *lta)
  \brief Applies an LTA matrix to the coords of the given surface.
  Automatically determines which direction the LTA goes by looking
  at the volume geometries of the LTA and the surface. The vol
  geometry of the surface is changed to that of the LTA destination
  (keeping in mind that the LTA might have been reversed). The 
  LTA itself is not changed. 
  See also:   MRISmatrixMultiply() and MRIStransform().
*/
int MRISltaMultiply(MRIS *surf, const LTA *lta)
{
  LTA *ltacopy;
  extern double vg_isEqual_Threshold;
  vg_isEqual_Threshold = 10e-4;

  // Make a copy of the LTA so source is not contaminated
  ltacopy = LTAcopy(lta,NULL);

  // Check whether the source and destination geometries are the same
  // since this will make it impossible to determine the direction
  // from the LTA. 
  if(vg_isEqual(&ltacopy->xforms[0].src, &ltacopy->xforms[0].dst)){
    // If they are the same, check whether the registration is the identity
    // in which case the direction is not important.
    int c,r,IsIdentity=1;
    double val;
    for(r=1; r<=4; r++){
      for(c=1; c<=4; c++){
        val = ltacopy->xforms[0].m_L->rptr[r][c];
        if(r==c && fabs(val-1.0) > 10e-4) IsIdentity = 0;
        if(r!=c && fabs(val)     > 10e-4) IsIdentity = 0;
      }
    }
    if(! IsIdentity){
      // Only print out a warning if if they are the same and the reg is not identity.
      printf("\nINFO: MRISltaMultiply(): LTA src and dst vg's are the same and reg is not identity.\n");
      printf("  Make sure you have the direction correct!\n\n");
    }
  }

  // Determine which direction the LTA goes by looking at which side
  // matches the surface volume geometry. Invert if necessary
  if(surf->vg.valid){
    if(!vg_isEqual(&ltacopy->xforms[0].src, &surf->vg)){
      if(!vg_isEqual(&ltacopy->xforms[0].dst, &surf->vg)){
	// If this fails a lot, try setting vg_isEqual_Threshold higher
	printf("vg surf ------------------------\n");
	fflush(stdout);	fflush(stderr);
	vg_print(&surf->vg);
	fflush(stdout);	fflush(stderr);
	printf("LTA ------------------------\n");
	LTAprint(stdout,lta);
	printf("------------------------\n");
	printf("ERROR: MRISltaMultiply(): LTA does not match surface vg\n");
	printf("   src diff no %d\n",vg_isNotEqualThresh(&ltacopy->xforms[0].src, &surf->vg,vg_isEqual_Threshold));
	printf("   dst diff no %d\n",vg_isNotEqualThresh(&ltacopy->xforms[0].dst, &surf->vg,vg_isEqual_Threshold));
	printf("   vg_isEqual_Threshold = %10.6e\n",vg_isEqual_Threshold);
	LTAfree(&ltacopy);
	return(1);
      }
      printf("MRISltaMultiply(): inverting LTA\n");
      LTAinvert(ltacopy,ltacopy);
    }
  }
  else {
    printf("WARNING: source surface vg is not valid, assuming direction of transform is right\n");
  }

  // Must be in regdat space to apply it to the surface
  if(ltacopy->type != REGISTER_DAT){
    printf("Changing type to REGISTER_DAT\n");
    LTAchangeType(ltacopy, REGISTER_DAT);
  }

  // In REGISTER_DAT format, the matrix actually goes from target/dst
  // to mov/src so have to invert the matrix.  Note can't use
  // LTAinvert() because it will also reverse the src and dst vol
  // geometries.
  MatrixInverse(ltacopy->xforms[0].m_L,ltacopy->xforms[0].m_L);

  // Print out the matrix that will be applied
  printf("MRISltaMultiply(): applying matrix\n");
  MatrixPrint(stdout,ltacopy->xforms[0].m_L);
  printf("-----------------------------------\n");

  // Finally, multiply the surf coords by the matrix
  MRISmatrixMultiply(surf,ltacopy->xforms[0].m_L);

  // Reverse the faces if the reg has neg determinant
  if(MatrixDeterminant(ltacopy->xforms[0].m_L) < 0){
    printf("Determinant of linear transform is negative, so reversing face order\n");
    MRISreverseFaceOrder(surf);
  }

  // Copy the volume geometry of the destination volume
  copyVolGeom(&(ltacopy->xforms[0].dst), &(surf->vg));

  LTAfree(&ltacopy);
  return(0);
}



/*-----------------------------------------------------
  MRIStransform

  Parameters: mris, mri, lta, mri_dst

  Returns value: int

  Description
  Apply a linear transform (possibly octree) to a surface.
  Note that the LT is given in MRI coordinates, so the
  surface must be transformed into that coordinate system
  before applying the linear transform
  See also MRISmatrixMultiply().
  ------------------------------------------------------*/
#include "gcamorph.h"
int MRIStransform(MRIS *mris, MRI *mri, TRANSFORM *transform, MRI *mri_dst)
{
  LTA *lta;
  int vno;
  VERTEX *v;
  double xw, yw, zw;
  MATRIX *m = 0;

  // for ras-to-ras transform
  MATRIX *RASFromSurfaceRAS = 0;
  MATRIX *surfaceRASFromRAS = 0;
  // for voxel-to-voxel transform
  MATRIX *voxelFromSurfaceRAS = 0;
  MATRIX *surfaceRASFromVoxel = 0;
  //
  MATRIX *surfaceRASFromSurfaceRAS = 0;
  LT *lt = 0;
  int srcPresent = 1;
  int dstPresent = 1;
  int error = NO_ERROR;
  char errMsg[256];
  int dstNotGiven = 0;
  int srcNotGiven = 0;

  if (transform->type == MORPH_3D_TYPE) {
    GCA_MORPH *gcam;
    double xs, ys, zs, xv, yv, zv;
    float xv2, yv2, zv2;
    VECTOR *v1, *v2;

    /*
      transform point from surface (ras) coords to the source (atlas) volume
      that was used for the gcam:
    */

    gcam = (GCA_MORPH *)(transform->xform);
    if (!mri_dst) {
      dstNotGiven = 1;
      mri_dst = MRIalloc(gcam->image.width, gcam->image.height, gcam->image.depth, MRI_UCHAR);
      useVolGeomToMRI(&gcam->image, mri_dst);
    }
    if (!mri) {
      srcNotGiven = 1;
      mri = MRIalloc(gcam->atlas.width, gcam->atlas.height, gcam->atlas.depth, MRI_UCHAR);
      useVolGeomToMRI(&gcam->atlas, mri);
    }
    if (gcam->type == GCAM_RAS) {
        GCAMrasToVox(gcam, mri_dst);
    }

    v1 = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(v1, 4) = 1.0;
    v2 = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(v2, 4) = 1.0;
    voxelFromSurfaceRAS = voxelFromSurfaceRAS_(mri);

    // now apply the transform
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }

      // transform vertex point to actual voxel point
      V3_X(v1) = v->x;
      V3_Y(v1) = v->y;
      V3_Z(v1) = v->z;
      MatrixMultiply(voxelFromSurfaceRAS, v1, v2);
      xv = V3_X(v2);
      yv = V3_Y(v2);
      zv = V3_Z(v2);
      GCAMsampleMorph(gcam, xv, yv, zv, &xv2, &yv2, &zv2);
      MRIvoxelToSurfaceRAS(mri_dst, xv2, yv2, zv2, &xs, &ys, &zs);
      v->x = xs;
      v->y = ys;
      v->z = zs;
    }
    mrisComputeSurfaceDimensions(mris);
    // save the volume information from dst
    getVolGeom(mri_dst, &mris->vg);
    VectorFree(&v1);
    VectorFree(&v2);
    MatrixFree(&voxelFromSurfaceRAS);
    if (dstNotGiven) {
      MRIfree(&mri_dst);
    }
    if (srcNotGiven) {
      MRIfree(&mri);
    }
  }
  else {
    lta = (LTA *)(transform->xform);
    if (!mri_dst) {
      dstNotGiven = 1;
    }

    // depend on the type of transform you have to handle differently
    //
    //       orig              ------>      RAS   c_(ras) != 0
    //        |                              |
    //        |                              | identity
    //        V                              V
    //    conformed vol        ------>      RAS   c_(ras) != 0
    //        |                              |
    //        | identity                     | surfaceRASFromRAS
    //        V                              V
    //    conformed vol        ------>   surfaceRAS  c_(ras) = 0
    //
    // given a volume transform you have to create a surfaceRAS transform
    //
    // Note that vertices are given by surfaceRAS coordinates
    //
    // RAS-to-RAS transform
    //                  orig                 dst
    //    surfaceRAS--->RAS --(ras-to-ras)-->RAS -->surfaceRAS
    //
    // VOX-to-Vox transform
    //
    //    surfaceRAS--->Vox---(vox-to-vox)-->Vox -->surfaceRAS
    //
    //
    if (lta->num_xforms > 1) {
      ErrorExit(ERROR_BADPARM, "we cannot handle multiple transforms\n");
    }
    if (lta->num_xforms == 0) {
      ErrorExit(ERROR_BADPARM, "transform does not have transform ;-) \n");
    }

    // if volumes are not given, then try to get them from transform
    lt = &lta->xforms[0];

    /////////////////////////////////////////////////////////////////////
    // The following rather fragile treatment about volume info is
    // that lta may or may not store the valid src and dst info
    // when the transform was created
    // Another problem is due to the convention of conformed volume
    ////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    // check the c_ras values for mri side
    if (lta->type == LINEAR_RAS_TO_RAS) {
      if (mri && lt->src.valid == 1) {
        if (!(FZERO(lt->src.c_r - lt->src.c_r) && FZERO(lt->src.c_a - lt->src.c_a) &&
              FZERO(lt->src.c_s - lt->src.c_s))) {
          fprintf(stderr, "WARNING:****************************************************\n");
          fprintf(stderr,
                  "WARNING: c_(ras) values are not equal for the "
                  "input volume and the transform src.\n");
          fprintf(stderr, "WARNING: The transformed surface position may be shifted.\n");
          fprintf(stderr, "WARNING:***************************************************\n");
        }
      }
      if (mri && mris->vg.valid == 1) {
        if (!(FZERO(mris->vg.c_r - mris->vg.c_r) && FZERO(mris->vg.c_a - mris->vg.c_a) &&
              FZERO(mris->vg.c_s - mris->vg.c_s))) {
          fprintf(stderr,
                  "WARNING:***************************"
                  "******************************\n");
          fprintf(stderr,
                  "WARNING: c_(ras) values are not equal for the input volume "
                  "and the surface stored volume.\n");
          fprintf(stderr, "WARNING: The transformed surface position may be shifted.\n");
          fprintf(stderr,
                  "WARNING:*******************************"
                  "**************************\n");
        }
      }
    }

    // if mri is not given, then use the one stored in the transform
    if (!mri && lt->src.valid == 1) {
      srcPresent = 0;
      fprintf(stderr, "INFO:try to get src info from transform.\n");
      mri = MRIallocHeader(lt->src.width, lt->src.height, lt->src.depth, MRI_UCHAR, 1);
      mri->x_r = lt->src.x_r;
      mri->y_r = lt->src.y_r;
      mri->z_r = lt->src.z_r;
      mri->c_r = lt->src.c_r;
      mri->x_a = lt->src.x_a;
      mri->y_a = lt->src.y_a;
      mri->z_a = lt->src.z_a;
      mri->c_a = lt->src.c_a;
      mri->x_s = lt->src.x_s;
      mri->y_s = lt->src.y_s;
      mri->z_s = lt->src.z_s;
      mri->c_s = lt->src.c_s;
      mri->xsize = lt->src.xsize;
      mri->ysize = lt->src.ysize;
      mri->zsize = lt->src.zsize;
      mri->ras_good_flag = 1;
      MRIreInitCache(mri);
    }
    // if mri is not given and transform does not
    // have it, get it from the surface
    else if (!mri && mris->vg.valid == 1) {
      fprintf(stderr, "INFO:try to get src info from the surface.\n");
      mri = MRIallocHeader(mris->vg.width, mris->vg.height, mris->vg.depth, MRI_UCHAR, 1);
      mri->x_r = mris->vg.x_r;
      mri->y_r = mris->vg.y_r;
      mri->z_r = mris->vg.z_r;
      mri->c_r = mris->vg.c_r;
      mri->x_a = mris->vg.x_a;
      mri->y_a = mris->vg.y_a;
      mri->z_a = mris->vg.z_a;
      mri->c_a = mris->vg.c_a;
      mri->x_s = mris->vg.x_s;
      mri->y_s = mris->vg.y_s;
      mri->z_s = mris->vg.z_s;
      mri->c_s = mris->vg.c_s;
      mri->xsize = mris->vg.xsize;
      mri->ysize = mris->vg.ysize;
      mri->zsize = mris->vg.zsize;
      mri->ras_good_flag = 1;
      MRIreInitCache(mri);
    }
    else if (!mri) {
      error = 1;
      strcpy(errMsg, "When mri == NULL, the transform must have the valid src info.\n");
      goto mristransform_cleanup;
    }
    ///////////////////////////////////////////////////////////////////////////
    // mri_dst side
    // Note: if mri_dst is not given, override the one stored in the transform
    if (!mri_dst && lt->dst.valid == 1) {
      dstPresent = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(stderr, "INFO: Get dst info from transform.\n");
      }
      lt = &lta->xforms[0];
      mri_dst = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR, 1);
      mri_dst->x_r = lt->dst.x_r;
      mri_dst->y_r = lt->dst.y_r;
      mri_dst->z_r = lt->dst.z_r;
      mri_dst->c_r = lt->dst.c_r;
      mri_dst->x_a = lt->dst.x_a;
      mri_dst->y_a = lt->dst.y_a;
      mri_dst->z_a = lt->dst.z_a;
      mri_dst->c_a = lt->dst.c_a;
      mri_dst->x_s = lt->dst.x_s;
      mri_dst->y_s = lt->dst.y_s;
      mri_dst->z_s = lt->dst.z_s;
      mri_dst->c_s = lt->dst.c_s;
      mri_dst->xsize = lt->dst.xsize;
      mri_dst->ysize = lt->dst.ysize;
      mri_dst->zsize = lt->dst.zsize;
      mri_dst->ras_good_flag = 1;
    }
    else if (!mri_dst) {
      fprintf(stderr, "WARNING:*********************************************************\n");
      fprintf(stderr, "WARNING: transform does not have valid destination volume.       \n");
      fprintf(stderr, "WARNING: The standard CORONAL volume with c_(ras) = 0 is assumed.\n");
      fprintf(stderr, "WARNING:*********************************************************\n");
      mri_dst = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR, 1);
      mri_dst->x_r = -1;
      mri_dst->y_r = 0;
      mri_dst->z_r = 0;
      mri_dst->c_r = 0;
      mri_dst->x_a = 0;
      mri_dst->y_a = 0;
      mri_dst->z_a = 1;
      mri_dst->c_a = 0;
      mri_dst->x_s = 0;
      mri_dst->y_s = -1;
      mri_dst->z_s = 0;
      mri_dst->c_s = 0;
      mri_dst->xsize = 1;
      mri_dst->ysize = 1;
      mri_dst->zsize = 1;
      mri_dst->ras_good_flag = 1;
    }
    // WATCH /////////////////////////////////////////////////////////
    // When dst is not given, our convention is that the
    // dst is conformed and thus we need
    // to modify mri_dst to be coronal conformed volume
    // (even though transform target was not,
    // as in the case of MNI average_305 being non-coronal, non-conformed
    //////////////////////////////////////////////////////////////////
    if (dstNotGiven) {
      // assuming mri is conformed
      mri_dst->width = mri->width;
      mri_dst->height = mri->height;
      mri_dst->depth = mri->depth;
      mri_dst->xsize = mri->xsize;
      mri_dst->ysize = mri->ysize;
      mri_dst->zsize = mri->zsize;
      mri_dst->x_r = -1;
      mri_dst->y_r = 0;
      mri_dst->z_r = 0;
      mri_dst->x_a = 0;
      mri_dst->y_a = 0;
      mri_dst->z_a = 1;
      mri_dst->x_s = 0;
      mri_dst->y_s = -1;
      mri_dst->z_s = 0;
      // this means that we only retain c_ras info
    }
    // you must reinitialise cache
    MRIreInitCache(mri_dst);

    //  MatrixPrint(stderr, mri_dst->i_to_r__);

    ////////////////////////////////////////////////////////////////
    // Now we can calculate
    if (lta->type == LINEAR_RAS_TO_RAS) {
      //  we follow the right hand side of the map
      //
      //    conformed -----> surfaceRAS (c_ras = 0)
      //       |              |
      //       V              V
      //    conformed -----> RAS (c_ras != 0)
      //       |              |
      //       |             xfm
      //       V              v
      //  conformed dst ---> RAS (c_ras != 0)
      //       |              |
      //       V              V
      //  conformed dst ---> surfaceRAS (c_ras = 0)
      //
      if (mris->useRealRAS) {
        surfaceRASFromSurfaceRAS = MatrixCopy(lta->xforms[0].m_L, NULL);
      }
      else {
        RASFromSurfaceRAS = RASFromSurfaceRAS_(mri,NULL);      // needs only c_(ras) info
        surfaceRASFromRAS = surfaceRASFromRAS_(mri_dst);  // need only c_(ras) info
        m = MatrixMultiply(lta->xforms[0].m_L, RASFromSurfaceRAS, NULL);
        surfaceRASFromSurfaceRAS = MatrixMultiply(surfaceRASFromRAS, m, NULL);
      }
    }
    else if (lta->type == LINEAR_VOX_TO_VOX) {
      if (mri->width != mri_dst->width || mri->height != mri_dst->height || mri->depth != mri_dst->depth) {
        fprintf(stderr, "WARNING:********************************************************\n");
        fprintf(stderr, "WARNING:voxel-to-voxel transform must have the same volume sizes.\n");
        fprintf(stderr,
                "WARNING:You gave src (%d, %dm, %d) vs. dst (%d, %d, %d).\n",
                mri->width,
                mri->height,
                mri->depth,
                mri_dst->width,
                mri_dst->height,
                mri_dst->depth);
        fprintf(stderr, "WARNING:The result of this transform is most likely wrong.\n");
        fprintf(stderr, "WARNING:********************************************************\n");
      }
      voxelFromSurfaceRAS = voxelFromSurfaceRAS_(mri);
      surfaceRASFromVoxel = surfaceRASFromVoxel_(mri_dst);
      //    MatrixPrint(stderr,  surfaceRASFromVoxel);
      m = MatrixMultiply(lta->xforms[0].m_L, voxelFromSurfaceRAS, NULL);
      surfaceRASFromSurfaceRAS = MatrixMultiply(surfaceRASFromVoxel, m, NULL);
    }
    // now apply the transform
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      // transform vertex point to actual voxel point
      TransformWithMatrix(surfaceRASFromSurfaceRAS, v->x, v->y, v->z, &xw, &yw, &zw);
      v->x = xw;
      v->y = yw;
      v->z = zw;
    }
    mrisComputeSurfaceDimensions(mris);
    // save the volume information from dst
    getVolGeom(mri_dst, &mris->vg);

  mristransform_cleanup:
    // free memory //////////////////////////////
    if (RASFromSurfaceRAS) {
      MatrixFree(&RASFromSurfaceRAS);
    }
    if (surfaceRASFromRAS) {
      MatrixFree(&surfaceRASFromRAS);
    }
    if (m) {
      MatrixFree(&m);
    }
    if (voxelFromSurfaceRAS) {
      MatrixFree(&voxelFromSurfaceRAS);
    }
    if (surfaceRASFromVoxel) {
      MatrixFree(&surfaceRASFromVoxel);
    }
    if (surfaceRASFromSurfaceRAS) {
      MatrixFree(&surfaceRASFromSurfaceRAS);
    }
    if (!srcPresent && mri) {
      MRIfree(&mri);
    }
    if (!dstPresent && mri_dst) {
      MRIfree(&mri_dst);
    }
  }

  if (error) {
    ErrorExit(ERROR_BADPARM, errMsg);
    return -1;  // just to satisfy compiler
  }
  else {
    return (NO_ERROR);
  }
}


/*! ----------------------------------------------------
  \fn int MRISclearCurvature(MRIS *mris)
  \brief sets v->curv=0 for unripped vertices.
  ------------------------------------------------------*/
int MRISclearCurvature(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = 0;
  }
  return (NO_ERROR);
}
void MRISclearCurvAndVal2(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX* v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = 0;
    v->val2 = 0;
  }
}


int MRISrectifyCurvature(MRIS *mris)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->curv = fabs(v->curv);
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}


int MRISnonmaxSuppress(MRIS *mris)
{
  double du, dv, up1, um1, vp1, vm1, src, dx, dy, dz, fp1, fm1, mag;
  VERTEX *v;
  int vno;
  float x, y, z, e1x, e1y, e1z, e2x, e2y, e2z, ux, uy, uz, vx, vy, vz;
  MRI_SP *mrisp, *mrisp_blur;
  double d_dist = D_DIST * mris->avg_vertex_dist;

  mrisp = MRIStoParameterization(mris, NULL, 1, 0);
  mrisp_blur = MRISPblur(mrisp, NULL, 20.0, 0);
  mrisComputeTangentPlanes(mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;
    x = v->x;
    y = v->y;
    z = v->z;
    src = MRISPfunctionVal(mrisp, mris->radius, x, y, z, 0);

    /* now compute gradient of template w.r.t. a change in vertex position */

    /*
      sample the curvature functions along the tangent plane axes and
      compute the derivates using them.
    */
    ux = 10 * e1x;
    uy = 10 * e1y;
    uz = 10 * e1z;
    vx = 10 * e2x;
    vy = 10 * e2y;
    vz = 10 * e2z;
    ux = e1x * d_dist;
    uy = e1y * d_dist;
    uz = e1z * d_dist;
    vx = e2x * d_dist;
    vy = e2y * d_dist;
    vz = e2z * d_dist;

    /* compute gradient using blurred image */
    up1 = MRISPfunctionVal(mrisp_blur, mris->radius, x + ux, y + uy, z + uz, 0);
    um1 = MRISPfunctionVal(mrisp_blur, mris->radius, x - ux, y - uy, z - uz, 0);
    vp1 = MRISPfunctionVal(mrisp_blur, mris->radius, x + vx, y + vy, z + vz, 0);
    vm1 = MRISPfunctionVal(mrisp_blur, mris->radius, x - vx, y - vy, z - vz, 0);
    du = (up1 - um1) / (2 * d_dist);
    dv = (vp1 - vm1) / (2 * d_dist);

    /* calculate curvature gradient */
    dx = (du * e1x + dv * e2x);
    dy = (du * e1y + dv * e2y);
    dz = (du * e1z + dv * e2z);
    mag = sqrt(dx * dx + dy * dy + dz * dz);
    if (FZERO(mag)) /* zero gradient */
    {
      v->curv = 0;
    }
    else {
      mag *= .1;
      dx = dx / mag;
      dy = dy / mag;
      dz = dz / mag;
    }
    fp1 = MRISPfunctionVal(mrisp, mris->radius, x + dx, y + dy, z + dz, 0);
    fm1 = MRISPfunctionVal(mrisp, mris->radius, x - dx, y - dy, z - dz, 0);

    if ((src >= fp1) && (src >= fm1)) /* local max */
    {
      v->curv = 1;
    }
    else if ((src <= fp1) && (src <= fm1)) /* local min */
    {
      v->curv = -1;
    }
    else {
      v->curv = 0;
    }
    v->curv = src;
  }
  MRISPfree(&mrisp);
  return (NO_ERROR);
}


int MRISscaleCurvatures(MRIS *mris, float min_curv, float max_curv)
{
  double old_min_curv, old_max_curv, mean, scale;
  int vno, vtotal;
  VERTEX *v;

  old_min_curv = 100000.0;
  old_max_curv = -100000.0;
  for (mean = 0.0, vtotal = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    vtotal++;
    mean += v->curv;
    if (v->curv > old_max_curv) {
      old_max_curv = v->curv;
    }
    if (v->curv < old_min_curv) {
      old_min_curv = v->curv;
    }
  }

  mean /= (double)vtotal;
  scale = (max_curv - min_curv) / (old_max_curv - old_min_curv);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = (v->curv - mean) * scale + mean;
  }

  return (NO_ERROR);
}


#define K_A 0.4f
static float kernel[] = {K_A, 0.25f, 0.25f - K_A / 2.0f};
int MRISsmoothCurvatures(MRIS *mris, int niterations)
{
  int vno, i, vn;
  double g, H, norm;

  for (i = 0; i < niterations; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
      VERTEX                * const vertex  = &mris->vertices         [vno];
      if (vertex->ripflag) {
        continue;
      }
      H = kernel[0] * vertex->H;
      g = kernel[1];
      for (vn = 0; vn < vertext->vnum; vn++) {
        H += g * mris->vertices[vertext->v[vn]].H;
      }
      g = kernel[2];
      for (; vn < vertext->v2num; vn++) {
        H += g * mris->vertices[vertext->v[vn]].H;
      }
      norm = kernel[0] + vertext->vnum * kernel[1] + (vertext->v2num - vertext->vnum) * kernel[2];
      vertex->d = H / norm;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const vertex = &mris->vertices[vno];
      if (vertex->ripflag) {
        continue;
      }
      vertex->H = vertex->d;
    }
  }

  return (NO_ERROR);
}


int MRISzeroMeanCurvature(MRIS *mris)
{
  double mean;
  int vno, vtotal;
  VERTEX *v;

  for (mean = 0.0f, vtotal = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    mean += v->curv;
    vtotal++;
  }

  mean /= (double)vtotal;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "curvature mean = %2.3f\n", mean);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv -= mean;
  }

  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}
/*-----------------------------------------------------
int MRISnormalizeCurvature(MRIS *mris, int which_norm)
Normalize the curvatures so they have unit standard deviation.
which_norm can be NORM_MEDIAN or NORM_MEAN. NORM_MEAN results
in the usual stddev computation. NORM_MEDIAN uses the median
instead of the mean for the stddev computation.
  ------------------------------------------------------*/
int MRISnormalizeCurvature(MRIS *mris, int which_norm)
{
  double mean, var, std, median;
  int vno, vtotal;
  VERTEX *v;

  if (which_norm == NORM_NONE) {
    return (NO_ERROR);
  }

  if (which_norm == NORM_MEDIAN) {
    float *curvs;
    curvs = (float *)calloc(mris->nvertices, sizeof(float));
    if (curvs == NULL) {
      ErrorExit(ERROR_NOMEMORY, "MRISnormalize(MEDIAN): couldn't alloc array");
    }

    MRISextractCurvatureVector(mris, curvs);
    qsort(curvs, mris->nvertices, sizeof(curvs[0]), mris_sort_compare_float);
    median = curvs[mris->nvertices / 2];
    free(curvs);
    for (var = 0.0f, vtotal = vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      std = (v->curv - median);
      vtotal++;
      var += std * std;
    }

    var /= (double)vtotal;
    std = sqrt(var);

    //    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "curvature median = %2.3f, std = %2.3f\n", median, std);

    // now normalize the curvatures so they have unit standard deviation
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = (v->curv - median) / std;
    }
  }
  else  // normalize mean
  {
    for (mean = 0.0f, vtotal = vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      mean += v->curv;
      vtotal++;
    }

    mean /= (double)vtotal;

    for (var = 0.0f, vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      std = v->curv - mean;
      var += std * std;
    }

    var /= (double)vtotal;
    std = sqrt(var);

    //    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "curvature mean = %2.3f, std = %2.3f\n", mean, std);

    /* now normalize the curvatures so they have unit standard deviation, but
       leave the mean alone */
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->curv = (v->curv - mean) / std /* + mean*/;
    }
  }
  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISnormalizeCurvatureVariance(MRIS *mris)
{
  double mean, var, std;
  int vno, vtotal;
  VERTEX *v;

  for (mean = 0.0f, vtotal = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    mean += v->curv;
    vtotal++;
  }

  mean /= (double)vtotal;

  for (var = 0.0f, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    std = v->curv - mean;
    var += std * std;
  }

  var /= (double)vtotal;
  std = sqrt(var);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "curvature mean = %2.3f, std = %2.3f\n", mean, std);
  }

  /* now normalize the curvatures so they have unit standard deviation, but
     leave the mean alone */
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = (v->curv - mean) / std + mean;
  }

  mrisComputeCurvatureMinMax(mris);
  return (NO_ERROR);
}



int mrisCountCompressed(MRIS *mris, double min_dist)
{
  int vno, n, num;
  double d, dx, dy, dz;

  for (vno = num = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      dx = vn->x - v->x;
      dy = vn->y - v->y;
      dz = vn->z - v->z;
      d = sqrt(dx * dx + dy * dy + dz * dz);
      if (d < min_dist) {
        num++;
      }
    }
  }
  return (num);
}


//=========================================================================
// Change the positions of the vertices
//=========================================================================

int MRISrestoreOldPositions(MRIS *mris)
{
  return (MRISrestoreVertexPositions(mris, TMP_VERTICES));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISstoreCurrentPositions(MRIS *mris)
{
  return (MRISsaveVertexPositions(mris, TMP_VERTICES));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISstoreMeanCurvature(MRIS *mris)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->H = v->curv;
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISstoreMetricProperties(MRIS *mris)
{
  int vno, nvertices, fno, n;
  FACE *f;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    v->origarea = v->area;
    if (v->dist) {
      
      if (!v->dist_orig)
        MRISmakeDistOrig(mris, vno);
       
      // Used to only go to vtotal, but that is v[nsizeCur]num, and the code can go to to v[nsizeMax]num
      int const vsize = mrisVertexVSize(mris,vno);
      for (n = 0; n < vsize; n++) {
        v->dist_orig[n] = v->dist[n];
      }
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    setFaceOrigArea(mris, fno, f->area);
    copyAnglesPerTriangle(f->orig_angle,f->angle);
  }
  mris->orig_area = mris->total_area;
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISrestoreMetricProperties(MRIS *mris)
{
  int vno, nvertices, fno, n;
  FACE *f;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    v->area = v->origarea;
    int const vsize = mrisVertexVSize(mris,vno);
    for (n = 0; n < vsize; n++) {
      v->dist[n] = v->dist_orig[n];
    }
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    f->area = getFaceOrigArea(mris, fno);
    copyAnglesPerTriangle(f->angle,f->orig_angle);
  }
  mris->total_area = mris->orig_area;
  return (NO_ERROR);
}


// Two ways of saving and restoring the VERTEX:xyz values
//
// The push/pop is preferable to using VERTEX members because it uses all the entries in a cache line
//
struct MRISsavedXYZ {
  int       nvertices;
  MRIS_XYZ* p;
};

MRISsavedXYZ* MRISsaveXYZ(MRIS *mris) {
  MRISsavedXYZ* saved = (MRISsavedXYZ*)malloc(sizeof(MRISsavedXYZ));
  int       const nvertices = saved->nvertices = mris->nvertices;
  MRIS_XYZ* const p         = saved->p         = (MRIS_XYZ*)malloc(nvertices*sizeof(MRIS_XYZ));
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    p[vno].x = v->x;
    p[vno].y = v->y;
    p[vno].z = v->z;
  }
  return saved;
}

void MRISloadXYZ (MRIS *mris, MRISsavedXYZ const * const pMRISsavedXYZ) {
  MRIS_XYZ* const p         = pMRISsavedXYZ->p;
  int       const nvertices = pMRISsavedXYZ->nvertices;
  cheapAssert(nvertices == mris->nvertices);
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    MRISsetXYZ(mris,vno,p[vno].x,p[vno].y,p[vno].z);
  }
}

void MRISfreeXYZ (MRIS *mris, MRISsavedXYZ ** ppMRISsavedXYZ) {
  MRISsavedXYZ* p = *ppMRISsavedXYZ; *ppMRISsavedXYZ = NULL;
  freeAndNULL(p->p);
  freeAndNULL(p);
}

void MRISpopXYZ (MRIS *mris, MRISsavedXYZ ** ppMRISsavedXYZ) {
  MRISloadXYZ(mris, *ppMRISsavedXYZ);
  MRISfreeXYZ(mris,  ppMRISsavedXYZ);
}


/*-----------------------------------------------------
    These functions abuse the ORIGINAL_VERTICES and other fields
    using them as temp storage.  This abuse is especially bad for the 
    orig_xyz because these are input the two different dist_orig algorithms
    so putting other values into there is especially bad.
  ------------------------------------------------------*/
int MRISsaveVertexPositions(MRIS *mris, int which)
{
  if (which == ORIGINAL_VERTICES) {

    MRISsetOriginalXYZfromXYZ(mris);

  } else {
  
    int const nvertices = mris->nvertices;

    int vno;
    for (vno = 0; vno < nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      switch (which) {
        case LAYERIV_VERTICES:
          v->l4x = v->x;
          v->l4y = v->y;
          v->l4z = v->z;
          break;
        case TARGET_VERTICES:
          v->targx = v->x;
          v->targy = v->y;
          v->targz = v->z;
          break;
        case WHITE_VERTICES:
          v->whitex = v->x;
          v->whitey = v->y;
          v->whitez = v->z;
          break;
        case PIAL_VERTICES:
          v->pialx = v->x;
          v->pialy = v->y;
          v->pialz = v->z;
          break;
        case INFLATED_VERTICES:
          v->infx = v->x;
          v->infy = v->y;
          v->infz = v->z;
          break;
        case FLATTENED_VERTICES:
          v->fx = v->x;
          v->fy = v->y;
          v->fz = v->z;
          break;
        case CANONICAL_VERTICES:
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;
          break;
        case TMP2_VERTICES:
          v->t2x = v->x;
          v->t2y = v->y;
          v->t2z = v->z;
          break;
        case TMP_VERTICES:
          v->tx = v->x;
          v->ty = v->y;
          v->tz = v->z;
          break;
        default:
          cheapAssert(false);
      }
    }
  }
  
  if (which == CANONICAL_VERTICES) {
    MRIScomputeCanonicalCoordinates(mris);
  }
  
  return (NO_ERROR);
}


int MRISrestoreVertexPositions(MRIS *mris, int which)
{
  if (which == ORIGINAL_VERTICES) {
    
    // Verify it is plausible to copy these values into the MRIS vertices xyz 
    // without changing the status
    //
    checkOrigXYZCompatible(mris->status,mris->origxyz_status);
  }
  
  int const nvertices = mris->nvertices;
  
  mris->dist_nsize = 0;
  
  int vno;
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    switch (which) {
      case TARGET_VERTICES:
        v->x = v->targx;
        v->y = v->targy;
        v->z = v->targz;
        break;
      case WHITE_VERTICES:
        v->x = v->whitex;
        v->y = v->whitey;
        v->z = v->whitez;
        break;
      case LAYERIV_VERTICES:
        v->x = v->l4x;
        v->y = v->l4y;
        v->z = v->l4z;
        break;
      case PIAL_VERTICES:
        v->x = v->pialx;
        v->y = v->pialy;
        v->z = v->pialz;
        break;
      case INFLATED_VERTICES:
        v->x = v->infx;
        v->y = v->infy;
        v->z = v->infz;
        break;
      case FLATTENED_VERTICES:
        v->x = v->fx;
        v->y = v->fy;
        v->z = v->fz;
        break;
      case CANONICAL_VERTICES:
        v->x = v->cx;
        v->y = v->cy;
        v->z = v->cz;
        break;
      case ORIGINAL_VERTICES:
        v->x = v->origx;
        v->y = v->origy;
        v->z = v->origz;
        break;
      case TMP2_VERTICES:
        v->x = v->t2x;
        v->y = v->t2y;
        v->z = v->t2z;
        break;
      default:
      case TMP_VERTICES:
        v->x = v->tx;
        v->y = v->ty;
        v->z = v->tz;
        break;
    }
  }
  if (mris->status == MRIS_SPHERICAL_PATCH) {
    for (vno = 0; vno < nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      switch (which) {
        case WHITE_VERTICES:
          v->cx = v->whitex;
          v->cy = v->whitey;
          v->cz = v->whitez;
          break;
        case PIAL_VERTICES:
          v->cx = v->pialx;
          v->cy = v->pialy;
          v->cz = v->pialz;
          break;
        case INFLATED_VERTICES:
          v->cx = v->infx;
          v->cy = v->infy;
          v->cz = v->infz;
          break;
        case FLATTENED_VERTICES:
          v->cx = v->fx;
          v->cy = v->fy;
          v->cz = v->fz;
          break;
        case CANONICAL_VERTICES:
          v->cx = v->cx;
          v->cy = v->cy;
          v->cz = v->cz;
          break;
        case ORIGINAL_VERTICES:
          v->cx = v->origx;
          v->cy = v->origy;
          v->cz = v->origz;
          break;
        case TMP2_VERTICES:
          v->cx = v->t2x;
          v->cy = v->t2y;
          v->cz = v->t2z;
          break;
        default:
        case TMP_VERTICES:
          v->cx = v->tx;
          v->cy = v->ty;
          v->cz = v->tz;
          break;
      }
    }
  }
  mrisComputeSurfaceDimensions(mris);
  return (NO_ERROR);
}


int MRISsaveNormals(MRIS *mris, int which)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    switch (which) {
      case TMP_VERTICES:
        v->tx = v->nx;
        v->ty = v->ny;
        v->tz = v->nz;
        break;
      case TMP2_VERTICES:
        v->t2x = v->nx;
        v->t2y = v->ny;
        v->t2z = v->nz;
        break;
      case PIAL_VERTICES:
        v->pnx = v->nx;
        v->pny = v->ny;
        v->pnz = v->nz;
        break;
      case WHITE_VERTICES:
        v->wnx = v->nx;
        v->wny = v->ny;
        v->wnz = v->nz;
        break;
      case INFLATED_VERTICES:
      case FLATTENED_VERTICES:
      case CANONICAL_VERTICES:
      case ORIGINAL_VERTICES:
      default:
        ErrorExit(ERROR_BADPARM, "MRISsaveNormals: unsupported which %d", which);
        break;
    }
  }
  return (NO_ERROR);
}


int MRISrestoreNormals(MRIS *mris, int which)
{
  int vno, nvertices;
  VERTEX *v;

  nvertices = mris->nvertices;
  for (vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    switch (which) {
      case TMP_VERTICES:
        v->nx = v->tx;
        v->ny = v->ty;
        v->nz = v->tz;
        break;
      case TMP2_VERTICES:
        v->nx = v->t2x;
        v->ny = v->t2y;
        v->nz = v->t2z;
        break;
      case PIAL_VERTICES:
        v->nx = v->pnx;
        v->ny = v->pny;
        v->nz = v->pnz;
        break;
      case WHITE_VERTICES:
      case INFLATED_VERTICES:
      case FLATTENED_VERTICES:
      case CANONICAL_VERTICES:
      case ORIGINAL_VERTICES:
      default:
        ErrorExit(ERROR_BADPARM, "MRISsaveNormals: unsupported which %d", which);
        break;
    }
  }
  return (NO_ERROR);
}


float FACES_commonEdgeLength_find(MRIS *apmris, FACE *apFACE_I, FACE *apFACE_J)
{
  //
  // PRECONDITIONS
  //  o The <apFACE>s should be triangles with 3 vertices each.
  //  o The FACES share a common edge.
  //  o The FACES are not the same.
  //
  // POSTCONDITIONS
  //  o Length of common edge is returned.
  //  o If common vertices != 2, then function ErrorExits.
  //

  static int calls = 0;
  const char *pch_function = "FACES_commonEdgeLength_find";
  VERTEX *pVERTEX_O = NULL;                 // Common vertex O
  VERTEX *pVERTEX_I = NULL;                 // Common vertex I
  static VECTOR *pVECTOR_O = NULL;          // Common vertex O cart. coords
  static VECTOR *pVECTOR_I = NULL;          // Common vertex I cart. coords
  static VECTOR *pVECTOR_edgeVoVi = NULL;   // Edge Vo->Vi
  static VECTOR *pv_commonVertices = NULL;  // Vector housing vertices that
  // are common between two
  // neighbouring faces.
  int commonVertices = 0;  // number of vertices in common
  // between two faces
  float f_edgeLength = 0.;  // Length of edge v->vI

  DebugEnterFunction(("%s", pch_function));
  if (!calls) {
    pv_commonVertices = VectorAlloc(3, MATRIX_REAL);
    pVECTOR_O = VectorAlloc(3, MATRIX_REAL);
    pVECTOR_I = VectorAlloc(3, MATRIX_REAL);
    pVECTOR_edgeVoVi = VectorAlloc(3, MATRIX_REAL);
  }
  commonVertices = VERTICES_commonInFaces_find(apFACE_I, apFACE_J, pv_commonVertices);
  if (commonVertices != 2)
    ErrorExit(-4, "%s: No common edge found! <commonVertices> = %d\n", pch_function, commonVertices);

  pVERTEX_O = &apmris->vertices[(int)VECTOR_ELT(pv_commonVertices, 1)];
  pVERTEX_I = &apmris->vertices[(int)VECTOR_ELT(pv_commonVertices, 2)];

  V3_LOAD(pVECTOR_O, pVERTEX_O->x, pVERTEX_O->y, pVERTEX_O->z);
  V3_LOAD(pVECTOR_I, pVERTEX_I->x, pVERTEX_I->y, pVERTEX_I->z);
  V3_SUBTRACT(pVECTOR_I, pVECTOR_O, pVECTOR_edgeVoVi);
  f_edgeLength = V3_LEN(pVECTOR_edgeVoVi);
  calls++;
  xDbg_PopStack();
  return (f_edgeLength);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeTangentPlanes(MRIS *mris)
{
  VECTOR *v_n, *v_e1, *v_e2, *v;
  int vno;
  VERTEX *vertex;

  v_n  = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v    = VectorAlloc(3, MATRIX_REAL);

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz);
    VECTOR_LOAD(v,   vertex->ny, vertex->nz, vertex->nx);
    V3_CROSS_PRODUCT(v_n, v, v_e1);
    if (VectorLen(v_e1) < 0.001) /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx);
      V3_CROSS_PRODUCT(v_n, v, v_e1);
    }

    if ((V3_LEN_IS_ZERO(v_e1)) && DIAG_VERBOSE_ON) /* happened to pick a parallel vector */
    {
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno);
    }
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2);
    V3_NORMALIZE(v_e1, v_e1);
    V3_NORMALIZE(v_e2, v_e2);
    if (V3_LEN(v_e1) < 0.5) {
      DiagBreak();
    }
    vertex->e1x = V3_X(v_e1);
    vertex->e2x = V3_X(v_e2);
    vertex->e1y = V3_Y(v_e1);
    vertex->e2y = V3_Y(v_e2);
    vertex->e1z = V3_Z(v_e1);
    vertex->e2z = V3_Z(v_e2);
  }

  VectorFree(&v);
  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  return (NO_ERROR);
}

int MRISstoreTangentPlanes(MRIS *mris, int which_vertices)
{
  VERTEX *v;
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->pe1x = v->e1x;
    v->pe1y = v->e1y;
    v->pe1z = v->e1z;
    v->pe2x = v->e2x;
    v->pe2y = v->e2y;
    v->pe2z = v->e2z;
  }
  return (NO_ERROR);
}


int mrisComputeOptimalPlane(MRIS *mris,
                            int *vertices,
                            int nvertices,
                            double *pxn,
                            double *pyn,
                            double *pzn,
                            double *pxc,
                            double *pyc,
                            double *pzc)
{
  int vno, n, vnum;
  MATRIX *M, *m_evectors;
  VERTEX *v;
  double a, b, c, xc, yc, zc, norm;
  float evalues[3];

  M = MatrixAlloc(3, 3, MATRIX_REAL);
  m_evectors = MatrixAlloc(3, 3, MATRIX_REAL);
  for (vnum = 0, xc = yc = zc = 0.0, n = 0; n < nvertices; n++) {
    vno = vertices[n];
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }

    xc += v->x;
    yc += v->y;
    zc += v->z;
    vnum++;
  }
  if (vnum < 4) {
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "mrisComputeOptimalPlane: must specify at least 4 non-ripped verices\n"));
  }

  for (xc = yc = zc = 0.0, n = 0; n < nvertices; n++) {
    vno = vertices[n];
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    *MATRIX_RELT(M, 1, 1) += SQR(v->x - xc);
    *MATRIX_RELT(M, 2, 2) += SQR(v->y - yc);
    *MATRIX_RELT(M, 3, 3) += SQR(v->z - zc);

    *MATRIX_RELT(M, 2, 1) += (v->x - xc) * (v->y - yc);
    *MATRIX_RELT(M, 1, 2) += (v->x - xc) * (v->y - yc);

    *MATRIX_RELT(M, 3, 1) += (v->x - xc) * (v->z - zc);
    *MATRIX_RELT(M, 1, 3) += (v->x - xc) * (v->z - zc);

    *MATRIX_RELT(M, 3, 2) += (v->y - yc) * (v->z - zc);
    *MATRIX_RELT(M, 2, 3) += (v->y - yc) * (v->z - zc);
  }
  MatrixEigenSystem(M, evalues, m_evectors);

  // evalues are distance squared to plane, so use smallest one, which is in 3rd col
  a = *MATRIX_RELT(m_evectors, 1, 3);
  b = *MATRIX_RELT(m_evectors, 2, 3);
  c = *MATRIX_RELT(m_evectors, 3, 3);
  norm = sqrt(a * a + b * b + c * c);
  a /= norm;
  b /= norm;
  c /= norm;

  *pxn = a;
  *pyn = b;
  *pzn = c;
  *pxc = xc;
  *pyc = yc;
  *pzc = zc;
  MatrixFree(&m_evectors);
  MatrixFree(&M);
  return (NO_ERROR);
}



/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute the folding of the surface.
  ------------------------------------------------------*/
double MRIScomputeFolding(MRIS *mris)
{
  int vno;
  VERTEX *vertex;
  double k1, k2, folding, area;

  folding = 0.0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    k1 = (double)vertex->k1;
    k2 = (double)vertex->k2;
    area = (double)vertex->area;
    folding += area * (k1 - k2) * (k1 - k2);
  }

  return (folding);
}


int MRIScomputeCanonicalCoordinates(MRIS *mris)
{
  float theta, phi, r, d, x, y, z;
  VERTEX *v;
  int vno;

  r = mris->radius = MRISaverageCanonicalRadius(mris);
  r = mris->radius = (float)nint(mris->radius);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    x = v->cx;
    y = v->cy;
    z = v->cz;
    theta = atan2f(y, x);
    if (theta < 0.0f) {
      theta = 2 * M_PI + theta; /* make it 0 --> 2*PI */
    }
    d = r * r - z * z;
    if (d < 0.0) {
      d = 0.0;
    }
    phi = atan2f(sqrt(d), z);
    v->theta = theta;
    v->phi = phi;
  }
  return (NO_ERROR);
}

int MRISvertexCoord2XYZ_float(VERTEX const * const v, int const which, float * const x, float * const y, float * const z)
{
#define CASE(WHICH, FIELD) \
  case WHICH: \
    cheapAssert(sizeof(*x) == sizeof(v->FIELD##x)); \
    *x = v->FIELD##x;  *y = v->FIELD##y;  *z = v->FIELD##z; \
    break;
        
  switch (which) {
    CASE(ORIGINAL_VERTICES,orig)
    CASE(TMP_VERTICES,t)
    CASE(CANONICAL_VERTICES,c)
    CASE(CURRENT_VERTICES,)
    CASE(INFLATED_VERTICES,inf)
    CASE(PIAL_VERTICES,pial)
    CASE(TMP2_VERTICES,t2)   
    CASE(FLATTENED_VERTICES,f)
    CASE(WHITE_VERTICES,white)
    CASE(VERTEX_NORMALS,n)
    CASE(PIAL_NORMALS,pn)
    CASE(WHITE_NORMALS,wn)
    default:
      ErrorExit(ERROR_UNSUPPORTED, "MRISvertexCoord2XYZ_float: unsupported which %d", which);
      break;
  }

#undef CASE
  
  return (NO_ERROR);
}

int MRISvertexCoord2XYZ_double(VERTEX const * const v, int const which, double * const x, double * const y, double * const z)
{
    float fx,fy,fz;
    auto result = MRISvertexCoord2XYZ_float(v, which, &fx, &fy, &fz);
    *x = fx; *y = fy; *z = fz;
    return result;
}

#define VERTEX_DIF(leg, v0, v1) leg[0] = v1->x - v0->x, leg[1] = v1->y - v0->y, leg[2] = v1->z - v0->z;

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
int mrisCalculateCanonicalFaceCentroid(MRIS *mris, int fno, float *px, float *py, float *pz)
{
  float x, y, z;
  VERTEX *v0, *v1, *v2;
  FACE *f;

  f = &mris->faces[fno];
  v0 = &mris->vertices[f->v[0]];
  v1 = &mris->vertices[f->v[1]];
  v2 = &mris->vertices[f->v[2]];

  /* first bisect v1->v2 line */

  x = (v1->cx + v2->cx) / 2.0f;
  y = (v1->cy + v2->cy) / 2.0f;
  z = (v1->cz + v2->cz) / 2.0f;

  /* now bisect v0->bisector line */
  *px = (v0->cx + x) / 2.0f;
  *py = (v0->cy + y) / 2.0f;
  *pz = (v0->cz + z) / 2.0f;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
static int mrisCalculateOriginalFaceCentroid(MRIS *mris, int fno, float *px, float *py, float *pz)
{
  float x, y, z;
  VERTEX *v0, *v1, *v2;
  FACE *f;

  f = &mris->faces[fno];
  v0 = &mris->vertices[f->v[0]];
  v1 = &mris->vertices[f->v[1]];
  v2 = &mris->vertices[f->v[2]];

  /* first bisect v1->v2 line */

  x = (v1->origx + v2->origx) / 2.0f;
  y = (v1->origy + v2->origy) / 2.0f;
  z = (v1->origz + v2->origz) / 2.0f;

  /* now bisect v0->bisector line */
  *px = (v0->origx + x) / 2.0f;
  *py = (v0->origy + y) / 2.0f;
  *pz = (v0->origz + z) / 2.0f;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
int mrisCalculateFaceCentroid(MRIS *mris, int fno, float *px, float *py, float *pz)
{
  float x, y, z;
  VERTEX *v0, *v1, *v2;
  FACE *f;

  f = &mris->faces[fno];
  v0 = &mris->vertices[f->v[0]];
  v1 = &mris->vertices[f->v[1]];
  v2 = &mris->vertices[f->v[2]];

  /* first bisect v1->v2 line */

  x = (v1->x + v2->x) / 2.0f;
  y = (v1->y + v2->y) / 2.0f;
  z = (v1->z + v2->z) / 2.0f;

  /* now bisect v0->bisector line */
  *px = (v0->x + x) / 2.0f;
  *py = (v0->y + y) / 2.0f;
  *pz = (v0->z + z) / 2.0f;
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Choose the face whose centroid is closest to the orig position
  of v.
  ------------------------------------------------------*/

int mrisChooseFace(MRIS *mris, MHT *mht, VERTEX *v)
{
  int fno, nfound, flist[1000], min_fno, i, j, total_found, *fptr;
  double dist, d;
  float dx, dy, dz, cx[1000], cy[1000], cz[1000], total_dist, max_dist;

  for (total_found = nfound = 0, dist = -.25; dist <= .25; dist += .25) {
    d = dist;
    fptr = &flist[total_found];
    nfound = mrisAllNormalDirectionCurrentTriangleIntersections(mris, v, mht, &d, fptr);
    if (nfound > 0) {
      for (i = 0; i < total_found; i++) {
        for (j = 0; j < nfound; j++) {
          if (flist[i] == fptr[j]) {
            fptr[j] = -1; /* was already found */
          }
        }
      }
      for (j = 0; j < nfound; j++) {
        if (fptr[j] >= 0) {
          if (total_found == 1000) {
            ErrorExit(ERROR_BADPARM, "Too many surface found");
          }
          flist[total_found++] = fptr[j];
        }
      }
    }
  }

  if (total_found <= 0) {
    return (-1);
  }

  /*
    use face which is furthest distance from negative faces.
  */
  max_dist = -10.0f;
  min_fno = -1;

  for (i = 0; i < total_found; i++) {
    fno = flist[i];
    mrisCalculateOriginalFaceCentroid(mris, fno, cx + i, cy + i, cz + i);
  }
  for (i = 0; i < total_found; i++) {
    fno = flist[i];
    if (mris->faces[fno].area < 0) {
      continue; /* don't map it to one with negative area */
    }

    for (total_dist = 0.0, j = 0; j < total_found; j++) {
      if (mris->faces[flist[j]].area > 0) {
        continue;
      }
      dx = cx[j] - cx[i];
      dy = cy[j] - cy[i];
      dz = cz[j] - cz[i];
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      total_dist += dist;
    }
    if (total_dist > max_dist) {
      max_dist = dist;
      min_fno = fno;
    }
  }
  return (min_fno);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Find the face in which v lies. If it lies in more than one
  face, return not found.
  ------------------------------------------------------*/
int mrisFindUnambiguousFace(MRIS *mris, MHT *mht, VERTEX *v, int *pnfound)
{
  int nfound, flist[1000], *fptr, total_found, i, j;
  double dist, d;

  for (total_found = nfound = 0, dist = -.25; dist <= .25; dist += .25) {
    d = dist;
    fptr = &flist[total_found];
    nfound = mrisAllNormalDirectionCurrentTriangleIntersections(mris, v, mht, &d, fptr);
    if (nfound > 0) {
      for (i = 0; i < total_found; i++) {
        for (j = 0; j < nfound; j++) {
          if (flist[i] == fptr[j]) {
            fptr[j] = -1; /* was already found */
          }
        }
      }
      for (j = 0; j < nfound; j++) {
        if (fptr[j] >= 0) {
          if (total_found == 1000) {
            ErrorExit(ERROR_BADPARM, "Too many ambiguous faces");
          }
          flist[total_found++] = fptr[j];
        }
      }
    }
  }

  *pnfound = total_found;
  return (total_found == 1 ? flist[0] : -1);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#include "stats.h"
int MRISsampleStatVolume(MRIS *mris, void *vsv, int time_point, int coords)
{
  VERTEX *v;
  int vno, xv, yv, zv, width, height, depth;
  double x, y, z, xt, yt, zt;
  STAT_VOLUME *sv = (STAT_VOLUME *)vsv;

  if (time_point >= sv->mri_pvals[0]->nframes)
    ErrorExit(ERROR_BADPARM,
              "MRISsampleStatVolume: time point (%d) out of bounds [%d, %d]\n",
              time_point,
              0,
              sv->mri_pvals[0]->nframes - 1);
  width = sv->mri_pvals[0]->width;
  height = sv->mri_pvals[0]->height;
  depth = sv->mri_pvals[0]->depth;
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == 47) {
      DiagBreak();
    }
    v = &mris->vertices[vno];
    x = (double)v->x;
    y = (double)v->y;
    z = (double)v->z;

    /* now convert them into talairach space */
    switch (coords) {
      case TALAIRACH_COORDS:
        MRIworldToTalairachVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt);
        break;
      case SPHERICAL_COORDS:
        x = (double)v->cx;
        y = (double)v->cy;
        z = (double)v->cz;
        MRIworldToVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt);
        break;
      default:
        MRIworldToVoxel(sv->mri_pvals[0], x, y, z, &xt, &yt, &zt);
        break;
    }
    xv = nint(xt);
    yv = nint(yt);
    zv = nint(zt);
    if (xv >= 0 && xv < width && yv >= 0 && yv < height && zv >= 0 && zv < depth) {
      v->val = MRIFseq_vox(sv->mri_pvals[0], xv, yv, zv, time_point);
    }
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScomputeFaceAreaStats(MRIS *mris, double *psigma, double *pmin, double *pmax)
{
  double total_area, mean, var, nf, sigma, min_area, max_area, area, area_scale;
  int fno;
  FACE *f;

  MRIScomputeMetricProperties(mris);

  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }

  min_area = 1000;
  max_area = -1;
  for (var = nf = total_area = 0.0, fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    area = area_scale * f->area;
    total_area += area;
    nf++;
    var += area * area;
    if (area > max_area) {
      max_area = area;
    }
    if (area < min_area) {
      min_area = area;
    }
  }
  mean = total_area / nf;
  *psigma = sigma = sqrt(var / nf - mean * mean);
  if (pmin) {
    *pmin = min_area;
  }
  if (pmax) {
    *pmax = max_area;
  }
  return (mean);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScomputeVertexSpacingStats(
    MRIS *mris, double *psigma, double *pmin, double *pmax, int *pvno, int *pvno2, int which_vertices)
{
  double total_dist, mean, var, nv, dist, sigma, min_dist, max_dist, dist_scale;
  int vno, n;

  MRIScomputeMetricProperties(mris);
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
  dist_scale = 1.0f;
  min_dist = 1000;
  max_dist = -1;
  for (var = nv = total_dist = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      nv++;
      switch (which_vertices) {
        default:
          ErrorExit(ERROR_BADPARM, "MRIScomputeVertexSpacingStats: unsupported vertex set %d", which_vertices);
          dist = 0;
          break;
        case ORIGINAL_VERTICES:
          dist = sqrt(SQR(vn->origx - v->origx) + SQR(vn->origy - v->origy) + SQR(vn->origz - v->origz));
          break;
        case CURRENT_VERTICES:
          dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
          break;
        case WHITE_VERTICES:
          dist = sqrt(SQR(vn->whitex - v->whitex) + SQR(vn->whitey - v->whitey) + SQR(vn->whitez - v->whitez));
          break;
        case PIAL_VERTICES:
          dist = sqrt(SQR(vn->pialx - v->pialx) + SQR(vn->pialy - v->pialy) + SQR(vn->pialz - v->pialz));
          break;
        case TARGET_VERTICES:
          dist = sqrt(SQR(vn->targx - v->targx) + SQR(vn->targy - v->targy) + SQR(vn->targz - v->targz));
          break;
      }

      dist *= dist_scale;
      if (dist > max_dist) {
        if (pvno) {
          *pvno = vno;
        }
        if (pvno2) {
          *pvno2 = vt->v[n];
        }
        max_dist = dist;
      }
      if (dist < min_dist) {
        min_dist = dist;
      }
      total_dist += dist;
      var += dist * dist;
    }
  }

  mean = total_dist / nv;
  
  if (psigma) {
    *psigma = sigma = sqrt(var / nv - mean * mean);
  }
  if (pmin) {
    *pmin = min_dist;
  }
  if (pmax) {
    *pmax = max_dist;
  }
  return (mean);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScomputeTotalVertexSpacingStats(
    MRIS *mris, double *psigma, double *pmin, double *pmax, int *pvno, int *pvno2)
{
  double total_dist, mean, var, nv, dist, sigma, min_dist, max_dist, dist_scale;
  int vno, n;

  MRIScomputeMetricProperties(mris);
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
  dist_scale = 1.0f;
  min_dist = 1000;
  max_dist = -1;
  for (var = nv = total_dist = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      nv++;
      dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
      dist *= dist_scale;
      if (dist > max_dist) {
        if (pvno) {
          *pvno = vno;
        }
        if (pvno2) {
          *pvno2 = vt->v[n];
        }
        max_dist = dist;
      }
      if (dist < min_dist) {
        min_dist = dist;
      }
      total_dist += dist;
      var += dist * dist;
    }
  }
  mean = total_dist / nv;
  if (psigma) {
    *psigma = sigma = sqrt(var / nv - mean * mean);
  }
  if (pmin) {
    *pmin = min_dist;
  }
  if (pmax) {
    *pmax = max_dist;
  }
  return (mean);
}



/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Expand the list of neighbors of each vertex, reallocating
  the v->v array to hold the expanded list.
  ------------------------------------------------------*/
static void MRISsetNeighborhoodSizeAndDistWkr(MRIS *mris, int nsize, bool alwaysDoDist);

static void MRISsetNeighborhoodSizeAndOptionallyDist(MRIS *mris, int nsize, bool alwaysDoDist) {
  static int count, limit;
  if (++count == limit) {
    fprintf(stdout, "%s:%d MRISsetNeighborhoodSizeAndDist count:%d\n",__FILE__,__LINE__,count);
  }
  mrisCheckVertexFaceTopology(mris);
  
  MRISsetNeighborhoodSizeAndDistWkr(mris, nsize, alwaysDoDist);
  
  if (!(mris->dist_alloced_flags & 1))  { cheapAssert(mris->dist_nsize == 0);     }
  else                                  { cheapAssert(mris->dist_nsize >= nsize); }
  mrisCheckVertexFaceTopology(mris);
}

void MRISsetNeighborhoodSize(MRIS *mris, int nsize) {
    MRISsetNeighborhoodSizeAndOptionallyDist(mris, nsize, false);
}

void MRISsetNeighborhoodSizeAndDist(MRIS *mris, int nsize) {
    MRISsetNeighborhoodSizeAndOptionallyDist(mris, nsize, true);
}


static void MRISsetNeighborhoodSizeAndDistWkr(MRIS *mris, int nsize, bool alwaysDoDist)
{
  cheapAssert(1 <= nsize && nsize < 4);
  
  /*
    now build a list of 2-connected neighbors. After this is done,
    reallocate the v->n list and arrange the 2-connected neighbors
    sequentially after it.
  */

  bool const max_nsize_grew = (nsize > mris->max_nsize);
  
  if (nsize <= mris->max_nsize) {
    cheapAssert(!max_nsize_grew);
    
    MRISresetNeighborhoodSize(mris, nsize);

  } else {
    cheapAssert(max_nsize_grew);

    alwaysDoDist |= (mris->dist_alloced_flags & 1);

    // rebuild all the neighbor caches
    //
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      int vlist[MAX_NEIGHBORS], hops[MAX_NEIGHBORS];
      MRISfindNeighborsAtVertex(mris, vno, nsize, MAX_NEIGHBORS, vlist, hops);
    }
    MRISresetNeighborhoodSize(mris, nsize);

    mris->max_nsize = nsize;
  }

  cheapAssert(mris->nsize == nsize);
  
  if (nsize <= mris->dist_nsize) {
    return;
  }
  	
  // Adjust the dist and dist_orig
  // unless they are already valid
  //
  bool useOldBehaviour = false;
  if (useOldBehaviour) {
    switch (copeWithLogicProblem("FREESURFER_fix_MRISsetNeighborhoodSizeAndDist","should not be creating dist_origs")) {
    case LogicProblemResponse_old:
      useOldBehaviour = true;
      break;
    case LogicProblemResponse_fix:
      useOldBehaviour = false; 
    }
  }
  
  bool shouldDoDist     = alwaysDoDist || max_nsize_grew || mris->dist_nsize < nsize ;
  bool shouldDoDistOrig = shouldDoDist && (mris->dist_alloced_flags & 2);

  if (useOldBehaviour) {
    if (shouldDoDist != alwaysDoDist) {
      copeWithLogicProblem(NULL,alwaysDoDist?"new not computing dist when old did":"old not computing dist");
      shouldDoDist = alwaysDoDist;
    }
    if (shouldDoDistOrig != alwaysDoDist) {
      copeWithLogicProblem(NULL,alwaysDoDist?"new not computing dist_orig when old did":"old not computing dist_orig");
      shouldDoDistOrig = alwaysDoDist;
    }
  }
  
  if (shouldDoDist)     mrisComputeVertexDistances(mris);
  if (shouldDoDistOrig) mrisComputeOriginalVertexDistances(mris);

  cheapAssert( (!(mris->dist_alloced_flags&1)) || (mris->dist_nsize      >= nsize) );
  cheapAssert( (!(mris->dist_alloced_flags&2)) || (mris->dist_orig_nsize >= nsize) );
}


int IsMRISselfIntersecting(MRIS *mris)
{
  MRIS_HASH_TABLE *mht;
  int fno;

  mht = MHTcreateFaceTable(mris);

  for (fno = 0; fno < mris->nfaces; fno++) {
    if (MHTdoesFaceIntersect(mht, mris, fno)) {
      MHTfree(&mht);
      return 1;
    }
  }

  MHTfree(&mht);
  return 0;
}


int face_barycentric_coords(
  double V0[3], double V1[3], double V2[3],
  double  cx,
  double  cy,
  double  cz,
  double *pl1,
  double *pl2,
  double *pl3)
{
  double point[3], proj[3], detT;
  double l1, l2, l3, x, y, x1, x2, x3, y1, y2, y3, e1[3], e2[3];
  
  int ret = 0;
  point[0] = cx;
  point[1] = cy;
  point[2] = cz;
  if (project_point_to_plane(point, V0, V1, V2, proj, e1, e2) < 0) {
    if (pl1) {
      *pl1 = *pl2 = *pl3 = 1.0 / 3;
    }
    return (-3);
  }
  x = DOT(proj, e1);
  y = DOT(proj, e2);
  x1 = DOT(V0, e1);
  x2 = DOT(V1, e1);
  x3 = DOT(V2, e1);
  y1 = DOT(V0, e2);
  y2 = DOT(V1, e2);
  y3 = DOT(V2, e2);
  detT = (x1 - x3) * (y2 - y3) - (y1 - y3) * (x2 - x3);
  if (DZERO(detT)) {
    if (pl1) {
      *pl1 = *pl2 = *pl3 = 1 / 3;
    }
    return (-2);
  }

  l1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
  l2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
  l3 = (1 - l1 - l2);

  if (pl1) {
    *pl1 = l1;
    *pl2 = l2;
    *pl3 = l3;
  }

  if (l1 > 0 && l1 < 1 && l2 > 0 && l2 < 1 && l3 > 0 && l3 < 1) {
    ret = 1;  // interior to the triangle
  }
  else {
    float l1d, l2d, l3d;

    l1d = MIN(l1 - 0, 1 - l1);
    l2d = MIN(l2 - 0, 1 - l2);
    l3d = MIN(l3 - 0, 1 - l3);
    if ((l1d < 0 && !FZERO(l1d)) || (l2d < 0 && !FZERO(l2d)) || (l3d < 0 && !FZERO(l3d))) {
      ret = -1;  // outside triangle
    }
    else {
      ret = 0;  // on  boundary of triangle
    }
  }
  return (ret);
}


int face_barycentric_coords(MRIS const *mris,
                            int fno,
                            int which_vertices,
                            double cx,
                            double cy,
                            double cz,
                            double *pl1,
                            double *pl2,
                            double *pl3)
{
  auto surface = SurfaceFromMRIS::Analysis::Surface(const_cast<MRIS*>(mris));
  return face_barycentric_coords_template(surface, fno, which_vertices,cx,cy,cz,pl1,pl2,pl3);   
}

int MRISsampleFaceCoords(MRIS *mris,
                         int fno,
                         double x,
                         double y,
                         double z,
                         int which_coords,
                         int which_barycentric,
                         float *px,
                         float *py,
                         float *pz)
{
  float xv, yv, zv;
  double lambda[3];
  int n, ret;
  FACE *face;
  VERTEX *v;

  face = &mris->faces[fno];

  xv = yv = zv = 0.0;  // to get rid of mac warnings
  ret = face_barycentric_coords(mris, fno, which_barycentric, x, y, z, &lambda[0], &lambda[1], &lambda[2]);
  if (ret < 0) {
    lambda[0] = lambda[1] = lambda[2] = 1.0 / 3.0;
  }

  *px = *py = *pz = 0;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[face->v[n]];
    MRISvertexCoord2XYZ_float(v, which_coords, &xv, &yv, &zv);
    *px += lambda[n] * xv;
    *py += lambda[n] * yv;
    *pz += lambda[n] * zv;
  }

  return (ret);
}

double MRISsampleFace(
    MRIS *mris, int fno, int which, double x, double y, double z, double val0, double val1, double val2)
{
  float xv, yv, zv;
  double lambda[3], val;
  int n, ret;
  FACE *face;
  VERTEX *v;

  face = &mris->faces[fno];

  xv = yv = zv = 0.0;  // to get rid of mac warnings
  ret = face_barycentric_coords(mris, fno, which, x, y, z, &lambda[0], &lambda[1], &lambda[2]);

  val = 0.0;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[face->v[n]];
    MRISvertexCoord2XYZ_float(v, which, &xv, &yv, &zv);
    switch (n) {
      default:
      case 0:
        val += lambda[0] * val0;
        break;
      case 1:
        val += lambda[1] * val1;
        break;
      case 2:
        val += lambda[2] * val2;
        break;
    }
  }

  return (val);
}
int MRISsampleFaceNormal(MRIS *mris, int fno, double x, double y, double z, float *px, float *py, float *pz)
{
  float d, dtotal, dx, dy, dz, xc, yc, zc;
  int n;
  FACE *face;
  VERTEX *v;

  face = &mris->faces[fno];

  xc = yc = zc = 0.0;  // to get rid of mac warnings
  for (dtotal = 0.0, n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[face->v[n]];
    MRISvertexCoord2XYZ_float(v, CANONICAL_VERTICES, &dx, &dy, &dz);
    dx -= x;
    dy -= y;
    dz -= z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < 0) {
      continue;
    }
    dtotal += d;
  }

  *px = *py = *pz = 0;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[face->v[n]];
    MRISvertexCoord2XYZ_float(v, CANONICAL_VERTICES, &dx, &dy, &dz);
    dx -= x;
    dy -= y;
    dz -= z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    d = 1 - (2 * d) / dtotal;

    if (d < 0) {
      continue;
    }
    *px += d * v->nx;
    *py += d * v->ny;
    *pz += d * v->nz;
  }

  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISfindClosestVertex(MRIS *mris, float x, float y, float z, float *dmin, int which_vertices)
{
  int vno, min_v = -1;
  VERTEX *v;
  float d, min_d, dx, dy, dz;

  min_d = 10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    switch (which_vertices)
    {
    case CURRENT_VERTICES:
      dx = v->x - x;
      dy = v->y - y;
      dz = v->z - z;
      break ;
    case CANONICAL_VERTICES:
      dx = v->cx - x;
      dy = v->cy - y;
      dz = v->cz - z;
      break ;
    default:
      dx = v->x - x;
      dy = v->y - y;
      dz = v->z - z;
      ErrorExit(ERROR_UNSUPPORTED, "MRISfindClosestVertex: unsupported vertex set %d\n", which_vertices);
      break ;
    }

    d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < min_d) {
      min_d = d;
      min_v = vno;
    }
  }
  if (dmin != NULL) {
    *dmin = min_d;
  }
  return (min_v);
}

/*!
  \fn double MRISfindMinDistanceVertexWithDotCheck(MRI_SURFACE *mris, double xs, double ys, double zs, MRI *mri, double dot_dir, int *pvtxno)
  \brief Finds the closest vertex to the given XYZ with the constraint that the vector from the vertex to the XYZ
  is (more-or-less) in the same (dot_dir +1) or opposite (dot_dir -1) direction as the surface normal. This is used
  to assure that the closest vertex is not on the other side of a mesh and is only closest because the truly closest
  vertex is in a big triangle. This goes through all the non-ripped vertices, so it can be very slow. If dot_dir=0,
  then it is the same as a brute force search to find the closest vertex. 
 */
double MRISfindMinDistanceVertexWithDotCheckXYZ(MRI_SURFACE *mris, double xs, double ys, double zs, MRI *mri, double dot_dir, int *pvtxno)
{
  int     vno, min_vno ;
  VERTEX  *v ;
  double  dist, dot, min_dist, dx, dy, dz ;

  min_vno = -1 ; min_dist = 1e10;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    v = &mris->vertices[vno] ;
    if (v->ripflag) continue ;
    dx = xs-v->x ; dy = ys-v->y ; dz = zs-v->z ;
    dot = v->nx*dx + v->ny*dy + v->nz*dz ;
    if (dot*dot_dir < 0) continue ;
    dist = sqrt(SQR(xs-v->x) + SQR(ys-v->y) + SQR(zs-v->z)) ;
    if (dist < min_dist){
      min_dist = dist ;
      min_vno = vno ;
    }
  }
  *pvtxno = min_vno ;
  return(min_dist);
}

/*!
  \fn double MRISfindMinDistanceVertexWithDotCheck(MRI_SURFACE *mris, int c, int r, int s, MRI *mri, double dot_dir, int *pvtxno)
  \brief Finds the closest vertex to the given CRS with a vector dot
  product constraint.  See MRISfindMinDistanceVertexWithDotCheckXYZ()
  for more info. This function originally appeared in mri_aparc2aseg
  as mrisFindMinDistanceVertexWithDotCheck().
 */
double MRISfindMinDistanceVertexWithDotCheck(MRI_SURFACE *mris, int c, int r, int s, MRI *mri, double dot_dir, int *pvtxno)
{
  double  xs, ys, zs, min_dist;

  MRIvoxelToSurfaceRAS(mri, c, r, s, &xs, &ys, &zs);
  min_dist = MRISfindMinDistanceVertexWithDotCheckXYZ(mris, xs, ys, zs, mri, dot_dir, pvtxno);
  return(min_dist);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISfindClosestOriginalVertex(MRIS *mris, float x, float y, float z)
{
  int vno, min_v = -1;
  VERTEX *v;
  float d, min_d, dx, dy, dz;

  min_d = 10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (vno == 91007 || vno == 91814) {
      DiagBreak();
    }
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    dx = v->origx - x;
    dy = v->origy - y;
    dz = v->origz - z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < min_d) {
      min_d = d;
      min_v = vno;
    }
  }

  return (min_v);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISfindClosestCanonicalVertex(MRIS *mris, float x, float y, float z)
{
  int vno, min_v = -1;
  VERTEX *v;
  float d, min_d, dx, dy, dz;

  min_d = 10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    dx = v->cx - x;
    dy = v->cy - y;
    dz = v->cz - z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < min_d) {
      min_d = d;
      min_v = vno;
    }
  }

  return (min_v);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISfindClosestWhiteVertex(MRIS *mris, float x, float y, float z)
{
  int vno, min_v = -1;
  VERTEX *v;
  float d, min_d, dx, dy, dz;

  min_d = 10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    dx = v->whitex - x;
    dy = v->whitey - y;
    dz = v->whitez - z;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d < min_d) {
      min_d = d;
      min_v = vno;
    }
  }

  return (min_v);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisFindAllOverlappingFaces(MRIS *mris, MHT *mht, int fno, int *flist)
{
  double x0, x1, y0, y1, z0, z1, x, y, z;
  int i, n, m, total_found, all_faces[1000000], nfaces;
  EDGE edge1, edge2;
  FACE *f1, *f2;
  VERTEX *v;

  f1 = &mris->faces[fno];
  if (fno == Gdiag_no) {
    DiagBreak();
  }
  x0 = y0 = z0 = 100000.0;
  x1 = y1 = z1 = -x0;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[f1->v[n]];
    x = v->x;
    y = v->y;
    z = v->z;
    x0 = MIN(x, x0);
    y0 = MIN(y, y0);
    z0 = MIN(z, z0);
    x1 = MAX(x, x1);
    y1 = MAX(y, y1);
    z1 = MAX(z, z1);
  }

  nfaces = total_found = 0;
  flist[total_found++] = fno;
  
  {
    MHBT* prev_bucket = NULL;
    for (x = x0; x <= x1; x += 0.5) {
      for (y = y0; y <= y1; y += 0.5) {
        for (z = z0; z <= z1; z += 0.5) {
          MHBT *bucket = MHTacqBucket(mht, x, y, z);
          if (!bucket) {
            continue;
          }
          if (bucket == prev_bucket) {
            MHTrelBucket(&bucket);
            continue;
          }
          prev_bucket = bucket;
          MHB *bin;
          for (bin = bucket->bins, i = 0; i < bucket->nused; i++, bin++) {
            f2 = &mris->faces[bin->fno];
            if (!f2->ripflag) /* only add it once */
            {
              if (nfaces == 1000000) {
                ErrorExit(ERROR_BADPARM, "Too many faces");
              }
              all_faces[nfaces++] = bin->fno;
              f2->ripflag = 1;
            }
          }
          MHTrelBucket(&bucket);
        }
      }
    }
  }
      
  for (i = 0; i < nfaces; i++) /* reset ripflag */
  {
    mris->faces[all_faces[i]].ripflag = 0;
  }

  for (i = 0; i < nfaces; i++) {
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (all_faces[i] < 0) /* duplicate */
    {
      continue;
    }
    f2 = &mris->faces[all_faces[i]];
    if (all_faces[i] == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      edge1.vno1 = f1->v[n];
      edge1.vno2 = f1->v[n == VERTICES_PER_FACE - 1 ? 0 : n + 1];
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (f2->ripflag) {
          continue;
        }
        edge2.vno1 = f2->v[m];
        edge2.vno2 = f2->v[m == VERTICES_PER_FACE - 1 ? 0 : m + 1];
        if (edge2.vno1 == edge1.vno1 || edge2.vno1 == edge1.vno2 || edge2.vno2 == edge1.vno1 ||
            edge2.vno2 == edge1.vno2) {
          continue;
        }
        if (edgesIntersect(mris, &edge1, &edge2)) {
          f2->ripflag = 1; /* use ripflag as a mark */
          if (total_found == MAX_INT_FACES)
            ErrorExit(ERROR_BADPARM,
                      "Too many intersected faces for face %d (%d, %d, %d)",
                      fno,
                      f1->v[0],
                      f1->v[1],
                      f1->v[2]);
          flist[total_found++] = all_faces[i];
        }
      }
    }
  }

  for (i = 0; i < total_found; i++) {
    f1 = &mris->faces[flist[i]];
    f1->ripflag = 0;
  }

  return (total_found);
}

int MRIScountNegativeFaces(MRIS *mris)
{
  int fno, neg;
  FACE *face;

  for (neg = fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) continue;
    if (face->area < 0) neg++;
  }
  return (neg);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Scale a surface anisotropically.
  ------------------------------------------------------*/
int MRISprintTessellationStats(MRIS *mris, FILE *fp)
{
  double mean, dsigma, dmin, dmax;
  int vno = 0, vno2 = 0;

  vno = vno2 = 0;
  mean = MRIScomputeVertexSpacingStats(mris, &dsigma, &dmin, &dmax, &vno, &vno2, CURRENT_VERTICES);
  fprintf(fp,
          "vertex spacing %2.2f +- %2.2f (%2.2f-->%2.2f) "
          "(max @ vno %d --> %d)\n",
          mean,
          dsigma,
          dmin,
          dmax,
          vno,
          vno2);
  
  mean = MRIScomputeFaceAreaStats(mris, &dsigma, &dmin, &dmax);
  fprintf(fp, "face area %2.2f +- %2.2f (%2.2f-->%2.2f)\n", mean, dsigma, dmin, dmax);

  if (dmax > 20) {
    int n;
    float dist;

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dist = sqrt(SQR(vn->x - v->x) + SQR(vn->y - v->y) + SQR(vn->z - v->z));
      if (dist > 20) {
        fprintf(stdout, "\t%d --> %d = %2.1f mm\n", vno, vt->v[n], dist);
      }
    }
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
// get the vertex number
int MRIStalairachToVertex(MRIS *mris, double xt, double yt, double zt)
{
  int vno;
  double xw, yw, zw;

  TransformWithMatrix(mris->TalSRASToSRAS_, xt, yt, zt, &xw, &yw, &zw);

  vno = MRISfindClosestOriginalVertex(mris, xw, yw, zw);

  return (vno);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScanonicalToVertex(MRIS *mris, double phi, double theta)
{
  int vno;
  double xw, yw, zw;

  MRIScanonicalToWorld(mris, phi, theta, &xw, &yw, &zw);
  vno = MRISfindClosestCanonicalVertex(mris, xw, yw, zw);
  return (vno);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
float MRISdistanceToSurface(MRIS *mris, MHT *mht, float x0, float y0, float z0, float nx, float ny, float nz)
{
  double dist, len;


  len = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= len;
  ny /= len;
  nz /= len;
  for (dist = 0.0f; dist < 128; dist += .25) {
    if (mrisDirectionTriangleIntersection(mris, x0, y0, z0, nx, ny, nz, mht, &dist, -1)) {
      return (dist);
    }
  }

  return (0.0);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRISrmsTPHeight(MRIS *mris)
{
  int vno, i, total_nbrs;
  double avg_height, dot, nx, ny, nz, x, y, z, d;

  if (mris->status == MRIS_PLANE) {
    return (NO_ERROR);
  }

  avg_height = 0.0;
  total_nbrs = 0;
  mrisComputeTangentPlanes(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX          const * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }

    if (vertext->vtotal <= 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    nx = vertex->nx;
    ny = vertex->ny;
    nz = vertex->nz;

    for (i = 0; i < vertext->vtotal; i++) {
      VERTEX const * const vnb = &mris->vertices[vertext->v[i]];
      if (vnb->ripflag) {
        continue;
      }
      x = vnb->x - vertex->x;
      y = vnb->y - vertex->y;
      z = vnb->z - vertex->z;
      d = (x * x + y * y + z * z);
      if (FZERO(d)) {
        continue;
      }
      /*
        calculate the projection of this vertex onto the surface normal
      */
      dot = nx * x + ny * y + nz * z; /* height above TpS */
      avg_height += dot * dot / d;
      total_nbrs++;
    }
  }

  return (sqrt(avg_height / (double)total_nbrs));
}
/*!
  \fn double mrisRmsValError(MRIS *mris, MRI *mri)
  \brief Samples mri at each vertex and computes the difference
  between the sample and v->val. The diff is squared and summed (SSE);
  that is then divided by the number of vertices hit and sqrt taken to
  give RMS. Similar to mrisComputeIntensityError() which returns the
  simple SSE. No changes are made to the input surface structure
  unless RmsValErrorRecord==1 in which case v->valbak takes the
  value of the sampled MRI value and v->val2bak takes the value
  of the error = sampled - target (ie, v->val2bak - v->val). #RMS
*/
double mrisRmsValError(MRIS *mris, MRI *mri)
{
  if (debugNonDeterminism) {
    mris_print_hash(stdout, mris, "mrisRmsValError", "\n");
  }
  int vno, n; // xv, yv, zv;
  double val, total, delta, x, y, z;
  VERTEX *v;
  extern int RmsValErrorRecord;

  for (total = 0.0, n = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) 
      continue;
    n++;
    // Sample mri at vertex
    MRISvertexToVoxel(mris, v, mri, &x, &y, &z);
    MRIsampleVolume(mri, x, y, z, &val);
    delta = (val - v->val);
    total += delta * delta;
    if(RmsValErrorRecord){
      v->valbak = val;
      v->val2bak = delta;
    }
  }
  if (debugNonDeterminism) {
    mris_print_hash(stdout, mris, "mrisRmsValError", "\n");
    printf("mrisRmsValError() total = %f, n=%d\n",total,n);
  }
  return (sqrt(total / (double)n));
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScomputeAnalyticDistanceError(MRIS *mris, int which, FILE *fp)
{
  int vno, n, vtotal, ndists;
  int const * pv;
  float d, xd, yd, zd, circumference = 0.0f, angle, odist;
  double pct_orig, pct, mean, mean_orig_error, mean_error, smean_error, smean_orig_error;
  VECTOR *v1, *v2;

  v1 = VectorAlloc(3, MATRIX_REAL);
  v2 = VectorAlloc(3, MATRIX_REAL);

  mean_orig_error = mean_error = pct_orig = pct = mean = 0.0;
  smean_orig_error = smean_error = 0.0;
  ndists = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vtotal = vt->vtotal;
    switch (which) {
      default: /* don't really know what to do in other cases */
      case MRIS_PLANE:
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX const * const vn = &mris->vertices[*pv++];
          if (vn->ripflag) {
            continue;
          }
          xd = v->origx - vn->origx;
          yd = v->origy - vn->origy;
          zd = v->origz - vn->origz;
          d = xd * xd + yd * yd + zd * zd;
          odist = sqrt(d);
          mean_orig_error += fabs(v->dist_orig[n] - odist);
          mean_error += fabs(v->dist[n] - odist);
          smean_orig_error += (v->dist_orig[n] - odist);
          smean_error += (v->dist[n] - odist);
          mean += odist;
          ndists++;
        }
        break;
      case MRIS_PARAMETERIZED_SPHERE:
      case MRIS_SPHERE:
        VECTOR_LOAD(v1, v->origx, v->origy, v->origz); /* radius vector */
        if (FZERO(circumference))                      /* only calculate once */
        {
          circumference = M_PI * 2.0 * V3_LEN(v1);
        }
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX const * const vn = &mris->vertices[*pv++];
          if (vn->ripflag) {
            continue;
          }
          VECTOR_LOAD(v2, vn->origx, vn->origy, vn->origz); /* radius vector */
          angle = fabs(Vector3Angle(v1, v2));
          d = circumference * angle / (2.0 * M_PI);
          odist = d;
          mean_orig_error += fabs(v->dist_orig[n] - odist);
          mean_error += fabs(v->dist[n] - odist);
          smean_orig_error += (v->dist_orig[n] - odist);
          smean_error += (v->dist[n] - odist);
          mean += fabs(odist);
          ndists++;
        }
        break;
    }
  }

  mean /= (double)ndists;
  mean_error /= (double)ndists;
  mean_orig_error /= (double)ndists;
  smean_orig_error /= (double)ndists;
  smean_error /= (double)ndists;
  pct = mean_error / mean;
  pct_orig = mean_orig_error / mean;
  fprintf(stdout,
          "mean orig = %2.3f mm (%%%2.2f), final = %2.3f mm (%%%2.2f)\n",
          mean_orig_error,
          100.0 * pct_orig,
          mean_error,
          100.0 * pct);
  fprintf(stdout, "signed mean orig error = %2.3f, final mean error = %2.3f\n", smean_orig_error, smean_error);
  if (fp) {
    char *cp;
    float measured_error, disturb_pct;

    cp = getenv("FS_DISTURB_DISTANCES");
    if (cp) {
      disturb_pct = atof(cp);
    }
    else {
      disturb_pct = 0.0;
    }
    measured_error = MRISpercentDistanceError(mris);
    fprintf(fp,
            "%2.3f  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f\n",
            100.0f * (float)MRISvalidVertices(mris) / (float)mris->nvertices,
            disturb_pct,
            100.0 * pct_orig,
            mean_orig_error,
            100.0 * pct,
            mean_error,
            measured_error);
  }
  VectorFree(&v1);
  VectorFree(&v2);
  return (pct);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRISstoreAnalyticDistances(MRIS *mris, int which)
{
  int vno, n, vtotal;
  int const * pv;
  float d, xd, yd, zd, circumference = 0.0f, angle, odist;
  double pct_orig, pct, mean, mean_orig_error, mean_error, smean_error, smean_orig_error;
  VECTOR *v1, *v2;

  v1 = VectorAlloc(3, MATRIX_REAL);
  v2 = VectorAlloc(3, MATRIX_REAL);

  mean_orig_error = mean_error = pct_orig = pct = mean = 0.0;
  smean_orig_error = smean_error = 0.0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vtotal = vt->vtotal;
    switch (which) {
      default: /* don't really know what to do in other cases */
      case MRIS_PLANE:
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX * const vn = &mris->vertices[*pv++];
          if (vn->ripflag) {
            continue;
          }
          xd = v->origx - vn->origx;
          yd = v->origy - vn->origy;
          zd = v->origz - vn->origz;
          d = xd * xd + yd * yd + zd * zd;
          odist = sqrt(d);
          v->dist_orig[n] = odist;
        }
        break;
      case MRIS_PARAMETERIZED_SPHERE:
      case MRIS_SPHERE:
        VECTOR_LOAD(v1, v->origx, v->origy, v->origz); /* radius vector */
        if (FZERO(circumference))                      /* only calculate once */
        {
          circumference = M_PI * 2.0 * V3_LEN(v1);
        }
        for (pv = vt->v, n = 0; n < vtotal; n++) {
          VERTEX * const vn = &mris->vertices[*pv++];
          if (vn->ripflag) {
            continue;
          }
          VECTOR_LOAD(v2, vn->origx, vn->origy, vn->origz); /* radius vector */
          angle = fabs(Vector3Angle(v1, v2));
          d = circumference * angle / (2.0 * M_PI);
          v->dist_orig[n] = d;
        }
        break;
    }
  }

  VectorFree(&v1);
  VectorFree(&v2);
  return (pct);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISdisturbOriginalDistances(MRIS *mris, double max_pct)
{
  int vno, n;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    for (n = 0; n < vt->vtotal; n++) {
      v->dist_orig[n] *= (1.0 + randomNumber(-max_pct / 100.0f, max_pct / 100.0f));
    }
  }

  return (NO_ERROR);
}


/*
  use MARS code from Mert and Thomas to compute the distance at each point on the surface
  to the boundary of a label.

  The distances will be returned in the vertex->val field.
*/

#include "MARS_DT_Boundary.h"
int MRISdistanceTransform(MRIS *mris, LABEL *area, int mode)
{
  int *vertices, *vertNbrs, vno, max_nbrs, j, index;
  double *cost, *vertDists;

  cost = (double *)calloc(mris->nvertices, sizeof(double));
  if (cost == NULL) {
    ErrorExit(ERROR_NOMEMORY, "MRISdistanceTransform: could not allocate %d cost array\n", mris->nvertices);
  }
  vertices = (int *)calloc(mris->nvertices, sizeof(int));
  if (vertices == NULL) {
    ErrorExit(ERROR_NOMEMORY, "MRISdistanceTransform: could not allocate %d vertex array\n", mris->nvertices);
  }

  for (vno = max_nbrs = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    //    if (v->ripflag)
    // continue ;
    if (vt->vnum > max_nbrs) {
      max_nbrs = vt->vnum;
    }
  }
  vertNbrs = (int *)calloc(mris->nvertices * max_nbrs, sizeof(int));
  if (vertNbrs == NULL)
    ErrorExit(
        ERROR_NOMEMORY, "MRISdistanceTransform: could not allocate nbrs (%d x %d) array\n", mris->nvertices, max_nbrs);
  vertDists = (double *)calloc(mris->nvertices * max_nbrs, sizeof(double));
  if (vertDists == NULL)
    ErrorExit(
        ERROR_NOMEMORY, "MRISdistanceTransform: could not allocate dist (%d x %d) array\n", mris->nvertices, max_nbrs);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    //    if (v->ripflag)
    //      continue ;
    for (j = 0; j < vt->vnum; j++) {
      //      if (mris->vertices[v->v[j]].ripflag)
      // continue ;
      index = index_2D_array(j, vno, max_nbrs);
      vertDists[index] = v->dist[j];
      vertNbrs[index] = vt->v[j] + 1;  // nbrs is 1-based
    }
  }


  for (vno = 0; vno < area->n_points; vno++)
    if (area->lv[vno].vno >= 0 && area->lv[vno].deleted == 0) {
      vertices[area->lv[vno].vno] = 1;
    }
  if (mode == DTRANS_MODE_INSIDE)  // mark exterior and compute distance from it
  {
    for (vno = 0; vno < mris->nvertices; vno++) {
      vertices[vno] = !vertices[vno];
    }
  }

  // compute outside distances
  MARS_DT_Boundary(vertices, mris->nvertices, max_nbrs, vertNbrs, vertDists, cost);
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].val = cost[vno];  // distance to inside points
  }

  if (mode == DTRANS_MODE_SIGNED || mode == DTRANS_MODE_UNSIGNED) {
    // compute distance to interior points
    for (vno = 0; vno < mris->nvertices; vno++) {
      vertices[vno] = !vertices[vno];  // 1=in exterior, 0=in interior
    }
    MARS_DT_Boundary(vertices, mris->nvertices, max_nbrs, vertNbrs, vertDists, cost);
    for (vno = 0; vno < mris->nvertices; vno++)
      if (cost[vno] > 0) {
        mris->vertices[vno].val = mode == DTRANS_MODE_SIGNED ? -cost[vno] : cost[vno];
      }
  }

  free(vertices);
  free(cost);
  free(vertDists);
  free(vertNbrs);
  return (NO_ERROR);
}
/*
  assumes v->val has the distances for label1 and v->val2 has the distances
  for label 2, computed using MRISdistanceTransform
*/
double MRIScomputeHausdorffDistance(MRIS *mris, int mode)
{
  double hdist = 0, d, dist;
  int vno, n, num;

  // assume mode is symmetric mean for now
  for (num = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      // look for zero crossings in val
      if (vn->val * v->val < 0)  // different signs - compute location of 0 between them
      {
        // compute distance of '0' point from v to vn
        if (vn->val < 0) {
          dist = -vn->val / (-vn->val + v->val);  // fraction of way
          dist = 1 - dist;                        // measure from v not vn
        }
        else {
          dist = -v->val / (-v->val + vn->val);  // fraction of way
        }
        d = dist * fabs(v->val2) + (1 - dist) * fabs(vn->val2);  // compute val2 at this point
        num++;
        hdist += d;
      }

      // look for zero crossings in val2
      if (vn->val2 * v->val2 < 0)  // different signs - compute location of 0 between them
      {
        // compute distance of '0' point from v to vn
        if (vn->val2 < 0) {
          dist = -vn->val2 / (-vn->val2 + v->val2);  // fraction of way
          dist = 1 - dist;                           // measure from v not vn
        }
        else {
          dist = -v->val2 / (-v->val2 + vn->val2);  // fraction of way
        }
        num++;
        d = dist * fabs(v->val) + (1 - dist) * fabs(vn->val);  // compute val2 at this point
        hdist += d;
      }
    }
  }
  if (num > 0) {
    hdist /= num;
  }

  return (hdist);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute the ratio of the two principal curvatures and
  store it in the vertex->curv variable.
  ------------------------------------------------------*/
int MRISuseCurvatureRatio(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float min_curv, max_curv, k1, k2, curv;

  /*  MRIScomputeSecondFundamentalForm(mris) ;*/

  min_curv = 10000.0f;
  max_curv = -min_curv;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (fabs(v->k1) > fabs(v->k2)) {
      k1 = v->k1;
      k2 = v->k2;
    }
    else {
      k1 = v->k2;
      k2 = v->k1;
    }

    if (!FZERO(k2)) {
      curv = fabs(k1 / k2);
      if (curv < min_curv) {
        min_curv = curv;
      }
      if (curv > max_curv) {
        max_curv = curv;
      }

      v->curv = curv;
    }
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (fabs(v->k1) > fabs(v->k2)) {
      k1 = v->k1;
      k2 = v->k2;
    }
    else {
      k1 = v->k2;
      k2 = v->k1;
    }

    if (FZERO(k2)) {
      if (FZERO(k1)) {
        curv = 0.0;
      }
      else {
        curv = k1 < 0 ? min_curv : max_curv;
      }
      v->curv = curv;
    }
  }

  mris->min_curv = min_curv;
  mris->max_curv = max_curv;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute the contrast of the two principal curvatures and
  store it in the vertex->curv variable.
  ------------------------------------------------------*/
int MRISuseCurvatureContrast(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float min_curv, max_curv, k1, k2, curv, min_k;

  /*  MRIScomputeSecondFundamentalForm(mris) ;*/

  min_k = 10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->k1 < min_k) {
      min_k = v->k1;
    }
    if (v->k2 < min_k) {
      min_k = v->k2;
    }
  }

  min_curv = 10000.0f;
  max_curv = -min_curv;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->k1 > v->k2) {
      k1 = v->k1;
      k2 = v->k2;
    }
    else {
      k1 = v->k2;
      k2 = v->k1;
    }
    k1 -= min_k;
    k2 -= min_k;
    curv = (k1 - k2) / (k1 + k2);

    if (curv < min_curv) {
      min_curv = curv;
    }
    if (curv > max_curv) {
      max_curv = curv;
    }

    v->curv = curv;
  }

  mris->min_curv = min_curv;
  mris->max_curv = max_curv;
  return (NO_ERROR);
}


#include "mrishash.h"

#define MAX_DIST 10
int MRISmeasureDistanceBetweenSurfaces(MRIS *mris, MRIS *mris2, int signed_dist)
{
  int vno;
  VERTEX *v1, *v2;
  MRIS_HASH_TABLE *mht;
  double dx, dy, dz;

  mht = MHTcreateVertexTable_Resolution(mris2, CURRENT_VERTICES, MAX_DIST);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v1 = &mris->vertices[vno];
    if (v1->ripflag) {
      continue;
    }
    float min_dist;
    v2 = MHTfindClosestVertex2(mht, mris2, mris, v1, &min_dist);
    if (v2 == NULL) {
      v1->curv = MAX_DIST;
      continue;
    }
    dx = v1->x - v2->x;
    dy = v1->y - v2->y;
    dz = v1->z - v2->z;
    v1->curv = sqrt(dx * dx + dy * dy + dz * dz);
    if (signed_dist) {
      double dot;

      dot = dx * v1->nx + dy * v1->ny + dz * v1->nz;
      if (dot < 0) {
        v1->curv *= -1;
      }
    }
  }

  MHTfree(&mht);
  return (NO_ERROR);
}



/*!
  \bf int StuffVertexCoords(MRIS *surf, int vertexno, double p[3])
  \brief Stuff vertex xyz coords into a double array[3]
*/
int StuffVertexCoords(MRIS *surf, int vertexno, double p[3])
{
  p[0] = surf->vertices[vertexno].x;
  p[1] = surf->vertices[vertexno].y;
  p[2] = surf->vertices[vertexno].z;
  return(0);
}
/*!
  \bf int StuffFaceCoords(MRIS *surf, int vertexno, double p[3])
  \brief Stuff xyz coords of a given face corner into a double array[3]
*/
int StuffFaceCoords(MRIS *surf, int faceno, int cornerno, double p[3])
{
  int vertexno = surf->faces[faceno].v[cornerno];
  StuffVertexCoords(surf,vertexno,p);
  return(0);
}
/*!
  \bf double MinDistToTriangleBF(double p1[3], double p2[3], double p3[3], double ptest[3], double pmin[3], double dL)
  \brief Compute the RMS distance from a given point (ptest) to the
  triangle defined by points (p1,p2,p3). The closest point is returned
  in pmin. Uses brute force (BF) by barycentrically tessellating the
  triangle into dL spaced points. It is possible to do this
  analytically, but this was easy.
*/
double MinDistToTriangleBF(double p1[3], double p2[3], double p3[3], double ptest[3], double pmin[3], double dL)
{
  double l1, l2, l3,d,dmin;
  double r[3], pr[3];
  int k;

  dmin = 10e10;
  for(l1=0; l1 < 1; l1 += dL){
    for(l2=0; l2 < 1; l2 += dL){
      l3 = 1.0 - (l1+l2);
      if(l3 < 0.0 || l3 > 1.0) continue;
      d = 0;
      for(k=0; k<3; k++) {
	r[k] = l1*p1[k] + l2*p2[k] + l3*p3[k];// location of barycentric point
	pr[k] = ptest[k]-r[k]; //dist to test point for this k
	d += (pr[k]*pr[k]); // accumulate
      }
      d = sqrt(d);
      // check for min
      if(d<dmin) {
	dmin=d;
	for(k=0; k<3; k++) pmin[k] = r[k];
      }
    }
  }
  return(dmin);
}

/*!
  \bf int MRISdistanceBetweenSurfacesExact(MRIS *surf1, MRIS *surf2)
  \brief Computes the distance from a vertex on surf1 to the closest
  point on surf2. It works by finding the vertex on surf2 that is
  closest to the given vertex on surf1, then finds the closest point
  on the neighboring faces of the surf2 vertex. The original purpose
  of this function was to help assess changes in surface position due
  to some processing where the number of vertices might have changed
  or their position might have shifted but the surface itself might
  not have moved. In that case, the closest vertex may be 1mm or so
  away eventhough the surfaces are nearly on top of each other. This
  function should not be used to compute thickness because it does not
  assure that the closest point is in the normal direction. Currently
  is not truly "exact" in that the faces are tessellated into 25
  points controlled by dL), and the distance is computed as the
  closest of those 25. An exact solution is possible, just harder 
  to program. The difference is put in to the curv field of surf1.
  The closest point is put into targ{xyz}
*/
int MRISdistanceBetweenSurfacesExact(MRIS *surf1, MRIS *surf2)
{
  int vno1, nthface, faceno2,vno2,k; 
  VERTEX_TOPOLOGY *vt2;
  MRIS_HASH_TABLE *mht;
  float dminv;
  double d,dmin,dL=0.2;
  double pv1[3], pf1[3], pf2[3], pf3[3], pmin[3], pmin0[3];
  VERTEX *v1, *v2;

  mht = MHTcreateVertexTable_Resolution(surf2, CURRENT_VERTICES, 10);

  for (vno1 = 0; vno1 < surf1->nvertices; vno1++) {
    v1 = &surf1->vertices[vno1];
    if (v1->ripflag)  continue;

    // Get the coords for this vertex on surf1
    StuffVertexCoords(surf1, vno1, pv1);

    // Get the closest vertex in surf2
    vno2 = MHTfindClosestVertexNo2(mht, surf2, surf1, v1, &dminv);
    v2 = &surf2->vertices[vno2];    
    vt2 = &surf2->vertices_topology[vno2];
    pmin[0] = v2->x;
    pmin[1] = v2->y;
    pmin[2] = v2->z;

    if(vno1 == Gdiag_no){
      printf("MRISdistanceBetweenSurfacesExact() %6d  %6d  %d a = [%7.3f,%7.3f,%7.3f]; b=[%7.3f,%7.3f,%7.3f]; \n",
	     vno1,vno2,vt2->num, v1->x,v1->y,v1->z,v2->x,v2->y,v2->z);
    }

    // Go through each face of this vertex
    dmin = dminv;
    for(nthface = 0; nthface < vt2->num; nthface++){
      faceno2 = vt2->f[nthface];
      StuffFaceCoords(surf2, faceno2, 0, pf1);
      StuffFaceCoords(surf2, faceno2, 1, pf2);
      StuffFaceCoords(surf2, faceno2, 2, pf3);
      // Find the closest point on the surface from vertex1 to this face
      d = MinDistToTriangleBF(pf1,pf2,pf3,pv1,pmin0,dL);
      if(d<dmin) {
	dmin = d;
	for(k=0; k<3; k++) pmin[k] = pmin0[k];
      }
    }
    //if(vno1 == Gdiag_no || vno2 == Gdiag_no) printf("%6d %6d  %6.4f %6.4f\n",vno1,vno2,dminv,dmin);

    v1->curv = dmin;
    v1->targx = pmin[0];
    v1->targy = pmin[1];
    v1->targz = pmin[2];
  }

  MHTfree(&mht);
  return (NO_ERROR);
}

short FACES_Hcurvature_determineSign(MRIS *apmris,
    	int			apFACE_O_fno,
    	int			apFACE_I_fno)
{
  FACE * const apFACE_O = &apmris->faces[apFACE_O_fno];
  FACE * const apFACE_I = &apmris->faces[apFACE_I_fno];
  
  FaceNormCacheEntry const * const apFNorm_O = getFaceNorm(apmris, apFACE_O_fno);
  FaceNormCacheEntry const * const apFNorm_I = getFaceNorm(apmris, apFACE_I_fno);

  //
  // PRECONDITIONS
  //  o Typically called from MRIS_Hcurvature_determineSign()
  //  o apFACE_I and apFACE_J are geometric neighbouring faces
  //    about a given vertex.
  //
  // POSTCONDITIONS
  //  o If the faces "diverge", i.e. have face normals that point "away"
  //    from each other, then the curvature sign return is -1.
  //  o If the faces "converge", i.e. have face normals that point "toward"
  //    each other, then the curvature sign return is +1
  //
  // HISTORY
  //  02 July 2007
  //  o Initial design and coding.
  //

  static int calls = 0;
  const char *pch_function = "FACES_Hcurvature_determineSign";
  int ret = 0;
  int vertexO = -1;
  int vertexI = -1;
  VERTEX *pVERTEX_O = NULL;
  VERTEX *pVERTEX_I = NULL;
  VECTOR *pv_O = NULL;               // Coords of 1st common vertex
  VECTOR *pv_normalO = NULL;         // Normal for face O
  VECTOR *pv_OnO = NULL;             // pv_O + normal
  VECTOR *pv_I = NULL;               // Coords of 2nd common vertex
  VECTOR *pv_normalI = NULL;         // Normal for face I
  VECTOR *pv_InI = NULL;             // pv_I + normal
  VECTOR *pv_connectOI = NULL;       // vector connecting normals
  VECTOR *pv_commonVertices = NULL;  // Vector housing vertices that
  // are common between two
  // neighbouring faces.
  int commonVertices = 0;  // number of vertices in common
  // between two faces
  float f_distNormal = 0;
  float f_distAntiNormal = 0;
  int sign = -1;

  pv_commonVertices = VectorAlloc(3, MATRIX_REAL);
  pv_I = VectorAlloc(3, MATRIX_REAL);
  pv_O = VectorAlloc(3, MATRIX_REAL);
  pv_normalI = VectorAlloc(3, MATRIX_REAL);
  pv_normalO = VectorAlloc(3, MATRIX_REAL);
  pv_InI = VectorAlloc(3, MATRIX_REAL);
  pv_OnO = VectorAlloc(3, MATRIX_REAL);
  pv_connectOI = VectorAlloc(3, MATRIX_REAL);

  commonVertices = VERTICES_commonInFaces_find(apFACE_O, apFACE_I, pv_commonVertices);
  if (commonVertices != 2)
    ErrorExit(-4, "%s: During call %d, the passed faces do not share an edge", pch_function, calls);

  vertexO = FACE_vertexIndexAtMask_find(apFACE_O, pv_commonVertices);
  vertexI = FACE_vertexIndexAtMask_find(apFACE_I, pv_commonVertices);
  pVERTEX_O = &apmris->vertices[vertexO];
  pVERTEX_I = &apmris->vertices[vertexI];
  V3_LOAD(pv_O, pVERTEX_O->x, pVERTEX_O->y, pVERTEX_O->z);
  V3_LOAD(pv_I, pVERTEX_I->x, pVERTEX_I->y, pVERTEX_I->z);
  V3_LOAD(pv_normalO, apFNorm_O->nx, apFNorm_O->ny, apFNorm_O->nz);
  V3_LOAD(pv_normalI, apFNorm_I->nx, apFNorm_I->ny, apFNorm_I->nz);

  for (sign = 1; sign >= -1; sign -= 2) {
    V3_SCALAR_MUL(pv_normalO, sign, pv_normalO);
    V3_SCALAR_MUL(pv_normalI, sign, pv_normalI);
    V3_ADD(pv_O, pv_normalO, pv_OnO);
    V3_ADD(pv_I, pv_normalI, pv_InI);
    V3_SUBTRACT(pv_OnO, pv_InI, pv_connectOI);

    if (sign == 1) {
      f_distNormal = V3_LEN(pv_connectOI);
    }
    else {
      f_distAntiNormal = V3_LEN(pv_connectOI);
    }
  }

  ret = (f_distNormal < f_distAntiNormal) ? +1 : -1;

  VectorFree(&pv_O);
  VectorFree(&pv_normalO);
  VectorFree(&pv_OnO);
  VectorFree(&pv_I);
  VectorFree(&pv_normalI);
  VectorFree(&pv_InI);
  VectorFree(&pv_connectOI);
  VectorFree(&pv_commonVertices);
  calls++;
  return ret;
}

int VERTEX_faceAngles_determine(MRIS *apmris, int avertex, VECTOR *apv_angle)
{
  //
  // PRECONDITIONS
  //  o <apmris> is a valid surface.
  //  o <apVERTEX> is a vertex to analyze.
  //
  // POSTCONDITIONS
  //  o The angle between each face in the <apex> normal
  //    is determined and returned in <apv_angle>.
  //  o The caller is responsible for clearing the memory allocated to
  //    <apv_angle>!
  //  o The number of faces processed (i.e. size of the <apv_angle
  //    vector) is returned in the function name.
  //
  // HISTORY
  //  30 July 2007
  //  o Initial design and coding.
  //

  const char *pch_function = "VERTEX_faceAngles_determine";
  int nfaces = -1;
  float f_angle = 0.;
  float f_lenApexNormal = 0.;
  float f_lenFaceNormal = 0.;
  float f_lenNormals = 0.;
  float f_acosArg = 0.;
  float f_dot = 0.;

  int face = 0;
  static int calls = 0;
  static VECTOR *pv_faceNormal = NULL;
  static VECTOR *pv_apexNormal = NULL;

  FACE *pFACE_side = NULL;

  if (!calls) {
    pv_faceNormal = VectorAlloc(3, MATRIX_REAL);
    pv_apexNormal = VectorAlloc(3, MATRIX_REAL);
  }
  VERTEX_TOPOLOGY const * const pVERTEXt_apex = &apmris->vertices_topology[avertex];
  VERTEX          const * const pVERTEX_apex  = &apmris->vertices         [avertex];
  nfaces = pVERTEXt_apex->num;
  VECTOR_ELT(pv_apexNormal, 1) = pVERTEX_apex->nx;
  VECTOR_ELT(pv_apexNormal, 2) = pVERTEX_apex->ny;
  VECTOR_ELT(pv_apexNormal, 3) = pVERTEX_apex->nz;
  f_lenApexNormal = V3_LEN(pv_apexNormal);

  for (face = 0; face < nfaces; face++) {
    pFACE_side = &apmris->faces[pVERTEXt_apex->f[face]];
    FaceNormCacheEntry const * const pFNorm_side = getFaceNorm(apmris, pVERTEXt_apex->f[face]);
    VECTOR_ELT(pv_faceNormal, 1) = pFNorm_side->nx;
    VECTOR_ELT(pv_faceNormal, 2) = pFNorm_side->ny;
    VECTOR_ELT(pv_faceNormal, 3) = pFNorm_side->nz;
    f_lenFaceNormal = V3_LEN(pv_faceNormal);
    f_lenNormals = f_lenApexNormal * f_lenFaceNormal;
    f_dot = V3_DOT(pv_apexNormal, pv_faceNormal);
    errno = 0;
    //  feclearexcept(FE_ALL_EXCEPT);
    f_acosArg = f_dot / f_lenNormals;
    // Check on the bounds of the acos argument. Without this bounds check,
    // it is quite possible to have 'nan' acos results, especially on 64-bit
    // builds.
    if (f_acosArg > 1.) {
      f_acosArg = 1.0;
    }
    if (f_acosArg < -1.) {
      f_acosArg = -1.0;
    }
    f_angle = acos(f_acosArg);
    if (errno) {
      f_angle = 0.;
      printf("%s: acos error - angle set to zero for vertex = %d, face = %d.\n", pch_function, avertex, face);
    }
    VECTOR_ELT(apv_angle, face + 1) = f_angle;
  }
  calls++;
  return face;
}

float FACES_angleNormal_find(MRIS *apmris, int apFACE_I_fno, int apFACE_J_fno)
{
  //FACE * const apFACE_I = &apmris->faces[apFACE_I_fno];
  //FACE * const apFACE_J = &apmris->faces[apFACE_J_fno];

  FaceNormCacheEntry const * const apFNorm_I = getFaceNorm(apmris, apFACE_I_fno); 
  FaceNormCacheEntry const * const apFNorm_J = getFaceNorm(apmris, apFACE_J_fno);
 
  //
  // PRECONDITIONS
  //  o The <apFACE>s should be triangles with 3 vertices each.
  //  o It is assumed (but not mandatory) that the FACES share
  //    at least one common vertex - or more often a common
  //    edge, i.e. two common vertices.
  //
  // POSTCONDITIONS
  //  o The angle between the normals on each FACE is returned.
  //

  static int calls = 0;                  // Used for vector allocation
  static VECTOR *pv_faceNormalI = NULL;  // Normal vector for face I
  static VECTOR *pv_faceNormalJ = NULL;  // Normal vector for face J
  static VECTOR *pv_crossIJ = NULL;      // Cross product of input vectors
  float f_faceNormalIlen = 0.;           // Length of face normal I
  float f_faceNormalJlen = 0.;           // Length of face normal J
  float f_faceNormalIJlen = 0.;          // Face normal length product
  float f_angleNormalIJ = 0.;            // Angle between face normals
  float f_acosArg = 0.;                  // Dot product arguments
  float f_dot = 0.;                      // Dot product
  short sign = 1;                        // Angle "sign"
  const char *pch_function = "FACES_angleNormal_find";

  DebugEnterFunction(("%s", pch_function));
  if (!calls) {
    pv_faceNormalI = VectorAlloc(3, MATRIX_REAL);
    pv_faceNormalJ = VectorAlloc(3, MATRIX_REAL);
    pv_crossIJ = VectorAlloc(3, MATRIX_REAL);
  }
  V3_LOAD(pv_faceNormalI, apFNorm_I->nx, apFNorm_I->ny, apFNorm_I->nz);
  V3_LOAD(pv_faceNormalJ, apFNorm_J->nx, apFNorm_J->ny, apFNorm_J->nz);
  f_faceNormalIlen = V3_LEN(pv_faceNormalI);
  f_faceNormalJlen = V3_LEN(pv_faceNormalJ);
  if (f_faceNormalIlen > 1.0001 || f_faceNormalJlen > 1.0001)
    ErrorExit(
        -4, "%s: face normal not unit length -- Ni: %f, Nj: %f\n", pch_function, f_faceNormalIlen, f_faceNormalJlen);
  f_faceNormalIJlen = f_faceNormalIlen * f_faceNormalJlen;
  f_dot = V3_DOT(pv_faceNormalI, pv_faceNormalJ);
  sign = FACES_Hcurvature_determineSign(apmris, apFACE_I_fno, apFACE_J_fno);
  f_acosArg = f_dot / f_faceNormalIJlen;
  // Check on the bounds of the acos argument. Without this bounds check,
  //  it is quite possible to have 'nan' acos results, especially on 64-bit
  //  builds.
  if (f_acosArg > 1.) {
    f_acosArg = 1.0;
  }
  if (f_acosArg < -1.) {
    f_acosArg = -1.0;
  }
  f_angleNormalIJ = acosf(f_acosArg) * sign;
  calls++;
  xDbg_PopStack();
  return f_angleNormalIJ;
}


int VERTEX_faceMinMaxAngles_determine(
    MRIS *apmris, int avertex, int *ap_minIndex, float *apf_minAngle, int *ap_maxIndex, float *apf_maxAngle)
{
  //
  // PRECONDITIONS
  //  o <apmris> is a valid surface.
  //  o <apVERTEX> is a vertex to analyze.
  //
  // POSTCONDITIONS
  //  o For the given <avertex>, the minimum and maximum
  //    face angles and their indices are returned in the
  //    argument pointers.
  //  o The number of faces is returned in the function name.
  //
  // HISTORY
  //  31 July 2007
  //  o Initial design and coding.
  //

  const char *pch_function = "VERTEX_faceMinMaxAngles_determine";
  int face = 0;                  // Face index counte
  int nfaces = 0;                // Number of faces at <avertex>
  float f_faceAngle = 0.;        // Actual face angle
  VECTOR *pv_faceAngles = NULL;  // vector containing angles
  // between each face normal
  // and apex normal

  // Determine the angles between each face and the vertex normal;
  //  find the min/max angles and indices
  VERTEX_TOPOLOGY const * const pVERTEXt = &apmris->vertices_topology[avertex];
  nfaces = pVERTEXt->num;
  pv_faceAngles = VectorAlloc(nfaces, MATRIX_REAL);
  nfaces = VERTEX_faceAngles_determine(apmris, avertex, pv_faceAngles);
  if (!nfaces) {
    ErrorExit(-4, "%s: error with determining face angles.", pch_function);
  }
  f_faceAngle = VECTOR_ELT(pv_faceAngles, 1);
  *apf_minAngle = f_faceAngle;
  *apf_maxAngle = f_faceAngle;
  for (face = 1; face < nfaces; face++) {
    f_faceAngle = VECTOR_ELT(pv_faceAngles, face + 1);  // base 1 index
    if (f_faceAngle < *apf_minAngle) {
      *apf_minAngle = f_faceAngle;
      *ap_minIndex = face;
    }
    if (f_faceAngle > *apf_maxAngle) {
      *apf_maxAngle = f_faceAngle;
      *ap_maxIndex = face;
    }
  }
  VectorFree(&pv_faceAngles);
  return face;
}

int MRIS_Hcurvature_determineSign(MRIS *apmris)
{
  //
  // NOTE
  //  This function is obsolete and should not be used! Mean curvature (H)
  //  determination is now done directly when processing faces and normal
  //  angles.
  //
  // PRECONDITIONS
  //  o <apmris> is a valid surface.
  //  o MRIS_facesAtVertices_reorder()
  //
  // POSTCONDITIONS
  //  o The face geometry at each vertex is examined, and the orientation
  //    of each face relative to its neighbor is determined as either
  //    converging (-) or diverging (+).
  //  o If the sum of each convergence/divergence is determined to be
  //    negative, the vertex is marked as diverging; otherwise it is
  //    marked as converging.
  //  o Convergence is indicated with pVERTEX->undefval=-1; divergence
  //    is marked with pVERTEX->undefval=1
  //
  // HISTORY
  //  02 July 2007
  //  o Initial design and coding.
  //

  int vertex = 0;
  int face = 0;
  int nfaces = 0;
  int ret = 1;
  int signSum = 0;

  for (vertex = 0; vertex < apmris->nvertices; vertex++) {
    MRIS_vertexProgress_print(apmris, vertex, "Determining H sign for vertex faces...");
    VERTEX_TOPOLOGY const * const pVERTEXt = &apmris->vertices_topology[vertex];
    VERTEX                * const pVERTEX  = &apmris->vertices         [vertex];
    nfaces = pVERTEXt->num;
    signSum = 0;
    for (face = 0; face < nfaces; face++) {
      signSum += FACES_Hcurvature_determineSign(apmris, pVERTEXt->f[face], pVERTEXt->f[(face + 1) % nfaces]);
    }
    pVERTEX->undefval = (signSum >= 0) ? 1 : -1;
  }
  return ret;
}

short MRIS_discreteKH_compute(MRIS *apmris)
{
  //
  // PRECONDITIONS
  //  o A valid SURFACE with computed triangle properties.
  //
  // POSTCONDITIONS
  //  o The discrete K and H curvatures at each vertex are
  //    computed.
  //

  const char *pch_function = "MRIS_discreteKH_compute";

  VECTOR *pv_geometricOrder = NULL;  // Geometrically ordered faces
  int vertex = 0;                    // Vertex index number
  FACE *pFACE_I = NULL;              // Face I with vertex apex
  FACE *pFACE_J = NULL;              // Face I+1 with vertex apex
  int face = 0;                      // face counter
  int faceI = 0;                     // face I index
  int faceJ = 0;                     // face J index
  int nfaces = 0;                    // total number of faces
  int *pFaceIndex = NULL;            // face index array at vertex
  int angleIndex = -1;               // angle index
  float f_faceAreaSum = 0.;          // area about vertex
  float f_angleDeficitSum = 0.;      // angle deficit about vertex
  float f_angleNormalIJ = 0.;        // angle between normals
  float f_angleNormalIJSum = 0.;     // sum angle between normals
  float f_edgeLength = 0.;           // Length of edge v->vI
  double f_K = 0.;                   // Gaussian curvature at vertex
  double f_H = 0.;                   // Mean curvature at vertex
  float f_Kmin = 0.;
  float f_Kmax = 0.;
  float f_Hmin = 0.;
  float f_Hmax = 0.;
  float f_Ktotal = 0.;

  DebugEnterFunction(("%s", pch_function));
  for (vertex = 0; vertex < apmris->nvertices; vertex++) {
    MRIS_vertexProgress_print(apmris, vertex, "Determining KH curvatures...");
    f_faceAreaSum = 0.;
    f_angleDeficitSum = 0.;
    f_angleNormalIJSum = 0.;
    
    VERTEX_TOPOLOGY const * const pVertext = &apmris->vertices_topology[vertex];
    VERTEX          const * const pVertex  = &apmris->vertices         [vertex];
    if (pVertex->marked)  // couldn't find geometrical packing of faces - can't process this vertex
    {
      apmris->vertices[vertex].K = 0;
      apmris->vertices[vertex].H = 0;
      continue;
    }

    nfaces = pVertext->num;
    pFaceIndex = pVertext->f;
    pv_geometricOrder = VectorAlloc(nfaces, MATRIX_REAL);
    for (face = 0; face < nfaces; face++) {
      faceI = face;
      faceJ = (face + 1) % nfaces;

      pFACE_I = &apmris->faces[pFaceIndex[faceI]];
      pFACE_J = &apmris->faces[pFaceIndex[faceJ]];

      f_angleNormalIJ = FACES_angleNormal_find  (apmris, pFaceIndex[faceI], pFaceIndex[faceJ]);
      f_edgeLength = FACES_commonEdgeLength_find(apmris, pFACE_I, pFACE_J);
      f_faceAreaSum += pFACE_I->area;
      angleIndex = FACE_vertexIndex_find(pFACE_I, vertex);
      if (angleIndex == -1)
        ErrorExit(-4, "%s:\n\tangleIndex lookup failure for vertex %d, face %d", pch_function, vertex, face);
      f_angleDeficitSum += pFACE_I->angle[angleIndex];
      f_angleNormalIJSum += f_angleNormalIJ * f_edgeLength;
    }
    VectorFree(&pv_geometricOrder);
    pv_geometricOrder = NULL;
    f_K = 3 / f_faceAreaSum * (2 * M_PI - f_angleDeficitSum);
    f_H = 0.75 / f_faceAreaSum * f_angleNormalIJSum;
    apmris->vertices[vertex].K = f_K;
    apmris->vertices[vertex].H = f_H;
    if (!vertex) {
      f_Kmin = f_Kmax = f_K;
      f_Hmin = f_Hmax = f_H;
    }
    if (f_K > f_Kmax) {
      f_Kmax = f_K;
    }
    if (f_K < f_Kmin) {
      f_Kmin = f_K;
    }
    if (f_H > f_Hmax) {
      f_Hmax = f_H;
    }
    if (f_H < f_Hmin) {
      f_Hmin = f_H;
    }
    f_Ktotal += f_K * pVertex->area;
  }
  apmris->Kmax = f_Kmax;
  apmris->Kmin = f_Kmin;
  apmris->Hmax = f_Hmax;
  apmris->Hmin = f_Hmin;
  apmris->Ktotal = f_Ktotal;
  xDbg_PopStack();
  return (NO_ERROR);
}

short MRIS_discretek1k2_compute(MRIS *apmris, short ab_signedPrinciples)
{
  //
  // PRECONDITIONS
  //  o A valid SURFACE with computed triangle properties.
  //  o A valid K and H at each vertex.
  //
  // POSTCONDITIONS
  //  o The discrete K and H curvatures at each vertex are
  //    computed.
  //  o If <ab_signedPrinciples> is true, then k1 and k2
  //    are assigned according to signed size, and not
  //    f_abs(..) size.
  //

  const char *pch_function = "MRIS_discretek1k2_compute";
  VERTEX *pVERTEX = NULL;
  float f_k1 = 0.;
  float f_k2 = 0.;
  float f_A = 0.;
  float f_B = 0.;
  float f_delta = 0.;
  float f_K = 0.;
  float f_H = 0.;
  int vertex = 0;
  int deltaViolations = 0;

  for (vertex = 0; vertex < apmris->nvertices; vertex++) {
    MRIS_vertexProgress_print(apmris, vertex, "Determining k1k2 curvatures...");
    pVERTEX = &apmris->vertices[vertex];
    f_K = pVERTEX->K;
    f_H = pVERTEX->H;
    f_delta = f_H * f_H - f_K;
    if (f_delta < 0) {
      deltaViolations++;
      f_delta = 0.;
    }
    if (f_delta < 0)
      ErrorExit(-4, "%s: f_delta = %f, vertex = %d, f_K = %f, f_H = %f\n", pch_function, f_delta, vertex, f_K, f_H);
    f_A = f_H + sqrt(f_delta);
    f_B = f_H - sqrt(f_delta);
    if (!ab_signedPrinciples) {
      f_k1 = fabs(f_A) >= fabs(f_B) ? f_A : f_B;
      f_k2 = fabs(f_A) <= fabs(f_B) ? f_A : f_B;
    }
    else {
      f_k1 = f_A >= f_B ? f_A : f_B;
      f_k2 = f_A <= f_B ? f_A : f_B;
    }
    pVERTEX->k1 = f_k1;
    pVERTEX->k2 = f_k2;
  }
  if (deltaViolations) {
    cprintd("deltaViolations", deltaViolations);
  }
  return (NO_ERROR);
}

short MRIScomputeSecondFundamentalFormDiscrete(MRIS *apmris, short ab_signedPrinciples)
{
  int retKH, retk1k2;

  retKH = 1;
  retk1k2 = 1;

  MRISclearMarks(apmris);   // not clear this is needed before the next line
  MRIScomputeTriangleProperties(apmris);

  // maybe this should have been here MRISclearMarks(apmris);
  MRIS_facesAtVertices_reorder(apmris);
  // this leaves marked set - don't know if this is required for the next step
  
  retKH = MRIS_discreteKH_compute(apmris);
  retk1k2 = MRIS_discretek1k2_compute(apmris, ab_signedPrinciples);
  return (retKH | retk1k2);
}


int MRIScomputeSurfaceNormals(MRIS *mris, int which, int navgs)
{
  int vno, n, i;
  double nx = 0.0, ny = 0.0, nz = 0.0, norm = 0.0;

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRISrestoreVertexPositions(mris, which);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v  = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    MRIScomputeNormal(mris, which, vno, &nx, &ny, &nz);
    switch (which) {
      default:
      case CURRENT_VERTICES:
        v->nx = nx;
        v->ny = ny;
        v->nz = nz;
        break;
      case WHITE_VERTICES:
        v->wnx = nx;
        v->wny = ny;
        v->wnz = nz;
        break;
      case PIAL_VERTICES:
        v->pnx = nx;
        v->pny = ny;
        v->pnz = nz;
        break;
      case ORIGINAL_VERTICES:
        v->onx = nx;
        v->ony = ny;
        v->onz = nz;
        break;
    }
  }
  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      switch (which) {
        default:
        case CURRENT_VERTICES:
          nx = v->nx;
          ny = v->ny;
          nz = v->nz;
          break;
        case WHITE_VERTICES:
          nx = v->wnx;
          ny = v->wny;
          nz = v->wnz;
          break;
        case PIAL_VERTICES:
          nx = v->pnx;
          ny = v->pny;
          nz = v->pnz;
          break;
        case ORIGINAL_VERTICES:
          nx = v->onx;
          ny = v->ony;
          nz = v->onz;
          break;
      }
      v->tdx = nx;
      v->tdy = ny;
      v->tdz = nz;
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        if (vn->ripflag) {
          continue;
        }
        switch (which) {
          default:
          case CURRENT_VERTICES:
            nx = vn->nx;
            ny = vn->ny;
            nz = vn->nz;
            break;
          case WHITE_VERTICES:
            nx = vn->wnx;
            ny = vn->wny;
            nz = vn->wnz;
            break;
          case PIAL_VERTICES:
            nx = vn->pnx;
            ny = vn->pny;
            nz = vn->pnz;
            break;
          case ORIGINAL_VERTICES:
            nx = vn->onx;
            ny = vn->ony;
            nz = vn->onz;
            break;
        }
        v->tdx += nx;
        v->tdy += ny;
        v->tdz += nz;
      }
      norm = sqrt(v->tdx * v->tdx + v->tdy * v->tdy + v->tdz * v->tdz);
      if (!FZERO(norm)) {
        v->tdx /= norm;
        v->tdy /= norm;
        v->tdz /= norm;
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      switch (which) {
        default:
        case CURRENT_VERTICES:
          v->nx = v->tdx;
          v->ny = v->tdy;
          v->nz = v->tdz;
          break;
        case WHITE_VERTICES:
          v->wnx = v->tdx;
          v->wny = v->tdy;
          v->wnz = v->tdz;
          break;
        case PIAL_VERTICES:
          v->pnx = v->tdx;
          v->pny = v->tdy;
          v->pnz = v->tdz;
          break;
        case ORIGINAL_VERTICES:
          v->onx = v->tdx;
          v->ony = v->tdy;
          v->onz = v->tdz;
          break;
      }
    }
  }

  MRISrestoreVertexPositions(mris, which);
  return (NO_ERROR);
}


int MRISsetCurvature(MRIS *mris, float val)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = val;
  }
  return (NO_ERROR);
}

int MRIStrinarizeCurvature(MRIS *mris, float trinarize_thresh)
{
  int vno;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->curv < -trinarize_thresh) {
      v->curv = -1;
    }
    else if (v->curv > trinarize_thresh) {
      v->curv = 1;
    }
    else {
      v->curv = 0;
    }
  }
  return (NO_ERROR);
}


int MRISscaleCurvature(MRIS *apmris, float af_scale)
{
  //
  // POSTCONDITIONS
  // o Each curvature value in apmris is scaled by:
  //
  //    (curv-f_mean)*<af_scale> + f_mean
  //
  //   where f_mean is the mean of all the surface curvatures
  //

  VERTEX *pvertex;
  int vno;
  int vtotal;
  double f_mean;

  for (f_mean = 0.0, vtotal = vno = 0; vno < apmris->nvertices; vno++) {
    pvertex = &apmris->vertices[vno];
    if (pvertex->ripflag) {
      continue;
    }
    vtotal++;
    f_mean += pvertex->curv;
  }
  f_mean /= (double)vtotal;

  for (vno = 0; vno < apmris->nvertices; vno++) {
    pvertex = &apmris->vertices[vno];
    if (pvertex->ripflag) {
      continue;
    }
    pvertex->curv = (pvertex->curv - f_mean) * af_scale + f_mean;
  }
  return (NO_ERROR);
}


int MRISbinarizeCurvature(MRIS *mris, float thresh, float low, float high, int use_abs)
{
  int vno;
  VERTEX *v;
  float val;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    val = (use_abs ? fabs(v->curv) : v->curv);
    if (val < thresh) {
      v->curv = low;
    }
    else {
      v->curv = v->curv < 0 ? -high : high;
    }
  }

  return (NO_ERROR);
}
int MRISthresholdCurvature(MRIS *mris, float thresh, int use_abs)
{
  int vno;
  VERTEX *v;
  float val;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    val = (use_abs ? fabs(v->curv) : v->curv);
    if (val < thresh) {
      v->curv = 0;
    }
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

void MRISnormalizeField(MRIS *mris, int distance_field, int which_norm)
{
  int n;
  VERTEX *v;
  float max_value;

  if (which_norm == NORM_NONE) {
    return;
  }

  if (distance_field || which_norm == NORM_MAX) {
    /* distance fields have their max value set to 1 */
    for (max_value = 0.0f, n = 0; n < mris->nvertices; n++) {
      v = &mris->vertices[n];
      if (v->ripflag) {
        continue;
      }
      if (v->curv > max_value) {
        max_value = v->curv;  // fabs(v->curv);
      }
    }
    if (!FZERO(max_value)) /* normalize the max value to 1.0f */
    {
      printf("normalizing max %2.1f to be in [0 1]\n", max_value);
      for (n = 0; n < mris->nvertices; n++) {
        v = &mris->vertices[n];
        if (v->ripflag) {
          continue;
        }
        v->curv = v->curv / max_value;
      }
    }
  }
  else /* gaussian */
  {
    MRISnormalizeCurvature(mris, which_norm);
  }
}


// remove outliers in the v->curv field
int MRIShistoThresholdCurvature(MRIS *mris, float thresh_pct)
{
  double min_curv, max_curv, curv_scale, total, thresh;
  int bin, zbin, nthresh = 0, vno;
  HISTOGRAM *h_curv;
  VERTEX *vertex;

  h_curv = HISTOalloc(1000);

  min_curv = 1000;
  max_curv = -1000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    if (vertex->curv > max_curv) {
      max_curv = vertex->curv;
    }
    if (vertex->curv < min_curv) {
      min_curv = vertex->curv;
    }
  }

  curv_scale = (h_curv->nbins - 1) / (max_curv - min_curv);
  h_curv->bin_size = 1.0 / curv_scale;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    bin = nint(curv_scale * (vertex->curv - min_curv));
    h_curv->counts[bin]++;
  }
  for (bin = 0; bin < h_curv->nbins; bin++) {
    h_curv->bins[bin] = (bin / curv_scale) + min_curv;
  }
  zbin = nint(curv_scale * (0.0 - min_curv));
  HISTOmakePDF(h_curv, h_curv);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(h_curv, "curv.plt");
  }

  zbin = HISTOfindHighestPeakInRegion(h_curv, 0, h_curv->nbins);

  for (total = 0, bin = zbin; bin < h_curv->nbins; bin++) {
    total += h_curv->counts[bin];
  }
  thresh = thresh_pct * total;
  for (total = 0, bin = zbin; bin < h_curv->nbins - 1; bin++) {
    total += h_curv->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  max_curv = h_curv->bins[bin];

  for (total = 0, bin = zbin; bin > 0; bin--) {
    total += h_curv->counts[bin];
  }
  thresh = thresh_pct * total;
  for (total = 0, bin = zbin; bin > 0; bin--) {
    total += h_curv->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  min_curv = h_curv->bins[bin];

  mris->min_curv = 10000.0f;
  mris->max_curv = -10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    if (vertex->curv > max_curv) {
      vertex->curv = max_curv;
      nthresh++;
    }
    if (vertex->curv < min_curv) {
      vertex->curv = min_curv;
      nthresh++;
    }
    if (vertex->curv < mris->min_curv) {
      mris->min_curv = vertex->curv;
    }
    if (vertex->curv > mris->max_curv) {
      mris->max_curv = vertex->curv;
    }
  }

  HISTOfree(&h_curv);
  if (nthresh > 0) fprintf(stderr, "%d vertices thresholded to be in [%2.2f %2.2f]\n", nthresh, min_curv, max_curv);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute the cortical thickness at each vertex by measuring the
  distance from it to the pial surface.

  This routine assumes that the white matter surface is stored in
  ORIGINAL_VERTICES, and that the current vertex positions reflect
  the pial surface.
  ------------------------------------------------------*/

int MRISfindClosestOrigVertices(MRIS *mris, int nbhd_size)
{
  int vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n, min_vno;
  float dx, dy, dz, dist, min_dist, nx, ny, nz, dot;

  memset(nbr_count, 0, 100 * sizeof(int));

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dx = v->x - v->origx;
    dy = v->y - v->origy;
    dz = v->z - v->origz;
    min_dist = sqrt(dx * dx + dy * dy + dz * dz);
    v->marked = 1;
    vtotal = 1;
    vlist[0] = vno;
    min_n = 0;
    min_vno = vno;
    for (ns = 1; ns <= nbhd_size; ns++) {
      vnum = 0; /* will be # of new neighbors added to list */
      for (i = 0; i < vtotal; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX          const * const vn  = &mris->vertices         [vlist[i]];
        if (vn->ripflag) {
          continue;
        }
        if (vn->marked && vn->marked < ns - 1) {
          continue;
        }
        for (n = 0; n < vnt->vnum; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]];
          if (vn2->ripflag || vn2->marked) /* already processed */
          {
            continue;
          }
          vlist[vtotal + vnum++] = vnt->v[n];
          vn2->marked = ns;
          dx = vn2->x - v->origx;
          dy = vn2->y - v->origy;
          dz = vn2->z - v->origz;
          dot = dx * nx + dy * ny + dz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dot = vn2->nx * nx + vn2->ny * ny + vn2->nz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dist = sqrt(dx * dx + dy * dy + dz * dz);
          if (dist < min_dist) {
            min_n = ns;
            min_dist = dist;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON) fprintf(stdout, "%d --> %d = %2.3f\n", vno, vnt->v[n], dist);
            min_vno = vnt->v[n];
          }
        }
      }
      vtotal += vnum;
    }

    nbr_count[min_n]++;
    for (n = 0; n < vtotal; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]];
      if (vn->ripflag) {
        continue;
      }
      vn->marked = 0;
    }
    v->curv = min_vno;  // BLETCH
  }

  for (n = 0; n <= nbhd_size; n++) {
    fprintf(stdout, "%d vertices at %d distance\n", nbr_count[n], n);
  }
  return (NO_ERROR);
}
/*
  find the closest pial vertex and use it's spherical coords to initialize
  the thickness minimization. Put the v->c[xyz] coords of the nearest pial vertex into v->[xyz] of each vertex.
*/
int MRISfindClosestPialVerticesCanonicalCoords(MRIS *mris, int nbhd_size)
{
  int vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n, min_vno;
  float dx, dy, dz, dist, min_dist, nx, ny, nz, dot;

  memset(nbr_count, 0, 100 * sizeof(int));

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    nx = v->wnx;
    ny = v->wny;
    nz = v->wnz;
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dx = v->pialx - v->whitex;
    dy = v->pialy - v->whitey;
    dz = v->pialz - v->whitez;
    min_dist = sqrt(dx * dx + dy * dy + dz * dz);
    v->marked = 1;
    vtotal = 1;
    vlist[0] = vno;
    min_n = 0;
    min_vno = vno;
    for (ns = 1; ns <= nbhd_size; ns++) {
      vnum = 0; /* will be # of new neighbors added to list */
      for (i = 0; i < vtotal; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX          const * const vn  = &mris->vertices         [vlist[i]];
        if (vn->ripflag) {
          continue;
        }
        if (vn->marked && vn->marked < ns - 1) {
          continue;
        }
        for (n = 0; n < vnt->vnum; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]];
          if (vn2->ripflag || vn2->marked) /* already processed */
          {
            continue;
          }
          vlist[vtotal + vnum++] = vnt->v[n];
          vn2->marked = ns;
          dx = vn2->pialx - v->whitex;
          dy = vn2->pialy - v->whitey;
          dz = vn2->pialz - v->whitez;
          dot = dx * nx + dy * ny + dz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dot = vn2->wnx * nx + vn2->wny * ny + vn2->wnz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dist = sqrt(dx * dx + dy * dy + dz * dz);
          if (dist < min_dist) {
            min_n = ns;
            min_dist = dist;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON) fprintf(stdout, "%d --> %d = %2.3f\n", vno, vnt->v[n], dist);
            min_vno = vnt->v[n];
          }
        }
      }
      vtotal += vnum;
    }

    nbr_count[min_n]++;
    for (n = 0; n < vtotal; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]];
      if (vn->ripflag) {
        continue;
      }
      vn->marked = 0;
    }
    v->curv = min_vno;                  // BLETCH
    v->x = mris->vertices[min_vno].cx;
    v->y = mris->vertices[min_vno].cy;
    v->z = mris->vertices[min_vno].cz;
  }

  for (n = 0; n <= nbhd_size; n++) {
    fprintf(stdout, "%d vertices at %d distance\n", nbr_count[n], n);
  }
  return (NO_ERROR);
}


int MRISmeasureCorticalThickness(MRIS *mris, int nbhd_size, float max_thick)
{
  int vno, n, vlist[100000], vtotal, ns, i, vnum, nbr_count[100], min_n, nwg_bad, ngw_bad;
  float dx, dy, dz, dist, min_dist, nx, ny, nz, dot;

  memset(nbr_count, 0, 100 * sizeof(int));
  nwg_bad = ngw_bad = 0;

  /* current vertex positions are gray matter, orig are white matter */
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (!(vno % 25000)) {
      fprintf(stdout, "%d of %d vertices processed\n", vno, mris->nvertices);
    }
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      v->curv = 0;
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dx = v->x - v->origx;
    dy = v->y - v->origy;
    dz = v->z - v->origz;
    min_dist = sqrt(dx * dx + dy * dy + dz * dz);
    v->marked = 1;
    vtotal = 1;
    vlist[0] = vno;
    min_n = 0;
    for (ns = 1; ns <= nbhd_size; ns++) {
      vnum = 0; /* will be # of new neighbors added to list */
      for (i = 0; i < vtotal; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX          const * const vn  = &mris->vertices         [vlist[i]];
        if (vn->ripflag) {
          continue;
        }
        if (vn->marked && vn->marked < ns - 1) {
          continue;
        }
        for (n = 0; n < vnt->vnum; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]];
          if (vn2->ripflag || vn2->marked) /* already processed */
          {
            continue;
          }
          vlist[vtotal + vnum++] = vnt->v[n];
          vn2->marked = ns;
          dx = vn2->x - v->origx;
          dy = vn2->y - v->origy;
          dz = vn2->z - v->origz;
          dot = dx * nx + dy * ny + dz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dot = vn2->nx * nx + vn2->ny * ny + vn2->nz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dist = sqrt(dx * dx + dy * dy + dz * dz);
	  if(Gdiag_no == vno) {
	    printf("vno=%d A %3d %6d (%g,%g,%g) (%g,%g,%g) %g %g\n",vno,n,vnt->v[n],v->origx,v->origy,v->origz,v->x,v->y,v->z,dist,min_dist);
	  }
          if (dist < min_dist) {
            min_n = ns;
            min_dist = dist;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON) fprintf(stdout, "%d --> %d = %2.3f\n", vno, vnt->v[n], dist);
          }
        }
      }
      vtotal += vnum;
    }

    nbr_count[min_n]++;
    for (n = 0; n < vtotal; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]];
      if (vn->ripflag) {
        continue;
      }
      vn->marked = 0;
    }
    if (min_dist > max_thick) {
      nwg_bad++;
      min_dist = max_thick;
    }
    if (min_dist < 0) {
      DiagBreak();
    }
    v->curv = min_dist;
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    if (!(vno % 25000)) {
      fprintf(stdout, "%d of %d vertices processed\n", vno, mris->nvertices);
    }
    VERTEX * const v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dx = v->x - v->origx;
    dy = v->y - v->origy;
    dz = v->z - v->origz;
    min_dist = sqrt(dx * dx + dy * dy + dz * dz);
    v->marked = 1;
    vtotal = 1;
    vlist[0] = vno;
    min_n = 0;
    for (ns = 1; ns <= nbhd_size; ns++) {
      vnum = 0; /* will be # of new neighbors added to list */
      for (i = 0; i < vtotal; i++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[i]];
        VERTEX          const * const vn  = &mris->vertices         [vlist[i]];
        if (vn->ripflag) {
          continue;
        }
        if (vn->marked && vn->marked < ns - 1) {
          continue;
        }
        for (n = 0; n < vnt->vnum; n++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n]];
          if (vn2->ripflag || vn2->marked) /* already processed */
          {
            continue;
          }
          vlist[vtotal + vnum++] = vnt->v[n];
          vn2->marked = ns;
          dx = v->x - vn2->origx;
          dy = v->y - vn2->origy;
          dz = v->z - vn2->origz;
          dot = dx * nx + dy * ny + dz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dot = vn2->nx * nx + vn2->ny * ny + vn2->nz * nz;
          if (dot < 0) /* must be outwards from surface */
          {
            continue;
          }
          dist = sqrt(dx * dx + dy * dy + dz * dz);
	  if(Gdiag_no == vno) {
	    printf("vno=%d B %3d %6d (%g,%g,%g) (%g,%g,%g) %g %g\n",vno,n,vnt->v[n],vn2->origx,vn2->origy,vn2->origz,v->x,v->y,v->z,dist,min_dist);
	  }
          if (dist < min_dist) {
            min_n = ns;
            min_dist = dist;
            if (min_n == nbhd_size && DIAG_VERBOSE_ON) fprintf(stdout, "%d --> %d = %2.3f\n", vno, vnt->v[n], dist);
          }
        }
      }
      vtotal += vnum;
    }

    nbr_count[min_n]++;
    for (n = 0; n < vtotal; n++) {
      VERTEX * const vn = &mris->vertices[vlist[n]];
      if (vn->ripflag) {
        continue;
      }
      vn->marked = 0;
    }
    if (DIAG_VERBOSE_ON && fabs(v->curv - min_dist) > 4.0)
      fprintf(stdout, "v %d, white->gray=%2.2f, gray->white=%2.2f\n", vno, v->curv, min_dist);
    if (min_dist > max_thick) {
      min_dist = max_thick;
      ngw_bad++;
    }
    if (min_dist < 0) {
      DiagBreak();
    }
    v->curv = (v->curv + min_dist) / 2;
    if (v->curv < 0) {
      DiagBreak();
    }
    if(Gdiag_no == vno) 
      printf("vno = %d, final measurment %g\n",vno,v->curv);
  }

  fprintf(stdout, "thickness calculation complete, %d:%d truncations.\n", nwg_bad, ngw_bad);
  for (n = 0; n <= nbhd_size; n++) {
    fprintf(stdout, "%d vertices at %d distance\n", nbr_count[n], n);
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisOrigNormalFace(MRIS *mris, int fac, int n, float norm[])
{
  int n0, n1;
  FACE *f;
  float v0[3], v1[3];
  VERTEX *v, *vn0, *vn1;

  n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  f = &mris->faces[fac];
  int const * pv = f->v;
  vn0 = &mris->vertices[pv[n0]];
  vn1 = &mris->vertices[pv[n1]];
  v = &mris->vertices[pv[n]];
  v0[0] = v->origx - vn0->origx;
  v0[1] = v->origy - vn0->origy;
  v0[2] = v->origz - vn0->origz;
  v1[0] = vn1->origx - v->origx;
  v1[1] = vn1->origy - v->origy;
  v1[2] = vn1->origz - v->origz;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1] * v0[2] + v0[1] * v1[2];
  norm[1] = v1[0] * v0[2] - v0[0] * v1[2];
  norm[2] = -v1[0] * v0[1] + v0[0] * v1[1];
  /*
    printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
    v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
  */
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisWhiteNormalFace(MRIS *mris, int fac, int n, float norm[])
{
  int n0, n1;
  FACE *f;
  float v0[3], v1[3];
  VERTEX *v, *vn0, *vn1;

  n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  f = &mris->faces[fac];
  int const * pv = f->v;
  vn0 = &mris->vertices[pv[n0]];
  vn1 = &mris->vertices[pv[n1]];
  v = &mris->vertices[pv[n]];
  v0[0] = v->whitex - vn0->whitex;
  v0[1] = v->whitey - vn0->whitey;
  v0[2] = v->whitez - vn0->whitez;
  v1[0] = vn1->whitex - v->whitex;
  v1[1] = vn1->whitey - v->whitey;
  v1[2] = vn1->whitez - v->whitez;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1] * v0[2] + v0[1] * v1[2];
  norm[1] = v1[0] * v0[2] - v0[0] * v1[2];
  norm[2] = -v1[0] * v0[1] + v0[0] * v1[1];
  /*
    printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
    v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
  */
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisPialNormalFace(MRIS *mris, int fac, int n, float norm[])
{
  int n0, n1;
  FACE *f;
  float v0[3], v1[3];
  VERTEX *v, *vn0, *vn1;

  n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  f = &mris->faces[fac];
  int const * pv = f->v;
  vn0 = &mris->vertices[pv[n0]];
  vn1 = &mris->vertices[pv[n1]];
  v = &mris->vertices[pv[n]];
  v0[0] = v->pialx - vn0->pialx;
  v0[1] = v->pialy - vn0->pialy;
  v0[2] = v->pialz - vn0->pialz;
  v1[0] = vn1->pialx - v->pialx;
  v1[1] = vn1->pialy - v->pialy;
  v1[2] = vn1->pialz - v->pialz;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1] * v0[2] + v0[1] * v1[2];
  norm[1] = v1[0] * v0[2] - v0[0] * v1[2];
  norm[2] = -v1[0] * v0[1] + v0[0] * v1[1];
  /*
    printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
    v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
  */
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeWhiteNormal(MRIS *mris, int vno, float norm[])
{
  float snorm[3];
  int n, num;

  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

  norm[0] = norm[1] = norm[2] = 0.0;
  for (num = n = 0; n < v->num; n++)
    if (!mris->faces[v->f[n]].ripflag) {
      num++;
      mrisWhiteNormalFace(mris, v->f[n], (int)v->n[n], snorm);
      norm[0] += snorm[0];
      norm[1] += snorm[1];
      norm[2] += snorm[2];
    }
  if (!num) {
    return (ERROR_BADPARM);
  }
  mrisNormalize(norm);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputePialNormal(MRIS *mris, int vno, float norm[])
{
  float snorm[3];
  int n, num;

  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

  norm[0] = norm[1] = norm[2] = 0.0;
  for (num = n = 0; n < v->num; n++)
    if (!mris->faces[v->f[n]].ripflag) {
      num++;
      mrisPialNormalFace(mris, v->f[n], (int)v->n[n], snorm);
      norm[0] += snorm[0];
      norm[1] += snorm[1];
      norm[2] += snorm[2];
    }
  if (!num) {
    return (ERROR_BADPARM);
  }
  mrisNormalize(norm);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeOrigNormal(MRIS *mris, int vno, float norm[])
{
  float snorm[3];
  int n, num;

  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

  norm[0] = norm[1] = norm[2] = 0.0;
  for (num = n = 0; n < v->num; n++)
    if (!mris->faces[v->f[n]].ripflag) {
      num++;
      mrisOrigNormalFace(mris, v->f[n], (int)v->n[n], snorm);
      norm[0] += snorm[0];
      norm[1] += snorm[1];
      norm[2] += snorm[2];
    }
  if (!num) {
    return (ERROR_BADPARM);
  }
  mrisNormalize(norm);
  return (NO_ERROR);
}


void computeVertexPseudoNormal(MRIS const *mris, int vno, float norm[3], int verbose)
{
  int n, n0, n1, n2;
  float v1[3], v2[3], alpha;
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];

  norm[0] = norm[1] = norm[2] = 0;

  for (n = 0; n < vt->num; n++) {
    int const fno = vt->f[n];
    FACE   const * const face = &mris->faces[fno];

    n0 = vt->n[n];
    n1 = (n0 == 2) ? 0 : n0 + 1;
    n2 = (n0 == 0) ? 2 : n0 - 1;

    if ((face->v[n0] != vno) || (face->v[n1] == vno) || (face->v[n2] == vno) || (face->v[n2] == face->v[n1])) {
      if (verbose >= VERBOSE_MODE_MEDIUM) {
        if (face->v[n0] != vno) {
          fprintf(WHICH_OUTPUT, "error for vno in face %d", vt->f[n]);
        }
        if (face->v[n1] == vno) {
          fprintf(WHICH_OUTPUT, "error for vn1 in face %d", vt->f[n]);
        }
        if (face->v[n2] == vno) {
          fprintf(WHICH_OUTPUT, "error for vn2 in face %d", vt->f[n]);
        }
        if (face->v[n2] == face->v[n1]) {
          fprintf(WHICH_OUTPUT, "error for vn in face %d", vt->f[n]);
        }

        fprintf(WHICH_OUTPUT, "face %d (%d,%d,%d) != (%d)\n", vt->f[n], face->v[n0], face->v[n1], face->v[n2], vno);

        if (verbose == VERBOSE_MODE_MEDIUM) {
          fprintf(stderr, "computeVertexPseudoNormal: SHOULD NOT HAPPEN\n");
        }

        //                      MRISwrite(mris,"rh.testdebug1");
        // MRISrestoreVertexPositions(mris,CANONICAL_VERTICES);
        // MRISwrite(mris,"rh.testdebug2");
        if (verbose == VERBOSE_MODE_HIGH) {
          ErrorExit(ERROR_BADPARM, "computeVertexPseudoNormal: SHOULD NOT HAPPEN\n");
        }
      }
    }

    v1[0] = mris->vertices[face->v[n1]].origx - v->origx;
    v1[1] = mris->vertices[face->v[n1]].origy - v->origy;
    v1[2] = mris->vertices[face->v[n1]].origz - v->origz;

    v2[0] = mris->vertices[face->v[n2]].origx - v->origx;
    v2[1] = mris->vertices[face->v[n2]].origy - v->origy;
    v2[2] = mris->vertices[face->v[n2]].origz - v->origz;

    alpha = MAX(0.0, MIN(1.0, F_DOT(v1, v2) / NORM3(v1) / NORM3(v2)));
    alpha = acos(alpha);

    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    norm[0] += alpha * fNorm->nx;
    norm[1] += alpha * fNorm->ny;
    norm[2] += alpha * fNorm->nz;
  }
}



int MRIScomputeNormal(MRIS *mris, int which, int vno, double *pnx, double *pny, double *pnz)
{
  float snorm[3], norm[3];
  int n, num;

  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

  norm[0] = norm[1] = norm[2] = 0.0;
  snorm[0] = snorm[1] = snorm[2] = 0.0;
  for (num = n = 0; n < vt->num; n++)
    if (!mris->faces[vt->f[n]].ripflag) {
      num++;
      switch (which) {
        case CURRENT_VERTICES:
          mrisNormalFace     (mris, vt->f[n], (int)vt->n[n], snorm);
          break;
        case ORIGINAL_VERTICES:
          mrisOrigNormalFace (mris, vt->f[n], (int)vt->n[n], snorm);
          break;
        case WHITE_VERTICES:
          mrisWhiteNormalFace(mris, vt->f[n], (int)vt->n[n], snorm);
          break;
        case PIAL_VERTICES:
          mrisPialNormalFace (mris, vt->f[n], (int)vt->n[n], snorm);
          break;
        default:
          ErrorExit(ERROR_BADPARM, "MRIScomputeNormal: which = %d not supported", which);
          break;
      }
      norm[0] += snorm[0];
      norm[1] += snorm[1];
      norm[2] += snorm[2];
    }
  if (!num) {
    return (ERROR_BADPARM);
  }

  mrisNormalize(norm);
  *pnx = norm[0];
  *pny = norm[1];
  *pnz = norm[2];

  return (NO_ERROR);
}


/*-----------------------------------------------------
    Support deferring the calculation of face norms until needed
    
    These values in computeDefectFaceNormal are the ones we are trying to avoid
    since calling that is an expensive part of mris_fix_topology
    
    The solution is to simply mark it as needing computation, and to compute it only 
    when getFaceNorm is called.
    
    To test this is working, there is an option for computing it immediately...
    
  ------------------------------------------------------*/

//#define CHECK_DEFERED_NORMS
template  <class Surface>
void computeDefectFaceNormal_calculate_template(
  Surface surface, int const fno, float* p_nx, float* p_ny, float* p_nz, float* p_orig_area)
{
  auto face = surface.faces(fno);
  
  auto v1 = face.v(0);
  auto v2 = face.v(1);
  auto v3 = face.v(2);

  /* compute the face normal on the original configuration */
  float a[3], b[3];

  a[0] = v2.origx() - v1.origx();
  b[0] = v3.origx() - v1.origx();
  a[1] = v2.origy() - v1.origy();
  b[1] = v3.origy() - v1.origy();
  a[2] = v2.origz() - v1.origz();
  b[2] = v3.origz() - v1.origz();

  float norm[3];
  F_CROSS(a, b, norm);

  float nx, ny, nz;
  nx = norm[0];
  ny = norm[1];
  nz = norm[2];

  /* normalize */
  float len = sqrtf(nx * nx + ny * ny + nz * nz);
  if (FZERO(len)) {
    // TO BE CHECKED
    //          fprintf(WHICH_OUTPUT,"face with a null normal (%f,%f,%f) - (%f,%f,%f) -
    //          (%f,%f,%f)",v1->origx,v1->origy,v1->origz,v2->origx,v2->origy,v2->origz,v3->origx,v3->origy,v3->origz);
    /* try another dot product */
    a[0] = 100.0 * (v3.origx() - v2.origx());
    b[0] = 100.0 * (v1.origx() - v2.origx());
    a[1] = 100.0 * (v3.origy() - v2.origy());
    b[1] = 100.0 * (v1.origy() - v2.origy());
    a[2] = 100.0 * (v3.origz() - v2.origz());
    b[2] = 100.0 * (v1.origz() - v2.origz());
    F_CROSS(a, b, norm);
    nx = norm[0];
    ny = norm[1];
    nz = norm[2];
    /* normalize */
    len = sqrt(nx * nx + ny * ny + nz * nz);
    if (FZERO(len)) {
      // fprintf(WHICH_OUTPUT,".");
      len = 1.0;
    }
    // fprintf(WHICH_OUTPUT,"\n");
  }
  *p_nx = nx / len;
  *p_ny = ny / len;
  *p_nz = nz / len;
  *p_orig_area = len / 2.0f;
}

static void computeDefectFaceNormal_calculate(
  MRIS const * const mris, int const fno, float* p_nx, float* p_ny, float* p_nz, float* p_orig_area)
{
  SurfaceFromMRIS::XYZPositionM::Surface surface(const_cast<MRIS*>(mris));
  computeDefectFaceNormal_calculate_template(surface, fno, p_nx, p_ny, p_nz, p_orig_area);
}

void computeDefectFaceNormal_calculate(
  MRIS_MP const * const mris, int const fno, float* p_nx, float* p_ny, float* p_nz, float* p_orig_area)
{
  SurfaceFromMRIS_MP::XYZPositionM::Surface surface(const_cast<MRIS_MP*>(mris));
  computeDefectFaceNormal_calculate_template(surface, fno, p_nx, p_ny, p_nz, p_orig_area);
}


void computeDefectFaceNormal(MRIS const * const mris, int const fno)
{
  float nx, ny, nz, orig_area;
  computeDefectFaceNormal_calculate(mris, fno, &nx, &ny, &nz, &orig_area);
  setFaceNorm    (mris, fno, nx, ny, nz);
  setFaceOrigArea(mris, fno, orig_area);
}



template <class SOME_MRIS>
void setFaceNorm_template(SOME_MRIS const * const mris, int fno, float nx, float ny, float nz) {
    FaceNormCacheEntry    * fNorm         = &mris->faceNormCacheEntries   [fno];
    FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];

    fNorm->nx = nx;
    fNorm->ny = ny;
    fNorm->nz = nz;
    fNormDeferred->deferred &= ~1;          // now known
}

void setFaceNorm(MRIS const * const mris, int fno, float nx, float ny, float nz) {
    setFaceNorm_template(mris, fno, nx, ny, nz);
}

void setFaceNorm(MRIS_MP const * const mris, int fno, float nx, float ny, float nz) {
    setFaceNorm_template(mris, fno, nx, ny, nz);
}


void setFaceOrigArea(MRIS const * const mris, int fno, float orig_area) {
    FaceNormCacheEntry    * fNorm         = &mris->faceNormCacheEntries   [fno];
    FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];
    fNorm->orig_area = orig_area;
    fNormDeferred->deferred &= ~2;            // now known
}

float getFaceOrigArea(MRIS const * const mris, int fno) {
    FaceNormCacheEntry    * fNorm         = &mris->faceNormCacheEntries   [fno];
    FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];
    if (fNormDeferred->deferred & 2) {      // must compute
        fprintf(stderr, "%s:%d NYI\n", __FILE__, __LINE__);
        *(int*)(-1) = 0;    // crash
    }
    return fNorm->orig_area;
}


template <class SOME_MRIS>
FaceNormCacheEntry const * getFaceNorm_template(SOME_MRIS const * const mris, int fno) {

    // volatile to stop the compiler from reordering the stores of the nx,ny,nz,orig with the stores of the deferred 
    //
    FaceNormCacheEntry    volatile * fNorm         = &mris->faceNormCacheEntries   [fno];
    FaceNormDeferredEntry volatile * fNormDeferred = &mris->faceNormDeferredEntries[fno];
    if (fNormDeferred->deferred) {
        float nx,ny,nz,orig_area;
	
	// Not locked so that can be done in parallel
	//
        computeDefectFaceNormal_calculate(mris, fno, &nx,&ny,&nz,&orig_area);
	
	// Just lock the update, since multiple threads would have got equally acceptable answers
	//
#ifdef HAVE_OMP
	#pragma omp critical
#endif
    	{
#ifdef CHECK_DEFERED_NORMS
            if (fNorm->deferred & 1) {      // must update
        	if (nx != fNorm->nx || ny != fNorm->ny || nz != fNorm->nz) {
                    fprintf(stderr, "%s:%d prediction of norm did not equal result\n",__FILE__, __LINE__);
        	}   
            }
#endif
            fNorm->nx = nx; fNorm->ny = ny; fNorm->nz = nz;
            if (fNormDeferred->deferred & 2) {      // must update
#ifdef CHECK_DEFERED_NORMS
        	if (orig_area != fNorm->orig_area) {
                    fprintf(stderr, "%s:%d prediction of norm did not equal result\n",__FILE__, __LINE__);
        	}
#endif
        	fNorm->orig_area = orig_area;
            }
            fNormDeferred->deferred = 0;
    	}
    }
    return (FaceNormCacheEntry*)fNorm;
}

FaceNormCacheEntry const * getFaceNorm(MRIS const * const mris, int fno) {
    return getFaceNorm_template(mris, fno);
}

FaceNormCacheEntry const * getFaceNorm(MRIS_MP const * const mris, int fno) {
    return getFaceNorm_template(mris, fno);
}


void mrisurf_deferSetFaceNorms(MRIS* mris) {
    static int use_parallel;
    static int once;
    if (!once) {
        once++;
        use_parallel = !!getenv("FREESURFER_deferSetFaceNorms_parallel");
    }
    if (!use_parallel) {
        // It looks like there is not enough work to go parallel...
#ifdef CHECK_DEFERED_NORMS
        int fno;
        for (fno = 0; fno < mris->nfaces; fno++) {
          FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];
          computeDefectFaceNormal_calculate(mris, fno, &fNorm->nx,&fNorm->ny,&fNorm->nz,&fNorm->orig_area);  // compute it now so can check later
          fNormDeferred->deferred = 3;  // compute them again later
        }
#else
        if (sizeof(mris->faceNormDeferredEntries[0]) != 1) *(int*)(-1) = 0;
        memset(&mris->faceNormDeferredEntries[0], 3, mris->nfaces);
#endif
    } else {
        int fno;
        ROMP_PF_begin  
#ifdef HAVE_OPENMP
        #pragma omp parallel for if_ROMP(assume_reproducible) 
#endif
        for (fno = 0; fno < mris->nfaces; fno++) {
          ROMP_PFLB_begin
          FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];
#ifdef CHECK_DEFERED_NORMS
          computeDefectFaceNormal_calculate(mris, fno, &fNorm->nx,&fNorm->ny,&fNorm->nz,&fNorm->orig_area);  // compute it now so can check later
#endif
          fNormDeferred->deferred = 3;  // compute them again later
          ROMP_PFLB_end
        }
        ROMP_PF_end
    }
}

static void recomputeOrUndeferFaceNorms(MRIS* mris, bool isRecompute) {
    int fno;
    ROMP_PF_begin  
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(assume_reproducible) 
#endif
    for (fno = 0; fno < mris->nfaces; fno++) {
      ROMP_PFLB_begin
      FaceNormDeferredEntry * fNormDeferred = &mris->faceNormDeferredEntries[fno];
#ifdef CHECK_DEFERED_NORMS
      if (isRecompute) fNormDeferred->deferred = 3;
      getFaceNorm(mris, fno);      
#else
      FaceNormCacheEntry * fNorm = &mris->faceNormCacheEntries[fno];
      if (!isRecompute && !fNormDeferred->deferred) continue;
      computeDefectFaceNormal_calculate(mris, fno, &fNorm->nx,&fNorm->ny,&fNorm->nz,&fNorm->orig_area);
      fNormDeferred->deferred = 0;
#endif
      ROMP_PFLB_end
    }
    ROMP_PF_end
}

void mrisurf_recomputeFaceNorms(MRIS* mris) {
    recomputeOrUndeferFaceNorms(mris, true);
}

void mrisurf_undeferSetFaceNorms(MRIS* mris) {
    recomputeOrUndeferFaceNorms(mris, false);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScountNegativeTriangles(MRIS *mris)
{
  int fno, negative;
  FACE *face;

  negative = 0;
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    if (face->area < 0.0f) {
      negative++;
    }
  }

  return (negative);
}


//==================================================================================================================
// More complicated properties
//
void mrisSetAvgInterVertexDist(MRIS *mris, double to) 
{
  double const * pc = &mris->avg_vertex_dist;
  double       * p  = (double*)pc;
  *p = to;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  same as  below, but tracks odx,ody,odz fields
  ------------------------------------------------------*/
int mrisTrackTotalDistanceNew(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float nc;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    nc = v->odx * v->nx + v->ody * v->ny + v->odz * v->nz;
    v->curv += nc;
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisTrackTotalDistance(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float nc;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    nc = v->dx * v->nx + v->dy * v->ny + v->dz * v->nz;
    v->curv += nc;
  }
  return (NO_ERROR);
}

int MRIScomputeAllDistances(MRIS *mris)
{
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "\ncomputing complete distance matrix...\n");
  }

  MRIScomputeMetricProperties(mris);

  int const nvalid = MRISvalidVertices(mris);

  LABEL* area = LabelAlloc(1, NULL, NULL);

  // this for loop can't be parallelized due to the use of MRISdistanceTransform and of marked
  //
  int done = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {

    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
    VERTEX          * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    if ((Gdiag & DIAG_SHOW) && (++done % (nvalid / 20)) == 0) {
      fprintf(stdout, "%%%1.0f done\n", 100.0f * (float)done / (float)nvalid);
      fflush(stdout);
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    area->lv[0].vno = vno;
    area->n_points = 1;
    MRISdistanceTransform(mris, area, DTRANS_MODE_UNSIGNED);        // reads v->vnum distances out of v->dist

    // make huge v, dist, dist_orig    
    //
    int const old_vtotal = vt->vtotal;
    vt->vtotal = nvalid - 1;                                        // set to something other then v->v#num
    vt->v = (int *)realloc(vt->v, (nvalid - 1)*sizeof(int));
    MRISmakeDist(mris, vno);
    if (!vt->v || !v->dist || !v->dist_orig) {
      ErrorExit(ERROR_NOMEMORY, "MRIScomputeAllDistances: could not allocate %d-sized arrays", nvalid - 1);
    }
    
    MRISclearMarks(mris);
    v->marked = 1;  // don't add self to list
    
    // keep the v[0..old_vtotal) so v#num stays right
    //
    int n;
    for (n = 0; n < old_vtotal; n++) {
      VERTEX * const vn = &mris->vertices[vt->v[n]];
      v->dist[n] = v->dist_orig[n] = vn->val;                       // v->dist will be restored by caller
      vn->marked = 1;
    }
    
    // fill in the rest
    // this trashes the values between vtotal and vnums[nsizeMax] which is bad
    // so note the loss
    //
    mris->max_nsize = mris->nsize;
    vt->nsizeMax = vt->nsizeCur;
    
    int vno2;
    for (vno2 = 0; vno2 < mris->nvertices; vno2++) {
      VERTEX * const vn = &mris->vertices[vno2];
      if (vn->ripflag || vn->marked) continue;

      v->dist[n] = v->dist_orig[n] = vn->val;                   // set to a different definition of 'distance'
      vt->v[n++] = vno2;
      vn->marked = 1;
      cheapAssert(n < nvalid);
    }
    
    if (vno == Gdiag_no) {
      char fname[STRLEN];
      sprintf(fname, "vno%d.mgz", vno);
      MRISwriteValues(mris, fname);

      for (n = vno2 = 0; vno2 < mris->nvertices; vno2++) {
        VERTEX const * const vn = &mris->vertices[vno2];
        if (vn->marked) {
          int m, found;
          n++;
          for (found = m = 0; m < vt->vtotal; m++)
            if (vt->v[m] == vno2) {
              found = 1;
              break;
            }
          if (found == 0 && vno2 != vno) DiagBreak();
        }
      }
    }
  }
  
  mrisCheckVertexFaceTopology(mris);
  
  LabelFree(&area);
  return (NO_ERROR);
}

double MRIStotalVariationDifference(MRIS *mris)
{
  int vno, nvertices;
  VERTEX *v;
  double var, delta;

  nvertices = mris->nvertices;
  for (var = 0.0, vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    delta = fabs(v->k1 - v->k2);
    delta *= delta;
    var += v->area * delta;
    if (!std::isfinite(var)) ErrorPrintf(ERROR_BADPARM, "curvature at vertex %d is not finite!\n", vno);
  }
  return (var);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIStotalVariation(MRIS *mris)
{
  int vno, nvertices;
  VERTEX *v;
  double var, delta;

  nvertices = mris->nvertices;
  for (var = 0.0, vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    delta = v->k1 * v->k1 + v->k2 * v->k2;
    var += v->area * delta;
    if (!std::isfinite(var)) ErrorPrintf(ERROR_BADPARM, "curvature at vertex %d is not finite!\n", vno);
  }
  return (var);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if this triangle contains another vertex within it.
  If so, it should not be added to the tessellation.
  ------------------------------------------------------*/
#if SPHERE_INTERSECTION
int containsAnotherVertexOnSphere(MRIS *mris, int vno0, int vno1, int vno2, int mode)
{
  int i, n, nvertices, vertices[1000];
  double sign[3], normal[3][3], orgn[3][3], tangent[3], test;

  // first compute the 3 normals, 3 original points and associated constants
  {
    VERTEX const * V0 = &mris->vertices[vno0];
    VERTEX const * V1 = &mris->vertices[vno1];
    VERTEX const * V2 = &mris->vertices[vno2];

    for (i = 0; i < 3; i++) {
      if (i) {
        // circular rotation
        VERTEX const * const v = V0;
        V0 = V1;
        V1 = V2;
        V2 = v;
      }

      // compute normal for edge V1<-->V2
      orgn[i][0] = (V1->cx + V2->cx) / 2.0;
      orgn[i][1] = (V1->cy + V2->cy) / 2.0;
      orgn[i][2] = (V1->cz + V2->cz) / 2.0;

      tangent[0] = V2->cx - V1->cx;
      tangent[1] = V2->cy - V1->cy;
      tangent[2] = V2->cz - V1->cz;
      // normal to edge in the planar basis
      F_CROSS(orgn[i], tangent, normal[i]);

      tangent[0] = V0->cx - orgn[i][0];
      tangent[1] = V0->cy - orgn[i][1];
      tangent[2] = V0->cz - orgn[i][2];

      sign[i] = F_DOT(tangent, normal[i]);
    }
  }
  
  // Next, list all the neighboring vertices
  memset(vertices, 0, 1000 * sizeof(int));
  for (nvertices = 0, i = 0; i < 3; i++) {
    int v0,v1,v2;
    switch (i) {
      case 0:
        v0 = vno0;
        v1 = vno1;
        v2 = vno2;
        break;
      case 1:
        v0 = vno1;
        v1 = vno0;
        v2 = vno2;
        break;
      default:
        v0 = vno2;
        v1 = vno0;
        v2 = vno1;
        break;
    }
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[v0];
    
    for (n = 0; n < vt->vnum; n++) {
      if (v1 == vt->v[n] || v2 == vt->v[n]) {
        continue;
      }
      VERTEX * const vn = &mris->vertices[vt->v[n]];
      if (vn->marked) {
        continue;
      }
      if (mode && vn->fixedval == 0) {
        continue;  // experimental
      }
      // add vn in the list
      vertices[nvertices++] = vt->v[n];
      vn->marked = 1;
    }
  }

  // unmark all vertices in the list
  for (n = 0; n < nvertices; n++) {
    mris->vertices[vertices[n]].marked = 0;
  }

  // then check intersection for all the neighboring vertices
  for (n = 0; n < nvertices; n++) {
    VERTEX const * const vn = &mris->vertices[vertices[n]];
    for (i = 0; i < 3; i++) {
      tangent[0] = vn->cx - orgn[i][0];
      tangent[1] = vn->cy - orgn[i][1];
      tangent[2] = vn->cz - orgn[i][2];

      test = F_DOT(tangent, normal[i]);

      if (sign[i] * test < 0) {
        break;
      }
    }
    if (i == 3) {
      return 1;
    }
  }

  return (0);
}
#else
int containsAnotherVertex(
    MRIS *mris, int vno0, int vno1, int vno2, double e0[3], double e1[3], double origin[3])
{
  int n, vno, i, n1;
  VERTEX *vn, *v0, *v1, *v2, *v;
  double cx, cy, cz, x0, y0, x1, y1, x2, y2, xn, yn, m, b;
  FILE *fp;

  v = &mris->vertices[vno0];

  if (((vno1 == 110718) || (vno1 == 110732) || (vno1 == 110748)) &&
      ((vno2 == 110718) || (vno2 == 110732) || (vno2 == 110748)) &&
      ((vno0 == 110718) || (vno0 == 110732) || (vno0 == 110748))) {
    fp = fopen("triangle.log", "w");
    for (n1 = 0; n1 < VERTICES_PER_FACE; n1++) {
      switch (n1) {
        default:
        case 0:
          v = &mris->vertices[vno0];
          break;
        case 1:
          v = &mris->vertices[vno1];
          break;
        case 2:
          v = &mris->vertices[vno2];
          break;
      }
      for (n = 0; n < v->vnum; n++) {
        vno = v->v[n];
        if (vno == vno0 || vno == vno1 || vno == vno2) {
          continue;
        }
        vn = &mris->vertices[vno];
        cx = vn->cx - origin[0];
        cy = vn->cy - origin[1];
        cz = vn->cz - origin[2];
        xn = cx * e0[0] + cy * e0[1] + cz * e0[2];
        yn = cx * e1[0] + cy * e1[1] + cz * e1[2];
        fprintf(fp, "# vertex %d\n", vno);
        fprintf(fp, "%f %f\n\n", xn, yn);
      }
    }
    fprintf(fp, "# vertices %d %d %d\n", vno0, vno1, vno2);
    v0 = &mris->vertices[vno0];
    v1 = &mris->vertices[vno1];
    v2 = &mris->vertices[vno2];

    cx = v0->cx - origin[0];
    cy = v0->cy - origin[1];
    cz = v0->cz - origin[2];
    x0 = cx * e0[0] + cy * e0[1] + cz * e0[2];
    y0 = cx * e1[0] + cy * e1[1] + cz * e1[2];

    cx = v1->cx - origin[0];
    cy = v1->cy - origin[1];
    cz = v1->cz - origin[2];
    x1 = cx * e0[0] + cy * e0[1] + cz * e0[2];
    y1 = cx * e1[0] + cy * e1[1] + cz * e1[2];

    cx = v2->cx - origin[0];
    cy = v2->cy - origin[1];
    cz = v2->cz - origin[2];
    x2 = cx * e0[0] + cy * e0[1] + cz * e0[2];
    y2 = cx * e1[0] + cy * e1[1] + cz * e1[2];

    fprintf(fp, "%f %f\n%f %f\n%f %f\n%f %f\n\n", x0, y0, x1, y1, x2, y2, x0, y0);
    fclose(fp);
  }

  for (n1 = 0; n1 < VERTICES_PER_FACE; n1++) {
    switch (n1) {
      default:
      case 0:
        v = &mris->vertices[vno0];
        break;
      case 1:
        v = &mris->vertices[vno1];
        break;
      case 2:
        v = &mris->vertices[vno2];
        break;
    }
    for (n = 0; n < v->vnum; n++) {
      vno = v->v[n];
      if (vno == vno0 || vno == vno1 || vno == vno2) {
        continue;
      }

      vn = &mris->vertices[v->v[n]];
      cx = vn->cx - origin[0];
      cy = vn->cy - origin[1];
      cz = vn->cz - origin[2];
      xn = cx * e0[0] + cy * e0[1] + cz * e0[2];
      yn = cx * e1[0] + cy * e1[1] + cz * e1[2];
      for (i = 0; i < VERTICES_PER_FACE; i++) {
        switch (i) {
          default:
          case 0:
            v0 = &mris->vertices[vno0];
            v1 = &mris->vertices[vno1];
            v2 = &mris->vertices[vno2];
            break;
          case 1:
            v0 = &mris->vertices[vno1];
            v1 = &mris->vertices[vno2];
            v2 = &mris->vertices[vno0];
            break;
          case 2:
            v0 = &mris->vertices[vno2];
            v1 = &mris->vertices[vno0];
            v2 = &mris->vertices[vno1];
            break;
        }
        cx = v0->cx - origin[0];
        cy = v0->cy - origin[1];
        cz = v0->cz - origin[2];
        x0 = cx * e0[0] + cy * e0[1] + cz * e0[2];
        y0 = cx * e1[0] + cy * e1[1] + cz * e1[2];

        cx = v1->cx - origin[0];
        cy = v1->cy - origin[1];
        cz = v1->cz - origin[2];
        x1 = cx * e0[0] + cy * e0[1] + cz * e0[2];
        y1 = cx * e1[0] + cy * e1[1] + cz * e1[2];

        cx = v2->cx - origin[0];
        cy = v2->cy - origin[1];
        cz = v2->cz - origin[2];
        x2 = cx * e0[0] + cy * e0[1] + cz * e0[2];
        y2 = cx * e1[0] + cy * e1[1] + cz * e1[2];

        if (FEQUAL(x1, x0)) /* vertical leg */
        {
          if (((x2 - x1) * (xn - x1)) < 0) /* on opposite sides of leg */
          {
            break;
          }
        }
        else {
          m = (y1 - y0) / (x1 - x0);
          b = y1 - m * x1;
          if ((y2 - (m * x2 + b)) * (yn - (m * xn + b)) < 0) {
            break;
          }
        }
      }
      if (i >= VERTICES_PER_FACE) /* all inside */
      {
        return (1);
      }
    }
  }
  return (0);
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScurvatureError(MRIS *mris, double Kd)
{
  int vno, nvertices, n;
  VERTEX *v;
  double Kerror, deltaK;

  nvertices = mris->nvertices;
  for (Kerror = 0.0, n = vno = 0; vno < nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    deltaK = (v->k1 - Kd);
    deltaK *= deltaK;
    Kerror += v->area * deltaK;
    n++;
    deltaK = (v->k2 - Kd);
    deltaK *= deltaK;
    Kerror += v->area * deltaK;
    n++;
  }
  /*  return(sqrt(Kerror/(double)n)) ;*/
  return (Kerror);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISsurfaceRASToTalairachVoxel(
    MRIS *mris, MRI *mri, double xw, double yw, double zw, double *pxv, double *pyv, double *pzv)
{
  //   conformed -->  SRAS
  //       |            |
  //       |            | SRASToRAS
  //       V            V
  //      src   --->   RAS
  //       |            |
  //       |            |  MNI xfm
  //       V            V
  //   talvol   --->  talRAS
  //        lta->dst stores this
  // we assume that talvol and src has the same direction cosine
  //   talvol -> talRAS =
  //
  // double  xt, yt, zt ;
  MRI *talVol = 0;
  VOL_GEOM talgeom;
  MATRIX *SRASToTalVol = 0;
  MATRIX *SRASToTalRAS = 0;  // not to be confused with SRASToTalSRAS_
  MATRIX *SRASToRAS = 0;
  MATRIX *talRASToTalVol = 0;

  LT *lt = &mris->lta->xforms[0];

  ////////////////////////////////////////////////////
  SRASToRAS = MatrixAlloc(4, 4, MATRIX_REAL);
  MatrixIdentity(4, SRASToRAS);
  *MATRIX_RELT(SRASToRAS, 1, 4) = lt->src.c_r;
  *MATRIX_RELT(SRASToRAS, 2, 4) = lt->src.c_a;
  *MATRIX_RELT(SRASToRAS, 3, 4) = lt->src.c_s;
  ///////////////////////////////////////////////////
  SRASToTalRAS = MatrixMultiply(lt->m_L, SRASToRAS, NULL);
  ///////////////////////////////////////////////////
  talgeom = lt->dst;
  talVol = MRIallocHeader(talgeom.width, talgeom.height, talgeom.depth, MRI_UCHAR, 1);
  useVolGeomToMRI(&talgeom, talVol);
  talRASToTalVol = extract_r_to_i(talVol);
  //////////////////////////////////////////////////
  // now combine to get SRAS->talVol
  SRASToTalVol = MatrixMultiply(talRASToTalVol, SRASToTalRAS, NULL);

  TransformWithMatrix(SRASToTalVol, xw, yw, zw, pxv, pyv, pzv);

  MatrixFree(&SRASToRAS);
  MatrixFree(&SRASToTalRAS);
  MatrixFree(&talRASToTalVol);
  MatrixFree(&SRASToTalVol);

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISvertexToVoxel(MRIS *mris, VERTEX *v, MRI *mri, double *pxv, double *pyv, double *pzv)
{
  double xw, yw, zw;

  xw = v->x;
  yw = v->y;
  zw = v->z;
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, pxv, pyv, pzv);
  return (NO_ERROR);
}

int MRISvertexNormalToVoxel(MRIS *mris, VERTEX *v, MRI *mri, double *pnx, double *pny, double *pnz)
{
  double xw, yw, zw;
  double xv0, yv0, zv0, xv1, yv1, zv1;

  xw = v->x;
  yw = v->y;
  zw = v->z;
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, &xv0, &yv0, &zv0);
  xw = v->x + v->nx;
  yw = v->y + v->ny;
  zw = v->z + v->nz;
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, &xv1, &yv1, &zv1);
  *pnx = xv1 - xv0;
  *pny = yv1 - yv0;
  *pnz = zv1 - zv0;
  return (NO_ERROR);
}
int MRISvertexNormalToVoxelScaled(MRIS *mris, VERTEX *v, MRI *mri, double *pnx, double *pny, double *pnz)
{
  double nx, ny, nz, len;

  MRISvertexNormalToVoxel(mris, v, mri, &nx, &ny, &nz);
  len = sqrt(nx * nx + ny * ny + nz * nz);
  if (FZERO(len)) return (ERROR_BADPARM);
  *pnx = nx / len;
  *pny = ny / len;
  *pnz = nz / len;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISvertexCoordToVoxel(MRIS *mris, VERTEX *v, MRI *mri, int coords, double *pxv, double *pyv, double *pzv)
{
  double xw = 0., yw = 0., zw = 0.;

  MRISgetCoords(v, coords, &xw, &yw, &zw);
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, pxv, pyv, pzv);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISorigVertexToVoxel(MRIS *mris, VERTEX *v, MRI *mri, double *pxv, double *pyv, double *pzv)
{
  double xw, yw, zw;

  xw = v->origx;
  yw = v->origy;
  zw = v->origz;
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, pxv, pyv, pzv);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISwhiteVertexToVoxel(MRIS *mris, VERTEX *v, MRI *mri, double *pxv, double *pyv, double *pzv)
{
  double xw, yw, zw;

  xw = v->whitex;
  yw = v->whitey;
  zw = v->whitez;
  MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, pxv, pyv, pzv);
  return (NO_ERROR);
}

/*
  on calling, the white, pial and spherical locations must all be loaded, where the current coordinates (v->[xyz])
  indicate the correspondence between white and pial, and the canonical spherical ones v->c[xyz] have the ?h.sphere in
  them.
*/
int MRISmeasureThicknessFromCorrespondence(MRIS *mris, MHT *mht, float max_thick)
{
  int vno, alloced;
  VERTEX *v;
  float xw, yw, zw, dx, dy, dz, xp, yp, zp;

  if (mht == NULL) {
    mht = MHTcreateFaceTable_Resolution(mris, CANONICAL_VERTICES, 1.0);  // to lookup closest face
    alloced = 1;
  }
  else
    alloced = 0;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      v->curv = 0;
      continue;
    }
    if (vno == Gdiag_no) DiagBreak();
    MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
    MRISsampleFaceCoordsCanonical(mht, mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);
    dx = xp - xw;
    dy = yp - yw;
    dz = zp - zw;
    v->curv = MIN(max_thick, sqrt(dx * dx + dy * dy + dz * dz));
  }

  if (alloced) MHTfree(&mht);

  return (NO_ERROR);
}


int mrisComputeNormalDotDistribution(MRIS *mris, HISTOGRAM *h_dot)
{
  int vno, bin, n, num;

  float bin_size, min_dot, max_dot, bin_val, dot, dx, dy, dz, nx, ny, nz, x, y, z;
  HISTOGRAM *h_raw;

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  MRIScomputeMetricProperties(mris);
  min_dot = 100000;
  max_dot = -100000;

  /* first compute min and max */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    x = v->x;
    y = v->y;
    z = v->z;
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - x;
      dy = vn->y - y;
      dz = vn->z - z;
      dot = dx * nx + dy * ny * dz * nz;
      if (dot < min_dot) {
        min_dot = dot;
      }
      if (dot > max_dot) {
        max_dot = dot;
      }
    }
  }

  /* add one bin at either end for almost zero probability events */
  bin_size = (max_dot - min_dot) / (h_dot->nbins - 2);
  h_dot->bin_size = bin_size;
  for (bin_val = min_dot - bin_size, bin = 0; bin < h_dot->nbins; bin++, bin_val += bin_size) {
    h_dot->bins[bin] = bin_val;
  }

  min_dot = h_dot->bins[0];

  /* now fill in distribution */
  for (num = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    x = v->x;
    y = v->y;
    z = v->z;
    for (n = 0; n < vt->vnum; n++) {
      num++;
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - x;
      dy = vn->y - y;
      dz = vn->z - z;
      dot = dx * nx + dy * ny * dz * nz;
      bin = (int)((dot - min_dot) / bin_size);
      if (bin == 0) {
        DiagBreak();
      }

      bin = MIN(h_dot->nbins, MAX(0, bin));
      h_dot->counts[bin]++;
    }
  }

  for (bin = 0; bin < h_dot->nbins; bin++) {
    if (h_dot->counts[bin] == 0) {
      h_dot->counts[bin] = 0.01;
    }
    h_dot->counts[bin] /= (float)num;
  }
  h_raw = HISTOcopy(h_dot, NULL);
  // to correct the bug in HISTOcopy..
  h_raw->bin_size = h_dot->bin_size;
  HISTOsmooth(h_raw, h_dot, 2.0);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(h_dot, "ndot.plt");
    HISTOplot(h_raw, "rdot.plt");
  }
  HISTOfree(&h_raw);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);
  return (NO_ERROR);
}
int mrisComputePrincipalCurvatureDistributions(MRIS *mris,
                                                      HISTOGRAM *h_k1,
                                                      HISTOGRAM *h_k2,
                                                      MRI *mri_k1_k2)
{
  int vno, bin, bink1, bink2, nvertices, x, y;
  VERTEX *v;
  float k1_bin_size, k2_bin_size, min_k1, max_k1, min_k2, max_k2, bin_val, norm;


  MRIScomputeSecondFundamentalForm(mris);
  min_k1 = min_k2 = 100000;
  max_k1 = max_k2 = -100000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (v->k1 < min_k1) {
      min_k1 = v->k1;
    }
    if (v->k2 < min_k2) {
      min_k2 = v->k2;
    }
    if (v->k1 > max_k1) {
      max_k1 = v->k1;
    }
    if (v->k2 > max_k2) {
      max_k2 = v->k2;
    }
  }

//    fprintf(stderr,"     k1: (min,max)=(%f,%f)\n",min_k1,max_k1);
//    fprintf(stderr,"     k2: (min,max)=(%f,%f)\n",min_k2,max_k2);

  min_k1 = MAX(-3, min_k1);
  max_k1 = MIN(3, max_k1);
  min_k2 = MAX(-3, min_k2);
  max_k2 = MIN(3, max_k2);

  k1_bin_size = (max_k1 - min_k1) / h_k1->nbins;
  k2_bin_size = (max_k2 - min_k2) / h_k2->nbins;

  mri_k1_k2->xsize = (max_k1 - min_k1) / (float)mri_k1_k2->width;
  mri_k1_k2->ysize = (max_k2 - min_k2) / (float)mri_k1_k2->height;
  mri_k1_k2->xstart = min_k1;
  mri_k1_k2->ystart = min_k2;

  h_k1->bin_size = k1_bin_size;
  h_k2->bin_size = k2_bin_size;

  for (bin_val = min_k1, bin = 0; bin < h_k1->nbins; bin++, bin_val += k1_bin_size) {
    h_k1->bins[bin] = bin_val;
  }
  for (bin_val = min_k2, bin = 0; bin < h_k2->nbins; bin++, bin_val += k2_bin_size) {
    h_k2->bins[bin] = bin_val;
  }

  nvertices = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    nvertices++;
    bin = MIN(h_k1->nbins - 1, MAX(0, (int)((v->k1 - min_k1) / k1_bin_size)));
    h_k1->counts[bin]++;
    bin = MIN(h_k2->nbins - 1, MAX(0, (int)((v->k2 - min_k2) / k2_bin_size)));
    h_k2->counts[bin]++;

    bink1 = MIN(mri_k1_k2->width - 1, MAX(0, (int)((v->k1 - mri_k1_k2->xstart) / mri_k1_k2->xsize)));

    bink2 = MIN(mri_k1_k2->height - 1, MAX(0, (int)((v->k2 - mri_k1_k2->ystart) / mri_k1_k2->ysize)));

    MRIFvox(mri_k1_k2, bink1, bink2, 0) += 1.0f;
  }

  for (bin = 0; bin < h_k1->nbins; bin++) {
    if (h_k1->counts[bin] == 0) {
      h_k1->counts[bin] = 0.01;
    }
    h_k1->counts[bin] /= (float)nvertices;
    ;
  }
  for (bin = 0; bin < h_k2->nbins; bin++) {
    if (h_k2->counts[bin] == 0) {
      h_k2->counts[bin] = 0.01;
    }
    h_k2->counts[bin] /= (float)nvertices;
    ;
  }

  for (x = 0; x < 100; x++)
    for (y = 0; y < 100; y++) {
      if (FZERO(MRIFvox(mri_k1_k2, x, y, 0))) {
        MRIFvox(mri_k1_k2, x, y, 0) = 0.1;
      }
    }
  for (norm = 0.0, x = 0; x < 100; x++)
    for (y = 0; y < 100; y++) {
      norm += MRIFvox(mri_k1_k2, x, y, 0);
    }

  for (x = 0; x < 100; x++)
    for (y = 0; y < 100; y++) {
      MRIFvox(mri_k1_k2, x, y, 0) = MRIFvox(mri_k1_k2, x, y, 0) / norm;
    }


  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(h_k1, "k1.plt");
    HISTOplot(h_k2, "k2.plt");
    MRIwrite(mri_k1_k2, "mri_k1_k2.mgh");
  }
  return (NO_ERROR);
}


//=====================================================================================================================
static int get_face_axes(
    MRIS *mris, FACE *face, float *pe1x, float *pe1y, float *pe1z, float *pe2x, float *pe2y, float *pe2z)
{
  double V0[3], V1[3], V2[3], e1[3], e2[3], tmp[3], norm, dot;

  V0[0] = mris->vertices[face->v[0]].x;
  V0[1] = mris->vertices[face->v[0]].y;
  V0[2] = mris->vertices[face->v[0]].z;
  V1[0] = mris->vertices[face->v[1]].x;
  V1[1] = mris->vertices[face->v[1]].y;
  V1[2] = mris->vertices[face->v[1]].z;
  V2[0] = mris->vertices[face->v[2]].x;
  V2[1] = mris->vertices[face->v[2]].y;
  V2[2] = mris->vertices[face->v[2]].z;

  // first basis vector
  SUB(e1, V1, V0);
  norm = VLEN(e1);
  if (FZERO(norm))  // degenerate case - not much to do
  {
    VERTEX *v = &mris->vertices[face->v[0]];
    *pe1x = v->e1x;
    *pe1y = v->e1y;
    *pe1z = v->e1z;
    *pe2x = v->e2x;
    *pe2y = v->e2y;
    *pe2z = v->e2z;
    return (NO_ERROR);
  }
  SCALAR_MUL(e1, 1.0 / norm, e1);

  // project 1st basis vector out of 2nd to orthonormalize it
  SUB(e2, V2, V0);
  norm = VLEN(e2);
  if (FZERO(norm))  // degenerate case - not much to do
  {
    VERTEX *v = &mris->vertices[face->v[0]];
    *pe1x = v->e1x;
    *pe1y = v->e1y;
    *pe1z = v->e1z;
    *pe2x = v->e2x;
    *pe2y = v->e2y;
    *pe2z = v->e2z;
    return (NO_ERROR);
  }
  SCALAR_MUL(e2, 1.0 / norm, e2);
  dot = DOT(e1, e2);
  SCALAR_MUL(tmp, dot, e1);
  SUB(e2, e2, tmp);
  norm = VLEN(e2);
  if (FZERO(norm))  // degenerate case - not much to do
  {
    VERTEX *v = &mris->vertices[face->v[0]];
    *pe1x = v->e1x;
    *pe1y = v->e1y;
    *pe1z = v->e1z;
    *pe2x = v->e2x;
    *pe2y = v->e2y;
    *pe2z = v->e2z;
    return (NO_ERROR);
  }
  SCALAR_MUL(e2, 1.0 / norm, e2);

  *pe1x = e1[0];
  *pe1y = e1[1];
  *pe1z = e1[2];
  *pe2x = e2[0];
  *pe2y = e2[1];
  *pe2z = e2[2];
  return (NO_ERROR);
}


//==================================================================================================================
// Deformity support
//
void INTEGRATION_PARMS_copy   (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src) {
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
    memcpy(dst, src, sizeof(*src)); // note: copies the const fp et. al.!
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
}

void INTEGRATION_PARMS_setFp  (INTEGRATION_PARMS* parms, FILE* file) {
  FILE* const* fpcp = &parms->fp;
  FILE*      * fpp  = (FILE**)fpcp;
  *fpp = file;
}

void INTEGRATION_PARMS_openFp (INTEGRATION_PARMS* parms, const char* name, const char* mode) {
  FILE* file = fopen(name, mode);
  if (!file) {
    fprintf(stderr, "%s:%d Error opening parms.fp using filename %s mode %s\n", __FILE__, __LINE__, name, mode);
    exit(1);
  }
  INTEGRATION_PARMS_setFp(parms, file);
}

void INTEGRATION_PARMS_closeFp(INTEGRATION_PARMS* parms) {
  if (!parms->fp) return;
  fclose(parms->fp);
  INTEGRATION_PARMS_setFp(parms, NULL);
}

void INTEGRATION_PARMS_copyFp (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src) {
  INTEGRATION_PARMS_setFp(dst, src->fp);
}



static double ashburnerTriangleEnergy(MRIS *mris, int fno, double lambda)
{
  static MATRIX *_m_x[_MAX_FS_THREADS], *_m_x_inv[_MAX_FS_THREADS], *_m_y[_MAX_FS_THREADS], *_m_J[_MAX_FS_THREADS],
      *_m_U[_MAX_FS_THREADS], *_m_V[_MAX_FS_THREADS], *_m_m[_MAX_FS_THREADS], *_m_evectors[_MAX_FS_THREADS];
  MATRIX *m_x, *m_x_inv, *m_y, *m_J, *m_U, *m_V, *m_m, *m_evectors = NULL;
  FACE *f;
  VERTEX *v1, *v2, *v3;
  float e1x, e1y, e1z, e2x, e2y, e2z, x11, x12, x13, x21, x22, x23, y11, y12, y13, y21, y22, y23, evalues[2], det, s11,
      s22;
  double energy /*, a, b, c, d*/;
#ifdef HAVE_OPENMP
  // thread ID
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif

  m_x = _m_x[tid];
  m_m = _m_m[tid];
  m_x_inv = _m_x_inv[tid];
  m_y = _m_y[tid];
  m_J = _m_J[tid];
  m_U = _m_U[tid];
  m_V = _m_V[tid];
  m_evectors = _m_evectors[tid];

  if (m_x == NULL) {
    m_x = MatrixAlloc(3, 3, MATRIX_REAL);
    m_m = MatrixAlloc(3, 3, MATRIX_REAL);
    m_x_inv = MatrixAlloc(3, 3, MATRIX_REAL);
    m_y = MatrixAlloc(3, 3, MATRIX_REAL);
    m_J = MatrixAlloc(2, 2, MATRIX_REAL);
    m_U = MatrixAlloc(2, 2, MATRIX_REAL);
    m_V = MatrixAlloc(2, 2, MATRIX_REAL);

    _m_x[tid] = m_x;
    _m_m[tid] = m_m;
    _m_x_inv[tid] = m_x_inv;
    _m_y[tid] = m_y;
    _m_J[tid] = m_J;
    _m_U[tid] = m_U;
    _m_V[tid] = m_V;
  }

  f = &mris->faces[fno];
  if (fno == Gdiag_no || f->area < 0) {
    DiagBreak();
  }
  get_face_axes(mris, f, &e1x, &e1y, &e1z, &e2x, &e2y, &e2z);
  // first index refers to coord # and second to which vertex
  v1 = &mris->vertices[f->v[0]];
  v2 = &mris->vertices[f->v[1]];
  v3 = &mris->vertices[f->v[2]];

  // original coords in tangent plane
  x11 = v1->cx * e1x + v1->cy * e1y + v1->cz * e1z;
  x21 = v1->cx * e2x + v1->cy * e2y + v1->cz * e2z;
  x12 = v2->cx * e1x + v2->cy * e1y + v2->cz * e1z;
  x22 = v2->cx * e2x + v2->cy * e2y + v2->cz * e2z;
  x13 = v3->cx * e1x + v3->cy * e1y + v3->cz * e1z;
  x23 = v3->cx * e2x + v3->cy * e2y + v3->cz * e2z;

  // mapped coords in tangent plane
  y11 = v1->x * e1x + v1->y * e1y + v1->z * e1z;
  y21 = v1->x * e2x + v1->y * e2y + v1->z * e2z;
  y12 = v2->x * e1x + v2->y * e1y + v2->z * e1z;
  y22 = v2->x * e2x + v2->y * e2y + v2->z * e2z;
  y13 = v3->x * e1x + v3->y * e1y + v3->z * e1z;
  y23 = v3->x * e2x + v3->y * e2y + v3->z * e2z;

  *MATRIX_RELT(m_x, 1, 1) = x11;
  *MATRIX_RELT(m_x, 1, 2) = x12;
  *MATRIX_RELT(m_x, 1, 3) = x13;
  *MATRIX_RELT(m_x, 2, 1) = x21;
  *MATRIX_RELT(m_x, 2, 2) = x22;
  *MATRIX_RELT(m_x, 2, 3) = x23;
  *MATRIX_RELT(m_x, 3, 1) = 1.0;
  *MATRIX_RELT(m_x, 3, 2) = 1.0;
  *MATRIX_RELT(m_x, 3, 3) = 1.0;

  *MATRIX_RELT(m_y, 1, 1) = y11;
  *MATRIX_RELT(m_y, 1, 2) = y12;
  *MATRIX_RELT(m_y, 1, 3) = y13;
  *MATRIX_RELT(m_y, 2, 1) = y21;
  *MATRIX_RELT(m_y, 2, 2) = y22;
  *MATRIX_RELT(m_y, 2, 3) = y23;
  *MATRIX_RELT(m_y, 3, 1) = 1.0;
  *MATRIX_RELT(m_y, 3, 2) = 1.0;
  *MATRIX_RELT(m_y, 3, 3) = 1.0;

#define SMALL 1e-8
  if (MatrixInverse(m_x, m_x_inv) == NULL) {
    return (log(SMALL) * log(SMALL));
  }
  MatrixMultiply(m_y, m_x_inv, m_m);
  MatrixCopyRegion(m_m, m_J, 1, 1, 2, 2, 1, 1);
  //  MatrixSVDEigenValues(m_J, evalues) ;
  det = MatrixDeterminant(m_J);
  if (det < 0) {
    return (log(SMALL) * log(SMALL));
  }
  m_evectors = MatrixEigenSystem(m_J, evalues, m_evectors);
  s11 = evalues[0];
  s22 = evalues[1];
  _m_evectors[tid] = m_evectors;

  if (!std::isfinite(s11)) return (log(SMALL) * log(SMALL));

  if (!std::isfinite(s22)) return (log(SMALL) * log(SMALL));

  if (s11 <= 0) s11 = SMALL;

  s11 = log(s11);
  s11 *= s11;

  if (s22 <= 0) {
    s22 = SMALL;
  }
  s22 = log(s22);
  s22 *= s22;

  energy = lambda * (1 + det) * (s11 + s22);  // log and square of snn already taken

  if (!std::isfinite(energy)) DiagBreak();

  return (energy / 2);
}


float mrisSampleAshburnerTriangleEnergy(
    MRIS * const mris, int const vno, INTEGRATION_PARMS * const parms, float cx, float cy, float cz)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];

  double sse, ox, oy, oz, sse_total;
  int n;

  if (v - mris->vertices == Gdiag_no) DiagBreak();

  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  ox = v->x;
  oy = v->y;
  oz = v->z;  //  save current location
  v->x = cx;
  v->y = cy;
  v->z = cz;  // change location to compute energy
  for (sse_total = 0.0, n = 0; n < vt->num; n++) {
    sse = ashburnerTriangleEnergy(mris, vt->f[n], parms->l_ashburner_lambda);
    if (sse < 0 || !std::isfinite(sse)) DiagBreak();

    sse_total += sse;
  }
  v->x = ox;
  v->y = oy;
  v->z = oz;  // restore original location of vertex
  sse_total /= vt->num;
  return (sse_total);
}


int mrisComputeCurvatureMinMax(MRIS *mris)
{
  int vno, found = 0;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (found == 0) {
      mris->max_curv = mris->min_curv = v->curv;
      found = 1;
    }
    else {
      if (v->curv > mris->max_curv) {
        mris->max_curv = v->curv;
      }
      if (v->curv < mris->min_curv) {
        mris->min_curv = v->curv;
      }
    }
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseCurvatureDifference(MRIS *mris)
{
  int vno;
  VERTEX *vertex;
  float kmin, kmax;

  kmin = 100000.0f;
  kmax = -100000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = fabs(vertex->k1 - vertex->k2);
    if (vertex->curv > kmax) {
      kmax = vertex->curv;
    }
    if (vertex->curv < kmin) {
      kmin = vertex->curv;
    }
  }

  mris->min_curv = kmin;
  mris->max_curv = kmax;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseCurvatureMax(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float kmin, kmax;

  kmin = 100000.0f;
  kmax = -100000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (fabs(v->k1) > fabs(v->k2)) {
      v->curv = v->k1;
    }
    else {
      v->curv = v->k2;
    }
    /*    v->curv = MAX(fabs(v->k1), fabs(v->k2)) ;*/
    if (v->curv > kmax) {
      kmax = v->curv;
    }
    if (v->curv < kmin) {
      kmin = v->curv;
    }
    if (v->curv < 0) {
      DiagBreak();
    }
  }

  fprintf(stdout, "kmin = %2.2f, kmax = %2.2f\n", kmin, kmax);
  mris->min_curv = kmin;
  mris->max_curv = kmax;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseCurvatureMin(MRIS *mris)
{
  int vno;
  VERTEX *v;
  float kmin, kmax;

  kmin = 100000.0f;
  kmax = -100000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (fabs(v->k1) > fabs(v->k2)) {
      v->curv = v->k2;
    }
    else {
      v->curv = v->k1;
    }
    /*    v->curv = MIN(v->k1, v->k2) ;*/
    if (v->curv > kmax) {
      kmax = v->curv;
    }
    if (v->curv < kmin) {
      kmin = v->curv;
    }
    if (v->curv < 0) {
      DiagBreak();
    }
  }

  mris->min_curv = kmin;
  mris->max_curv = kmax;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseCurvatureStretch(MRIS *mris)
{
  int vno, n;
  float kmin, kmax, dist, dist_orig, curv, dist_scale, max_stretch, stretch;

  dist_scale = sqrt(mris->orig_area / mris->total_area);
  kmin = 100000.0f;
  kmax = -100000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    max_stretch = 0.0f;
    for (curv = 0.0f, n = 0; n < vt->vtotal; n++) {
      dist = dist_scale * v->dist[n];
      dist_orig = v->dist_orig[n];
      stretch = dist - dist_orig;
      if (stretch > max_stretch) {
        max_stretch = stretch;
      }
    }
    v->curv = max_stretch;
    if (v->curv > kmax) {
      kmax = v->curv;
    }
    if (v->curv < kmin) {
      kmin = v->curv;
    }
    if (v->curv < 0) {
      DiagBreak();
    }
  }

  mris->min_curv = kmin;
  mris->max_curv = kmax;
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseNegCurvature(MRIS *mris)
{
  int vno, fno;
  FACE *f;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    v->curv = 0;
    for (fno = 0; fno < vt->num; fno++) {
      f = &mris->faces[vt->f[fno]];
      if (f->area < 0.0f) {
        v->curv = 1.0f;
      }
    }
  }

  mris->min_curv = 0;
  mris->max_curv = 1;
  return (NO_ERROR);
}


/*-----------------------------------------------------
  MRISaverageGradients() - spatially smooths gradients (dx,dy,dz)
  using num_avgs nearest-neighbor averages. See also
  MRISaverageGradientsFast()
  ------------------------------------------------------*/
typedef struct MRISaverageGradients_Data {
  float dx,dy,dz;
} MRISaverageGradients_Data;

typedef struct MRISaverageGradients_Control {    
  int numNeighbors;
  int firstNeighbor;
} MRISaverageGradients_Control;

static int                              MRISaverageGradients_neighbors_size;
static int*	                        MRISaverageGradients_neighbors;
static int  	                        MRISaverageGradients_controls_and_data_size;
static MRISaverageGradients_Control*    MRISaverageGradients_controls;
static MRISaverageGradients_Data*       MRISaverageGradients_datas_init;

static int*     MRISaverageGradients_new2old_indexMap;
static int*     MRISaverageGradients_old2new_indexMap;


// optimizing the program
//  There are two things that can be done to optimize the code
//  (a) put those that have the same numNeighbors together, to avoid the branch mispredicts
//  (b) keep the locality of reference by reordering only nearby indices
// both are important...

static bool MRISaverageGradients_index_le(int lhs, int rhs) {
 
    MRISaverageGradients_Control* controls = MRISaverageGradients_controls;
    
    // L2 caches are 256KB = 64KB fp numbers = about 20K vertices
    // so, for locality of reference purposes, we keep about 16768 = 2^14 together
    //
    int const lhsLocality = lhs>>14;
    int const rhsLocality = rhs>>14;
    
    if (       lhsLocality != rhsLocality)
    	return lhsLocality <= rhsLocality;

    if (       controls[lhs].numNeighbors != controls[rhs].numNeighbors)
    	return controls[lhs].numNeighbors <= controls[rhs].numNeighbors;

    return lhs <= rhs;
}

#define SORT_NAME           sort_MRISaverageGradients_index
#define SORT_NAME_partition sort_MRISaverageGradients_index_partition
#define SORT_NAME_isSorted  sort_MRISaverageGradients_index_isSorted
#define SORT_NAME_small     sort_MRISaverageGradients_index_small
   
#define SORT_ELEMENT    int
#define SORT_LE         MRISaverageGradients_index_le
#include "sort_definition.h"

static void MRISaverageGradients_optimize() {
  typedef MRISaverageGradients_Data    Data;
  typedef MRISaverageGradients_Control Control;

  static bool once, optimize;
  if (!once) { once = true;
    optimize = !getenv("FREESURFER_MRISaverageGradients_optimize_suppress");
    if (!optimize) {
      fprintf(stdout, "FREESURFER_MRISaverageGradients_optimize_suppress so not optimizing\n");
    }
  }
  
  if (!optimize) {
    free(MRISaverageGradients_new2old_indexMap);
    free(MRISaverageGradients_old2new_indexMap);
    return;
  }
  
  int                                 const controls_and_data_size  = MRISaverageGradients_controls_and_data_size;
  const MRISaverageGradients_Control* const controls                = MRISaverageGradients_controls;
  const MRISaverageGradients_Data*    const datas_init              = MRISaverageGradients_datas_init;
  const int*                          const neighbors               = MRISaverageGradients_neighbors;
  int                                 const neighbors_size          = MRISaverageGradients_neighbors_size;

  int* new2old_indexMap = MRISaverageGradients_new2old_indexMap;
  int* old2new_indexMap = MRISaverageGradients_old2new_indexMap;

  // sort the old indices to get the order to process the old indices in
  new2old_indexMap = (int*)realloc(new2old_indexMap, controls_and_data_size*sizeof(int));
  int ni;
  for (ni = 0; ni < controls_and_data_size; ni++) {
    new2old_indexMap[ni] = ni;
  }
  
  sort_MRISaverageGradients_index(new2old_indexMap, controls_and_data_size, true);

  // compute the inverse
  old2new_indexMap = (int*)realloc(old2new_indexMap, controls_and_data_size*sizeof(int));
  for (ni = 0; ni < controls_and_data_size; ni++) {
    int oi = new2old_indexMap[ni];
    old2new_indexMap[oi] = ni;
  }
  
  // permute the program
  int*	   const new_neighbors  = (int    *)malloc(neighbors_size         * sizeof(int));
  Data*	   const new_datas_init = (Data   *)malloc(controls_and_data_size * sizeof(Data));
  Control* const new_controls   = (Control*)malloc(controls_and_data_size * sizeof(Control));
  int            new_neighborsSize = 0;
  for (ni = 0; ni < controls_and_data_size; ni++) {
    int const old_index         = new2old_indexMap[ni];
    int const old_firstNeighbor = controls[old_index].firstNeighbor;
    int const numNeighbors      = controls[old_index].numNeighbors;
    new_controls[ni].numNeighbors  = numNeighbors;
    new_controls[ni].firstNeighbor = new_neighborsSize;
    int j;
    for (j = 0; j < numNeighbors; j++) {
      new_neighbors[new_neighborsSize++] = old2new_indexMap[neighbors[old_firstNeighbor + j]];	
    }
    new_datas_init[ni] = datas_init[old_index];
  }
    
  free(MRISaverageGradients_datas_init); MRISaverageGradients_datas_init = new_datas_init;
  free(MRISaverageGradients_controls);   MRISaverageGradients_controls   = new_controls;
  free(MRISaverageGradients_neighbors);  MRISaverageGradients_neighbors  = new_neighbors;

  MRISaverageGradients_new2old_indexMap = new2old_indexMap;
  MRISaverageGradients_old2new_indexMap = old2new_indexMap;
}


static void MRISaverageGradients_write(
  const char*                           filename,
  int                                   num_avgs,
  int                                   neighbors_size,
  const int*                            neighbors,
  int                                   index_to_vno_size,
  const MRISaverageGradients_Control*   controls, 
  const MRISaverageGradients_Data*      datas_inp) 
{
  fprintf(stdout, "%s:%d writing %s\n",__FILE__, __LINE__, filename);
  int i;
  FILE* f = fopen(filename, "w");
  fprintf(f, "num_avgs:%d\n", num_avgs);
  fprintf(f, "neighbors_size:%d\n", neighbors_size);
  for (i = 0; i < neighbors_size; i++) fprintf(f, "  %5d %5d\n", i, neighbors[i]);
  fprintf(f, "controls and data size:%d\n", index_to_vno_size);
  for (i = 0; i < index_to_vno_size; i++) {
    fprintf(f, "  %5d %5d %5d %g %g %g\n", 
	          i, 
		  controls[i].numNeighbors, controls[i].firstNeighbor, 
		  datas_inp[i].dx, datas_inp[i].dy, datas_inp[i].dz);
  }
  fclose(f);
  fprintf(stdout, "%s:%d finished writing %s\n",__FILE__, __LINE__, filename);
}

int MRISaverageGradients(MRIS *mris, int num_avgs)
{
  int i, vno;
  float sigma;
  VERTEX *v;
  MRI_SP *mrisp, *mrisp_blur;
  const char *UFSS;

  // Must explicity "setenv USE_FAST_SURF_SMOOTHER 0" to turn off fast
  UFSS = getenv("USE_FAST_SURF_SMOOTHER");
  if (!UFSS) {
    UFSS = "1";
  }
#ifdef HAVE_OPENMP
  UFSS = "0";  // MRISaverageGradientsFast can't use openmp
#endif
  if (strcmp(UFSS, "0")) {
    i = MRISaverageGradientsFast(mris, num_avgs);
    return (i);
  }
  if (Gdiag_no > 0 && DIAG_VERBOSE_ON) {
    printf("MRISaverageGradients()\n");
  }

  // Below is not used if using Fast
  if (num_avgs <= 0) {
    return (NO_ERROR);
  }

  if (Gdiag_no >= 0) {
    v = &mris->vertices[Gdiag_no];
    fprintf(stdout, "before averaging %d times dot = %2.2f ", num_avgs, v->dx * v->nx + v->dy * v->ny + v->dz * v->nz);
  }
  // use convolution
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE) {
    sigma = sqrt((float)num_avgs) / 4.0;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0);
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1);
    MRISgradientFromParameterization(mrisp_blur, mris);
    MRISPfree(&mrisp);
    MRISPfree(&mrisp_blur);
  }
  else
  {
    // This code shows the num_avgs typically is repeated 1-4 times as 1024, 256, 64, 16, 4, and then 1
    if (0) {
      static int num_avgs_being_counted = 0;
      static int count = 0;
      static int limit = 1;
      int changing = (num_avgs_being_counted != num_avgs);
      if (changing || count == limit) {
        if (count > 0) {
          fprintf(stderr,"%s:%d num_avgs:%d count:%d\n",__FILE__,__LINE__,num_avgs_being_counted,count);
          fprintf(stdout,"%s:%d num_avgs:%d count:%d\n",__FILE__,__LINE__,num_avgs_being_counted,count);
        }
        if (!changing) limit *= 2;
        else {
          num_avgs_being_counted = num_avgs;
          count = 0;
          limit = 1;
        }
      }
      count++;
    }

    // Since the num_avgs is so high there are far more efficient ways to iterate.
    // The following implementation is one of them.

    // The following code supports three ways of doing the two algorithms
    // - only the old algorithm
    // - only the new algorithm
    // - both algorithms, and compare them
    
    static int doNew = 1;
    static int doOld = 0;
    
    if (doOld && doNew && num_avgs > 10) {
        fprintf(stdout, "%s:%d reducing num_avgs:%d to 10\n", __FILE__, __LINE__, num_avgs);
        num_avgs = 10;
    }
     
    // This is the new algorithm
    int *vno_to_index = NULL;
    int *index_to_vno = NULL;
    int  index_to_vno_size = 0;
    
    typedef MRISaverageGradients_Data Data;
    Data *datas_inp = NULL;
    Data *datas_out = NULL;

    typedef MRISaverageGradients_Control Control;
    Control *controls = NULL;

    int  neighbors_capacity = 0;
    int  neighbors_size     = 0;
    int *neighbors          = NULL;
    
    if (doNew) {
    
      // Malloc the needed storage
      vno_to_index = (int*)    malloc(sizeof(int)     * mris->nvertices) ;
      index_to_vno = (int*)    malloc(sizeof(int)     * mris->nvertices) ;
      datas_inp    = (Data*)   malloc(sizeof(Data)    * mris->nvertices) ;
      datas_out    = (Data*)   malloc(sizeof(Data)    * mris->nvertices) ;
      controls     = (Control*)malloc(sizeof(Control) * mris->nvertices) ;
    
      // Find all the vertices that will change
      // Ripped ones, and those whose neighbors are all ripped,
      // will not change and hence are not entered
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX          const * const v  = &mris->vertices         [vno];
        
        int num = 0;
        if (!v->ripflag) {
          num = 1;
          // At this stage the neighbors with higher vno's will not have 
          // been entered into the vno_to_index table
          int  vnum = vt->vnum;
          int const *pnb  = vt->v;
          int  vnb;
          for (vnb = 0; vnb < vnum; vnb++) {
            int neighborVno = *pnb++;
            VERTEX* vn = &mris->vertices[neighborVno]; // neighboring vertex pointer
            if (vn->ripflag) continue;
            num++;
          }
        }
        
        if (num <= 1) {
          vno_to_index[vno] = -1;
        } else {
          vno_to_index[vno] = index_to_vno_size;
          index_to_vno[index_to_vno_size] = vno;
          Data *d = datas_inp + index_to_vno_size;
          d->dx = v->dx;
          d->dy = v->dy;
          d->dz = v->dz;
          Control *c = controls + index_to_vno_size;
          c->numNeighbors = num - 1;
          index_to_vno_size++;
        }
      }
      
      // For each vertex that will change, fill in the index of the vertices it must average with.
      // There could be many of them, so this array vector grows as needed.
      neighbors_capacity = 6 * index_to_vno_size;
      neighbors_size     = 0;
      neighbors          = (int*)malloc(sizeof(int) * neighbors_capacity);
      
      int index;
      for (index = 0; index < index_to_vno_size; index++) {
        const int vno = index_to_vno[index];
        Control *const c = controls + index;
        c->firstNeighbor = neighbors_size;
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        int const vnum = vt->vnum;
        int const *pnb = vt->v;
        int vnb;
        for (vnb = 0; vnb < vnum; vnb++) {
          int const neighborVno = *pnb++;
          int const neighborIndex = vno_to_index[neighborVno];
          if (neighborIndex < 0) continue;  // ripflag must be set
          if (neighbors_size == neighbors_capacity) {
            int  const new_neighbors_capacity = neighbors_capacity + neighbors_capacity/2;
            int* const new_neighbors = (int*)malloc(sizeof(int) * new_neighbors_capacity);
            int i; 
            for (i = 0; i < neighbors_size; i++) new_neighbors[i] = neighbors[i];
            free(neighbors);
            neighbors = new_neighbors ;
            neighbors_capacity = new_neighbors_capacity ;
          }
          neighbors[neighbors_size++] = neighborIndex ;
        }      
      }
      
      // Print out for later investigation
      //
      if (0) {
        static int once;
	if (!once++) {
          MRISaverageGradients_write("tesselation_averaging_data.txt",
            num_avgs, neighbors_size, neighbors, index_to_vno_size, controls, datas_inp);
	}
      }

      // Optimize enough iterations to payback the cost
      //
      freeAndNULL(MRISaverageGradients_new2old_indexMap);
      freeAndNULL(MRISaverageGradients_old2new_indexMap);

      if (num_avgs > 150) {

        MRISaverageGradients_neighbors_size         = neighbors_size;
        MRISaverageGradients_neighbors              = neighbors;
        MRISaverageGradients_controls_and_data_size = index_to_vno_size;
        MRISaverageGradients_controls               = controls;
        MRISaverageGradients_datas_init             = datas_inp;

        MRISaverageGradients_optimize();

        neighbors_size    = MRISaverageGradients_neighbors_size;
        neighbors         = MRISaverageGradients_neighbors;
        index_to_vno_size = MRISaverageGradients_controls_and_data_size;
        controls          = MRISaverageGradients_controls;
        datas_inp         = MRISaverageGradients_datas_init;
      }
      
      if (0) {
        static int once;
	if (!once++) {
          MRISaverageGradients_write("tesselation_averaging_data.optimized.txt", 
            num_avgs, neighbors_size, neighbors, index_to_vno_size, controls, datas_inp);
	}
      }

      // Choose chunk size assuming should evenly distribute the number of neighbors to process
      //
      typedef struct Chunk { int indexLo; } Chunk;
#define chunksCapacity 17  // 4 * omp_get_max_threads() + 1 // must be at least 3
      Chunk chunks[chunksCapacity];
      size_t chunksSize = 0;
      { int const neighborsPerChunk = neighbors_size/(chunksCapacity-1);
        int chunkSumNumNeighbors = 0;
        chunks[chunksSize++].indexLo = 0;
        int index;
        for (index = 0; index < index_to_vno_size; index++ ) {
          Control *c = controls  + index;
          chunkSumNumNeighbors += c->numNeighbors;
          if (chunkSumNumNeighbors < neighborsPerChunk) continue;
          chunks[chunksSize++].indexLo = index;
          chunkSumNumNeighbors = 0;
          if (chunksSize == chunksCapacity-1) break;
        }
        chunks[chunksSize++].indexLo = index_to_vno_size;
      }
#undef chunksCapacity

      // Do all the iterations
      for (i = 0 ; i < num_avgs ; i++) {
        
        unsigned int chunksIndex;
        ROMP_PF_begin
#ifdef HAVE_OPENMP
        #pragma omp parallel for if_ROMP(assume_reproducible)
#endif
        for (chunksIndex = 0; chunksIndex < chunksSize-1 ; chunksIndex++ ) {
          int
            indexLo = chunks[chunksIndex  ].indexLo, 
            indexHi = chunks[chunksIndex+1].indexLo;
          int index;
          for (index = indexLo; index < indexHi; index++ ) {
	    ROMP_PFLB_begin

            Data *d = datas_inp + index;
            Control *c = controls  + index;
            float dx = d->dx, dy = d->dy, dz = d->dz;
            
            int *nearby = neighbors + c->firstNeighbor;
            int n;
            for (n = 0; n < c->numNeighbors; n++) {
              int neighborIndex = nearby[n];
              Data *neighborData = datas_inp + neighborIndex;
              dx += neighborData->dx; dy += neighborData->dy; dz += neighborData->dz;
            }

            d = datas_out + index;
            float inv_num = 1.0f/(c->numNeighbors + 1);
            
            d->dx = dx*inv_num ; d->dy = dy*inv_num ; d->dz = dz*inv_num; 
          }	  
	  ROMP_PFLB_end
        }
	ROMP_PF_end
  
        // swap the output and the input going into the next round
        Data* tmp = datas_inp ; datas_inp = datas_out ; datas_out = tmp;
      }
    }
    
    // Only perform the old algorithm when needed
    if (doOld) {
    
      // This is the old algorithm
      for (i = 0 ; i < num_avgs ; i++) {
  
        ROMP_PF_begin
#ifdef HAVE_OPENMP
        #pragma omp parallel for if_ROMP(experimental)
#endif
  
        for (vno = 0; vno < mris->nvertices; vno++) {
	  ROMP_PFLB_begin
	  
          float dx, dy, dz, num;
          int vnb, vnum;
  
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
          VERTEX                * const v  = &mris->vertices         [vno];
          if (v->ripflag) ROMP_PF_continue;
  
          dx  = v->dx;
          dy  = v->dy;
          dz  = v->dz;
          int const * pnb = vt->v;
          // vnum = vt->v2num ? vt->v2num : vt->vnum;
          vnum = vt->vnum;
          for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
            VERTEX const * const vn = &mris->vertices[*pnb++]; // neighboring vertex pointer
            if (vn->ripflag) continue;
  
            num++;
            dx += vn->dx;
            dy += vn->dy;
            dz += vn->dz;
          }
          num++;
          v->tdx = dx / num;
          v->tdy = dy / num;
          v->tdz = dz / num;
	  ROMP_PFLB_end
        }
	ROMP_PF_end
	
	ROMP_PF_begin
#ifdef HAVE_OPENMP
	#pragma omp parallel for if_ROMP(experimental) schedule(static, 1)
#endif
        for (vno = 0; vno < mris->nvertices; vno++) {
	  ROMP_PFLB_begin
	  
          VERTEX * const v = &mris->vertices[vno];
          if (v->ripflag) ROMP_PF_continue;
  
          v->dx = v->tdx;
          v->dy = v->tdy;
          v->dz = v->tdz;
	  
	  ROMP_PFLB_end
        }
	ROMP_PF_end
	
      }  // end of one of the num_avgs iterations

    }  // if (doOld) 

    if (doOld && doNew) {
      fprintf(stdout, "%s:%d comparing old and new algorithms\n", __FILE__, __LINE__);
      int errors = 0;
      for (vno = 0 ; vno < mris->nvertices ; vno++) {
        VERTEX *v = &mris->vertices[vno];
        int old_index = vno_to_index[vno];
        if (old_index == -1) {
          // check has not changed - but this is not currently possible
          // but could be enabled by making all these nodes go into the high end of the initial datas_inp vector
        } else {
          int new_index = MRISaverageGradients_old2new_indexMap ? MRISaverageGradients_old2new_indexMap[old_index] : old_index;
          Data *d = datas_inp + new_index;
          if (!closeEnough(d->dx,v->dx) || !closeEnough(d->dy,v->dy) || !closeEnough(d->dz,v->dz)) {
            errors++;
            if (errors < 10) {
              fprintf(stdout,"%s:%d vno:%d d->dx:%g v->dx:%g d->dy:%g v->dy:%g d->dz:%g v->dz:%g - datas_inp\n",__FILE__,__LINE__,
                      vno, d->dx, v->dx, d->dy, v->dy, d->dz, v->dz);
            }
          }
          
        }
      }
      fprintf(stderr,"%s:%d errors:%d\n",__FILE__,__LINE__,errors);
      if (errors > 0) exit (103);
    }

    // Put the results where they belong
    // Note: this code does not set tdx, etc. - I believe they are temporaries   
    if (doNew) {
      int errors = 0;
      int new_index;
      for (new_index = 0 ; new_index < index_to_vno_size; new_index++) {

        int old_index = MRISaverageGradients_new2old_indexMap ? MRISaverageGradients_new2old_indexMap[new_index] : new_index;

        int vno   = index_to_vno[old_index];
        VERTEX *v = &mris->vertices[vno];
        Data   *d = datas_inp + new_index;
        
        v->dx = d->dx;
        v->dy = d->dy;
        v->dz = d->dz;
      }

      if (errors > 0) {
        fprintf(stderr,"%s:%d errors:%d\n",__FILE__,__LINE__,errors);
        exit (104);
      }
      
      free(datas_out); 
      free(datas_inp);
      free(neighbors); 
      free(controls);
      
      free(index_to_vno);
      free(vno_to_index);
    }
  }

  if (Gdiag_no >= 0) {
    float dot;
    v = &mris->vertices[Gdiag_no] ;
    dot = v->nx * v->dx + v->ny * v->dy + v->nz * v->dz;
    fprintf(stdout, " after dot = %2.2f D = (%2.2f, %2.2f, %2.2f), N = (%2.1f, %2.1f, %2.1f)\n",
            dot,v->dx, v->dy, v->dz, v->nx, v->ny, v->nz);
    if (fabs(dot) > 50) DiagBreak();
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  MRISaverageGradientsFast() - spatially smooths gradients (dx,dy,dz)
  using num_avgs nearest-neighbor averages. This is a faster version
  of MRISaverageGradients(). Should be about 40% faster. It also gives
  identical results as can be verified with MRISaverageGradientsFastCheck().
  ------------------------------------------------------*/
int MRISaverageGradientsFast(MRIS *mris, int num_avgs)
{
  int vno, vnb, *pnb, num, nthavg;

  float const **pdx, **pdy, **pdz;
  float sumdx, sumdy, sumdz;
  float const **pdx0, **pdy0, **pdz0;
  float *tdx, *tdy, *tdz, *tdx0, *tdy0, *tdz0;
  int *nNbrs, *nNbrs0, *rip, *rip0, nNbrsMax;

  if (Gdiag_no > 0 && DIAG_VERBOSE_ON) {
    printf("MRISaverageGradientsFast()\n");
  }

  // Alloc arrays. If there are ripped vertices, then only rip
  // needs nvertices elements
  nNbrsMax = 12;  // Should measure this, but overalloc does not hurt
  pdx = (float const **)calloc(mris->nvertices * nNbrsMax, sizeof(float *));
  pdy = (float const **)calloc(mris->nvertices * nNbrsMax, sizeof(float *));
  pdz = (float const **)calloc(mris->nvertices * nNbrsMax, sizeof(float *));
  tdx = (float        *)calloc(mris->nvertices,            sizeof(float));
  tdy = (float        *)calloc(mris->nvertices,            sizeof(float));
  tdz = (float        *)calloc(mris->nvertices,            sizeof(float));
  nNbrs = (int        *)calloc(mris->nvertices,            sizeof(int));
  rip = (int          *)calloc(mris->nvertices,            sizeof(int));

  pdx0 = pdx;
  pdy0 = pdy;
  pdz0 = pdz;
  tdx0 = tdx;
  tdy0 = tdy;
  tdz0 = tdz;
  rip0 = rip;
  nNbrs0 = nNbrs;

  // Set up pointers
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      rip[vno] = 1;
      continue;
    }
    rip[vno] = 0;
    *pdx++ = &(v->dx);
    *pdy++ = &(v->dy);
    *pdz++ = &(v->dz);
    pnb = vt->v;
    num = 1;
    for (vnb = 0; vnb < vt->vnum; vnb++) {
      VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
      if (vn->ripflag) {
        continue;
      }
      *pdx++ = &(vn->dx);
      *pdy++ = &(vn->dy);
      *pdz++ = &(vn->dz);
      num++;
    }
    *nNbrs++ = num;
  }

  // Loop through the iterations
  for (nthavg = 0; nthavg < num_avgs; nthavg++) {
    // Init pointers for this iteration
    pdx = pdx0;
    pdy = pdy0;
    pdz = pdz0;
    rip = rip0;
    nNbrs = nNbrs0;
    tdx = tdx0;
    tdy = tdy0;
    tdz = tdz0;
    // Loop through vertices, average nearest neighbors
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (*rip++) {
        continue;
      }
      sumdx = *(*pdx++);
      sumdy = *(*pdy++);
      sumdz = *(*pdz++);
      for (vnb = 0; vnb < (*nNbrs) - 1; vnb++) {
        sumdx += *(*pdx++);
        sumdy += *(*pdy++);
        sumdz += *(*pdz++);
      }
      *tdx++ = sumdx / (*nNbrs);
      *tdy++ = sumdy / (*nNbrs);
      *tdz++ = sumdz / (*nNbrs);
      nNbrs++;
    }
    // Load up for the next round
    rip = rip0;
    tdx = tdx0;
    tdy = tdy0;
    tdz = tdz0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (*rip++) {
        continue;
      }
      VERTEX * const v = &mris->vertices[vno];
      v->dx = *tdx++;
      v->dy = *tdy++;
      v->dz = *tdz++;
    }
  }

  rip = rip0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (*rip++) {
      continue;
    }
    VERTEX * const v = &mris->vertices[vno];
    v->tdx = v->dx;
    v->tdy = v->dy;
    v->tdz = v->dz;
  }

  free(pdx0);
  free(pdy0);
  free(pdz0);
  free(tdx0);
  free(tdy0);
  free(tdz0);
  free(rip0);
  free(nNbrs0);

  return (NO_ERROR);
}
/*--------------------------------------------------------*/
int MRISaverageGradientsFastCheck(int num_avgs)
{
  char tmpstr[2000], *UFSS;
  MRIS *mrisA, *mrisB;
  int k, msec, nerrs;
  Timer mytimer;
  float e;

  // Make sure to turn off override (restored later)
  UFSS = getenv("USE_FAST_SURF_SMOOTHER");
  setenv("USE_FAST_SURF_SMOOTHER", "0", 1);

  printf("MRISaverageGradientsFastCheck() num_avgs = %d\n", num_avgs);

  sprintf(tmpstr, "%s/subjects/fsaverage/surf/lh.white", getenv("FREESURFER_HOME"));
  printf("Reading surface %s\n", tmpstr);
  mrisA = MRISread(tmpstr);
  mrisB = MRISclone(mrisA);

  printf("Init\n");
  for (k = 0; k < mrisA->nvertices; k++) {
    mrisA->vertices[k].dx = drand48();
    mrisA->vertices[k].dy = drand48();
    mrisA->vertices[k].dz = drand48();
    mrisA->vertices[k].tdx = 0;
    mrisA->vertices[k].tdy = 0;
    mrisA->vertices[k].tdz = 0;
    mrisB->vertices[k].dx = mrisA->vertices[k].dx;
    mrisB->vertices[k].dy = mrisA->vertices[k].dy;
    mrisB->vertices[k].dz = mrisA->vertices[k].dz;
    mrisB->vertices[k].tdx = 0;
    mrisB->vertices[k].tdy = 0;
    mrisB->vertices[k].tdz = 0;
  }
  // Make sure at least 1 is ripped
  printf("ripping vertex 100\n");
  mrisA->vertices[100].ripflag = 1;
  mrisB->vertices[100].ripflag = 1;

  printf("Running Original Version\n");
  mytimer.reset();
  MRISaverageGradients(mrisA, num_avgs);
  msec = mytimer.milliseconds();
  printf("Original %6.2f min\n", msec / (1000 * 60.0));

  printf("Running Fast Version\n");
  mytimer.reset();
  MRISaverageGradientsFast(mrisB, num_avgs);
  msec = mytimer.milliseconds();
  printf("Fast %6.2f min\n", msec / (1000 * 60.0));

  printf("Checking for differences\n");
  nerrs = 0;
  for (k = 0; k < mrisA->nvertices; k++) {
    e = fabs(mrisA->vertices[k].dx - mrisB->vertices[k].dx);
    if (e > .0000001) {
      printf("ERROR: dx v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].dx, mrisB->vertices[k].dx);
      nerrs++;
    }
    e = fabs(mrisA->vertices[k].dy - mrisB->vertices[k].dy);
    if (e > .0000001) {
      printf("ERROR: dy v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].dy, mrisB->vertices[k].dy);
      nerrs++;
    }
    e = fabs(mrisA->vertices[k].dz - mrisB->vertices[k].dz);
    if (e > .0000001) {
      printf("ERROR: dz v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].dz, mrisB->vertices[k].dz);
      nerrs++;
    }
    e = fabs(mrisA->vertices[k].tdx - mrisB->vertices[k].tdx);
    if (e > .0000001) {
      printf("ERROR: tdx v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].tdx, mrisB->vertices[k].tdx);
      nerrs++;
    }
    e = fabs(mrisA->vertices[k].tdy - mrisB->vertices[k].tdy);
    if (e > .0000001) {
      printf("ERROR: tdy v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].tdy, mrisB->vertices[k].tdy);
      nerrs++;
    }
    e = fabs(mrisA->vertices[k].tdz - mrisB->vertices[k].tdz);
    if (e > .0000001) {
      printf("ERROR: tdz v=%d, e=%g, %g %g\n", k, e, mrisA->vertices[k].tdz, mrisB->vertices[k].tdz);
      nerrs++;
    }
  }
  printf("Found %d differences\n", nerrs);

  MRISfree(&mrisA);
  MRISfree(&mrisB);
  setenv("USE_FAST_SURF_SMOOTHER", UFSS, 1);
  return (nerrs);
}


int mrisAverageSignedGradients(MRIS *mris, int num_avgs)
{
  int i, vno;
  float sigma;

  MRI_SP *mrisp, *mrisp_blur;

  if (num_avgs <= 0) {
    return (NO_ERROR);
  }

  if (Gdiag_no >= 0) {
    VERTEX const * const v = &mris->vertices[Gdiag_no];
    fprintf(stdout, "before averaging dot = %2.2f ", v->dx * v->nx + v->dy * v->ny + v->dz * v->nz);

    {
      int vno;
      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * const vn = &mris->vertices[vno];
        if (vn->x > 1 && vn->dx > .1) DiagBreak();
      }
    }
    {
      VERTEX const *vmax = NULL;
      double mx = 0;
      int vno;
      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * const vn = &mris->vertices[vno];
        if (vn->x > 1 && vn->dx > .1) DiagBreak();
        if (vn->dx > mx && vn->x > 1) {
          mx = vn->dx;
          vmax = vn;
        }
      }
      if (mx > 0) DiagBreak();
    }
  }
  if (0 && mris->status == MRIS_PARAMETERIZED_SPHERE) /* use convolution */
  {
    sigma = sqrt((float)num_avgs) / 4.0;
    mrisp = MRISgradientToParameterization(mris, NULL, 1.0);
    mrisp_blur = MRISPblur(mrisp, NULL, sigma, -1);
    MRISgradientFromParameterization(mrisp_blur, mris);
    MRISPfree(&mrisp);
    MRISPfree(&mrisp_blur);
  }
  else
    for (i = 0; i < num_avgs; i++) {
      ROMP_PF_begin
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(experimental)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
	
        double dx, dy, dz, dot;
        int vnum, num, vnb;

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX                * const v  = &mris->vertices         [vno];
        if (v->ripflag) {
          continue;
        }
        dx = v->dx;
        dy = v->dy;
        dz = v->dz;
        int const * pnb = vt->v;
        /*      vnum = v->v2num ? v->v2num : v->vnum ;*/
        vnum = vt->vnum;
        for (num = 0, vnb = 0; vnb < vnum; vnb++) {
          VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
          if (vn->ripflag) {
            continue;
          }
          dot = vn->dx * v->dx + vn->dy * v->dy + vn->dz * v->dz;
          if (dot < 0) {
            continue; /* pointing in opposite directions */
          }

          num++;
          dx += vn->dx;
          dy += vn->dy;
          dz += vn->dz;
        }
        num++;
        v->tdx = dx / (float)num;
        v->tdy = dy / (float)num;
        v->tdz = dz / (float)num;
	ROMP_PFLB_end
      }
      ROMP_PF_end
      
      ROMP_PF_begin
#ifdef HAVE_OPENMP
      #pragma omp parallel for if_ROMP(experimental)
#endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
	
        VERTEX *v = &mris->vertices[vno];
        if (v->ripflag) ROMP_PF_continue;

        v->dx = v->tdx;
        v->dy = v->tdy;
        v->dz = v->tdz;
	ROMP_PFLB_end
      }
      ROMP_PF_end
    }
  if (Gdiag_no >= 0) {
    float dot;
    VERTEX const * const v = &mris->vertices[Gdiag_no];
    dot = v->nx * v->dx + v->ny * v->dy + v->nz * v->dz;
    fprintf(stdout, " after dot = %2.2f (%2.3f, %2.3f, %2.3f)\n", dot, v->dx, v->dy, v->dz);
    if (fabs(dot) > 50) {
      DiagBreak();
    }
    {
      VERTEX const *vmax = NULL;
      double mx = 0;
      int vno;
      for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * const vn = &mris->vertices[vno];
        if (vn->x > 1 && vn->dx > .1) DiagBreak();
        if (vn->dx > mx && vn->x > 1) {
          mx = vn->dx;
          vmax = vn;
        }
      }
      if (mx > 0) DiagBreak();
    }
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISsmoothSurfaceNormals(MRIS *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float nx, ny, nz, num, len;

  for (i = 0; i < navgs; i++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag) {
        continue;
      }
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (num = 0.0f, vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) {
          continue;
        }
        num++;
        nx += vn->nx;
        ny += vn->ny;
        nz += vn->nz;
      }
      num++; /* account for central vertex */
      v->tdx = nx / num;
      v->tdy = ny / num;
      v->tdz = nz / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      len = sqrt(v->tdx * v->tdx + v->tdy * v->tdy + v->tdz * v->tdz);
      if (FZERO(len)) {
        len = 1;
      }
      v->nx = v->tdx / len;
      v->ny = v->tdy / len;
      v->nz = v->tdz / len;
    }
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeBoundaryNormals(MRIS *mris)
{
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScomputeMeanCurvature(MRIS *mris)
{
  VECTOR *v_n, *v_e1, *v_e2, *v_i;
  int vno, i, N;
  float rsq, z, H, u, v, Hmin, Hmax;

  mrisComputeTangentPlanes(mris);
  v_n = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_i = VectorAlloc(3, MATRIX_REAL);

  Hmin = 10000.0f;
  Hmax = -Hmin;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    VECTOR_LOAD(v_n,  vertex->nx,  vertex->ny,  vertex->nz);
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);

    H = 0.0;
    N = 0;
    for (i = 0; i < vertext->vnum; i++) /* for each neighbor */
    {
      VERTEX const * const vnb = &mris->vertices[vertext->v[i]];
      if (vnb->ripflag) {
        continue;
      }
      VECTOR_LOAD(v_i, vnb->x - vertex->x, vnb->y - vertex->y, vnb->z - vertex->z);

      /* calculate projection onto tangent plane */
      u = V3_DOT(v_i, v_e1);
      v = V3_DOT(v_i, v_e2);
      rsq = u * u + v * v;
      z = V3_DOT(v_i, v_n); /* height above tangent plane */
      if (!FZERO(rsq)) {
        H += z / rsq;
        N++;
        if ((fabs(z / rsq) > 5.0) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
          fprintf(stdout, "%d --> %d: curvature = %2.1f\n", vno, vertext->v[i], z / rsq);
      }
    }
    if (N > 0) {
      vertex->H = 0.5f * H / (float)N;
    }
    else {
      vertex->H = 0;
    }
    if (vertex->H > Hmax) {
      Hmax = vertex->H;
    }
    if (vertex->H < Hmin) {
      Hmin = vertex->H;
    }
  }

  mris->min_curv = Hmin;
  mris->max_curv = Hmax;

  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "mean curvature range [%2.3f --> %2.3f]\n", Hmin, Hmax);
  }
  VectorFree(&v_i);
  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  return (NO_ERROR);
}

int MRIScomputeSecondFundamentalForm(MRIS *mris)
{
  return (MRIScomputeSecondFundamentalFormThresholded(mris, -1));
}

int MRIScomputeSecondFundamentalFormThresholded(MRIS *mris, double pct_thresh)
{
  double min_k1, min_k2, max_k1, max_k2, k1_scale, k2_scale, total, thresh, orig_rsq_thresh;
  int bin, zbin1, zbin2, nthresh = 0;
  int vno, i, n, vmax, nbad = 0, niter;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse, *m_eigen, *m_Q;
  VECTOR *v_c, *v_z, *v_n, *v_e1, *v_e2, *v_yi;
  float k1, k2, evalues[3], a11, a12, a21, a22, cond_no, kmax, kmin, rsq, k;
  double ui, vi, total_area = 0.0, max_error, vmean, vsigma, rsq_thresh;
  FILE *fp = NULL;
  HISTOGRAM *h_k1, *h_k2;

  if (mris->status == MRIS_PLANE) {
    return (NO_ERROR);
  }

  if (pct_thresh >= 0) {
    vmean = MRIScomputeTotalVertexSpacingStats(mris, &vsigma, NULL, NULL, NULL, NULL);
    rsq_thresh = MIN(vmean * .25, 1.0);
    rsq_thresh *= rsq_thresh;
  }
  else {
    rsq_thresh = 0.0;
  }

  orig_rsq_thresh = rsq_thresh;

  mrisComputeTangentPlanes(mris);

  v_c  = VectorAlloc(3, MATRIX_REAL);
  v_n  = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_yi = VectorAlloc(3, MATRIX_REAL);
  m_Q     = MatrixAlloc(2, 2, MATRIX_REAL); /* the quadratic form */
  m_eigen = MatrixAlloc(2, 2, MATRIX_REAL);

  mris->Kmin = mris->Hmin = 10000.0f;
  mris->Kmax = mris->Hmax = -10000.0f;
  mris->Ktotal = 0.0f;
  vmax = -1;
  max_error = -1.0;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    fp = fopen("curv.dat", "w");
  }
  
  for (vno = 0; vno < mris->nvertices; vno++) {
  
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];

    if (vertex->ripflag) {
      continue;
    }

    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz);
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);

    if (vertext->vtotal <= 0) {
      continue;
    }

    m_U = MatrixAlloc(vertext->vtotal, 3, MATRIX_REAL);
    v_z = VectorAlloc(vertext->vtotal,    MATRIX_REAL);

    if (vno == Gdiag_no) {
      DiagBreak();
    }

    /* fit a quadratic form to the surface at this vertex */
    rsq_thresh = orig_rsq_thresh;
    niter = 0;
    do {
      kmin = 10000.0f;
      kmax = -kmin;
      
      int kmaxI = -1;
      for (n = i = 0; i < vertext->vtotal; i++) {
        VERTEX const * const vnb = &mris->vertices[vertext->v[i]];

        if (vnb->ripflag) {
          continue;
        }
        /*
          calculate the projection of this vertex
          onto the local tangent plane
        */
        VECTOR_LOAD(v_yi, vnb->x - vertex->x, vnb->y - vertex->y, vnb->z - vertex->z);
        ui = V3_DOT(v_yi, v_e1);
        vi = V3_DOT(v_yi, v_e2);

        *MATRIX_RELT(m_U, n + 1, 1) = ui * ui;
        *MATRIX_RELT(m_U, n + 1, 2) = 2 * ui * vi;
        *MATRIX_RELT(m_U, n + 1, 3) = vi * vi;
        VECTOR_ELT(v_z, n + 1) = V3_DOT(v_n, v_yi); /* height above TpS */
        rsq = ui * ui + vi * vi;
        if (!FZERO(rsq) && rsq > rsq_thresh) {

          k = VECTOR_ELT(v_z, n + 1) / rsq;
          if (k > kmax) {
            kmax = k;
            kmaxI = i;
          }
          if (k < kmin) {
            kmin = k;
          }
          n++;
        }
      }
      
      rsq_thresh *= 0.25;
      if (n < 4) {
        DiagBreak();
      }
      if (niter++ > 100) {
        break;
      }
    } while (n < 4);

    MatrixBuffer m_Ut_buffer, m_tmp2_buffer;
    
    m_Ut   = MatrixAlloc2(m_U->cols, m_U->rows, MATRIX_REAL, &m_Ut_buffer);
    m_Ut   = MatrixTranspose(m_U, m_Ut);        /* Ut cols x rows */
    
    m_tmp2 = MatrixAlloc2(m_Ut->rows, m_U->cols, MATRIX_REAL, &m_tmp2_buffer);
    m_tmp2 = MatrixMultiply (m_Ut, m_U, m_tmp2);  /* Ut U  cols x rows * rows x cols = cols x cols */
    cond_no = MatrixConditionNumber(m_tmp2);

    MatrixBuffer m_inverse_buffer;
    m_inverse = MatrixAlloc2(m_tmp2->cols, m_tmp2->rows, MATRIX_REAL, &m_inverse_buffer);
    m_inverse = MatrixSVDInverse(m_tmp2, m_inverse); /* (Ut U)^-1 */
    if (!m_inverse) /* singular matrix - must be planar?? */
    {
      nbad++;
      evalues[0] = evalues[1] = 0.0;
      MatrixIdentity(m_eigen->rows, m_eigen);
    }
    else {
      MatrixBuffer m_tmp1_buffer;
      m_tmp1 = MatrixAlloc2(m_Ut->rows, v_z->cols, MATRIX_REAL, &m_tmp1_buffer); 
      m_tmp1 = MatrixMultiply(m_Ut, v_z, m_tmp1); /* Ut z */
      MatrixMultiply(m_inverse, m_tmp1, v_c);   /* (Ut U)^-1 Ut z */

      /* now build Hessian matrix */
      *MATRIX_RELT(m_Q, 1, 1) = 2 * VECTOR_ELT(v_c, 1);
      *MATRIX_RELT(m_Q, 1, 2) = *MATRIX_RELT(m_Q, 2, 1) = 2 * VECTOR_ELT(v_c, 2);
      *MATRIX_RELT(m_Q, 2, 2) = 2 * VECTOR_ELT(v_c, 3);

      if (cond_no >= ILL_CONDITIONED) {

        vertex->k1 = k1 = kmax;
        vertex->k2 = k2 = kmin;

        vertex->K = k1 * k2;
        vertex->H = (k1 + k2) / 2;
        MatrixFree(&m_Ut);
        MatrixFree(&m_tmp2);
        MatrixFree(&m_U);
        VectorFree(&v_z);
        MatrixFree(&m_tmp1);
        MatrixFree(&m_inverse);
        continue;
      }

      /* the columns of m_eigen will be the eigenvectors of m_Q */
      if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL) {

        nbad++;
        MatrixSVDEigenValues(m_Q, evalues);
        vertex->k1 = k1 = evalues[0];
        vertex->k2 = k2 = evalues[1];
        vertex->K = k1 * k2;
        vertex->H = (k1 + k2) / 2;
        MatrixFree(&m_Ut);
        MatrixFree(&m_tmp2);
        MatrixFree(&m_U);
        VectorFree(&v_z);
        MatrixFree(&m_tmp1);
        MatrixFree(&m_inverse);

        continue;
      }

      MatrixFree(&m_tmp1);
      MatrixFree(&m_inverse);
    }
    k1 = evalues[0];
    k2 = evalues[1];
    vertex->k1 = k1;
    vertex->k2 = k2;
    vertex->K = k1 * k2;
    vertex->H = (k1 + k2) / 2;
    if (fp) fprintf(fp, "%d %f %f %f %f\n", vno, k1, k2, vertex->K, vertex->H);
    if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
      fprintf(
          stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n", vno, vertex->k1, vertex->k2, vertex->K, vertex->H);
    if (vertex->K < mris->Kmin) {
      mris->Kmin = vertex->K;
    }
    if (vertex->H < mris->Hmin) {
      mris->Hmin = vertex->H;
    }
    if (vertex->K > mris->Kmax) {
      mris->Kmax = vertex->K;
    }
    if (vertex->H > mris->Hmax) {
      mris->Hmax = vertex->H;
    }
    mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area;
    total_area += (double)vertex->area;

    /* now update the basis vectors to be the principal directions */
    a11 = *MATRIX_RELT(m_eigen, 1, 1);
    a12 = *MATRIX_RELT(m_eigen, 1, 2);
    a21 = *MATRIX_RELT(m_eigen, 2, 1);
    a22 = *MATRIX_RELT(m_eigen, 2, 2);
    if (V3_LEN(v_e1) < 0.5) {
      DiagBreak();
    }

    vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21;
    vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21;
    vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21;
    vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22;
    vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22;
    vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22;
    if (SQR(vertex->e1x) + SQR(vertex->e1y) + SQR(vertex->e1z) < 0.5) {
      DiagBreak();
    }

    MatrixFree(&m_Ut);
    MatrixFree(&m_tmp2);
    MatrixFree(&m_U);
    VectorFree(&v_z);
  }

  if (fp) {
    fclose(fp);
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "max H error=%2.3f at %d\n", max_error, vmax);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "total area = %2.3f\n", total_area);
  }

  if (Gdiag & DIAG_SHOW && (nbad > 0)) {
    fprintf(stdout, "%d ill-conditioned points\n", nbad);
  }
  MatrixFree(&m_eigen);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  VectorFree(&v_c);
  VectorFree(&v_n);
  VectorFree(&v_yi);
  MatrixFree(&m_Q);

  if (pct_thresh < 0) {
    return (NO_ERROR);
  }

  // now remove outliers
  h_k1 = HISTOalloc(1000);
  h_k2 = HISTOalloc(1000);

  min_k1 = min_k2 = 1000;
  max_k1 = max_k2 = -1000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    if (vertex->k1 > max_k1) {
      max_k1 = vertex->k1;
    }
    if (vertex->k2 > max_k2) {
      max_k2 = vertex->k2;
    }
    if (vertex->k1 < min_k1) {
      min_k1 = vertex->k1;
    }
    if (vertex->k2 < min_k2) {
      min_k2 = vertex->k2;
    }
  }

  k1_scale = (h_k1->nbins - 1) / (max_k1 - min_k1);
  k2_scale = (h_k2->nbins - 1) / (max_k2 - min_k2);
  h_k1->bin_size = 1.0 / k1_scale;
  h_k2->bin_size = 1.0 / k2_scale;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    bin = nint(k1_scale * (vertex->k1 - min_k1));
    h_k1->counts[bin]++;
    bin = nint(k2_scale * (vertex->k2 - min_k2));
    h_k2->counts[bin]++;
  }
  for (bin = 0; bin < h_k1->nbins; bin++) {
    h_k1->bins[bin] = (bin / k1_scale) + min_k1;
  }
  for (bin = 0; bin < h_k2->nbins; bin++) {
    h_k2->bins[bin] = (bin / k2_scale) + min_k2;
  }
  zbin1 = nint(k1_scale * (0.0 - min_k1));
  zbin2 = nint(k2_scale * (0.0 - min_k2));
  HISTOmakePDF(h_k1, h_k1);
  HISTOmakePDF(h_k2, h_k2);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(h_k1, "k1.plt");
    HISTOplot(h_k2, "k2.plt");
  }

  zbin1 = HISTOfindHighestPeakInRegion(h_k1, 0, h_k1->nbins);
  zbin2 = HISTOfindHighestPeakInRegion(h_k2, 0, h_k2->nbins);

  for (total = 0, bin = zbin1; bin < h_k1->nbins; bin++) {
    total += h_k1->counts[bin];
  }
  thresh = pct_thresh * total;
  for (total = 0, bin = zbin1; bin < h_k1->nbins - 1; bin++) {
    total += h_k1->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  max_k1 = h_k1->bins[bin];

  for (total = 0, bin = zbin1; bin > 0; bin--) {
    total += h_k1->counts[bin];
  }
  thresh = pct_thresh * total;
  for (total = 0, bin = zbin1; bin > 0; bin--) {
    total += h_k1->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  min_k1 = h_k1->bins[bin];

  for (total = 0, bin = zbin2; bin < h_k2->nbins; bin++) {
    total += h_k2->counts[bin];
  }
  thresh = pct_thresh * total;
  for (total = 0, bin = zbin2; bin < h_k2->nbins - 1; bin++) {
    total += h_k2->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  max_k2 = h_k2->bins[bin];

  for (total = 0, bin = zbin2; bin > 0; bin--) {
    total += h_k2->counts[bin];
  }
  thresh = pct_thresh * total;
  for (total = 0, bin = zbin2; bin > 0; bin--) {
    total += h_k2->counts[bin];
    if (total > thresh) {
      break;
    }
  }
  min_k2 = h_k2->bins[bin];

  mris->Kmin = mris->Hmin = 10000.0f;
  mris->Kmax = mris->Hmax = -10000.0f;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const vertex = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vertex->ripflag) {
      continue;
    }
    if (vertex->k1 > max_k1) {
      vertex->k1 = max_k1;
      nthresh++;
    }
    if (vertex->k2 > max_k2) {
      vertex->k2 = max_k2;
      nthresh++;
    }
    if (vertex->k1 < min_k1) {
      vertex->k1 = min_k1;
      nthresh++;
    }
    if (vertex->k2 < min_k2) {
      vertex->k2 = min_k2;
      nthresh++;
    }
    vertex->K = vertex->k1 * vertex->k2;
    vertex->H = (vertex->k1 + vertex->k2) / 2.0;
    if (vertex->K < mris->Kmin) {
      mris->Kmin = vertex->K;
    }
    if (vertex->H < mris->Hmin) {
      mris->Hmin = vertex->H;
    }
    if (vertex->K > mris->Kmax) {
      mris->Kmax = vertex->K;
    }
    if (vertex->H > mris->Hmax) {
      mris->Hmax = vertex->H;
    }
  }

  HISTOfree(&h_k1);
  HISTOfree(&h_k2);
  if (nthresh > 0)
    fprintf(stderr,
            "%d vertices thresholded to be in k1 ~ "
            "[%2.2f %2.2f], k2 ~ [%2.2f %2.2f]\n",
            nthresh,
            min_k1,
            max_k1,
            min_k2,
            max_k2);

  return (NO_ERROR);
}

int MRIScomputeSecondFundamentalFormAtVertex(MRIS *mris, int vno, int *vertices, int vnum)
{
  int i, n, nbad = 0;
  VERTEX *vertex, *vnb;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse;
  VECTOR *v_z;
  static MATRIX *m_Q, *m_eigen;
  static VECTOR *v_c = NULL, *v_n, *v_e1, *v_e2, *v_yi;
  float k1, k2, evalues[3], a11, a12, a21, a22, cond_no, rsq, k, kmin, kmax;
  double ui, vi;

  if (mris->status == MRIS_PLANE) {
    return (NO_ERROR);
  }

  if (v_c == NULL) {
    v_c = VectorAlloc(3, MATRIX_REAL);
    v_n = VectorAlloc(3, MATRIX_REAL);
    v_e1 = VectorAlloc(3, MATRIX_REAL);
    v_e2 = VectorAlloc(3, MATRIX_REAL);
    v_yi = VectorAlloc(3, MATRIX_REAL);
    m_Q = MatrixAlloc(2, 2, MATRIX_REAL); /* the quadratic form */
    m_eigen = MatrixAlloc(2, 2, MATRIX_REAL);
  }

  vertex = &mris->vertices[vno];
  if (vertex->ripflag) {
    return (ERROR_BADPARM);
  }

  if (vno == 142915) {
    DiagBreak();
  }
  VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz);
  VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
  VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);

  if (vnum <= 0) {
    return (ERROR_BADPARM);
  }

  m_U = MatrixAlloc(vnum, 3, MATRIX_REAL);
  v_z = VectorAlloc(vnum, MATRIX_REAL);

  if (vno == Gdiag_no) {
    DiagBreak();
  }

  /* fit a quadratic form to the surface at this vertex */
  kmin = 10000.0f;
  kmax = -kmin;
  for (n = i = 0; i < vnum; i++) {
    vnb = &mris->vertices[vertices[i]];
    if (vnb->ripflag) {
      continue;
    }
    /*
      calculate the projection of this vertex onto the local tangent plane
    */
    VECTOR_LOAD(v_yi, vnb->x - vertex->x, vnb->y - vertex->y, vnb->z - vertex->z);
    ui = V3_DOT(v_yi, v_e1);
    vi = V3_DOT(v_yi, v_e2);
    *MATRIX_RELT(m_U, n + 1, 1) = ui * ui;
    *MATRIX_RELT(m_U, n + 1, 2) = 2 * ui * vi;
    *MATRIX_RELT(m_U, n + 1, 3) = vi * vi;
    VECTOR_ELT(v_z, n + 1) = V3_DOT(v_n, v_yi); /* height above TpS */
    rsq = ui * ui + vi * vi;
    if (!FZERO(rsq)) {
      k = VECTOR_ELT(v_z, n + 1) / rsq;
      if (k > kmax) {
        kmax = k;
      }
      if (k < kmin) {
        kmin = k;
      }
    }
    n++;
  }

  m_Ut = MatrixTranspose(m_U, NULL);        /* Ut */
  m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL); /* Ut U */
  cond_no = MatrixConditionNumber(m_tmp2);
  m_inverse = MatrixSVDInverse(m_tmp2, NULL); /* (Ut U)^-1 */
  if (!m_inverse) /* singular matrix - must be planar?? */
  {
    nbad++;
    evalues[0] = evalues[1] = 0.0;
  }
  else {
    m_tmp1 = MatrixMultiply(m_Ut, v_z, NULL); /* Ut z */
    MatrixMultiply(m_inverse, m_tmp1, v_c);   /* (Ut U)^-1 Ut z */

    /* now build Hessian matrix */
    *MATRIX_RELT(m_Q, 1, 1) = 2 * VECTOR_ELT(v_c, 1);
    *MATRIX_RELT(m_Q, 1, 2) = *MATRIX_RELT(m_Q, 2, 1) = 2 * VECTOR_ELT(v_c, 2);
    *MATRIX_RELT(m_Q, 2, 2) = 2 * VECTOR_ELT(v_c, 3);

    if (cond_no >= ILL_CONDITIONED) {
      vertex->k1 = k1 = kmax;
      vertex->k2 = k2 = kmin;
      // vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
      if (k1 * k2 < 0) {
        vertex->K = mris->Kmin;
        vertex->H = mris->Hmin;
      }
      else {
        vertex->K = mris->Kmax;
        vertex->H = mris->Hmax;
      }
      MatrixFree(&m_Ut);
      MatrixFree(&m_tmp2);
      MatrixFree(&m_U);
      VectorFree(&v_z);
      MatrixFree(&m_tmp1);
      MatrixFree(&m_inverse);
      return (ERROR_BADPARM);
    }

    /* the columns of m_eigen will be the eigenvectors of m_Q */
    if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL) {
      nbad++;
      MatrixSVDEigenValues(m_Q, evalues);
      vertex->k1 = k1 = evalues[0];
      vertex->k2 = k2 = evalues[1];
      vertex->K = k1 * k2;
      vertex->H = (k1 + k2) / 2;
      MatrixFree(&m_Ut);
      MatrixFree(&m_tmp2);
      MatrixFree(&m_U);
      VectorFree(&v_z);
      MatrixFree(&m_tmp1);
      MatrixFree(&m_inverse);
      return (ERROR_BADPARM);
    }

    MatrixFree(&m_tmp1);
    MatrixFree(&m_inverse);
  }
  k1 = evalues[0];
  k2 = evalues[1];
  vertex->k1 = k1;
  vertex->k2 = k2;
  vertex->K = k1 * k2;
  vertex->H = (k1 + k2) / 2;
  if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
    fprintf(stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n", vno, vertex->k1, vertex->k2, vertex->K, vertex->H);
  if (vertex->K < mris->Kmin) {
    mris->Kmin = vertex->K;
  }
  if (vertex->H < mris->Hmin) {
    mris->Hmin = vertex->H;
  }
  if (vertex->K > mris->Kmax) {
    mris->Kmax = vertex->K;
  }
  if (vertex->H > mris->Hmax) {
    mris->Hmax = vertex->H;
  }
  mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area;

  /* now update the basis vectors to be the principal directions */
  a11 = *MATRIX_RELT(m_eigen, 1, 1);
  a12 = *MATRIX_RELT(m_eigen, 1, 2);
  a21 = *MATRIX_RELT(m_eigen, 2, 1);
  a22 = *MATRIX_RELT(m_eigen, 2, 2);
  if (V3_LEN(v_e1) < 0.5) {
    DiagBreak();
  }
  vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21;
  vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21;
  vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21;
  vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22;
  vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22;
  vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22;
  if (SQR(vertex->e1x) + SQR(vertex->e1y) + SQR(vertex->e1z) < 0.5) {
    DiagBreak();
  }

  MatrixFree(&m_Ut);
  MatrixFree(&m_tmp2);
  MatrixFree(&m_U);
  VectorFree(&v_z);

  if (Gdiag & DIAG_SHOW && (nbad > 0)) {
    fprintf(stdout, "%d ill-conditioned points\n", nbad);
  }
  return (NO_ERROR);
}

// PerThreadMRIDistance is used to avoid locks and cache-bouncing below
//
PerThreadMRIDistance* makePerThreadMRIDistance(MRI const * const mri_distance) {
  PerThreadMRIDistance* ptd = (PerThreadMRIDistance*)malloc(sizeof(PerThreadMRIDistance));
  ptd->mri_distance 	= mri_distance;
  ptd->heightTimesDepth = mri_distance->height*mri_distance->depth;
  ptd->depth        	= mri_distance->depth;
  int size          	= mri_distance->width*mri_distance->height*mri_distance->depth;
  float* elts = ptd->elts = (float*)malloc(size*sizeof(float));
  int i;
  for (i = 0; i < size; i++) elts[i] = NPY;
  return ptd;
}

float* perThreadMRIDistanceElt(PerThreadMRIDistance* ptd, int i, int j, int k) {
  return &ptd->elts[i*ptd->heightTimesDepth + j*ptd->depth + k];
}

void freePerThreadMRIDistance(PerThreadMRIDistance** ptdp) {
  PerThreadMRIDistance* ptd = *ptdp;
  *ptdp = NULL;
  if (!ptd) return;
  freeAndNULL(ptd->elts);
  freeAndNULL(ptd);
}


void updateDistanceElt(volatile float* f, float distance, bool lockNeeded) {
    
#ifdef HAVE_OPENMP
  if (lockNeeded) {
    if (fabsf(distance) > fabsf(*f)) return;  	// avoid the lock if possible
    #pragma omp critical
    updateDistanceElt(f, distance, false);
    return;
  } 
#endif

  // There was a reproducibility problem here.
  // If the smallest positive distance and the smallest negative distance is the same fabs()
  // then this code would randomly choose between them.
  //
  // The solution is to have the positive be the preferred of two equal values.
  //
  if (fabsf(distance) < fabsf(*f)) {
    *f = distance;	            
  } else if (fabsf(distance) == fabsf(*f)) {
    // They are equal.  Prefer the positive.
    if (distance > 0) *f = distance;
  }
}


void updateDistanceEltFromSignArgAndSquareLockNeeded(volatile float* f, float distanceSignArg, float distanceSquared) {
  // Avoid calculating the sqrt unless definitely needed
  //
  // The obvious test has problems if sqrtf(distanceSquared) == f and distanceSign is positive, because this would reject that solution
  // when the above code would prefer it.  Hence the 1.01f margin of error.
  //
  // if (squaref(*f) < distanceSquared) return;
  //
  if (squaref(*f)*1.01f < distanceSquared) return;
  
  updateDistanceElt(f, SIGN(distanceSignArg)*sqrtf(distanceSquared), true);
}

void updateDistanceEltFromSignArgAndSquareNoLockNeeded(
    volatile float* f, 
    float distanceSignArg, 
    float distanceSquared,
    float sqrtfDistanceSquared) {
  // This function is performance-critical for mris_fix_topology, hence hand-optimized
  //
  // Avoid waiting for the sqrt unless definitely needed
  //
  // The obvious test has problems if sqrtf(distanceSquared) == f and distanceSign is positive, because this would reject that solution
  // when the above code would prefer it.  Hence the 1.01f margin of error.
  //
  // if (squaref(*f) < distanceSquared) return;
  //
  if (squaref(*f)*1.01f < distanceSquared) return;          // don't wait for the sqrtf to complete if it is not needed
  
  float positiveDistance = sqrtfDistanceSquared;
  if (positiveDistance > fabsf(*f)) return;

  float distance = SIGN(distanceSignArg)*positiveDistance;
  
  if (positiveDistance < fabsf(*f)) {
    *f = distance;	            
  } else if (positiveDistance == fabsf(*f)) {
    // They are equal.  Prefer the positive.
    if (distance > 0) *f = distance;
  }
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseAreaErrors(MRIS *mris)
{
  int vno, fi, n;
  float area, orig_area;
  FACE *face;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    area = orig_area = 0.0f;

    /* use average area of all faces this vertex is part of -
       this is not really correct, but should be good enough for
       visualization purposes.
    */
    for (n = fi = 0; fi < vertext->num; fi++) {
      face = &mris->faces[vertext->f[fi]];
      area += face->area;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, vertext->f[fi]);
      orig_area += fNorm->orig_area;
    }
    vertex->curv = (area - orig_area) / (float)n;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISuseGaussianCurvature(MRIS *mris)
{
  int vno;
  VERTEX *vertex;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = vertex->K;
  }

  mris->min_curv = mris->Kmin;
  mris->max_curv = mris->Kmax;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  int MRISuseMeanCurvature(MRIS *mris)
  Set vertex->curv = vertex->H ;
  ------------------------------------------------------*/
int MRISuseMeanCurvature(MRIS *mris)
{
  int vno;
  VERTEX *vertex;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = vertex->H;
  }

  mris->min_curv = mris->Hmin;
  mris->max_curv = mris->Hmax;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISusePrincipalCurvature(MRIS *mris)
{
  int vno;
  VERTEX *vertex;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = vertex->k1;
  }

  mris->min_curv = mris->Hmin;
  mris->max_curv = mris->Hmax;
  return (NO_ERROR);
}

int MRISuseK1Curvature(MRIS *mris)
{
  int vno;
  VERTEX *vertex;

  float f_min = mris->vertices[0].curv;
  float f_max = f_min;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = vertex->k1;
    if (vertex->curv < f_min) {
      f_min = vertex->curv;
    }
    if (vertex->curv > f_max) {
      f_max = vertex->curv;
    }
  }

  mris->min_curv = f_min;
  mris->max_curv = f_max;
  return (NO_ERROR);
}

int MRISuseK2Curvature(MRIS *mris)
{
  int vno;
  VERTEX *vertex;

  float f_min = mris->vertices[0].curv;
  float f_max = f_min;

  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    vertex->curv = vertex->k2;
    if (vertex->curv < f_min) {
      f_min = vertex->curv;
    }
    if (vertex->curv > f_max) {
      f_max = vertex->curv;
    }
  }

  mris->min_curv = f_min;
  mris->max_curv = f_max;
  return (NO_ERROR);
}

int MRISusePrincipalCurvatureFunction(MRIS *pmris, float (*f)(float k1, float k2))
{
  //
  // PRECONDITIONS
  //  o The principal curvatures k1 and k2 for each vertex point on the
  //    surface have been defined.
  //
  // POSTCONDITIONS
  //  o Each vertex 'curv' value is replaced with the result from the
  //    the (*f)(k1, k2) function.
  //  o Surface min and max values are set appropriately.
  //

  int vno;
  VERTEX *pvertex;
  float f_k1;
  float f_k2;

  float f_min = pmris->vertices[0].curv;
  float f_max = f_min;

  for (vno = 0; vno < pmris->nvertices; vno++) {
    pvertex = &pmris->vertices[vno];
    if (pvertex->ripflag) {
      continue;
    }
    f_k1 = pvertex->k1;
    f_k2 = pvertex->k2;
    pvertex->curv = (*f)(f_k1, f_k2);
    if (pvertex->curv < f_min) {
      f_min = pvertex->curv;
    }
    if (pvertex->curv > f_max) {
      f_max = pvertex->curv;
    }
  }

  pmris->min_curv = f_min;
  pmris->max_curv = f_max;
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute the folding (FI) and intrinsic curvature (ICI)
  indices of the surface.

  Drury and Van Essen:  ICI = 56  (lh) and 54  (rh)
  FI = 500 (lh) and 520 (rh)
  ------------------------------------------------------*/
int MRIScomputeCurvatureIndices(MRIS *mris, double *pici, double *pfi)
{
  int vno;
  VERTEX *vertex;
  double k1, k2, ici, fi, area, Kmax, Kmin;

  ici = fi = 0.0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    k1 = (double)vertex->k1;
    k2 = (double)vertex->k2;
    area = (double)vertex->area;
    if (vertex->K > 0) {
      ici += area * (double)vertex->K;
    }
    Kmax = (double)fabs(k1);
    Kmin = (double)fabs(k2);
    fi += area * Kmax * (Kmax - Kmin);
  }

  *pfi = fi / (4.0 * M_PI);
  *pici = ici / (4.0 * M_PI);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRISmaxRadius(MRIS *mris)
{
  double radius;
  int vno, n;
  VERTEX *vertex;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0, r;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = (double)vertex->x;
    y = (double)vertex->y;
    z = (double)vertex->z;
    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
  }
  mris->xhi = xhi;
  mris->xlo = xlo;
  mris->yhi = yhi;
  mris->ylo = ylo;
  mris->zhi = zhi;
  mris->zlo = zlo;
  x0 = (xlo + xhi) / 2.0f;
  y0 = (ylo + yhi) / 2.0f;
  z0 = (zlo + zhi) / 2.0f;
  for (radius = 0.0, n = vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    n++;
    x = (double)vertex->x - x0;
    y = (double)vertex->y - y0;
    z = (double)vertex->z - z0;
    r = sqrt(x * x + y * y + z * z);
    if (r > radius) {
      radius = r;
    }
  }

  return (radius);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRISaverageRadius(MRIS *mris)
{
  double radius;
  int vno, n;
  VERTEX *vertex;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = (double)vertex->x;
    y = (double)vertex->y;
    z = (double)vertex->z;
    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
    if (zlo < -1000) {
      DiagBreak();
    }
  }
  x0 = (xlo + xhi) / 2.0f;
  y0 = (ylo + yhi) / 2.0f;
  z0 = (zlo + zhi) / 2.0f;
  for (radius = 0.0, n = vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    n++;
    x = (double)vertex->x - x0;
    y = (double)vertex->y - y0;
    z = (double)vertex->z - z0;
    radius += sqrt(x * x + y * y + z * z);
  }

  radius /= (double)n;
  return (radius);
}


/* --------------------------------------------------------------------
   MRISrescaleMetricProperties() - rescale metric properties (area,
   dist) of group surfaces so that they match that of the average of
   the input set. Does not change the vertex xyz. Requires that
   TBD() has been run. Returns the distance scaling factor. 
   Has no effect if surface is not a group surface.
   -------------------------------------------------------------------- */
double MRISrescaleMetricProperties(MRIS *surf)
{
  int VtxNo, nthNNbr, nNNbrs, NbrVtxNo;
  double scale;

  if (surf->group_avg_surface_area == 0) {
    return (0.0);  // not a group surface
  }

  scale = sqrt((double)surf->group_avg_surface_area / surf->total_area);
  printf("scale = %lf\n", scale);

  for (VtxNo = 0; VtxNo < surf->nvertices; VtxNo++) {
    VERTEX * const vtx1 = &surf->vertices[VtxNo];
    if (vtx1->ripflag) {
      continue;
    }
    vtx1->area *= (scale * scale);
    nNNbrs = surf->vertices_topology[VtxNo].vnum;
    for (nthNNbr = 0; nthNNbr < nNNbrs; nthNNbr++) {
      NbrVtxNo = surf->vertices_topology[VtxNo].v[nthNNbr];
      VERTEX const * const vtx2 = &surf->vertices[NbrVtxNo];
      if (vtx2->ripflag) {
        continue;
      }
      vtx1->dist[nthNNbr] *= scale;
    }
  }
  surf->total_area *= (scale * scale);
  mrisSetAvgInterVertexDist(surf, surf->avg_vertex_dist * scale);
  surf->std_vertex_dist *= scale;
  surf->avg_vertex_area *= (scale * scale);
  return (scale);
}


int MRISnormalize(MRIS *mris, int dof, int which)
{
  int vno;
  VERTEX *v;
  float fdof = (float)dof, mean;

  if (dof <= 0) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRISnormalize: invalid dof %d", dof));

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    switch (which) {
      case VERTEX_AREA:
        v->origarea /= fdof;
        mean = v->origarea;
        break;
      case VERTEX_CURV:
        v->curv /= fdof;
        mean = v->curv;
        break;
      default:
      case VERTEX_LOGODDS:
      case VERTEX_VAL:
        v->val /= fdof;
        mean = v->val;
        break;
    }
    v->val2 = v->val2 / fdof - mean * mean;
  }
  return (NO_ERROR);
}


double MRISaverageCanonicalRadius(MRIS *mris)
{
  double radius;
  int vno, n;
  VERTEX *vertex;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;
  for (vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    x = (double)vertex->cx;
    y = (double)vertex->cy;
    z = (double)vertex->cz;
    if (x > xhi) {
      xhi = x;
    }
    if (x < xlo) {
      xlo = x;
    }
    if (y > yhi) {
      yhi = y;
    }
    if (y < ylo) {
      ylo = y;
    }
    if (z > zhi) {
      zhi = z;
    }
    if (z < zlo) {
      zlo = z;
    }
  }
  x0 = (xlo + xhi) / 2.0f;
  y0 = (ylo + yhi) / 2.0f;
  z0 = (zlo + zhi) / 2.0f;
  for (radius = 0.0, n = vno = 0; vno < mris->nvertices; vno++) {
    vertex = &mris->vertices[vno];
    n++;
    x = (double)vertex->cx - x0;
    y = (double)vertex->cy - y0;
    z = (double)vertex->cz - z0;
    radius += sqrt(x * x + y * y + z * z);
  }

  radius /= (double)n;
  return (radius);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISupdateSurface(MRIS *mris)
{
  MRIScomputeMetricProperties(mris);
  MRIScomputeSecondFundamentalForm(mris);
  return (NO_ERROR);
}


#if MATRIX_ALLOCATION

MATRIX *VoxelFromSRASmatrix = NULL;

int mriSurfaceRASToVoxel(
    double xr, double yr, double zr, 
    double *xv, double *yv, double *zv)
{
  static int once;
  static bool try_both_way = false;
  static bool try_new_way  = true;
  static bool try_old_way  = false;
  if (!once) {
    once++; 
    try_both_way = getenv("FREESURFER_mriSurfaceRASToVoxel_both");
    if  (try_both_way) 
        try_new_way = try_old_way = true;
    else if (getenv("FREESURFER_mriSurfaceRASToVoxel_old"))
        try_new_way = false, try_old_way = true;
  }

  bool const new_way = 
    try_new_way &&
      VoxelFromSRASmatrix->type == MATRIX_REAL &&
      VoxelFromSRASmatrix->rows == 4 &&
      VoxelFromSRASmatrix->cols == 4;

  bool const old_way = 
    try_old_way ||
    !new_way;

  float nx = 0, ny = 0, nz = 0;
  if (new_way) {
    float fxr = xr, fyr = yr, fzr = zr;     // cvt to float to match old way
    nx = 
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 1,1)*fxr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 1,2)*fyr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 1,3)*fzr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 1,4) ;
    ny = 
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 2,1)*fxr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 2,2)*fyr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 2,3)*fzr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 2,4) ;
    nz = 
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 3,1)*fxr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 3,2)*fyr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 3,3)*fzr +
      (double)*MATRIX_RELT(VoxelFromSRASmatrix, 3,4) ;
  }
  
  if (!old_way) {
    *xv = nx;
    *yv = ny;
    *zv = nz;
    return (NO_ERROR);
  }

  VECTOR* sr = VectorAlloc(4, MATRIX_REAL);
  V4_LOAD(sr, xr, yr, zr, 1.);

  VECTOR* vv = MatrixMultiply(VoxelFromSRASmatrix, sr, NULL);

  float ox = V3_X(vv);
  float oy = V3_Y(vv);
  float oz = V3_Z(vv);

  VectorFree(&sr);
  VectorFree(&vv);

  if (new_way) {
    if (closeEnough(ox,nx) || closeEnough(oy,ny) || closeEnough(oz, nz)) {
      fprintf(stderr, "%s:%d not same answer\n", __FILE__, __LINE__);
      *(int*)-1 = 0;
    }
  }

  *xv = ox;
  *yv = oy;
  *zv = oz;

  return (NO_ERROR);
}

#endif


//////////////////////////////////////////////////////////////
// Find inside and outside voxels of a surface mris
//////////////////////////////////////////////////////////////
MRI *MRISbinarizeVolume(MRIS *mris, MRI_REGION *region, float resolution, float distance_from_surface)
{
  int k, i, j, p, width, height, depth;
  float x0, x1, x2, y0, y1, y2, z0, z1, z2;
  float x, y, z;

  int fno, fn1;
  int delta;
  FACE *face;
  int imin, imax, jmin, jmax, kmin, kmax;
  float distance, sign;
  float n_f[3], n_e0[3], n_e1[3], n_e2[3], n_v0[3], n_v1[3], n_v2[3];
  float vec[3], vec0[3], vec1[3], vec2[3], e0[3], e1[3], e2[3], n0[3], n1[3], n2[3];
  float val, valu, val0, val1, val2;
  MRI *mri_distance;

  /* resolution*/
  if (resolution < 1.0f) {
    resolution = 1.0f;
  }
  /* distance */
  if (distance_from_surface < 2.0f) {
    distance_from_surface = 2.0f;
  }
  /* look at approximately +/- dist mm */
  delta = (int)(distance_from_surface / resolution);
  if (delta < 2) {
    delta = 2;
  }

  // find the region of interest in this coordinate system

  /* allocate the volume */
  width = ceil(resolution * (region->dx));
  height = ceil(resolution * (region->dy));
  depth = ceil(resolution * (region->dz));
  mri_distance = MRIalloc(width, height, depth, MRI_FLOAT);

  fprintf(WHICH_OUTPUT, "mri volume size : %d by %d by %d (resolution = %d)\n", width, height, depth, (int)resolution);

  mri_distance->xstart = region->x;
  mri_distance->xsize = resolution;

  mri_distance->ystart = region->y;
  mri_distance->ysize = resolution;

  mri_distance->zstart = region->z;
  mri_distance->zsize = resolution;

  /* optimize by listing the concerned faces */

  for (k = 0; k < mri_distance->depth; k++)
    for (j = 0; j < mri_distance->height; j++)
      for (i = 0; i < mri_distance->width; i++) {
        MRIFvox(mri_distance, i, j, k) = NPY;
      }

  ROMP_PF_begin  
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) 
#endif
  for (p = 0; p < mris->nfaces; p++) {
    ROMP_PFLB_begin
    computeDefectFaceNormal(mris, p);
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  /* find the distance to each surface voxels */
  for (p = 0; p < mris->nfaces; p++) {
    fno = p;  // tp->faces[p];
    face = &mris->faces[fno];

    // calculate three vertices
    x0 = mris->vertices[face->v[0]].x;
    y0 = mris->vertices[face->v[0]].y;
    z0 = mris->vertices[face->v[0]].z;
    x1 = mris->vertices[face->v[1]].x;
    y1 = mris->vertices[face->v[1]].y;
    z1 = mris->vertices[face->v[1]].z;
    x2 = mris->vertices[face->v[2]].x;
    y2 = mris->vertices[face->v[2]].y;
    z2 = mris->vertices[face->v[2]].z;

    /* find the bounding box */
    imin = iVOL(mri_distance, MIN3(x0, x1, x2)) - delta;
    imax = iVOL(mri_distance, MAX3(x0, x1, x2)) + delta;

    jmin = jVOL(mri_distance, MIN3(y0, y1, y2)) - delta;
    jmax = jVOL(mri_distance, MAX3(y0, y1, y2)) + delta;

    kmin = kVOL(mri_distance, MIN3(z0, z1, z2)) - delta;
    kmax = kVOL(mri_distance, MAX3(z0, z1, z2)) + delta;


    imin = MAX(imin, 0);
    imax = MIN(imax, mri_distance->width - 1);

    jmin = MAX(jmin, 0);
    jmax = MIN(jmax, mri_distance->height - 1);

    kmin = MAX(kmin, 0);
    kmax = MIN(kmax, mri_distance->depth - 1);

    /* generating the pseudo-normals for edges and vertices */
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    n_f[0] = fNorm->nx;
    n_f[1] = fNorm->ny;
    n_f[2] = fNorm->nz;

    /* edge0: x0 <--> x1 */
    e0[0] = x1 - x0;
    e0[1] = y1 - y0;
    e0[2] = z1 - z0;
    F_CROSS(n_f, e0, n0);
    fn1 = findOtherEdgeFace(mris, fno, face->v[0], face->v[1]);
    FaceNormCacheEntry const * fNorm0 = getFaceNorm(mris, fn1);
    n_e0[0] = fNorm->nx + fNorm0->nx;
    n_e0[1] = fNorm->ny + fNorm0->ny;
    n_e0[2] = fNorm->nz + fNorm0->nz;

    /* edge1: x1 <--> x2 */
    e1[0] = x2 - x1;
    e1[1] = y2 - y1;
    e1[2] = z2 - z1;
    F_CROSS(n_f, e1, n1);
    fn1 = findOtherEdgeFace(mris, fno, face->v[1], face->v[2]);
    FaceNormCacheEntry const * fNorm1 = getFaceNorm(mris, fn1);
    n_e1[0] = fNorm->nx + fNorm1->nx;
    n_e1[1] = fNorm->ny + fNorm1->ny;
    n_e1[2] = fNorm->nz + fNorm1->nz;

    /* edge2: x2 <--> x0 */
    e2[0] = x0 - x2;
    e2[1] = y0 - y2;
    e2[2] = z0 - z2;
    F_CROSS(n_f, e2, n2);
    fn1 = findOtherEdgeFace(mris, fno, face->v[2], face->v[0]);
    FaceNormCacheEntry const * fNorm2 = getFaceNorm(mris, fn1);
    n_e2[0] = fNorm->nx + fNorm2->nx;
    n_e2[1] = fNorm->ny + fNorm2->ny;
    n_e2[2] = fNorm->nz + fNorm2->nz;

    /* vertex pseudo-normals */
    computeVertexPseudoNormal(mris, face->v[0], n_v0, 0);
    computeVertexPseudoNormal(mris, face->v[1], n_v1, 0);
    computeVertexPseudoNormal(mris, face->v[2], n_v2, 0);

/* finding distance to surface */
    for (k = kmin; k <= kmax; k++)
      for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
          x = xSURF(mri_distance, i);
          y = ySURF(mri_distance, j);
          z = zSURF(mri_distance, k);

          vec0[0] = x - x0;
          vec0[1] = y - y0;
          vec0[2] = z - z0;
          vec1[0] = x - x1;
          vec1[1] = y - y1;
          vec1[2] = z - z1;
          vec2[0] = x - x2;
          vec2[1] = y - y2;
          vec2[2] = z - z2;
          vec[0] = (vec0[0] + vec1[0] + vec2[0]) / 3.0;
          vec[1] = (vec0[1] + vec1[1] + vec2[1]) / 3.0;
          vec[2] = (vec0[2] + vec1[2] + vec2[2]) / 3.0;

          /* compute distance to face */
          /* where is the point */
          val0 = F_DOT(vec0, n0);
          val1 = F_DOT(vec1, n1);
          val2 = F_DOT(vec2, n2);

          if ((val0 >= 0) && (val1 >= 0) && (val2 >= 0)) {
            /* the projection of the vertex is inside */
            val = F_DOT(n_f, vec);
            valu = 1;
            sign = val;
            distance = val; /* n_f is already normalized */
          }
          else {
            distance = NPY;
            sign = 0;
            valu = 0;

            if (val0 <= 0) {
              /* compute distance to edge0 */
              val = F_DOT(vec0, e0);
              if (val < 0) {
                /* closer to x0 */
                sign = F_DOT(n_v0, vec0);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
              }
              else if (val < SQR3(e0)) {
                /* closer to edge0 */
                sign = F_DOT(n_e0, vec0);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec0) - SQR(val) / SQR3(e0))));
              }
              else {
                /* closer to x1 */
                sign = F_DOT(n_v1, vec1);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
              }
            };
            if (val1 <= 0) {
              val = F_DOT(vec1, e1);
              if (val < 0) {
                /* closer to x1 */
                sign = F_DOT(n_v1, vec1);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
              }
              else if (val < SQR3(e1)) {
                /* closer to edge1 */
                sign = F_DOT(n_e1, vec1);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec1) - SQR(val) / SQR3(e1))));
              }
              else {
                /* closer to x2 */
                sign = F_DOT(n_v2, vec2);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
              }
            };
            if (val2 <= 0) {
              val = F_DOT(vec2, e2);
              if (val < 0) {
                /* closer to x2 */
                sign = F_DOT(n_v2, vec2);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
              }
              else if (val < SQR3(e2)) {
                /* closer to edge2 */
                sign = F_DOT(n_e2, vec2);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec2) - SQR(val) / SQR3(e2))));
              }
              else {
                /* closer to x0 */
                sign = F_DOT(n_v0, vec0);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
              }
            };
          }

          /* update distance map */
          if (fabs(distance) < fabs(MRIFvox(mri_distance, i, j, k))) {
            MRIFvox(mri_distance, i, j, k) = distance;
          }
        }
  }


  return mri_distance;
}

/* This function finds the orientation changes in the triangles
   constituting the tessellation. Writes the orientation changes in
   curv. Note that this function can easily be modified to check the
   orientation changes for any tessellations (not only triangulations) */
int MRISmarkOrientationChanges(MRIS *mris)
{
  int face1, face2, vn0, p, vn1, vn2, d1, d2, count;

  /* to avoid compiler warnings */
  face1 = face2 = 0;
  d1 = d2 = 0;

  /* erase curvature for diag purposes */
  MRISclearCurvature(mris);

  for (count = 0, vn0 = 0; vn0 < mris->nvertices; vn0++) {
    int nfaces, n;
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vn0];
    VERTEX                * const v  = &mris->vertices         [vn0];
    if (v->ripflag) {
      continue;
    }

    /* go through neighbors */
    if (vt->vnum == 0) {
      //                        continue;
      fprintf(WHICH_OUTPUT, "vertex %d without neighbors ! \n", vn0);
      v->curv = -5;
    }
    if (vt->num == 0) {
      fprintf(WHICH_OUTPUT, "vertex %d without faces ! \n", vn0);
    }

    for (p = 0; p < vt->vnum; p++) {
      vn1 = vt->v[p];
      if (mris->vertices[vn1].ripflag) fprintf(WHICH_OUTPUT, "vertex %d connected to ripped vertex %d !\n ", vn0, vn1);

      /* find the neighboring faces of vn0 and vn1 */
      nfaces = 0;
      for (n = 0; n < vt->vnum; n++) {
        vn2 = vt->v[n];

        if (vn2 == vn1) {
          continue;
        }

        if (isFace(mris, vn0, vn1, vn2)) {
          nfaces++;
          if (nfaces == 1) {
            face1 = findFace(mris, vn0, vn1, vn2);
          }
          else {
            face2 = findFace(mris, vn0, vn1, vn2);
          }
        }
      }

      if (nfaces == 1) {
        fprintf(WHICH_OUTPUT, "edge (%d<-->%d) has one single face %d !\n", vn0, vn1, face1);
        v->curv = -2;
      }

      if (nfaces == 0) {
        fprintf(WHICH_OUTPUT, "edge (%d<-->%d) has no face !\n", vn0, vn1);
        v->curv = -3;
      }

      if (nfaces > 2) {
        fprintf(WHICH_OUTPUT,
                "edge (%d<-->%d) has too many faces "
                "(at least : %d %d) !\n",
                vn0,
                vn1,
                face1,
                face2);
        v->curv = -4;
      }
      /* check the orientation for face 1 and face 2 */
      if (nfaces == 2 && (computeOrientation(mris, face1, vn0, vn1) * computeOrientation(mris, face2, vn0, vn1) > 0)) {
        v->curv = -1;
        fprintf(WHICH_OUTPUT,
                "\ndefective orientation at "
                "vertex %d(%d) with faces %d and %d\n",
                vn0,
                vn1,
                face1,
                face2);
        count++;
      }
    }
  }
  fprintf(WHICH_OUTPUT,
          "\n%2.3f %% of the vertices (%d vertices) "
          "exhibit an orientation change\n",
          100. * count / mris->nvertices,
          count);

  return count;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#if (!SPHERE_INTERSECTION)
static int mrisComputeCanonicalEdgeBasis(
    MRIS *mris, EDGE *edge1, EDGE *edge2, double origin[3], double e0[3], double e1[3])
{
  VERTEX *v0, *v1, *v2, *v3;
  double len, normal[3];
  float fx, fy, fz;

  v0 = &mris->vertices[edge1->vno1];
  v1 = &mris->vertices[edge1->vno2];
  v2 = &mris->vertices[edge2->vno1];
  v3 = &mris->vertices[edge2->vno2];
  fx = (v0->cx + v1->cx + v2->cx + v3->cx) / 4;
  fy = (v0->cy + v1->cy + v2->cy + v3->cy) / 4;
  fz = (v0->cz + v1->cz + v2->cz + v3->cz) / 4;
  normal[0] = origin[0] = (double)fx;
  normal[1] = origin[1] = (double)fy;
  normal[2] = origin[2] = (double)fz;
  len = 1.0f / VLEN(normal);
  SCALAR_MUL(normal, len, normal);

  /* pick any non-parallel vector and cross it with normal */
  e1[0] = normal[1];
  e1[1] = normal[2];
  e1[2] = normal[0];
  CROSS(e0, normal, e1);
  if ((VZERO(e0))) /* happened to pick parallel vector */
  {
    e1[0] = normal[1];
    e1[1] = -normal[2];
    e1[2] = normal[0];
    CROSS(e0, normal, e1);
  }
  CROSS(e1, e0, normal);
  len = 1.0f / VLEN(e0);
  SCALAR_MUL(e0, len, e0);
  len = 1.0f / VLEN(e1);
  SCALAR_MUL(e1, len, e1);
  len = DOT(e0, e1);
  return (NO_ERROR);
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisDumpDefectiveEdge(MRIS *mris, int vno1, int vno2)
{
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisCheckSurface(MRIS *mris)
{
  int    fno, vno, n, nfaces, m, vno2, nbad, flist[100];

  FACE   *f ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    for (n = 0; n < vt->num; n++) 
      if (vt->f[n] >= mris->nfaces)
	DiagBreak() ;
    for (n = 0; n < vt->vtotal; n++) 
      if (vt->v[n] >= mris->nvertices)
	DiagBreak() ;
  }

  for (fno = 0; fno < mris->nfaces; fno++) 
  {
    f = &mris->faces[fno];
    if (f->ripflag)
      continue ;
    for (m = 0; m < VERTICES_PER_FACE; m++) 
    {
      VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[f->v[m]];
      for (n = 0; n < vt->num && fno != vt->f[n]; n++) 
      {
        ;
      }
      if (n == vt->num) /* face has vertex, but vertex doesn't have face */
      {
        printf("mrisCheckSurface: %s: face[%d].v[%d] = %d, "
                  "but face %d not in vertex %d "
                  "face list\n",
                  mris->fname.data(),
                  fno,
                  m,
                  f->v[m],
                  fno,
                  f->v[m]);
	DiagBreak() ;
      }
    }
  }

  return(NO_ERROR) ;

  /*  fprintf(stdout, "\n") ;*/
  for (nbad = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[f->v[m]];
    VERTEX          const * const v  = &mris->vertices         [f->v[m]];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < vt->vnum; n++) {
      vno2 = vt->v[n];
      if (vno2 < vno) {
        continue;
      }
      for (nfaces = m = 0; m < vt->vnum; m++) {
        if (vt->v[m] == vno2) {
          continue;
        }
        if (vertexNeighbor(mris, vno2, vt->v[m]) && isFace(mris, vno, vno2, vt->v[m])) {
          flist[nfaces] = findFace(mris, vno, vno2, vt->v[m]);
          nfaces++;
        }
      }
      if (((mris->patch == 0) && (nfaces != 2)) || (mris->patch && (nfaces < 1 || nfaces > 2))) {
        int i;

        nbad++;
        fprintf(stdout, "edge %d <--> %d has %d face(s)! ", vno, vno2, nfaces);
        fprintf(stdout, "(");
        for (i = 0; i < nfaces; i++) {
          fprintf(stdout, "%d%s", flist[i], i < nfaces - 1 ? "," : "");
        }
        fprintf(stdout, ")\n");
        fprintf(stdout,
                "%d %d %d !!!\n",
                mris->faces[flist[0]].v[0],
                mris->faces[flist[0]].v[1],
                mris->faces[flist[0]].v[2]);
        mrisDumpDefectiveEdge(mris, vno, vno2);
        DiagBreak();
      }
    }
  }
  fprintf(stdout, "%d defective edges\n", nbad);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISmarkNegativeVertices(MRIS *mris, int mark)
{
  int fno, n;
  FACE *f;

  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->area < 0)
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mris->vertices[f->v[n]].marked = mark;
      }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISripNegativeVertices(MRIS *mris)
{
  int fno, n;
  FACE *f;

  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->area < 0) {
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mris->vertices[f->v[n]].ripflag = 1;
      }
      f->ripflag = 1;
    }
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRIScomputeAverageCurvature(MRIS *mris, double *psigma)
{
  double mean, var, total, total_sq, nv, d;
  int vno;
  VERTEX *v;

  for (vno = 0, nv = total_sq = total = 0.0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    d = (double)v->curv;
    total_sq += d * d;
    total += d;
    nv += 1.0;
  }
  if (nv) {
    mean = total / nv;
    var = total_sq / nv - (mean * mean);
  }
  else {
    var = mean = 0.0;
  }
  if (psigma) {
    *psigma = sqrt(var);
  }
  return (mean);
}


int MRISmaskLabel(MRIS *mris, LABEL *area)
{
  int i;
  VERTEX *v;

  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->curv = v->stat = v->val = v->imag_val = v->val2 = v->valbak = v->val2bak = 0.0;
  }
  return (NO_ERROR);
}

int MRISmaskNotLabel(MRIS *mris, LABEL *area)
{
  int i, vno;
  VERTEX *v;

  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->marked = 1;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked) {
      continue;
    }
    v->curv = v->stat = v->val = v->imag_val = v->val2 = v->valbak = v->val2bak = 0.0;
  }
  for (i = 0; i < area->n_points; i++) {
    v = &mris->vertices[area->lv[i].vno];
    v->marked = 0;
  }
  return (NO_ERROR);
}

int MRISsubsampleDist(MRIS *mris, float spacing)
{
  int k, m, n, sub_num;

  sub_num = 0;
  for (k = 0; k < mris->nvertices; k++) {
    VERTEX * const v = &mris->vertices[k];
    v->d = 10000;
    v->val = 0;
  }
  for (k = 0; k < mris->nvertices; k++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
    VERTEX                * const v  = &mris->vertices         [k];
    for (m = 0; m < vt->vnum; m++) {
      if (mris->vertices[vt->v[m]].d + 1 < v->d) {
        v->d = mris->vertices[vt->v[m]].d + 1;
      }
    }
    if (v->d >= spacing) {
      v->d = 0;
      v->val = 1;
      sub_num++;
    }
    for (m = 0; m < vt->vnum; m++) {
      if (mris->vertices[vt->v[m]].d > v->d + 1) {
        mris->vertices[vt->v[m]].d = v->d + 1;
      }
    }
  }
  for (k = mris->nvertices - 1; k >= 0; k--) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
    VERTEX                * const v  = &mris->vertices         [k];
    for (m = 0; m < vt->vnum; m++) {
      if (mris->vertices[vt->v[m]].d + 1 < v->d) {
        v->d = mris->vertices[vt->v[m]].d + 1;
      }
      if (mris->vertices[vt->v[m]].d > v->d + 1) {
        mris->vertices[vt->v[m]].d = v->d + 1;
      }
    }
  }

  if (spacing == 2)
    for (k = 0; k < mris->nvertices; k++)
      if (mris->vertices[k].d > 0) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
        VERTEX                * const v  = &mris->vertices         [k];
        n = 0;
        for (m = 0; m < vt->vnum; m++) {
          if (mris->vertices[vt->v[m]].d == 0) {
            n++;
          }
        }
        if (n <= 2) {
          v->d = 0;
          v->val = 1;
          v->fixedval = TRUE;
          sub_num++;
        }
        for (m = 0; m < vt->vnum; m++) {
          if (mris->vertices[vt->v[m]].d > v->d + 1) {
            mris->vertices[vt->v[m]].d = v->d + 1;
          }
        }
      }

  return (sub_num);
}


/*-----------------------------------------------------------------------
  MRISarN() - computes spatial autocorrelation function at each vertex
  by averaging the ARs within the neighborhood of a vertex. Note: does
  not try to take into account different distances between
  neighbors. N is the number of hops.  arN will have N frames where
  frame 0 is always 1, frame 1 is the average AR between the vertex
  and the vertices 1 hop away, etc.
  -----------------------------------------------------------------------*/
MRI *MRISarN(MRIS *surf, MRI *src, MRI *mask, MRI *arN, int N)
{
  int **crslut, nvox, vtx;

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox) {
    printf("ERROR: MRISarN: Surf/Src dimension mismatch.\n");
    return (NULL);
  }

  if (arN == NULL) {
    arN = MRIcloneBySpace(src, MRI_FLOAT, N);
    if (arN == NULL) {
      printf("ERROR: could not alloc\n");
      return (NULL);
    }
  }

  // Build LUT to map from col,row,slice to vertex
  crslut = MRIScrsLUT(surf, src);

  ROMP_PF_begin
#ifdef _OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    ROMP_PFLB_begin
    
    int nnbrs, frame, nbrvtx, nthnbr, c, r, s;
    int cnbr, rnbr, snbr, nnbrs_actual;
    double valvtx, valnbr, arsum, sumsqvtx, vtxvar, sumsqnbr, sumsqx, nbrvar;
    SURFHOPLIST *shl;
    int nthhop;

    if (surf->vertices[vtx].ripflag) continue;
    c = crslut[0][vtx];
    r = crslut[1][vtx];
    s = crslut[2][vtx];
    if (mask)
      if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;

    // Compute variance for this vertex
    sumsqvtx = 0;
    for (frame = 0; frame < src->nframes; frame++) {
      valvtx = MRIFseq_vox(src, c, r, s, frame);
      sumsqvtx += (valvtx * valvtx);
    }
    if (sumsqvtx == 0) continue;  // exclude voxels with 0 variance
    vtxvar = sumsqvtx / src->nframes;

    // Zero hop is always 1
    MRIFseq_vox(arN, c, r, s, 0) = 1;

    // Create structure to manage the multiple hops for this vertex
    shl = SetSurfHopList(vtx, surf, N);

    // loop through hops
    for (nthhop = 1; nthhop < N; nthhop++) {
      nnbrs = shl->nperhop[nthhop];
      nnbrs_actual = 0;
      arsum = 0;
      // loop through the neighbors nthhop links away
      for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
        nbrvtx = shl->vtxlist[nthhop][nthnbr];
        if (surf->vertices[nbrvtx].ripflag) continue;
        cnbr = crslut[0][nbrvtx];
        rnbr = crslut[1][nbrvtx];
        snbr = crslut[2][nbrvtx];
        if (mask)
          if (MRIgetVoxVal(mask, cnbr, rnbr, snbr, 0) < 0.5) continue;
        sumsqnbr = 0;
        sumsqx = 0;
        for (frame = 0; frame < src->nframes; frame++) {
          valvtx = MRIFseq_vox(src, c, r, s, frame);
          valnbr = MRIFseq_vox(src, cnbr, rnbr, snbr, frame);
          sumsqnbr += (valnbr * valnbr);
          sumsqx += (valvtx * valnbr);
        }
        if (sumsqnbr == 0) continue;
        nbrvar = sumsqnbr / src->nframes;
        arsum += (sumsqx / src->nframes) / sqrt(vtxvar * nbrvar);
        nnbrs_actual++;
      } /* end loop over hop neighborhood */

      if (nnbrs_actual != 0) MRIFseq_vox(arN, c, r, s, nthhop) = (arsum / nnbrs_actual);
    } /* end loop over hop */

    SurfHopListFree(&shl);
    
    ROMP_PFLB_end
  } /* end loop over vertex */
  ROMP_PF_end
  
  MRIScrsLUTFree(crslut);

  return (arN);
}

/*-----------------------------------------------------------------------
  MRISsmoothKernel() -
  kernel = ACF^2
  -----------------------------------------------------------------------*/
MRI *MRISsmoothKernel(MRIS *surf, MRI *src, MRI *mask, MRI *mrikern, MATRIX *globkern, SURFHOPLIST ***pshl, MRI *out)
{
  int vtx, **crslut, nvox;
  double *kern;
  int n, nhops;
  SURFHOPLIST **shl;
  Timer mytimer;
  int msecTime;

  if (mrikern && globkern) {
    printf("ERROR: MRISsmoothKernel(): cannot spec both mrikern and globkern\n");
    return (NULL);
  }

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox) {
    printf("ERROR: MRISsmoothKernel(): Surf/Src dimension mismatch.\n");
    return (NULL);
  }

  if (out == NULL) {
    out = MRIcloneBySpace(src, MRI_FLOAT, src->nframes);
    if (out == NULL) {
      printf("ERROR: MRISsmoothKernel(): could not alloc\n");
      return (NULL);
    }
  }
  else {
    if (out == src) {
      printf("ERROR: MRISsmoothKernel(): cannot run in-place\n");
      return (NULL);
    }
    // Zero the output
  }

  nhops = 0;
  if (mrikern) nhops = mrikern->nframes;
  if (globkern) nhops = globkern->rows;
  printf("nhops = %d\n", nhops);
  kern = (double *)calloc(nhops, sizeof(double));
  if (globkern)
    for (n = 0; n < nhops; n++) kern[n] = globkern->rptr[n + 1][1];

  // Build LUT to map from col,row,slice to vertex
  crslut = MRIScrsLUT(surf, src);

  if (*pshl == NULL) {
    mytimer.reset();
    printf("Allocating shl %d\n", surf->nvertices);
    shl = (SURFHOPLIST **)calloc(sizeof(SURFHOPLIST *), surf->nvertices);
    ROMP_PF_begin
#ifdef _OPENMP
    #pragma omp parallel for if_ROMP(experimental)
#endif
    for (vtx = 0; vtx < surf->nvertices; vtx++) {
      ROMP_PFLB_begin
      
      // Create structure to manage the multiple hops for this vertex
      shl[vtx] = SetSurfHopList(vtx, surf, nhops);
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    *pshl = shl;
    msecTime = mytimer.milliseconds();
    printf("Done allocating shl %d, %g sec\n", surf->nvertices, msecTime / 1000.0);
  }
  else
    shl = *pshl;

  printf("Starting loop over %d vertices\n", surf->nvertices);
  mytimer.reset();
  ROMP_PF_begin
#ifdef _OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    ROMP_PFLB_begin
    
    int nnbrs, frame, nbrvtx, nthnbr, c, r, s;
    int cnbr, rnbr, snbr, nnbrs_actual;
    double vtxval = 0, *vkern, ksum, kvsum;
    int nthhop;

    // if(vtx%10000 == 0) printf("%4.1f ",(100.0*vtx)/surf->nvertices);

    if (surf->vertices[vtx].ripflag) continue;
    c = crslut[0][vtx];
    r = crslut[1][vtx];
    s = crslut[2][vtx];
    if (mask)
      if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) continue;

    if (mrikern) {
      vkern = (double *)calloc(nhops, sizeof(double));
      for (nthhop = 0; nthhop < nhops; nthhop++) vkern[nthhop] = MRIgetVoxVal(mrikern, vtx, 0, 0, nthhop);
    }
    else
      vkern = kern;

    for (frame = 0; frame < src->nframes; frame++) {
      // loop through hops and neighbors
      if (frame == 0) ksum = 0;
      kvsum = 0;
      nnbrs_actual = 0;
      for (nthhop = 0; nthhop < nhops; nthhop++) {
        nnbrs = shl[vtx]->nperhop[nthhop];
        // loop through the neighbors nthhop links away
        for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
          nbrvtx = shl[vtx]->vtxlist[nthhop][nthnbr];
          if (surf->vertices[nbrvtx].ripflag) continue;
          cnbr = crslut[0][nbrvtx];
          rnbr = crslut[1][nbrvtx];
          snbr = crslut[2][nbrvtx];
          if (mask)
            if (MRIgetVoxVal(mask, cnbr, rnbr, snbr, 0) < 0.5) continue;
          vtxval = MRIFseq_vox(src, cnbr, rnbr, snbr, frame);
          kvsum += (vtxval * kern[nthhop]);
          // fabs() wont make diff if low pass filter
          if (frame == 0) ksum += fabs(kern[nthhop]);
          // if(vtx==1031) printf("%d %d %d %g %g %g\n",nnbrs_actual,vtx,nbrvtx,vtxval,kern[nthhop],kvsum);
          nnbrs_actual++;
        } /* end loop over hop neighborhood */
      }   /* end loop over hop */
      // normalize - might not be a good idea for high-pass filtering
      if (nnbrs_actual != 0) MRIFseq_vox(out, c, r, s, frame) = (kvsum / ksum);
      // if(vtx==1031)  printf("%5d %d %d %d %d %g %g
      // %g\n",vtx,c,r,s,nnbrs_actual,kvsum,ksum,MRIFseq_vox(out,c,r,s,frame));
    }  // end loop over frame
    // SurfHopListFree(&shl[vtx]);

    if (mrikern) free(vkern);
    ROMP_PFLB_end
  } /* end loop over vertex */
  ROMP_PF_end
  
  printf("\n");
  msecTime = mytimer.milliseconds();
  printf("Finished loop over %d vertices, %d hops, t=%g sec\n", surf->nvertices, nhops, msecTime / 1000.0);

  printf("vtxval %g\n", MRIFseq_vox(out, 1031, 0, 0, 0));

  free(kern);
  MRIScrsLUTFree(crslut);

  return (out);
}

/*-----------------------------------------------------------------------
  MRISar1() - computes spatial AR1 at each vertex by averaging the AR1s
  within the neighborhood of a vertex. Note: does not try to take into
  account different distances between neighbors. Note: theoretical AR1
  for one smoothing step is 4/7=0.57.
  -----------------------------------------------------------------------*/
MRI *MRISar1(MRIS *surf, MRI *src, MRI *mask, MRI *ar1)
{
  int nnbrs, frame, vtx, nbrvtx, nthnbr, **crslut, c, r, s, nvox;
  int cnbr, rnbr, snbr, nnbrs_actual;
  double valvtx, valnbr, ar1sum, sumsqvtx, vtxvar, sumsqnbr, sumsqx, nbrvar;

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox) {
    printf("ERROR: MRISar1: Surf/Src dimension mismatch.\n");
    return (NULL);
  }

  if (ar1 == NULL) {
    ar1 = MRIcloneBySpace(src, MRI_FLOAT, 1);
    if (ar1 == NULL) {
      printf("ERROR: could not alloc\n");
      return (NULL);
    }
  }

  // Build LUT to map from col,row,slice to vertex
  crslut = MRIScrsLUT(surf, src);

  for (vtx = 0; vtx < surf->nvertices; vtx++) {
    if (surf->vertices[vtx].ripflag) {
      continue;
    }

    c = crslut[0][vtx];
    r = crslut[1][vtx];
    s = crslut[2][vtx];
    if (mask)
      if (MRIgetVoxVal(mask, c, r, s, 0) < 0.5) {
        continue;
      }

    nnbrs = surf->vertices_topology[vtx].vnum;
    sumsqvtx = 0;
    for (frame = 0; frame < src->nframes; frame++) {
      valvtx = MRIFseq_vox(src, c, r, s, frame);
      sumsqvtx += (valvtx * valvtx);
    }
    if (sumsqvtx == 0) {
      continue;  // exclude voxels with 0 variance
    }
    vtxvar = sumsqvtx / src->nframes;

    nnbrs_actual = 0;
    ar1sum = 0;
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      nbrvtx = surf->vertices_topology[vtx].v[nthnbr];
      if (surf->vertices[nbrvtx].ripflag) {
        continue;
      }
      cnbr = crslut[0][nbrvtx];
      rnbr = crslut[1][nbrvtx];
      snbr = crslut[2][nbrvtx];
      if (mask)
        if (MRIgetVoxVal(mask, cnbr, rnbr, snbr, 0) < 0.5) {
          continue;
        }
      sumsqnbr = 0;
      sumsqx = 0;
      for (frame = 0; frame < src->nframes; frame++) {
        valvtx = MRIFseq_vox(src, c, r, s, frame);
        valnbr = MRIFseq_vox(src, cnbr, rnbr, snbr, frame);
        sumsqnbr += (valnbr * valnbr);
        sumsqx += (valvtx * valnbr);
      }
      if (sumsqnbr == 0) {
        continue;
      }
      nbrvar = sumsqnbr / src->nframes;
      ar1sum += (sumsqx / src->nframes) / sqrt(vtxvar * nbrvar);
      nnbrs_actual++;
    } /* end loop over neighbor */

    if (nnbrs_actual != 0)  // bug fix
    {
      MRIFseq_vox(ar1, c, r, s, 0) = (ar1sum / nnbrs_actual);
    }
  } /* end loop over vertex */

  MRIScrsLUTFree(crslut);
  return (ar1);
}
/*---------------------------------------------------------------
  MRIScrsLUT() - constructs a lookup table to quickly map from
  a vertex to the col, row, slice in an MRI struct, where the
  MRI struct has ncols*nrows*nslices = nvertices.
  See also MRIScrsLUTFee().
  ---------------------------------------------------------------*/
int **MRIScrsLUT(MRIS *surf, MRI *src)
{
  int **crslut, c, r, s, nvox, vtx;
  crslut = (int **)calloc(3, sizeof(int *));

  nvox = src->width * src->height * src->depth;
  if (surf->nvertices != nvox) {
    printf("ERROR: MRIScrsLUT: surf/src dimension mismatch.\n");
    return (NULL);
  }

  crslut[0] = (int *)calloc(nvox, sizeof(int));
  crslut[1] = (int *)calloc(nvox, sizeof(int));
  crslut[2] = (int *)calloc(nvox, sizeof(int));
  vtx = 0;
  for (s = 0; s < src->depth; s++) {
    for (r = 0; r < src->height; r++) {
      for (c = 0; c < src->width; c++) {
        // Do columns fastest
        crslut[0][vtx] = c;
        crslut[1][vtx] = r;
        crslut[2][vtx] = s;
        vtx++;
      }
    }
  }
  return (crslut);
}
/*----------------------------------------------------------
  MRIScrsLUTFree() - frees memory allocated by MRIScrsLUT.
  ----------------------------------------------------------*/
int MRIScrsLUTFree(int **crslut)
{
  free(crslut[0]);
  free(crslut[1]);
  free(crslut[2]);
  free(crslut);
  return (0);
}

/*-------------------------------------------------------------------
  MRISextendedNeighbors() - read everything! Finds the set of
  "extended" neighbors of a given target vertex and a distance
  threshold. Call with CurVtxNo=TargVtxNo. There are other
  recursive calls where CurVtxNo!=TargVtxNothat.

  The distance metric is the dot product between the vectors
  pointing from the origin to the target and current vertices. For
  any given target vertex, the dot product must be greater than
  (NOT less than) the DotProdThresh to be included.

  XNbrVtxNo is the (already allocated) list of vertex numbers of
  vertices that are within threshold.

  XNbrDotProd is the (already allocated) list of dot products of
  vertices that are within threshold.

  nXNbrs is a pointer to the total number of vertices that are
  within threshold.

  nXNbrsMax is the maximum number allowed (ie, the allocation
  lengths of XNbrVtxNo and XNbrDotProd.

  NOTE: IMPORTANT!
  1. This really only works on the spherical surface.
  2. Assuming a spherical surface, the distance along the
  sphere between the two vertices can be computed as
  costheta = dotprod/(Radius^2);
  theta = acos(costheta);
  d = Radius * theta;
  This also allows you to work backwards to get a
  dot product threshold from a given distance threshold.
  3. Modifies vertex->val2bak. It is important that this
  be set to -1 for all vertices prior to the first
  call or before calling with the same target vertex
  again.
  4. It is expected that this will be run for each vertex
  on a surface, so care has been taken to keep the
  computations light (eg, not allocating XNbrnVtxNo
  and XNbrnDotProd, using a dot product threshold
  instead of a distance threshold, not resetting val2bak).
  5. The extended neighborhood of a vertex includes itself.
  -------------------------------------------------------------------*/
int MRISextendedNeighbors(MRIS *SphSurf,
                          int TargVtxNo,
                          int CurVtxNo,
                          double DotProdThresh,
                          int *XNbrVtxNo,
                          double *XNbrDotProd,
                          int *nXNbrs,
                          int nXNbrsMax,
                          int DistType)
{
  static int ncalls = 0;
  VERTEX *vtarg, *vcur;
  int nNNbrs, n, NbrVtxNo, err;
  double DotProd, dx, dy, dz;

  // Get the current vertex
  vcur = &SphSurf->vertices[CurVtxNo];

  // Return if this vertex has been hit
  if ((int)vcur->val2bak == TargVtxNo) {
    return (0);
  }

  // Return if this vertex is ripped
  if (vcur->ripflag) {
    return (0);
  }

  // Keep track of the number of recursive calls
  if (CurVtxNo == TargVtxNo) {
    *nXNbrs = 0;
    ncalls = 0;
  }
  ncalls++;

  // Get the target vertex
  vtarg = &SphSurf->vertices[TargVtxNo];

  if (DistType == 1) {
    // DotProduct - dist along sphere
    // Compute the dot product between the two
    DotProd = (vtarg->x * vcur->x) + (vtarg->y * vcur->y) + (vtarg->z * vcur->z);
    DotProd = fabs(DotProd);
    // printf("c %d %d %d %g %d\n",ncalls,TargVtxNo,CurVtxNo,DotProd,*nXNbrs);
    if (DotProd <= DotProdThresh) {
      return (0);
    }
  }
  else {
    // Cartesian - dist squared along sphere (so thresh should be squared)
    dx = vtarg->x - vcur->x;
    dy = vtarg->y - vcur->y;
    dz = vtarg->z - vcur->z;
    DotProd = dx * dx + dy * dy + dz * dz;
    DotProd = fabs(DotProd);
    // printf("c %d %d %d %g %d\n",ncalls,TargVtxNo,CurVtxNo,DotProd,*nXNbrs);
    if (DotProd >= DotProdThresh) {
      return (0);
    }
  }

  // Check whether another neigbor can be added
  if (*nXNbrs >= nXNbrsMax - 1) {
    return (1);
  }

  // OK, add this vertex as an extended neighbor
  XNbrVtxNo[*nXNbrs] = CurVtxNo;
  XNbrDotProd[*nXNbrs] = DotProd;
  (*nXNbrs)++;
  vcur->val2bak = TargVtxNo;  // record a hit

  // Now, loop over the current nearest neighbors
  nNNbrs = SphSurf->vertices_topology[CurVtxNo].vnum;
  for (n = 0; n < nNNbrs; n++) {
    NbrVtxNo = SphSurf->vertices_topology[CurVtxNo].v[n];
    err = MRISextendedNeighbors(
        SphSurf, TargVtxNo, NbrVtxNo, DotProdThresh, XNbrVtxNo, XNbrDotProd, nXNbrs, nXNbrsMax, DistType);
    if (err) {
      return (err);
    }
  }

  return (0);
}


int mrisMarkIntersections(MRIS *mris, int FillHoles)
{
  MRIS_HASH_TABLE *mht;
  FACE *f;
  int fno, n, num = 0;

  mht = MHTcreateFaceTable(mris);

  MRISclearMarks(mris);
  for (num = fno = 0; fno < mris->nfaces; fno++) {
    if (MHTdoesFaceIntersect(mht, mris, fno)) {
      num++;
      f = &mris->faces[fno];
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mris->vertices[f->v[n]].marked = 1;
      }
    }
  }
  if(FillHoles){
    printf("mrisMarkIntersections(): Filling any holes\n");
    for(n=0; n < mris->nvertices; n++){
      mris->vertices[n].val = !mris->vertices[n].marked; // invert
    }
    SURFCLUSTERSUM *SurfClustList;
    int nClusters;
    SurfClustList = sclustMapSurfClusters(mris,0.5,-1,1,0,&nClusters,NULL,NULL);
    printf("  Found %d holes\n",nClusters);
    if(nClusters > 0){
      for(n=0; n < mris->nvertices; n++){
	if(mris->vertices[n].undefval > 1) mris->vertices[n].marked = 1;
      }
    }
    free(SurfClustList);
  }

  MHTfree(&mht);
  return (num);
}


#define WM_VAL 1
#define GM_VAL 2
#define CSF_VAL 3

MRI *MRIcomputeLaminarVolumeFractions(MRIS *mris, double resolution, MRI *mri_src, MRI *mri_fractions)
{
  int width, height, depth, nvox;
  MRI *mri_layers, *mri_interior_pial, *mri_tmp;
  MATRIX *m_vox2vox;
  static MRI *mri_interior_wm = NULL;

  MRISsaveVertexPositions(mris, TMP2_VERTICES);

  nvox = nint(MAX(MAX((mri_src->xsize / resolution), (mri_src->ysize / resolution)), (mri_src->zsize / resolution)));
  if (Gdiag & DIAG_VERBOSE_ON)
    printf("packing each voxel with %d voxels to compute laminar fractions\n", (int)pow(nvox, 3));
  width = (int)ceil(mri_src->width * nvox);
  height = (int)ceil(mri_src->height * nvox);
  depth = (int)ceil(mri_src->depth * nvox);
  mri_layers = MRIupsampleN(mri_src, NULL, nvox);

  mri_interior_pial = MRIclone(mri_layers, NULL);
  if (mri_interior_wm == NULL) {
    MRISrestoreVertexPositions(mris, WHITE_VERTICES);
    MRIScomputeMetricProperties(mris);
    MRISfillInterior(mris, resolution, mri_layers);
    mri_interior_wm = MRIcopy(mri_layers, NULL);
  }
  else
    MRIcopy(mri_interior_wm, mri_layers);

  MRISrestoreVertexPositions(mris, PIAL_VERTICES);
  MRIScomputeMetricProperties(mris);
  MRISfillInterior(mris, resolution, mri_interior_pial);

  mri_tmp = MRInot(mri_layers, NULL);
  MRIand(mri_interior_pial, mri_tmp, mri_tmp, 1);
  MRIreplaceValuesOnly(mri_tmp, mri_tmp, 1, GM_VAL);
  MRIcopyLabel(mri_tmp, mri_layers, GM_VAL);
  MRIfree(&mri_tmp);
  MRIreplaceValuesOnly(mri_layers, mri_layers, 0, CSF_VAL);
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_layers, mri_src);
  if (mri_fractions == NULL) {
    mri_fractions = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT, 3);
    MRIcopyHeader(mri_src, mri_fractions);
  }
  MRIcomputeVolumeFractions(mri_src, m_vox2vox, mri_layers, mri_fractions);
  MatrixFree(&m_vox2vox);
  MRIfree(&mri_layers);
  MRIfree(&mri_interior_pial);
  MRISrestoreVertexPositions(mris, TMP2_VERTICES);
  MRIScomputeMetricProperties(mris);
  return (mri_fractions);
}


int MRISscaleVertexCoordinates(MRIS *mris, double scale)
{
  return MRISscale(mris, scale);
}



MRIS *MRISscaleBrain(MRIS *mris_src, MRIS *mris_dst, float scale)
{
  cheapAssert(mris_src == mris_dst);
  
  MRISanisotropicScale(mris_src, scale,scale,scale);

  return mris_dst;
}


void mrisFindMiddleOfGray(MRIS *mris) {
  int     vno ;
  VERTEX  *v ;
  float   nx, ny, nz, thickness ;

  MRISaverageCurvatures(mris, 3) ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris);
  
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    nx = v->nx ;
    ny = v->ny ;
    nz = v->nz ;
    thickness = 0.5 * v->curv ;
    MRISsetXYZ(mris,vno,
      v->origx + thickness * nx,
      v->origy + thickness * ny,
      v->origz + thickness * nz);
  }
}


// Cloning should be the union of a surface and an empty surface
// but that is NYI
//
MRIS* MRISunion(MRIS const * mris, MRIS const * mris2) {
    int vno,vno2,vno3;
    MRIS * const mris3 = MRISalloc(mris->nvertices+mris2->nvertices,
                                   mris->nfaces+mris2->nfaces);
    copyVolGeom(&mris->vg,&mris3->vg);
    for (vno=0,vno3=0; vno < mris->nvertices; vno++, vno3++)
    {
      MRISsetXYZ(mris3,vno3,
        mris->vertices[vno].x,
        mris->vertices[vno].y,
        mris->vertices[vno].z);
    }
    for (vno2=0; vno2 < mris2->nvertices; vno2++, vno3++)
    {
      MRISsetXYZ(mris3,vno3,
        mris2->vertices[vno2].x,
        mris2->vertices[vno2].y,
        mris2->vertices[vno2].z);
    }
    int fno,fno2,fno3;
    for (fno=0,fno3=0; fno < mris->nfaces; fno++, fno3++)
    {
      mris3->faces[fno3].v[0] = mris->faces[fno].v[0];
      mris3->faces[fno3].v[1] = mris->faces[fno].v[1];
      mris3->faces[fno3].v[2] = mris->faces[fno].v[2];
    }
    int offset = mris->nvertices;
    for (fno2=0; fno2 < mris2->nfaces; fno2++, fno3++)
    {
      mris3->faces[fno3].v[0] = mris2->faces[fno2].v[0] + offset;
      mris3->faces[fno3].v[1] = mris2->faces[fno2].v[1] + offset;
      mris3->faces[fno3].v[2] = mris2->faces[fno2].v[2] + offset;
    }
    mrisCheckVertexFaceTopology(mris3);
    return mris3;
}    

MRIS* MRISclone(MRIS const * mris_src)
{
  // Cloning could be a copy the input data and recompute the derived,
  // but it is quicker to copy the derived data also which means the derived data must be written, 
  // which means placing this where the writing functions are visible,
  // which means not in mrisurf_base.c
  //
  mrisCheckVertexFaceTopology(mris_src);

  MRI_SURFACE *mris_dst;
  int vno, fno, n;
  FACE *fsrc, *fdst;

  mris_dst = MRISalloc(mris_src->nvertices, mris_src->nfaces);

  mris_dst->type = mris_src->type;
  mris_dst->status = mris_src->status;
  mris_dst->origxyz_status = mris_src->origxyz_status;
  
  mris_dst->nsize                   = mris_src->nsize;
  mris_dst->max_nsize               = mris_src->max_nsize;
  mris_dst->vtotalsMightBeTooBig    = mris_src->vtotalsMightBeTooBig;
  mris_dst->nsizeMaxClock           = mris_src->nsizeMaxClock;
  
  mris_dst->hemisphere = mris_src->hemisphere;
  mris_dst->xctr = mris_src->xctr;
  mris_dst->yctr = mris_src->yctr;
  mris_dst->zctr = mris_src->zctr;
  mris_dst->xlo = mris_src->xlo;
  mris_dst->ylo = mris_src->ylo;
  mris_dst->zlo = mris_src->zlo;
  mris_dst->xhi = mris_src->xhi;
  mris_dst->yhi = mris_src->yhi;
  mris_dst->zhi = mris_src->zhi;
  mris_dst->min_curv = mris_src->min_curv;
  mris_dst->max_curv = mris_src->max_curv;
  mris_dst->total_area = mris_src->total_area;
  mris_dst->orig_area = mris_src->orig_area;

  mris_dst->radius = mris_src->radius;  // to be checked

  mris_dst->avg_vertex_area = mris_src->avg_vertex_area;
  mrisSetAvgInterVertexDist(mris_dst, mris_src->avg_vertex_dist);
  mris_dst->std_vertex_dist = mris_src->std_vertex_dist;

  // just copy the pointer ///////////////////////////////////
  mris_dst->lta = mris_src->lta;
  mris_dst->SRASToTalSRAS_ = mris_src->SRASToTalSRAS_;
  mris_dst->TalSRASToSRAS_ = mris_src->TalSRASToSRAS_;
  mris_dst->free_transform = 0;  // mark not to try to free them
  //                             // BUG - THE mris_src may still free them!  reference counting needed.
  /////////////////////////////////////////////////////////////
  
  if (mris_src->v_frontal_pole)
    mris_dst->v_frontal_pole = &mris_dst->vertices[mris_src->v_frontal_pole - mris_src->vertices];
  if (mris_src->v_occipital_pole)
    mris_dst->v_occipital_pole = &mris_dst->vertices[mris_src->v_occipital_pole - mris_src->vertices];
  if (mris_src->v_temporal_pole)
    mris_dst->v_temporal_pole = &mris_dst->vertices[mris_src->v_temporal_pole - mris_src->vertices];
    
  for (vno = 0; vno < mris_src->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    VERTEX_TOPOLOGY const * const vsrct = &mris_src->vertices_topology[vno];
    VERTEX          const * const vsrc  = &mris_src->vertices         [vno];
    VERTEX_TOPOLOGY       * const vdstt = &mris_dst->vertices_topology[vno];
    VERTEX                * const vdst  = &mris_dst->vertices         [vno];
    
    vdst->ripflag = vsrc->ripflag;
        //
        // Even in ripped vertices, the topology below must be maintained
        
    vdst->x = vsrc->x;
    vdst->y = vsrc->y;
    vdst->z = vsrc->z;
    vdst->nx = vsrc->nx;
    vdst->ny = vsrc->ny;
    vdst->nz = vsrc->nz;
    vdst->cx = vsrc->cx;
    vdst->cy = vsrc->cy;
    vdst->cz = vsrc->cz;
    vdst->curv = vsrc->curv;
    vdstt->num = vsrct->num;


    if (vdstt->num) {
      vdstt->f = (int *)calloc(vdstt->num, sizeof(int));
      vdstt->n = (uchar *)calloc(vdstt->num, sizeof(uchar));
      if (!vdstt->f) ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d faces", vdstt->num);
      if (!vdstt->n) ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d num", vdstt->num);
      for (n = 0; n < vdstt->num; n++) {
        vdstt->n[n] = vsrct->n[n];
        vdstt->f[n] = vsrct->f[n];
      }
    }
    
    modVnum(mris_dst,vno,vsrct->vnum,true);
    vdstt->v2num    = vsrct->v2num;
    vdstt->v3num    = vsrct->v3num;
    vdstt->nsizeMax = vsrct->nsizeMax;
    vdstt->nsizeMaxClock = vsrct->nsizeMaxClock;

    MRIS_setNsizeCur(mris_dst,vno,vsrct->nsizeCur);
    
    {
      int vSize = mrisVertexVSize(mris_src, vno);

      if (vSize) {
        vdstt->v = (int *)calloc(vSize, sizeof(int));
        if (!vdstt->v) ErrorExit(ERROR_NO_MEMORY, "MRISclone: could not allocate %d nbrs", vSize);
        for (n = 0; n < vSize; n++) {
          vdstt->v[n] = vsrct->v[n];
        }
      }

      if (vsrc->dist) {
        MRISmakeDist(mris_dst, vno);
        for (n = 0; n < vSize; n++) {
          vdst->dist[n] = vsrc->dist[n];
        }
      }
      
      if (vsrc->dist_orig) {
        MRISmakeDistOrig(mris_dst, vno);
        for (n = 0; n < vSize; n++) {
          vdst->dist_orig[n] = vsrc->dist_orig[n];
        }
      }
    }
    
    vdst->border   = vsrc->border;
    vdst->area     = vsrc->area;
    vdst->origarea = vsrc->origarea;
  }

  for (fno = 0; fno < mris_src->nfaces; fno++) {
    fsrc = &mris_src->faces[fno];
    fdst = &mris_dst->faces[fno];
    memmove(fdst, fsrc, sizeof(FACE));
    
    // The pointer fields need special handling
    if (fsrc->norm) fdst->norm = DMatrixCopy(fsrc->norm, NULL);
    int i;
    for (i = 0; i < 3; i++) {
      if (fsrc->gradNorm[i]) fdst->gradNorm[i] = DMatrixCopy(fsrc->gradNorm[i], NULL);
    }
  }

  mrisCheckVertexFaceTopology(mris_dst);

  // copy geometry info
  copyVolGeom(&mris_src->vg, &mris_dst->vg);

  return (mris_dst);
}



int MRISprintVertexStats(MRIS *mris, int vno, FILE *fp, int which_vertices)
{
  double mn, d;
  int n, num;
  float x0, y0, z0, x, y, z, dx, dy, dz;

  if (vno < 0) {
    return (NO_ERROR);
  }
  x0 = y0 = z0 = x = y = z = 0;
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
  VERTEX          * const v  = &mris->vertices         [vno];
  MRISvertexCoord2XYZ_float(v, which_vertices, &x0, &y0, &z0);
  printf("vertex %d spacing for %s surface\n",
         vno,
         which_vertices == ORIGINAL_VERTICES
             ? "orig"
             : which_vertices == CURRENT_VERTICES ? "current" : which_vertices == WHITE_VERTICES ? "white" : "unknown");
  for (mn = 0.0, num = n = 0; n < vt->vnum; n++) {
    VERTEX * const vn = &mris->vertices[vt->v[n]];
    if (vn->ripflag) {
      printf("nbr %d = %d is ripped\n", n, vt->v[n]);
      continue;
    }
    num++;
    MRISvertexCoord2XYZ_float(vn, which_vertices, &x, &y, &z);
    dx = x - x0;
    dy = y - y0;
    dz = z - z0;
    d = sqrt(dx * dx + dy * dy + dz * dz);
    printf("\tvn %d: %d = %2.3f mm distant\n", n, vt->v[n], d);
    mn += d;
  }
  printf("\tmean = %2.3f\n", mn / (double)num);

  return (NO_ERROR);
}


/*!
  \fn int MRISprintVertexInfo(FILE *fp, MRIS *surf, int vertexno)
  \brief Prints out info about a vertex. Neighbor info if formated:
  (1) Ordinal Neighbor number
  (2) Vertex index number of vertex neighbor
  (3) Distance to the vertex neighbor (mm)
  (4) "Area" of the neighbor vertex (mm2)
  (5) Face index number of nth face neighbor
  (6) Area of face 
  See also int MRISprintVertexStats()
 */
int MRISprintVertexInfo(FILE *fp, MRIS *surf, int vertexno)
{
  int n, vnno, faceno;
  VERTEX *v = &(surf->vertices[vertexno]);
  VERTEX *vn;
  VERTEX_TOPOLOGY *vt = &(surf->vertices_topology[vertexno]);
  double dx, dy, dz, dist;
  FACE *face;

  fprintf(fp,"vertexno %d\n",vertexno);
  fprintf(fp,"xyz  %7.4f %7.4f %7.4f\n",v->x,v->y,v->z);
  fprintf(fp,"nxyz %7.4f %7.4f %7.4f\n",v->nx,v->ny,v->nz);
  fprintf(fp,"vertex_area %lf\n",v->area);
  fprintf(fp,"neighbors %d\n",vt->num);
  for(n=0; n < vt->num; n++){
    vnno = vt->v[n];
    vn = &(surf->vertices[vnno]);
    dx = v->x - vn->x;
    dy = v->y - vn->y;
    dz = v->z - vn->z;
    dist = sqrt(dx*dx+dy*dy+dz*dz);
    faceno = vt->f[n];
    face = &(surf->faces[faceno]);
    // nbr nbrno vno dist varea faceno farea
    fprintf(fp,"nbr %d  %5d %6.4lf %6.4lf %5d %6.4f\n",n,vnno,dist,vn->area,faceno,face->area);
  }

  if(surf->edges){
    int edgeno;
    for(edgeno = 0; edgeno < surf->nedges; edgeno++){
      MRI_EDGE *e = &(surf->edges[edgeno]);
      int evno;
      for(evno = 0; evno < 2; evno ++){
	if(e->vtxno[evno] == vertexno){
	  printf("Edge %6d %d  %4d %4d %7.3lf %7.3lf\n",edgeno,evno,e->vtxno[0],e->vtxno[1],e->len,e->angle);
	}
      }
    }
  }

  return(0);
}


int MRISprintSurfQualityStats(FILE *fp, MRIS *surf)
{
  double *estats, *hstats, *astats, *cstats;
  astats = MRIStriangleAreaStats(surf, NULL, NULL); // trangle area
  cstats = MRIScornerStats(surf, 1, NULL, NULL); // corner angle (deg)
  estats = MRISedgeStats(surf, 0, NULL, NULL); // edge length
  hstats = MRISedgeStats(surf, 2, NULL, NULL); // hinge angle (deg)

  // mean, stddev, min, max
  fprintf(fp,"%7.5f %7.5f %7.5f %7.5f  ",astats[1],astats[2],astats[3],astats[4]);
  fprintf(fp,"%7.5f %7.5f %7.5f %7.5f  ",cstats[1],cstats[2],cstats[3],cstats[4]);
  fprintf(fp,"%7.5f %7.5f %7.5f %7.5f  ",estats[1],estats[2],estats[3],estats[4]);
  fprintf(fp,"%7.5f %7.5f %7.5f %7.5f  ",hstats[1],hstats[2],hstats[3],hstats[4]);
  fprintf(fp,"\n");
  fflush(fp);

  free(astats);
  free(cstats);
  free(estats);
  free(hstats);

  return(0);
}

int MRISprettyPrintSurfQualityStats(FILE *fp, MRIS *surf)
{
  double *estats, *hstats, *astats, *cstats;
  astats = MRIStriangleAreaStats(surf, NULL, NULL); // trangle area
  cstats = MRIScornerStats(surf, 1, NULL, NULL); // corner angle (deg)
  estats = MRISedgeStats(surf, 0, NULL, NULL); // edge length
  hstats = MRISedgeStats(surf, 2, NULL, NULL); // hinge angle (deg)

  // mean, stddev, min, max
  fprintf(fp,"Area   %7d %8.5f %8.5f %8.6f %8.4f\n",(int)astats[0],astats[1],astats[2],astats[3],astats[4]);
  fprintf(fp,"Corner %7d %8.5f %8.5f %8.6f %8.4f\n",(int)cstats[0],cstats[1],cstats[2],cstats[3],cstats[4]);
  fprintf(fp,"Edge   %7d %8.5f %8.5f %8.6f %8.4f\n",(int)estats[0],estats[1],estats[2],estats[3],estats[4]);
  fprintf(fp,"Hinge  %7d %8.5f %8.5f %8.6f %8.4f\n",(int)hstats[0],hstats[1],hstats[2],hstats[3],hstats[4]);
  if(0){
    // Number of vertices that belong to a face that intersects another face
    mrisMarkIntersections(surf,0);
    int n, nintersections=0;
    for(n=0; n < surf->nvertices; n++){
      if(surf->vertices[n].marked) nintersections++;
    }
    fprintf(fp,"Intersections %d \n",nintersections);
  }
  fflush(fp);

  free(astats);
  free(cstats);
  free(estats);
  free(hstats);

  return(0);
}

/*!
  \fn int MRISfindExpansionRegions(MRI_SURFACE *mris)
  \brief Sets v->curv=0 unless the given vertex has more than 25% of
  the neighbors whose v->d value is greater than mean+2*std (mean and
  std are the global mean and stddev distances). The dist, eg, is the
  distance along the normal to point of the max gradient. The v->val
  is the max gradient. Moved from mris_make_surfaces. This is essentially
  masking out the curv field of vertices with long distances.
  Hidden parameters: 0.25 and mean+2*std
 */
int MRISfindExpansionRegions(MRI_SURFACE *mris)
{
  int    vno, num, n, num_long, total ;
  float  d, dsq, mean, std, dist ;

  // Compute the mean and stddev of the distance to max gradient
  d = dsq = 0.0f ;
  for (total = num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->ripflag || v->val <= 0)
    {
      continue ;
    }
    num++ ;
    dist = fabs(v->d) ;
    d += dist ;
    dsq += (dist*dist) ;
  }
  mean = d / num ;
  std = sqrt(dsq/num - mean*mean) ;
  printf("mean absolute distance = %2.2f +- %2.2f\n", mean, std); fflush(stdout);

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    v->curv = 0 ;

    if (v->ripflag || v->val <= 0) continue ;

    if (fabs(v->d) < mean+2*std) continue ;

    // Only gets here if distance is not too big

    // Count number of neighbors with big distances
    for (num_long = num = 1, n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->val <= 0 || v->ripflag) continue ;
      if (fabs(vn->d) >= mean+2*std)
        num_long++ ;
      num++ ;
    }

    // If the number of big dist neighbors is greater than 25%
    // Set v->curv = fabs(v->d) and increment total
    if ((float)num_long / (float)num > 0.25)
    {
      v->curv = fabs(v->d) ;
      total++ ; // not used for anything except diagnostic
    }

  } // end loop over vertices

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "%d vertices more than 2 sigmas from mean.\n", total) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRISwriteCurvature(mris, "long") ;

  return(NO_ERROR) ;
}

// consolidated from two identical copies in mris_flatten and mris_sphere
int MRISscaleUp(MRI_SURFACE *mris)
{
  int     vno, n, max_v, max_n ;
  float   ratio, max_ratio ;

  max_ratio = 0.0f ;
  max_v = max_n = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < vt->vnum ; n++)
    {
      if (FZERO(v->dist[n]))   /* would require infinite scaling */
        continue ;
      ratio = v->dist_orig[n] / v->dist[n] ;
      if (ratio > max_ratio)
      {
        max_v = vno ;
        max_n = n ;
        max_ratio = ratio ;
      }
    }
  }

  fprintf(stderr, "max @ (%d, %d), scaling brain by %2.3f\n",
          max_v, max_n, max_ratio) ;
#if 0
  MRISscaleBrain(mris, mris, max_ratio) ;
#else
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < vt->vnum ; n++)
      v->dist_orig[n] /= max_ratio ;
  }
#endif
  return(NO_ERROR) ;
}
