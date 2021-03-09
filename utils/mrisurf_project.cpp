#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
#define COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND
/*
 *
 */
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
#include "mrisurf_project.h"
#include "mrisurf_base.h"

/* project onto the sphere of radius DEFAULT_RADIUS */
void mrisSphericalProjectXYZ(float xs, float ys, float zs, float *xd, float *yd, float *zd)
{
  double dist, lambda;

  dist = sqrt(SQR(xs) + SQR(ys) + SQR(zs));
  lambda = DEFAULT_RADIUS / dist;

  /* making sure things are stable : double projection */
  *xd = xs * lambda;
  *yd = ys * lambda;
  *zd = zs * lambda;

  xs = *xd;
  ys = *yd;
  zs = *zd;
  dist = sqrt(SQR(xs) + SQR(ys) + SQR(zs));
  lambda = DEFAULT_RADIUS / dist;

  *xd = xs * lambda;
  *yd = ys * lambda;
  *zd = zs * lambda;
}


void mrisSphericalProjection(MRIS *mris)
{
  int n;
  VERTEX *v;
  //    fprintf(stderr,"spherical projection\n");
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->ripflag) {
      continue;
    }

    /*
      if(n == 88 )
      fprintf(stderr,"bf sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"bf sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"bf sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
    */

    mrisSphericalProjectXYZ(v->x, v->y, v->z, &v->x, &v->y, &v->z);
    v->cx = v->x;
    v->cy = v->y;
    v->cz = v->z;

    /*
      if(n == 88 )
      fprintf(stderr,"af sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"af sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"af sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);

      mrisSphericalProjectXYZ(v->x,v->y,v->z,&v->x,&v->y,&v->z);

      if(n == 88 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"af 2 sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
    */
  }
}


static MRIS* MRISprojectOntoSphereWkr(MRIS* mris, double r)
{
  if (FZERO(r)) {
    r = DEFAULT_RADIUS;
  }

  if ((mris->status != MRIS_SPHERE) && (mris->status != MRIS_PARAMETERIZED_SPHERE)) {
    MRIScenter(mris, mris);
  }

  mris->radius = r;

  MRISfreeDistsButNotOrig(mris);

  bool const showTotalDist = ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON);
  
  double total_dist = 0.0; 
  for (int vno = 0; vno < mris->nvertices; vno++) {
  
    auto v = &mris->vertices[vno];
    if (v->ripflag) /* shouldn't happen */
      continue;

    double const x = v->x, x2 = x*x;
    double const y = v->y, y2 = y*y;
    double const z = v->z, z2 = z*z;

    double const dist = sqrt(x2 + y2 + z2);
    double const d = FZERO(dist) ? 0 : (1 - r / dist);
    
    double const dx = d * x;
    double const dy = d * y;
    double const dz = d * z;
    
    v->x = x - dx;
    v->y = y - dy;
    v->z = z - dz;

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z)) {
      DiagBreak();
    }

    if (showTotalDist) total_dist += sqrt((double)(dx*dx + dy*dy + dz*dz));
  }

  if (showTotalDist) {
    fprintf(stdout, "sphere_project: total dist = %f\n", total_dist);
  }

  mris->status = (mris->status == MRIS_PARAMETERIZED_SPHERE) ? MRIS_PARAMETERIZED_SPHERE : MRIS_SPHERE;

  MRISupdateEllipsoidSurface(mris);

  return mris;
}


MRIS* MRISprojectOntoSphere(MRIS* mris, double r)
{
  return MRISprojectOntoSphereWkr(mris, r);
}

MRIS_MP* MRISprojectOntoSphere(MRIS_MP* mris, double r)
{
  cheapAssert(false);
  return mris;
}


MRIS* MRISprojectOntoSphere(MRIS* mris_src, MRIS* mris_dst, double r)
{
  cheapAssert(mris_src == mris_dst);
  auto result = MRISprojectOntoSphere(mris_src, r);
  cheapAssert(result == mris_dst);
  return result;
}



void mrisAssignFaces(MRI_SURFACE *mris, MHT *mht, int which_vertices)
{
  int vno;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin 
    
    int fno;
    VERTEX *v;
    double fdist;
    FACE *face;

    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    project_point_onto_sphere(v->x, v->y, v->z, mris->radius, &v->x, &v->y, &v->z);
    MHTfindClosestFaceGeneric(mht, mris, v->x, v->y, v->z, 8, 8, 1, &face, &fno, &fdist);
    if (fno < 0) MHTfindClosestFaceGeneric(mht, mris, v->x, v->y, v->z, 1000, -1, -1, &face, &fno, &fdist);

    v->fno = fno;
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Perform a projection onto a cylinder moving each
  point on the cortical surface to the closest cylindrical
  coordinate.
  ------------------------------------------------------*/
int MRISprojectOntoCylinder(MRI_SURFACE *mris, float radius)
{
  VERTEX *v;
  int k;
  float x, y, z, x2, z2, dx, dz;
  float d;

  MRIScenter(mris, mris);

  for (k = 0; k < mris->nvertices; k++) {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;

    x2 = x * x;
    z2 = z * z;

    d = (-1.0 + (float)radius / sqrt(x2 + z2));
    if (!std::isfinite(d)) {
      ErrorPrintf(ERROR_BADPARM, "point (%2.2f,%2.2f,%2.2f) cannot be projected on cylinder", x, y, z);
    }
    dx = d * x;
    dz = d * z;
    v->x = x + dx;
    v->z = z + dz;

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z)) {
      DiagBreak();
    }
  }
  MRISupdateSurface(mris);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Perform a projection onto an ellipsoid moving each
  point on the cortical surface to the closest ellipsoidal
  coordinate.
  ------------------------------------------------------*/
extern double sqrt(double);

MRI_SURFACE *MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, float a, float b, float c)
{
  VERTEX *v;
  int k;
  float x, y, z, x2, y2, z2, dx, dy, dz, a2, b2, c2, a4, b4, c4, a6, b6, c6;
  float f, g, h, d, dist, avgdist = 0.0f;

  if (FZERO(a)) {
    a = DEFAULT_A;
    b = DEFAULT_B;
    c = DEFAULT_C;
  }

  if (!mris_dst) {
    mris_dst = MRISclone(mris_src);
  }

  MRIScenter(mris_dst, mris_dst);

  mris_dst->a = a;
  mris_dst->b = b;
  mris_dst->c = c;

  /*  printf("ellipsoid_project(%f,%f,%f)\n",a,b,c);*/
  a2 = a * a;
  b2 = b * b;
  c2 = c * c;
  a4 = a2 * a2;
  b4 = b2 * b2;
  c4 = c2 * c2;
  a6 = a2 * a4;
  b6 = b2 * b4;
  c6 = c2 * c4;

#if 0
  /* rescale brain so that it is contained within the ellipsoid */
  xscale = mris_dst->xhi / a ;
  yscale = mris_dst->yhi / b ;
  zscale = mris_dst->zhi / c ;
  if ((xscale > yscale) && (xscale > zscale))
  {
    scale = 1.0f / xscale ;
  }
  else if (yscale > zscale)
  {
    scale = 1.0f / yscale ;
  }
  else
  {
    scale = 1.0f / zscale ;
  }

  MRISscaleBrain(mris_dst, mris_dst, scale) ;
#endif

  for (k = 0; k < mris_dst->nvertices; k++) {
    v = &mris_dst->vertices[k];
    /*
      printf("%6d: before: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
    */
    x = v->x;
    y = v->y;
    z = v->z;
#if 0
    if ((fabs(x) > a) || (fabs(y) > b) || (fabs(z) > c))
    {
      return(MRISradialProjectOntoEllipsoid(mris_src, mris_dst, a, b, c)) ;
    }
#endif

    x2 = x * x;
    y2 = y * y;
    z2 = z * z;
    f = x2 / a6 + y2 / b6 + z2 / c6;
    g = 2 * (x2 / a4 + y2 / b4 + z2 / c4);
    h = x2 / a2 + y2 / b2 + z2 / c2 - 1;
    d = (-g + (float)sqrt((double)(g * g - 4 * f * h))) / (2 * f);
    if (!std::isfinite(d)) {
      ErrorPrintf(ERROR_BADPARM,
                  "point (%2.2f,%2.2f,%2.2f) cannot be projected on ell "
                  "(%2.0f,%2.0f,%2.0f...\n",
                  x,
                  y,
                  z,
                  a,
                  b,
                  c);

      return (MRISradialProjectOntoEllipsoid(mris_src, mris_dst, a, b, c));
    }
    dx = d * x / a2;
    dy = d * y / b2;
    dz = d * z / c2;
    v->x = x + dx;
    v->y = y + dy;
    v->z = z + dz;

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z)) {
      DiagBreak();
    }

    if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
      dist = (float)sqrt((double)(dx * dx + dy * dy + dz * dz));
      avgdist += dist;
    }
    /*
      printf("%6d: after: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
    */
  }
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
    fprintf(stdout, "ellipsoid_project: avgdist = %f\n", avgdist / mris_dst->nvertices);
  MRISupdateEllipsoidSurface(mris_dst);
  if (FZERO(a - b) && FZERO(b - c)) {
    mris_dst->status = MRIS_SPHERE;
  }
  else {
    mris_dst->status = MRIS_ELLIPSOID;
  }
  return (mris_dst);
}

/*
  this one projects along the line from the origin to the ellipsoidal
  surface - not orthographic unless the ellipsoid is a sphere.
*/
MRI_SURFACE *MRISradialProjectOntoEllipsoid(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, float a, float b, float c)
{
  int vno;
  VERTEX *vsrc, *vdst;
  float x0, y0, z0, x1, y1, z1, denom, asq_bsq, asq_csq, bsq_csq, x1sq, y1sq, z1sq, abc;

  if (FZERO(a)) {
    a = DEFAULT_A;
    b = DEFAULT_B;
    c = DEFAULT_C;
  }

  if (!mris_dst) {
    mris_dst = MRISclone(mris_src);
  }

  x0 = mris_dst->xctr;
  y0 = mris_dst->yctr;
  z0 = mris_dst->zctr;
  asq_bsq = a * a * b * b;
  bsq_csq = b * b * c * c;
  asq_csq = a * a * c * c;
  abc = a * b * c;

  for (vno = 0; vno < mris_src->nvertices; vno++) {
    vsrc = &mris_src->vertices[vno];
    vdst = &mris_dst->vertices[vno];
    x1 = (vsrc->x - x0);
    y1 = (vsrc->y - y0);
    z1 = (vsrc->z - z0);
    x1sq = x1 * x1;
    y1sq = y1 * y1;
    z1sq = z1 * z1;

    /* right out of mathematica (almost) */
    denom = sqrt(bsq_csq * x1sq + asq_csq * y1sq + asq_bsq * z1sq);

    vdst->x = abc * x1 / denom /* + x0 */;
    vdst->y = abc * y1 / denom /* + y0 */;
    vdst->z = abc * z1 / denom /* + z0 */;
  }

  x0 = y0 = z0 = 0; /* set center of ellipsoid at origin */
#if 0
  if (mris_dst->v_temporal_pole)
  {
    mris_dst->v_temporal_pole->x = x0 ;
    mris_dst->v_temporal_pole->y = y0 ;
    mris_dst->v_temporal_pole->z = -c+z0 ;
    mris_dst->v_temporal_pole->tethered = TETHERED_TEMPORAL_POLE ;
  }
  if (mris_dst->v_frontal_pole)
  {
    mris_dst->v_frontal_pole->x = x0 ;
    mris_dst->v_frontal_pole->y = b+y0 ;
    mris_dst->v_frontal_pole->z = z0 ;
    mris_dst->v_frontal_pole->tethered = TETHERED_FRONTAL_POLE ;
  }
  if (mris_dst->v_occipital_pole)
  {
    mris_dst->v_occipital_pole->x = x0 ;
    mris_dst->v_occipital_pole->y = -b+y0 ;
    mris_dst->v_occipital_pole->z = z0 ;
    mris_dst->v_occipital_pole->tethered = TETHERED_OCCIPITAL_POLE ;
  }
#endif

  MRISupdateEllipsoidSurface(mris_dst);
  if (FZERO(a - b) && FZERO(b - c)) {
    mris_dst->status = MRIS_SPHERE;
  }
  else {
    mris_dst->status = MRIS_ELLIPSOID;
  }
  return (mris_dst);
}


