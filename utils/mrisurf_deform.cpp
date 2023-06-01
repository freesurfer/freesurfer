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
#include "mrisurf_deform.h"

#include "mrisurf_project.h"
#include "mrisurf_sseTerms.h"
#include "mrisurf_compute_dxyz.h"
#include "mrisurf_io.h"
#include "mrisurf_MRIS_MP.h"

#include "mrisurf_base.h"
#include "mrisutils.h"


#define MAX_VOXELS          mrisurf_sse_MAX_VOXELS
#define MAX_DISPLACEMENT    mrisurf_sse_MAX_DISPLACEMENT 
#define DISPLACEMENT_DELTA  mrisurf_sse_DISPLACEMENT_DELTA


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISremoveTopologicalDefects(MRI_SURFACE *mris, float curv_thresh)
{
  VECTOR *v_n, *v_e1, *v_e2, *v_i;
  int vno, i;
  float rsq, z, u, v;

  mrisComputeTangentPlanes(mris);
  v_n  = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_i  = VectorAlloc(3, MATRIX_REAL);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX          const * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    VECTOR_LOAD(v_n,  vertex->nx,  vertex->ny,  vertex->nz);
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);

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
        if (fabs(z / rsq) > curv_thresh) {
          mrisRemoveLink(mris, vno, vertext->v[i--]);
        }
      }
    }
  }

  /*  MRISsetRipInFacesWithRippedVertices(mris) ;*/
  VectorFree(&v_i);
  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  return (NO_ERROR);
}

int MRISsmoothOnSphere(MRIS *mris, int niters)
{
  int n, p;
  float x, y, z;

  while (niters--) {
    for (n = 0; n < mris->nvertices; n++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[n];
      VERTEX                * const v  = &mris->vertices         [n];

      x = y = z = 0.0f;

      for (p = 0; p < vt->vnum; p++) {
        VERTEX * const vp = &mris->vertices[vt->v[p]];
        x += vp->x;
        y += vp->y;
        z += vp->z;
      }
      if (vt->vnum == 0) {
        v->tx = v->x;
        v->ty = v->y;
        v->tz = v->z;
        DiagBreak();
      }
      else {
        v->tx = x / vt->vnum;
        v->ty = y / vt->vnum;
        v->tz = z / vt->vnum;
      }
      if (!std::isfinite(v->tx)) {
        DiagBreak();
      }
    }

    for (n = 0; n < mris->nvertices; n++) {
      VERTEX * const v  = &mris->vertices[n];
      mrisSphericalProjectXYZ(v->tx, v->ty, v->tz, &v->x, &v->y, &v->z);
      if (!std::isfinite(v->x)) {
        DiagBreak();
      }
    }
  }

  return NO_ERROR;
}


static void mris_project_point_into_face(
  MRIS *mris, FACE *face, int which, double x, double y, double z, double *px, double *py, double *pz)
{
  double point[3], V0[3], V1[3], V2[3], proj[3];

  point[0] = x;
  point[1] = y;
  point[2] = z;
  switch (which) {
    default:
      ErrorExit(ERROR_BADPARM, "mris_project_point_into_face: which %d not supported", which);
      break;
    case FLATTENED_VERTICES:
      V0[0] = mris->vertices[face->v[0]].fx;
      V0[1] = mris->vertices[face->v[0]].fy;
      V0[2] = mris->vertices[face->v[0]].fz;
      V1[0] = mris->vertices[face->v[1]].fx;
      V1[1] = mris->vertices[face->v[1]].fy;
      V1[2] = mris->vertices[face->v[1]].fz;
      V2[0] = mris->vertices[face->v[2]].fx;
      V2[1] = mris->vertices[face->v[2]].fy;
      V2[2] = mris->vertices[face->v[2]].fz;
      break;
    case CURRENT_VERTICES:
      V0[0] = mris->vertices[face->v[0]].x;
      V0[1] = mris->vertices[face->v[0]].y;
      V0[2] = mris->vertices[face->v[0]].z;
      V1[0] = mris->vertices[face->v[1]].x;
      V1[1] = mris->vertices[face->v[1]].y;
      V1[2] = mris->vertices[face->v[1]].z;
      V2[0] = mris->vertices[face->v[2]].x;
      V2[1] = mris->vertices[face->v[2]].y;
      V2[2] = mris->vertices[face->v[2]].z;
      break;
    case PIAL_VERTICES:
      V0[0] = mris->vertices[face->v[0]].pialx;
      V0[1] = mris->vertices[face->v[0]].pialy;
      V0[2] = mris->vertices[face->v[0]].pialz;
      V1[0] = mris->vertices[face->v[1]].pialx;
      V1[1] = mris->vertices[face->v[1]].pialy;
      V1[2] = mris->vertices[face->v[1]].pialz;
      V2[0] = mris->vertices[face->v[2]].pialx;
      V2[1] = mris->vertices[face->v[2]].pialy;
      V2[2] = mris->vertices[face->v[2]].pialz;
      break;
    case CANONICAL_VERTICES:
      V0[0] = mris->vertices[face->v[0]].cx;
      V0[1] = mris->vertices[face->v[0]].cy;
      V0[2] = mris->vertices[face->v[0]].cz;
      V1[0] = mris->vertices[face->v[1]].cx;
      V1[1] = mris->vertices[face->v[1]].cy;
      V1[2] = mris->vertices[face->v[1]].cz;
      V2[0] = mris->vertices[face->v[2]].cx;
      V2[1] = mris->vertices[face->v[2]].cy;
      V2[2] = mris->vertices[face->v[2]].cz;
      break;
  }

  project_point_to_plane(point, V0, V1, V2, proj, NULL, NULL);
  *px = proj[0];
  *py = proj[1];
  *pz = proj[2];
}

static void mrisProjectOntoSurface(MRI_SURFACE *mris, int which_vertices)
{
  int vno, fno;
  double px, py, pz, fdist;
  FACE *face;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    MHTfindClosestFaceGeneric((MHT *)(mris->mht), mris, v->x, v->y, v->z, 8.0, 8, 1, &face, &fno, &fdist);
    if (face == NULL || fdist > 0.5)  // current coord doesn't project easily into a nearby face
    {
      double lambda[3], xv, yv, zv, norm;
      int n;

      // no projection in this call - it will find a face
      MHTfindClosestFaceGeneric((MHT *)(mris->mht), mris, v->x, v->y, v->z, 1000, -1, -1, &face, &fno, &fdist);

      // crop barycentric coords to [0 1] so that it is in face (on border) and recompute coords
      face_barycentric_coords(mris, fno, which_vertices, v->x, v->y, v->z, &lambda[0], &lambda[1], &lambda[2]);

      for (n = 0; n < VERTICES_PER_FACE; n++) {
        if (lambda[n] < 0) {
          lambda[n] = 0;
        }
      }

      norm = VLEN(lambda);
      if (FZERO(norm)) {
        lambda[0] = 1.0;  // arbitrary choice to handle pathological case
      }
      else {
        SCALAR_MUL(lambda, 1.0 / norm, lambda);
      }
      lambda[2] = 1 - (lambda[0] + lambda[1]);
      if (lambda[2] < 0) {
        DiagBreak();
        lambda[2] = 0;
        lambda[1] = 1 - lambda[0];
        if (lambda[1] < 0) {
          lambda[1] = 0;
          lambda[0] = 1;
        }
      }
      px = py = pz = 0.0;
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        VERTEX *vn;
        vn = &mris->vertices[face->v[n]];
        MRISvertexCoord2XYZ_double(vn, which_vertices, &xv, &yv, &zv);
        px += lambda[n] * xv;
        py += lambda[n] * yv;
        pz += lambda[n] * zv;
      }
      if (vno == Gdiag_no) {
        printf("v %d:  lambda = (%2.2f %2.2f %2.2f)\n", vno, lambda[0], lambda[1], lambda[2]);
      }
    }
    else {
      mris_project_point_into_face(mris, face, which_vertices, v->x, v->y, v->z, &px, &py, &pz);
    }

    if (vno == Gdiag_no)
      printf("v %d: (%2.2f %2.2f %2.2f) projected to (%2.2f %2.2f %2.2f), delta = (%2.2f %2.2f %2.2f)\n",
             vno,
             v->x,
             v->y,
             v->z,
             px,
             py,
             pz,
             v->x - px,
             v->y - py,
             v->z - pz);
    v->x = px;
    v->y = py;
    v->z = pz;
    if (face_barycentric_coords(mris, fno, which_vertices, v->x, v->y, v->z, NULL, NULL, NULL) < 0) {
      DiagBreak();
    }
  }
}



template <class MRIS_Representation>
class mrisSurfaceProjector {
public:
  bool canDo(int mris_status);
  void project(MRIS_Representation* representation);
};


template<>
class mrisSurfaceProjector<MRIS_MP> {
public:
  
  bool canDo(int mris_status) { 
    switch (mris_status) {
    //case MRIS_PARAMETERIZED_SPHERE: return true;
    case MRIS_SPHERE:               return true;
    default:;
    }
    static int shown = 0;
    int shownMask = (1 << mris_status);
    if (!(shown&shownMask)) {
      shown |= shownMask;
      fs::debug() << "need to process surface status " << mris_status;
    }
    return false;
  }
  
  void project(MRIS_MP* mris_mp) {
    switch (mris_mp->status) {
      case MRIS_SPHERE:   MRISprojectOntoSphere(mris_mp, mris_mp->radius);  break;
      default: break;
    }
    cheapAssert(false);
  }
  
};


template<>
class mrisSurfaceProjector<MRIS> {
public:
  bool canDo(int mris_status) { return true; }
  void project(MRIS* mris)
  {
    switch (mris->status) {
    case MRIS_PLANE:
      MRISflattenPatch(mris);
      break;
    case MRIS_SPHERICAL_PATCH:
      mrisSphericalProjection(mris);
      break;
    case MRIS_PARAMETERIZED_SPHERE:
      MRISprojectOntoSphere(mris, mris->radius);
      break;
    case MRIS_SPHERE:
      MRISprojectOntoSphere(mris, mris->radius);
      break;
    case MRIS_ELLIPSOID:
      MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      break;
    //se PROJECT_PLANE:
    //mrisOrientPlane(mris)
    //break ;
    case MRIS_RIGID_BODY:
      break;
    case MRIS_PIAL_SURFACE:
      mrisProjectOntoSurface(mris, PIAL_VERTICES);
      break;
    default:
      // seen to happen     cheapAssert(!"was just fall thru before - Bevin looking for cases");
      break;
    }
  }
};

int mrisProjectSurface(MRIS* mris) 
{
  mrisSurfaceProjector<MRIS>().project(mris);
  return (NO_ERROR);
}

int mrisProjectSurface(MRIS_MP* mris) 
{
  mrisSurfaceProjector<MRIS_MP>().project(mris);
  return (NO_ERROR);
}

static bool mrisProjectSurface_CanDo(int mris_status)
{
  return mrisSurfaceProjector<MRIS_MP>().canDo(mris_status);
}

int mrisComputePlaneTerm(MRIS* mris, double l_plane, double l_spacing)
{
  int vno, n, vnum;
  MATRIX *M, *m_evectors;
  double dx, dy, dz, norm, dist, a, b, c, d, xc, yc, zc, dxt, dyt, dzt;
  float evalues[3];
  d = 0.0f;

  if (FZERO(l_plane) && FZERO(l_spacing)) {
    return (NO_ERROR);
  }

  M = MatrixAlloc(3, 3, MATRIX_REAL);
  m_evectors = MatrixAlloc(3, 3, MATRIX_REAL);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }

    vnum = vt->vnum;  // try with just closest nbrs
    if (vt->vnum < 4) {
      vnum = vt->v2num;
    }
    vnum = vt->vtotal;

    if (vnum < 4) {
      continue;  // can't estimate it
    }
    for (xc = yc = zc = 0.0, n = 0; n < vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      xc += vn->x;
      yc += vn->y;
      zc += vn->z;
    }
    xc /= vnum;
    yc /= vnum;
    zc /= vnum;
    MatrixClear(M);
    for (dxt = dyt = dzt = 0.0, n = 0; n < vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      *MATRIX_RELT(M, 1, 1) += SQR(vn->x - xc);
      *MATRIX_RELT(M, 2, 2) += SQR(vn->y - yc);
      *MATRIX_RELT(M, 3, 3) += SQR(vn->z - zc);

      *MATRIX_RELT(M, 2, 1) += (vn->x - xc) * (vn->y - yc);
      *MATRIX_RELT(M, 1, 2) += (vn->x - xc) * (vn->y - yc);

      *MATRIX_RELT(M, 3, 1) += (vn->x - xc) * (vn->z - zc);
      *MATRIX_RELT(M, 1, 3) += (vn->x - xc) * (vn->z - zc);

      *MATRIX_RELT(M, 3, 2) += (vn->y - yc) * (vn->z - zc);
      *MATRIX_RELT(M, 2, 3) += (vn->y - yc) * (vn->z - zc);

      dx = v->x - vn->x;
      dy = v->y - vn->y;
      dz = v->z - vn->z;
      norm = sqrt(dx * dx + dy * dy + dz * dz);
      if (DZERO(norm) || norm > .5) {
        continue;
      }
      //      norm *= norm;   // 1/r^2
      dx /= norm;
      dy /= norm;
      dz /= norm;
      norm = sqrt(dx * dx + dy * dy + dz * dz);
      dxt += l_spacing * dx;
      dyt += l_spacing * dy;
      dzt += l_spacing * dz;
    }
    dxt /= vnum;
    dyt /= vnum;
    dzt /= vnum;
    MatrixEigenSystem(M, evalues, m_evectors);

    // evalues are distance squared to plane, so use smallest one, which is in 3rd col
    a = *MATRIX_RELT(m_evectors, 1, 3);
    b = *MATRIX_RELT(m_evectors, 2, 3);
    c = *MATRIX_RELT(m_evectors, 3, 3);
    norm = sqrt(a * a + b * b + c * c);
    a /= norm;
    b /= norm;
    c /= norm;
    dist = (a * (v->x - xc) + b * (v->y - yc) + c * (v->z - zc));
    if (vno == Gdiag_no)
      printf("vno %d, best fitting plane = (%2.1f %2.1f %2.1f %2.1f), dist = %2.2f, moving by (%2.2f, %2.2f, %2.2f)\n",
             vno,
             a,
             b,
             c,
             d,
             dist,
             dist * a,
             dist * b,
             dist * c);
    dist *= l_plane;
    dxt -= dist * a;
    dyt -= dist * b;
    dzt -= dist * c;
    norm = sqrt(dxt * dxt + dyt * dyt + dzt * dzt);
    if (norm > 1) {
      dxt /= norm;
      dyt /= norm;
      dzt /= norm;
    }
    v->dx += dxt;
    v->dy += dyt;
    v->dz += dzt;
  }

  MatrixFree(&m_evectors);
  MatrixFree(&M);
  return (NO_ERROR);
}
// fit and return point/normal representation of best-fitting (LS) plane

int mrisLogStatus(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, FILE *fp, float dt, float old_sse)
{
  if (!(Gdiag & DIAG_SHOW)) {
    return (NO_ERROR);
  }

  int const doSSE = (mris->dist_alloced_flags == 3);
  int const negative = MRIScountNegativeTriangles(mris);

  int fyi = 0;
  if (parms->flags & IP_USE_MULTIFRAMES) {
    if (FZERO(parms->l_corr)) {
      /* just for your information */
      /* check if we can load curvature information */
      int n;
      for (n = 0; n < parms->nfields; n++)
        if (parms->fields[n].field == CURVATURE_CORR_FRAME) {
          parms->frame_no = parms->fields[n].frame * IMAGES_PER_SURFACE;
          fyi = 1;
          parms->l_corr = 1.0f;
          break;
        }
    }
  }

  if (parms->l_thick_min || parms->l_thick_parallel || parms->l_thick_normal || parms->l_ashburner_triangle ||
      parms->l_tspring || parms->l_thick_spring) {

    float const sse = !doSSE ? 0.0f : MRIScomputeSSE(mris, parms);

    float pct_change;
    if (old_sse > 0) {
      pct_change = 100 * (old_sse - sse) / (old_sse);
    } else {
      pct_change = 0.0;
    }
    
    fprintf(fp,
            "%3.3d: dt: %2.4f, sse: %2.1f  "
            "neg: %d (%%%2.3f:%%%2.2f), avgs: %d, %2.2f%%\n",
            parms->t,
            dt,
            sse,
            negative,
            100.0 * mris->neg_area / (mris->neg_area + mris->total_area),
            100.0 * mris->neg_orig_area / (mris->orig_area),
            parms->n_averages,
            pct_change);
  }
  else {
    float area_rms = 0, angle_rms = 0, curv_rms = 0, dist_rms = 0, corr_rms = 0;
    float sse = !doSSE ? 0.0f : mrisComputeError(mris, parms, &area_rms, &angle_rms, &curv_rms, &dist_rms, &corr_rms);

    if (fyi) {
      parms->l_corr = 0.0f;
    }

#if 0
    sse /= (float)MRISvalidVertices(mris) ;
    sse = sqrt(sse) ;
#endif

    if (mris->status == MRIS_SPHERICAL_PATCH) {
      return NO_ERROR;
    }

    if (FZERO(parms->l_corr) && FZERO(parms->l_pcorr) && ((parms->flags & IP_USE_MULTIFRAMES) == 0)) {
      fprintf(fp,
              "%3.3d: dt: %2.2f, sse: %2.1f (%2.3f, %2.1f, %2.3f), "
              "neg: %d (%%%2.3f:%%%2.2f), avgs: %d\n",
              parms->t,
              dt,
              sse,
              area_rms,
              (float)DEGREES(angle_rms),
              dist_rms,
              negative,
              100.0 * mris->neg_area / (mris->neg_area + mris->total_area),
              100.0 * mris->neg_orig_area / (mris->orig_area),
              parms->n_averages);
      if (dist_rms > 20) {
        DiagBreak();
      }
    
    } else if (parms->flags & IP_USE_MULTIFRAMES) {

      float nv = (float)MRISvalidVertices(mris);

      fprintf(fp,
              "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
              "neg: %d (%%%2.2f:%%%2.2f), avgs: %d\n",
              parms->t,
              dt,
              sse,
              area_rms,
              (float)DEGREES(angle_rms),
              dist_rms,
              corr_rms,
              negative,
              100.0 * mris->neg_area / (mris->neg_area + mris->total_area),
              100.0 * mris->neg_orig_area / (mris->orig_area),
              parms->n_averages);
      int n;
      for (n = 0; n < parms->nfields; n++) {
        if (FZERO(parms->fields[n].l_corr + parms->fields[n].l_pcorr)) {
          continue;
        }
        fprintf(stdout, "  (%d: %2.3f : %2.3f)", n, parms->fields[n].sse, sqrt(parms->fields[n].sse / nv));
      }
      fprintf(stdout, "\n");

    } else {
        fprintf(fp,
                "%3.3d: dt: %2.3f, sse: %2.1f (%2.3f, %2.1f, %2.3f, %2.3f), "
                "neg: %d (%%%2.2f:%%%2.2f), avgs: %d\n",
                parms->t,
                dt,
                sse,
                area_rms,
                (float)DEGREES(angle_rms),
                dist_rms,
                corr_rms,
                negative,
                100.0 * mris->neg_area / (mris->neg_area + mris->total_area),
                100.0 * mris->neg_orig_area / (mris->orig_area),
                parms->n_averages);
    }
  }

  fflush(fp);
  return (NO_ERROR);
}


int MRISmarkedSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];

    if (vertex->ripflag || vertex->marked == 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = vertex->x;
    y = vertex->y;
    z = vertex->z;

    sx = sy = sz = 0.0;

    n = 0;
    for (m = 0; m < vertext->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]];
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n > 0) {
      sx = sx / n;
      sy = sy / n;
      sz = sz / n;
    }
    sx = l_spring * sx; /* move in normal direction */
    sy = l_spring * sy;
    sz = l_spring * sz;

    vertex->dx += sx;
    vertex->dy += sy;
    vertex->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d marked spring term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}

int MRISnormalSpringTermWithGaussianCurvature(MRI_SURFACE *mris, double gaussian_norm, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, scale, nc, nx, ny, nz;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = vertex->x;
    y = vertex->y;
    z = vertex->z;
    nx = vertex->nx;
    ny = vertex->ny;
    nz = vertex->nz;

    sx = sy = sz = 0.0;

    n = 0;
    for (m = 0; m < vertext->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]];
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n > 0) {
      sx = sx / n;
      sy = sy / n;
      sz = sz / n;
    }
    nc = sx * nx + sy * ny + sz * nz; /* projection onto normal */
    sx = l_spring * nc * nx;          /* move in normal direction */
    sy = l_spring * nc * ny;
    sz = l_spring * nc * nz;
    scale = pow(fabs(vertex->K), gaussian_norm);
    if (!std::isfinite(scale)) {
      scale = 0;
    };
    if (scale > 1) {
      scale = 1;
    }
    scale *= l_spring;
    sx *= scale; /* move in normal direction */
    sy *= scale;
    sz *= scale;

    vertex->dx += sx;
    vertex->dy += sy;
    vertex->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d Gaussian normal term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}
int MRISspringTermWithGaussianCurvature(MRI_SURFACE *mris, double gaussian_norm, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, scale;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = vertex->x;
    y = vertex->y;
    z = vertex->z;

    sx = sy = sz = 0.0;
    n = 0;
    for (m = 0; m < vertext->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]];
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n > 0) {
      sx = sx / n;
      sy = sy / n;
      sz = sz / n;
    }
    scale = pow(fabs(vertex->K), gaussian_norm);
    if (!std::isfinite(scale)) {
      scale = 0;
    };
    if (scale > 1) {
      scale = 1;
    }
    scale *= l_spring;
    sx *= scale; /* move in normal direction */
    sy *= scale;
    sz *= scale;

    vertex->dx += sx;
    vertex->dy += sy;
    vertex->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d Gaussian normal term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}

int MRISnormalTermWithGaussianCurvature(MRI_SURFACE *mris, double lambda)
{
  int vno;
  VERTEX *v;

  if (FZERO(lambda)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    v->dx = -v->nx * v->K * lambda;
    v->dy = -v->ny * v->K * lambda;
    v->dz = -v->nz * v->K * lambda;
  }

  return (NO_ERROR);
}


/* preserve the topology of the deformation field
   this constraint is too strong - useless!
*/
int mrisApplyGradientPositiveAreaPreserving(MRI_SURFACE *mris, double dt)
{
  int vno, nvertices, n;
  float x, y, z, dx, dy, dz;
  float orig_area, area;

  int last_step = 5;
  float epsilon[6] = {1.0, 0.8, 0.6, 0.4, 0.2, 0.0}, eps;
  int step;

  int count;
  float neg_area;

  nvertices = mris->nvertices;
  MRISstoreCurrentPositions(mris);

  neg_area = 0.0;
  count = 0;
  for (n = 0; n < mris->nfaces; n++) {
    orig_area = mrisComputeArea(mris, n, 0);
    if (orig_area <= 0) {
      count++;
      neg_area += orig_area;
    }
  }

  fprintf(stderr, "before : neg = %d - area = %f \n", count, neg_area);

  for (vno = 0; vno < nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z)) {
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
    }
    if (!std::isfinite(v->dx) || !std::isfinite(v->dy) || !std::isfinite(v->dz)) {
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
    }

    dx = dt * v->dx;
    dy = dt * v->dy;
    dz = dt * v->dz;
    x = v->x;
    y = v->y;
    z = v->z;

    step = 0;
    /* preserve triangle area */
    for (n = 0; n < vt->num && step < last_step; n++) {
      mrisSphericalProjectXYZ(x, y, z, &v->x, &v->y, &v->z);
      orig_area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
      if (orig_area <= 0) {
        continue;
      }
      while (step < last_step) {
        eps = epsilon[step];
        mrisSphericalProjectXYZ(x + eps * dx, y + eps * dy, z + eps * dz, &v->x, &v->y, &v->z);
        area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
        if (area > 0) {
          break; /* we can stop here */
        }
        step++;
      }
    }

    eps = epsilon[step];
    mrisSphericalProjectXYZ(x + eps * dx, y + eps * dy, z + eps * dz, &v->x, &v->y, &v->z);
  }

  neg_area = 0.0;
  count = 0;
  for (n = 0; n < mris->nfaces; n++) {
    orig_area = mrisComputeArea(mris, n, 0);
    if (orig_area <= 0) {
      count++;
      neg_area += orig_area;
    }
  }
  fprintf(stderr, "after : neg = %d - area = %f \n", count, neg_area);

  return (NO_ERROR);
}

/* to be implemented
   constrain the gradient triangle by triangle
   but do not ensure that the final gradient is topology preserving
*/
int mrisApplyGradientPositiveAreaMaximizing(MRI_SURFACE *mris, double dt)
{
  return mrisApplyGradientPositiveAreaPreserving(mris, dt);
}

int MRISapplyGradient(MRIS* mris, double dt)
{
  MRISstoreCurrentPositions(mris);
  if (mris->status == MRIS_RIGID_BODY) {
    MRISrotate(mris, mris, dt * mris->alpha, dt * mris->beta, dt * mris->gamma);
  } else {
    MRIStranslate_along_vertex_dxdydz(mris, mris, dt);
  }
  return (NO_ERROR);
}


struct MRIScomputeSSE_asThoughGradientApplied_ctx::Impl {
  Impl() : inited(false), dx(nullptr), dy(nullptr), dz(nullptr) {}
  ~Impl() 
  {
    freeAndNULL(dx); freeAndNULL(dy); freeAndNULL(dz);
  }
  
  bool inited;
  MRIS_MP curr, orig;
  float *dx, *dy, *dz;
  
  void init(MRIS* mris) {
    cheapAssert(!inited);
    
    MRISMP_ctr(&orig);
    MRISMP_ctr(&curr);

    MRISmemalignNFloats(mris->nvertices, &dx, &dy, &dz);
    MRISMP_load(&orig, mris, true, dx,dy,dz);               // needs to load the outputs because those are inputs to ProjectSurface
      
    inited = true;
  }
};


MRIScomputeSSE_asThoughGradientApplied_ctx::MRIScomputeSSE_asThoughGradientApplied_ctx() 
  : _impl(new Impl)
{
}
  
MRIScomputeSSE_asThoughGradientApplied_ctx::~MRIScomputeSSE_asThoughGradientApplied_ctx() 
{
  delete _impl;
}

double MRIScomputeSSE_asThoughGradientApplied(
  MRIS*              mris, 
  double             delta_t, 
  INTEGRATION_PARMS* parms,
  MRIScomputeSSE_asThoughGradientApplied_ctx& ctx)
{
  bool const canUseNewBehaviour = mrisProjectSurface_CanDo(mris->status) && MRIScomputeSSE_canDo((MRIS_MP*)nullptr,parms);

  bool useOldBehaviour = !!getenv("FREESURFER_OLD_MRIScomputeSSE_asThoughGradientApplied");
  bool useNewBehaviour = !!getenv("FREESURFER_NEW_MRIScomputeSSE_asThoughGradientApplied") || !useOldBehaviour;

  useNewBehaviour = false;
  
  if (!canUseNewBehaviour) {
    useNewBehaviour = false;
    useOldBehaviour = true;
  }
  
  double new_result = 0.0;
  if (useNewBehaviour) {
    auto & ctxImpl = *ctx._impl;
    
    if (!ctxImpl.inited) ctxImpl.init(mris);
    
    MRISMP_copy(&ctxImpl.curr, &ctxImpl.orig, 
      false,  // needs to copy the outputs because those are inputs to ProjectSurface
      true);  // no need to copy v_x[*] etc. because they are obtained by the translate_along_vertex_dxdydxz below

    MRIStranslate_along_vertex_dxdydz(&ctxImpl.orig, &ctxImpl.curr, delta_t);
    mrisProjectSurface(&ctxImpl.curr);
    MRIScomputeMetricProperties(&ctxImpl.curr);
    new_result = MRIScomputeSSE(&ctxImpl.curr, parms);
  }

  double old_result = 0.0;
  if (useOldBehaviour) {
    MRISapplyGradient(mris, delta_t);
    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    old_result = MRIScomputeSSE(mris, parms);

    MRISrestoreOldPositions(mris);
  }
  
  if (useOldBehaviour && useNewBehaviour && (old_result != new_result)) {
    cheapAssert(false);
  }
  
  return useOldBehaviour ? old_result : new_result;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description: apply gradient to the sphere, so that the
  topology of the border is preserved
  (check mrisFindOptimalDefectMapping)

  ------------------------------------------------------*/

typedef struct
{
  EDGE *inside_edges, *border_edges;
  int n_inside_edges, n_border_edges;
} EDGE_LIST_INFO;

#define DEBUG_PRESERVING_GRADIENT 0
#if DEBUG_PRESERVING_GRADIENT
/* for debugging purposes */
static int ver1 = -1, ver2 = -1, ver3 = -1, ver4 = -1;
#endif

int mrisApplyTopologyPreservingGradient(MRI_SURFACE *mris, double dt, int which_gradient)
{
  int vno, nvertices, n, m;
  EDGE e1, *e2;
  FACE *face;
#if DEBUG_PRESERVING_GRADIENT
  EDGE e3, e4;
#endif
  double x, y, z, dx, dy, dz;
  float orig_area, area;
  int v1, v2, v3;

#if 1
  int last_step = 5;
  double epsilon[6] = {1.0, 0.8, 0.6, 0.4, 0.2, 0.0};
#else
  int last_step = 2;
  double epsilon[3] = {1.0, 0.5, 0.0};
#endif
  int step;

  int intersect, ninside, nborder;
  EDGE *inside, *border;

  EDGE_LIST_INFO *eli = (EDGE_LIST_INFO *)mris->vp;
  ninside = eli->n_inside_edges;
  inside = eli->inside_edges;
  nborder = eli->n_border_edges;
  border = eli->border_edges;

  nvertices = mris->nvertices;
  MRISstoreCurrentPositions(mris);

  /* just making sure */
  for (vno = 0; vno < nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->cx = v->x;
    v->cy = v->y;
    v->cz = v->z;
  }

  for (vno = 0; vno < nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    if (!std::isfinite(v->x) || !std::isfinite(v->y) || !std::isfinite(v->z)) {
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
    }

    x = v->x;
    y = v->y;
    z = v->z;

    if (which_gradient) {
      if (!std::isfinite(v->odx) || !std::isfinite(v->ody) || !std::isfinite(v->odz)) {
        ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
      }

      dx = v->odx;
      dy = v->ody;
      dz = v->odz;
    }
    else {
      if (!std::isfinite(v->dx) || !std::isfinite(v->dy) || !std::isfinite(v->dz)) {
        ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
      }

      dx = dt * v->dx;
      dy = dt * v->dy;
      dz = dt * v->dz;
    }

#if DEBUG_PRESERVING_GRADIENT
    if (vno == ver1 || vno == ver2 || vno == ver3 || vno == ver4) {
      fprintf(stderr, "\nbf %d : %f %f %f - %f %f %f  ", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);

      /* test */
      e3.vno1 = ver1;
      e3.vno2 = ver2;
      e4.vno1 = ver3;
      e4.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "\nXXX intersection %d !\n", vno);
      }
      /* test */
      e3.vno1 = ver1;
      e3.vno2 = ver2;
      e4.vno2 = ver3;
      e4.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "YYY intersection %d !\n", vno);
      }
      /* test */
      e3.vno2 = ver1;
      e3.vno1 = ver2;
      e4.vno2 = ver3;
      e4.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "ZZZ intersection %d !\n", vno);
      }
      /* test */
      e3.vno2 = ver1;
      e3.vno1 = ver2;
      e4.vno1 = ver3;
      e4.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "UUU intersection %d !\n", vno);
      }
      /* test */
      e4.vno1 = ver1;
      e4.vno2 = ver2;
      e3.vno1 = ver3;
      e3.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "AAA intersection %d !\n", vno);
      }
      /* test */
      e4.vno1 = ver1;
      e4.vno2 = ver2;
      e3.vno2 = ver3;
      e3.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "BBB intersection %d !\n", vno);
      }
      /* test */
      e4.vno2 = ver1;
      e4.vno1 = ver2;
      e3.vno2 = ver3;
      e3.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "CCC intersection %d !\n", vno);
      }
      /* test */
      e4.vno2 = ver1;
      e4.vno1 = ver2;
      e3.vno1 = ver3;
      e3.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "DDD intersection %d !\n", vno);
      }
    }
#endif

    switch (v->flags) {
      case VERTEX_CHULL:
        step = 0;
        /* preserve triangle area */
        for (n = 0; n < vt->num && step < last_step; n++) {
          mrisSphericalProjectXYZ(x, y, z, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;
          orig_area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
          if ((orig_area <= 0)) {
            //&&(parms->verbose>=VERBOSE_MODE_MEDIUM)) {

            fprintf(stderr, "negative area : this should never happen!\n");
            fprintf(stderr,
                    "face %d (%d,%d,%d) at vertex %d\n",
                    vt->f[n],
                    mris->faces[vt->f[n]].v[0],
                    mris->faces[vt->f[n]].v[1],
                    mris->faces[vt->f[n]].v[2],
                    vno);
            v1 = mris->faces[vt->f[n]].v[0];
            v2 = mris->faces[vt->f[n]].v[1];
            v3 = mris->faces[vt->f[n]].v[2];
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v1,
                    mris->vertices[v1].x,
                    mris->vertices[v1].y,
                    mris->vertices[v1].z);
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v2,
                    mris->vertices[v2].x,
                    mris->vertices[v2].y,
                    mris->vertices[v2].z);
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v3,
                    mris->vertices[v3].x,
                    mris->vertices[v3].y,
                    mris->vertices[v3].z);
            // if(parms->verbose==VERBOSE_MODE_HIGH)
            ErrorExit(ERROR_BADPARM, "mrisApplyTopologyPreservingGradient:SHOULD NOT HAPPEN\n");
          }
          while (step < last_step) {
            mrisSphericalProjectXYZ(
                x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
            v->cx = v->x;
            v->cy = v->y;
            v->cz = v->z;
            area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
            if (area > 0) {
              break; /* we can stop here */
            }
            step++;
          }
        }

/* apply gradient */
#if DEBUG_PRESERVING_GRADIENT
        // step=5;
        if (step != 0 && step != 5) {
          fprintf(stderr, "%d+", step);
        }
        else {
          fprintf(stderr, ".");
        }
#endif
        mrisSphericalProjectXYZ(
            x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
        v->cx = v->x;
        v->cy = v->y;
        v->cz = v->z;
        break;

      case VERTEX_BORDER:
        step = 0;
        /* preserve triangle area for border/chull vertices */
        for (n = 0; n < vt->num && step < last_step; n++) {
          face = &mris->faces[vt->f[n]];
          if (mris->vertices[face->v[0]].flags == VERTEX_INTERIOR ||
              mris->vertices[face->v[1]].flags == VERTEX_INTERIOR ||
              mris->vertices[face->v[2]].flags == VERTEX_INTERIOR) {
            continue;
          }
          /* test : could be removed */
          mrisSphericalProjectXYZ(x, y, z, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;
          orig_area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
          if ((orig_area <= 0)) {
            //&&(parms->verbose>=VERBOSE_MODE_MEDIUM)) {
            fprintf(stderr, "negative area : should not happen!\n");
            fprintf(stderr,
                    "face %d (%d,%d,%d) at vertex %d\n",
                    vt->f[n],
                    mris->faces[vt->f[n]].v[0],
                    mris->faces[vt->f[n]].v[1],
                    mris->faces[vt->f[n]].v[2],
                    vno);
            v1 = mris->faces[vt->f[n]].v[0];
            v2 = mris->faces[vt->f[n]].v[1];
            v3 = mris->faces[vt->f[n]].v[2];
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v1,
                    mris->vertices[v1].x,
                    mris->vertices[v1].y,
                    mris->vertices[v1].z);
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v2,
                    mris->vertices[v2].x,
                    mris->vertices[v2].y,
                    mris->vertices[v2].z);
            fprintf(stderr,
                    "cur: vertex %d (%f,%f,%f)\n",
                    v3,
                    mris->vertices[v3].x,
                    mris->vertices[v3].y,
                    mris->vertices[v3].z);
            // if(parms->verbose==VERBOSE_MODE_HIGH)
            ErrorExit(ERROR_BADPARM, "mrisApplyTopologyPreservingGradient:SHOULD NOT HAPPEN\n");
          }
          /* end of test */
          while (step < last_step) {
            mrisSphericalProjectXYZ(
                x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
            v->cx = v->x;
            v->cy = v->y;
            v->cz = v->z;
            area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
            if (area > 0) {
              break; /* we can stop here */
            }
            step++;
          }
        }

        /* check intersection/inversion with all the border edges
           For EVERY time step, all edges should not intersect anything
           This is because the border does not have to be convex!!
        */
        /* test : could be removed */
        for (n = 0; n < vt->vnum; n++) {
          VERTEX const * const vn = &mris->vertices[vt->v[n]];
#if DEBUG_PRESERVING_GRADIENT
// if(vno==412)
// fprintf(stderr,"\n%d(%d) and %d(%d)\n",vno,v->flags,v->v[n],vn->flags);
#endif
          if (vn->flags != VERTEX_BORDER) {
            continue;
          }
          e1.vno1 = vno;
          e1.vno2 = vt->v[n];
          mrisSphericalProjectXYZ(x, y, z, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;
          for (m = 0; m < ninside; m++) {
            e2 = &inside[m];
            /* intersection */
            if (edgesIntersect(mris, &e1, e2)) {
              {
                VERTEX const * v;   // NOTE HIDES OUTER DEFINITION
                
                // if(parms->verbose>=VERBOSE_MODE_MEDIUM){
                fprintf(stderr, "edge intersection : should not happen\n");
                fprintf(stderr, "edge %d-%d with edge %d %d \n", e1.vno1, e1.vno2, e2->vno1, e2->vno2);
                vno = e1.vno1;
                v = &mris->vertices[vno];
                fprintf(stderr, "%d : %f %f %f  - %f %f %f \n", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);
                vno = e1.vno2;
                v = &mris->vertices[vno];
                fprintf(stderr, "%d : %f %f %f - %f %f %f \n", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);
                vno = e2->vno1;
                v = &mris->vertices[vno];
                fprintf(stderr, "%d : %f %f %f - %f %f %f \n", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);
                vno = e2->vno2;
                v = &mris->vertices[vno];
                fprintf(stderr, "%d : %f %f %f - %f %f %f \n", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);
                //                                      if(parms->verbose==VERBOSE_MODE_HIGH)
                ErrorExit(ERROR_BADPARM, "mrisApplyTopologyPreservingGradient:SHOULD NOT HAPPEN\n");
              }
              break;
            }
          }
        }
        /* end of test */
        while (step < last_step) {
          intersect = 0;
          /* new coordinates */
          mrisSphericalProjectXYZ(
              x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;

          /* check for all edges */
          for (n = 0; (n < vt->vnum) && (!intersect); n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            if (vn->flags != VERTEX_BORDER) {
              continue;
            }
            e1.vno1 = vno;
            e1.vno2 = vt->v[n];

            for (m = 0; m < ninside; m++) {
              e2 = &inside[m];
              /* intersection */
              if (edgesIntersect(mris, &e1, e2)) {
                intersect = 1;
                break;
              }
            }
          }
          if (!intersect) {
            break;
          }
          step++;
        }

/* apply gradient */
#if DEBUG_PRESERVING_GRADIENT
        // step=5;
        if (step != 0 && step != 5) {
          fprintf(stderr, "%d_", step);
        }
        else {
          fprintf(stderr, ":");
        }
#endif
        mrisSphericalProjectXYZ(
            x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
        v->cx = v->x;
        v->cy = v->y;
        v->cz = v->z;
        break;

      case VERTEX_INTERIOR: /* preserve edge-edge intersection and angle inversion */
        step = 0;
        /* check intersection/inversion with all the border edges
           For EVERY time step, all edges should not intersect anything
           This is because the border does not have to be convex!!
        */
        /* test : could be removed */
        for (n = 0; n < vt->vnum; n++) {
          VERTEX const * const vn = &mris->vertices[vt->v[n]];
#if DEBUG_PRESERVING_GRADIENT
// if(vno==412)
// fprintf(stderr,"\n%d(%d) and %d(%d)\n",vno,v->flags,v->v[n],vn->flags);
#endif
          if (vn->flags != VERTEX_INTERIOR) {
            continue;
          }
          e1.vno1 = vno;
          e1.vno2 = vt->v[n];
          mrisSphericalProjectXYZ(x, y, z, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;
          for (m = 0; m < nborder; m++) {
            e2 = &border[m];
            /* intersection */
            if (edgesIntersect(mris, &e1, e2)) {
              {
                // if(parms->verbose>=VERBOSE_MODE_MEDIUM){
                fprintf(stderr, "Error; edge intersection : should not happen\n");
                fprintf(stderr, "edge %d-%d with edge %d %d \n", e1.vno1, e1.vno2, e2->vno1, e2->vno2);
                // if(parms->verbose==VERBOSE_MODE_HIGH)
                ErrorExit(ERROR_BADPARM, "mrisApplyTopologyPreservingGradient:SHOULD NOT HAPPEN\n");
              }
              break;
            }
          }
        }
        /* end of test */

        while (step < last_step) {
          intersect = 0;
          /* new coordinates */
          mrisSphericalProjectXYZ(
              x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
          v->cx = v->x;
          v->cy = v->y;
          v->cz = v->z;

          /* check for all edges */
          for (n = 0; n < vt->vnum && (!intersect); n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            if (vn->flags != VERTEX_INTERIOR) {
              continue;
            }
            e1.vno1 = vno;
            e1.vno2 = vt->v[n];

            for (m = 0; m < nborder; m++) {
              e2 = &border[m];
              /* intersection */
              if (edgesIntersect(mris, &e1, e2)) {
                intersect = 1;
                break;
              }
            }
          }
          if (!intersect) {
            break;
          }
          step++;
        }

/* apply gradient */
#if DEBUG_PRESERVING_GRADIENT
        // step=5;
        if (step != 0 && step != 5) {
          fprintf(stderr, "%d-", step);
        }
        else {
          fprintf(stderr, ",");
        }
#endif
        mrisSphericalProjectXYZ(
            x + epsilon[step] * dx, y + epsilon[step] * dy, z + epsilon[step] * dz, &v->x, &v->y, &v->z);
        v->cx = v->x;
        v->cy = v->y;
        v->cz = v->z;
        break;
    }
#if DEBUG_PRESERVING_GRADIENT
    if (vno == ver1 || vno == ver2 || vno == ver3 || vno == ver4) {
      fprintf(stderr, "\naf %d : %f %f %f - %f %f %f  ", vno, v->x, v->y, v->z, v->cx, v->cy, v->cz);

      /* test */
      e4.vno1 = ver1;
      e4.vno2 = ver2;
      e3.vno1 = ver3;
      e3.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "\nAAA intersection %d !\n", vno);
      }
      /* test */
      e4.vno1 = ver1;
      e4.vno2 = ver2;
      e3.vno2 = ver3;
      e3.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "BBB intersection %d !\n", vno);
      }
      /* test */
      e4.vno2 = ver1;
      e4.vno1 = ver2;
      e3.vno2 = ver3;
      e3.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "CCC intersection %d !\n", vno);
      }
      /* test */
      e4.vno2 = ver1;
      e4.vno1 = ver2;
      e3.vno1 = ver3;
      e3.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "DDD intersection %d !\n", vno);
      }
      /* test */
      e3.vno1 = ver1;
      e3.vno2 = ver2;
      e4.vno1 = ver3;
      e4.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "XXX intersection %d !\n", vno);
      }
      /* test */
      e3.vno1 = ver1;
      e3.vno2 = ver2;
      e4.vno2 = ver3;
      e4.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "YYY intersection %d !\n", vno);
      }
      /* test */
      e3.vno2 = ver1;
      e3.vno1 = ver2;
      e4.vno2 = ver3;
      e4.vno1 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "ZZZ intersection %d !\n", vno);
      }
      /* test */
      e3.vno2 = ver1;
      e3.vno1 = ver2;
      e4.vno1 = ver3;
      e4.vno2 = ver4;
      if (edgesIntersect(mris, &e3, &e4)) {
        fprintf(stderr, "UUU intersection %d !\n", vno);
      }
    }
#endif
  }
  //    fprintf(stderr,"\n\n\n");
  return (NO_ERROR);
}


#define MAX_TOTAL_MOVEMENT 5.0

#define DEBUG_V 33100

int mrisComputePosteriorTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  static const double MAX_MOVEMENT = 0.2;

  MRI *mri = parms->mri_brain;
  static MRI *mri_white_dist;
  FILE *fp;
  double dist, xs, ys, zs, dn, best_dn, best_ll, ll;
  float vmin, vmax, val, wdist;
  HISTOGRAM *hin, *hout;
  int x, y, z, vno, n;
  VECTOR *v1, *v2;
  MATRIX *m_vox2vox;
  MHT *mht;

  VOXEL_LIST **vl, **vl2;

  if (FZERO(parms->l_map)) return (NO_ERROR);

  vl = vlst_alloc(mris, MAX_VOXELS);
  vl2 = vlst_alloc(mris, MAX_VOXELS);

  hin = parms->hgm;
  hout = parms->hout;  // created in mrisComputeNegativeLogPosterior

  if (mri_white_dist == NULL) {
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, WHITE_VERTICES);
    mri_white_dist = MRIScomputeDistanceToSurface(mris, NULL, 0.5);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
  }
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);

  m_vox2vox = MRIgetVoxelToVoxelXform(mri, mri_white_dist);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_white_dist, "wd.mgz");

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }
  MRISclearMarks(mris);

  // build histogram estimates of PDFs of interior and exterior of ribbon
  for (x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (FZERO(val)) continue;
        V3_X(v1) = x;
        V3_Y(v1) = y;
        V3_Z(v1) = z;
        MatrixMultiply(m_vox2vox, v1, v2);

        if (MRIindexNotInVolume(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2))) continue;
        wdist = MRIgetVoxVal(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2), 0);
        if (wdist < 0) continue;

        // add this voxel to the list of voxels of the vertex it is closest to
        MRIvoxelToSurfaceRAS(mri, x, y, z, &xs, &ys, &zs);
        MHTfindClosestVertexGeneric(mht, xs, ys, zs, 10, 4, &vno, &dist);
        if (vno < 0) continue;
        VERTEX * v = &mris->vertices[vno];
        if (vno == Gdiag_no) DiagBreak();
        v->marked++;
        VLSTadd(vl[vno], x, y, z, xs, ys, zs);
      }

  // find vertices that don't have enough data and pool across nbrs
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) DiagBreak();
    vlst_add_to_list(vl[vno], vl2[vno]);
    if (v->ripflag || vlst_enough_data(mris, vno, vl[vno], 0.0)) continue;
    for (n = 0; n < vt->vnum; n++) {
      if (vl2[vno]->nvox + vl[vt->v[n]]->nvox >= vl2[vno]->max_vox) break;
      vlst_add_to_list(vl[vt->v[n]], vl2[vno]);
    }
    v->marked = vl[vno]->nvox;
  }

  vlst_free(mris, &vl);
  vl = vl2;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    if (v->ripflag || v->marked == 0) continue;
    if (vlst_enough_data(mris, vno, vl[vno], 0.0) == 0) {
      v->marked = 0;
      continue;
    }

    if (vno == Gdiag_no) {
      char fname[STRLEN];
      sprintf(fname, "vno%d.%d.l.dat", vno, parms->t);
      fp = fopen(fname, "w");
    }
    else
      fp = NULL;
    best_ll = -1e10;
    best_dn = 0;
    for (dn = -MAX_DISPLACEMENT; dn <= MAX_DISPLACEMENT; dn += DISPLACEMENT_DELTA) {
      ll = -vlst_loglikelihood(mris, mri, vno, dn, vl[vno], hin, hout);
      if (devIsnan(ll)) DiagBreak();
      if (fp) fprintf(fp, "%f %f\n", dn, ll);
      if (ll > best_ll) {
        best_dn = dn;
        best_ll = ll;
      }
    }
    if (fp) fclose(fp);

#if 1
    if (fabs(best_dn) > MAX_MOVEMENT) best_dn = MAX_MOVEMENT * best_dn / fabs(best_dn);
#endif
    if (vno == Gdiag_no) {
      int i;
      char fname[STRLEN];
      double dx, dy, dz, dist, dot, pin, pout;

      sprintf(fname, "vno%d.%d.vox.l.dat", vno, parms->t);
      fp = fopen(fname, "w");

      for (i = 0; i < vl[vno]->nvox; i++) {
        val = MRIgetVoxVal(mri, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], 0);
        dx = vl[vno]->xd[i] - v->x;
        dy = vl[vno]->yd[i] - v->y;
        dz = vl[vno]->zd[i] - v->z;
        dist = dx * dx + dy * dy + dz * dz;
        pin = HISTOgetCount(hin, val);
        pout = HISTOgetCount(hout, val);
        dot = dx * v->nx + dy * v->ny + dz * v->nz;
        if (dot < 0) dist *= -1;
        fprintf(fp,
                "%d %d %d %d %d %f %f %f %f\n",
                vno,
                i,
                vl[vno]->xi[i],
                vl[vno]->yi[i],
                vl[vno]->zi[i],
                val,
                dist,
                pin,
                pout);
      }
      fclose(fp);

      printf("l_map: vno %d, best displacement %2.3f, ll = %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",
             vno,
             best_dn,
             best_ll,
             best_dn * v->nx * parms->l_map,
             best_dn * v->ny * parms->l_map,
             best_dn * v->nz * parms->l_map);
      DiagBreak();
    }
    v->dx += best_dn * v->nx * parms->l_map;
    v->dy += best_dn * v->ny * parms->l_map;
    v->dz += best_dn * v->nz * parms->l_map;
    v->d = best_dn;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    int num;
    double dn = 0.0;

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    
    if (vno == Gdiag_no) DiagBreak();
    if (v->marked > 0)  // already estimated
      continue;
    for (n = num = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->marked == 0) continue;
      num++;
      dn += vn->d;
    }
    if (num > 0) {
      dn /= num;
      v->dx += dn * v->nx * parms->l_map;
      v->dy += dn * v->ny * parms->l_map;
      v->dz += dn * v->nz * parms->l_map;
      v->d = dn;
      if (vno == Gdiag_no)
        printf("l_map: vno %d, soap bubble displacement %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",
               vno,
               dn,
               dn * v->nx * parms->l_map,
               dn * v->ny * parms->l_map,
               dn * v->nz * parms->l_map);
    }
  }

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN], path[STRLEN];

    FileNamePath(mris->fname, path);
    int req = snprintf(fname,
		       STRLEN,
		       "%s/%s.%d.dist.mgz",
		       path,
		       mris->hemisphere == LEFT_HEMISPHERE ? "lh" : mris->hemisphere == BOTH_HEMISPHERES ? "both" : "rh",
		       parms->t);
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    MRISwriteD(mris, fname);
    DiagBreak();
  }
  MHTfree(&mht);
  //  HISTOfree(&hin) ; HISTOfree(&hout) ;
  VectorFree(&v1);
  VectorFree(&v2);
  MatrixFree(&m_vox2vox);
  vlst_free(mris, &vl);
  return (NO_ERROR);
}


int MRISrestoreExtraGradients(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;

  if (!mris->dx2) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->dx = mris->dx2[vno];
    v->dy = mris->dy2[vno];
    v->dz = mris->dz2[vno];
  }

  return (NO_ERROR);
}


int mrisComputePositioningGradients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  MRI *mri_brain = parms->mri_brain;
  int avgs;

  avgs = parms->n_averages;

  MHT* mht_v_orig = NULL;
  if (!FZERO(parms->l_surf_repulse)) {
    mht_v_orig = MHTcreateVertexTable(mris, ORIGINAL_VERTICES);
  }

  MHT *mht_v_current = nullptr, *mht_f_current = nullptr;
  if (!FZERO(parms->l_repulse)) {
    mht_v_current = MHTcreateVertexTable(mris, CURRENT_VERTICES);
    mht_f_current = MHTcreateFaceTable(mris);
  }

  mrisComputeTargetLocationTerm(mris, parms->l_location, parms);
  //MRISpointSetLocationError(mris, parms->l_targetpointset, parms->TargetPointSet,1);
  parms->TargetPointSet->CostAndGrad(parms->l_targetpointset,1);
  mrisComputeIntensityTerm(mris, parms->l_intensity, mri_brain, mri_brain, parms->sigma, parms);
  mrisComputeDuraTerm(mris, parms->l_dura, parms->mri_dura, parms->dura_thresh);
  mrisComputeIntensityGradientTerm(mris, parms->l_grad, mri_brain, mri_brain);
  mrisComputeSurfaceRepulsionTerm(mris, parms->l_surf_repulse, mht_v_orig);

  mrisComputeLaplacianTerm(mris, parms->l_lap);
  MRISaverageGradients(mris, avgs);

  /* smoothness terms */
  mrisComputeSpringTerm(mris, parms->l_spring);
  mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm);
  mrisComputeRepulsiveTerm(mris, parms->l_repulse, mht_v_current, mht_f_current);
  mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth, parms);
  mrisComputeThicknessMinimizationTerm(mris, parms->l_thick_min, parms);
  mrisComputeThicknessParallelTerm(mris, parms->l_thick_parallel, parms);
  mrisComputeNormalSpringTerm(mris, parms->l_nspring);
  mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv);
  /*    mrisComputeAverageNormalTerm(mris, avgs, parms->l_nspring) ;*/
  /*    mrisComputeCurvatureTerm(mris, parms->l_curv) ;*/
  mrisComputeNonlinearSpringTerm(mris, parms->l_nlspring, parms);
  mrisComputeTangentialSpringTerm(mris, parms->l_tspring);
  mrisComputeNonlinearTangentialSpringTerm(mris, parms->l_nltspring, parms->min_dist);

  if (mht_v_orig) {
    MHTfree(&mht_v_orig);
  }
  if (mht_v_current) {
    MHTfree(&mht_v_current);
  }
  if (mht_f_current) {
    MHTfree(&mht_f_current);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

double mrisRmsValError_mef(MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5, float weight30, float weight5)
{
  int vno, n, max_vno;  //, xv, yv, zv ;
  double val30, val5, total, delta, x, y, z, error, max_del;
  VERTEX *v;

  max_del = 0;
  max_vno = 0;
  for (total = 0.0, n = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    n++;
    MRISvertexToVoxel(mris, v, mri_30, &x, &y, &z);
    //      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    MRIsampleVolume(mri_30, x, y, z, &val30);
    MRIsampleVolume(mri_5, x, y, z, &val5);
    delta = (val30 - v->val);
    total += delta * delta * weight30;
    error = delta * delta;
    if (error > v->val2bak && v->val2bak > 0) {
      DiagBreak();
    }
    if (error - v->val2bak > max_del) {
      max_del = error - v->val2bak;
      max_vno = vno;
    }
    v->val2bak = error;
    delta = (val5 - v->valbak);
    total += delta * delta * weight5;
  }
  return (sqrt(total / (double)n));
}


double mrisComputeSSE_MEF(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, MRI *mri30, MRI *mri5, double weight30, double weight5, MHT *mht)
{
  double sse, l_intensity, rms;

  l_intensity = parms->l_intensity;
  parms->l_intensity = 0;

  sse = MRIScomputeSSE(mris, parms);
  parms->l_intensity = l_intensity;
  rms = mrisRmsValError_mef(mris, mri30, mri5, weight30, weight5);
  sse += l_intensity * rms * mris->nvertices;
  if (!FZERO(parms->l_surf_repulse))
    sse += parms->l_surf_repulse * mrisComputeSurfaceRepulsionEnergy(mris, parms->l_surf_repulse, mht);

  return (sse);
}


static float neg_area_ratios[] = {
    //    1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-3, 1e-2, 5e-2,1e-1
    1e-6,
    1e-5,
    1e-3,
    1e-2,
    1e-1};
#define MAX_PASSES (sizeof(neg_area_ratios) / sizeof(neg_area_ratios[0]))
int mrisRemoveNegativeArea(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int base_averages, float min_area_pct, int max_passes)
{
  unsigned int npasses;
  int total_steps, done, steps, n_averages, old_averages, niter;
  float pct_neg, ratio;
  const char *snum, *sdenom;
  float l_area, l_parea, l_corr, l_spring, l_dist, l_nlarea, *pnum, *pdenom, cmod;
  double tol;

  if (Gdiag & DIAG_WRITE && parms->fp == NULL) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.%s.out", 
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name); 
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w");
    }
    else {
      INTEGRATION_PARMS_openFp(parms, fname, "a");
    }
    if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  if (Gdiag & DIAG_SHOW) {
    mrisLogIntegrationParms(stderr, mris, parms);
  }
  pct_neg = 100.0 * mris->neg_area / (mris->neg_area + mris->total_area);
  if (pct_neg <= min_area_pct) {
    return (0); /* no steps */
  }

  tol = parms->tol;
#if 0
  parms->tol = 1e-2 ;
#endif
  niter = parms->niterations;
  old_averages = parms->n_averages;
  /*  parms->integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  /* parms->niterations = 25 ;*/
  /*  base_averages = 1024 ;*/
  l_area = parms->l_area;
  l_parea = parms->l_parea;
  l_nlarea = parms->l_nlarea;
  l_spring = parms->l_spring;
  l_dist = parms->l_dist;
  l_corr = parms->l_corr;
  parms->l_area = parms->l_parea = parms->l_dist = parms->l_corr = parms->l_spring = parms->l_nlarea = 0.0;

  /* there is one negative area removing term (area, nlarea, parea, spring),
     and one term we are seaking to retain (corr, dist).
  */
  cmod = 1.0f;
  if (!FZERO(l_corr)) {
    sdenom = "corr";
    pdenom = &parms->l_corr; /*cmod = 10.0f ;*/
  }
  else {
    sdenom = "dist";
    pdenom = &parms->l_dist;
  }

  if (!FZERO(l_area)) {
    snum = "area";
    pnum = &parms->l_area;
  }
  else if (!FZERO(l_nlarea)) {
    snum = "nlarea";
    pnum = &parms->l_nlarea;
  }
  else if (!FZERO(l_parea))
#if 0
  {
    snum = "parea" ;
    pnum = &parms->l_parea  ;
  }
#else
  {
    snum = "area";
    pnum = &parms->l_area;
  }
#endif
    else
    {
      snum = "spring";
      pnum = &parms->l_spring;
    }

  for (total_steps = npasses = 0; npasses < MAX_PASSES; npasses++) {
    *pnum = 1.0;
    *pdenom = cmod * (npasses >= MAX_PASSES ? neg_area_ratios[MAX_PASSES - 1] : neg_area_ratios[npasses]);

    ratio = *pnum / *pdenom;

    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "%s/%s = %2.3f\n", snum, sdenom, ratio);
    }
    if (Gdiag & DIAG_WRITE) {
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio);
    }
    if (ratio > 10000) {
      DiagBreak();
    }
    for (done = 0, n_averages = base_averages; !done; n_averages /= 4) {
      parms->n_averages = n_averages;
      steps = MRISintegrate(mris, parms, n_averages);
      parms->start_t += steps;
      total_steps += steps;
      pct_neg = 100.0 * mris->neg_area / (mris->neg_area + mris->total_area);
#if 0
      if (pct_neg < min_area_pct)
      {
        break ;
      }
#endif
      done = n_averages == 0;
      /* finished integrating at smallest scale */
    }
#if 0
    if (pct_neg < min_area_pct)
    {
      break ;
    }
#endif
  }
#if 0
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
#endif
  parms->n_averages = old_averages; /* hack, but no time to clean up now */
  parms->l_area = l_area;
  parms->l_parea = l_parea;
  parms->l_spring = l_spring;
  parms->l_dist = l_dist;
  parms->l_corr = l_corr;
  parms->niterations = niter;
  parms->tol = tol;
  return (total_steps);
}


int MRISflattenPatchRandomly(MRI_SURFACE *mris)
{
  float extent;
  int vno;

  extent = sqrt(mris->total_area);
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].x = randomNumber(-extent, extent);
    mris->vertices[vno].y = randomNumber(-extent, extent);
    mris->vertices[vno].z = 0;
    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout,
              "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",
              vno,
              mris->vertices[vno].x,
              mris->vertices[vno].y,
              mris->vertices[vno].z);
  }
  mris->status = MRIS_PLANE;
  MRIScomputeMetricProperties(mris);
  mrisOrientPlane(mris);
  return (NO_ERROR);
}


int MRISflattenPatch(MRI_SURFACE *mris)
{
  float x, y, z, d, d1, d2;
  float nx, ny, nz;
  VERTEX *v;
  int k, an;

  if (mris->status == MRIS_PLANE)  // already flattened
  {
    for (k = 0; k < mris->nvertices; k++)
      if (!mris->vertices[k].ripflag) {
        v = &mris->vertices[k];
        v->z = 0;
      }
    return (NO_ERROR);
  }

  x = y = z = nx = ny = nz = 0;
  an = 0;

  /* calculate average normal and vertex position */
  MRIScomputeNormals(mris);
  for (k = 0; k < mris->nvertices; k++)
    if (!mris->vertices[k].ripflag) {
      v = &mris->vertices[k];
      x += v->x;
      y += v->y;
      z += v->z;
#if 0
      if (!FZERO(v->nx))
        fprintf(stdout, "vertex %d, normal = (%2.3f, %2.3f, %2.3f)\n",
                k, v->nx, v->ny, v->nz) ;
#endif
      nx += v->nx;
      ny += v->ny;
      nz += v->nz;
      an++;
    }
  x /= an;
  y /= an;
  z /= an;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "flatten: avg p = {%2.1f, %2.1f, %2.1f}\n", x, y, z);
    fprintf(stdout, "flatten: sum n = {%2.2f, %2.2f, %2.2f}\n", nx, ny, nz);
  }
#if 0
  /* or override with direct front,back */
  if (project==POSTERIOR)
  {
    nx = nz = 0.0;
    ny = -1.0;
  }
  if (project==ANTERIOR)
  {
    nx = nz = 0.0;
    ny = 1.0;
  }
#endif

  /* make the average normal unit length */
  d = sqrt(nx * nx + ny * ny + nz * nz);
  nx /= d;
  ny /= d;
  nz /= d;
  d = sqrt(nx * nx + ny * ny);
  if (!FZERO(d)) /* not already in a plane */
  {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      fprintf(stdout, "flatten: norm n = {%2.2f, %2.2f, %2.2f}\n", nx, ny, nz);
    }
    for (k = 0; k < mris->nvertices; k++)
      if (!mris->vertices[k].ripflag) {
        v = &mris->vertices[k];
        v->x -= x;
        v->y -= y;
        v->z -= z;
        d1 = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
        transform(&v->x, &v->y, &v->z, nx, ny, nz, d);
        d2 = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
        if (fabs(d1 - d2) > 0.0001) {
          printf("flatten: d1=%f, d2=%f\n", d1, d2);
        }
        transform(&v->nx, &v->ny, &v->nz, nx, ny, nz, d);
      }

/* print transform matrix in tmp dir */
#if 0
    sprintf(fname,"%s/surfer.mat",dir);
    fp = fopen(fname,"w");
    if (fp==NULL)
    {
      ErrorPrintf(ERROR_NOFILE, "flatten: can't create file %s\n",fname);
    }
    else
    {
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",
              nx*nz/d,  -nx,  ny/d,  0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",
              d,       nz,   0.0, 0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",
              -ny*nz/d,  ny,   nx/d, 0.0);
      fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",
              0.0,      0.0,  0.0, 1.0);
      fclose(fp);
      if (Gdiag & DIAG_SHOW)
      {
        printf("flatten: file %s written\n",fname);
      }
    }
#endif

    transform(&nx, &ny, &nz, nx, ny, nz, d);
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "flatten: transformed n = {%2.1f, %2.1f, %2.1f}\n", nx, ny, nz);
  }
  for (k = 0; k < mris->nvertices; k++) {
    mris->vertices[k].z = 0;
    if (k == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout,
              "vertex %d flattened @ (%2.2f, %2.2f, %2.2f)\n",
              k,
              mris->vertices[k].x,
              mris->vertices[k].y,
              mris->vertices[k].z);
  }

  mris->status = MRIS_PLANE;
  MRIScomputeMetricProperties(mris);
  if (Gdiag & DIAG_SHOW && Gdiag_no >= 0) {
    int n;
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];
    VERTEX          const * const v  = &mris->vertices         [Gdiag_no];
    fprintf(stdout,
            "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, "
            "oa %2.3f, nz=%2.3f, vnum=%d\n",
            Gdiag_no,
            v->x,
            v->y,
            v->z,
            v->area,
            v->origarea,
            v->nz,
            vt->vnum);
    for (n = 0; n < vt->vnum; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      VERTEX          const * const vn  = &mris->vertices         [vt->v[n]];
      fprintf(stdout,
              "%d @ (%2.1f, %2.1f, %2.1f): area %2.3f, oa %2.3f, nz=%2.3f, "
              "vnum=%d, d=%2.2f, od=%2.2f\n",
              vt->v[n],
              vn->x,
              vn->y,
              vn->z,
              vn->area,
              v->origarea,
              vn->nz,
              vnt->vnum,
              v->dist[n],
              v->dist_orig[n]);
    }
    fprintf(stdout, "\n");
  }

  return (NO_ERROR);
}


#define EXPLODE_ITER 1
int mrisTearStressedRegions(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int    vno, nrip = 0 ;
  double max_stress = 0 ;
  static int ncalls = 0 ;
  ncalls++ ;

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX *v ;

    v = &mris->vertices[vno] ;
    if (v->ripflag)
      nrip++ ;
    v->oripflag = v->ripflag ;
  }
//  printf("starting with %d of %d ripped\n", nrip, mris->nvertices) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (!v->oripflag)
    {
      double x, y, z, nx, ny, nz, sd, sx, sy, sz, dx, dy, dz, stress ;
      int    n, m ;

      x = v->tx;  y = v->ty; z = v->tz;
      nx = v->nx;  ny = v->ny; nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m = 0 ; m < vt->vnum ; m++)
	if (!mris->vertices[vt->v[m]].oripflag)
          {
	    double d ;
            sx += dx = mris->vertices[vt->v[m]].tx - x;
            sy += dy = mris->vertices[vt->v[m]].ty - y;
            sz += dz = mris->vertices[vt->v[m]].tz - z;
            d = sqrt(dx*dx+dy*dy+dz*dz);
	    sd += d ;
#if 0
	    d = sqrt(d) ;
	    if (d> max_stress)
	      max_stress = d ;
	    if (d > parms->stressthresh && parms->explode_flag && ncalls>EXPLODE_ITER)
	    {
	      nrip++ ;
	      mrisRemoveLink(mris, vno, vt->v[m]);
	      mris->patch = 1 ;
	    }
#endif	      
            n++;
          }
        if (n>0)
        {
          sx = sx/n;
          sy = sy/n;
          sz = sz/n;
          sd = sd/n;
          stress = sd;
#if 1
	  if (stress > max_stress)
	    max_stress = stress ;
          if (parms->explode_flag && stress>=parms->stressthresh)
          {
            nrip++;
	    mris->patch = 1 ;
//            v->ripflag = TRUE;
	    for (m = 0 ; m < vt->vnum ; m++)
	      mrisRemoveLink(mris, vno, vt->v[m]);
          }
#endif
	}
    }
  }

  MRISremoveRipped(mris) ;
  if (parms->explode_flag)
    printf("max stress = %2.2f, nrip = %d, nv = %d, nf = %d\n", max_stress, nrip, mris->nvertices, mris->nfaces) ;
  if (parms->explode_flag && ncalls == EXPLODE_ITER && parms->stressthresh < 0)
  {
    parms->stressthresh = max_stress ;
    printf("!!!!!!!!! setting thresh flag to %2.2f\n", parms->stressthresh) ;
  }
  if (parms->explode_flag && ncalls >= EXPLODE_ITER)
  {
    int eno, nvertices, nfaces, nedges ;

    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
    if (nedges < 300000)
      DiagBreak() ;
  }
  return(NO_ERROR) ;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisSmoothBoundaryNormals(MRI_SURFACE *mris, int niter)
{
#if 0
  int iter,k,m,n;
  vertex_type *v;
  float sumx,sumy,r;

  for (iter=0; iter<niter; iter++)
  {
    for (k=0; k<mris->nvertices; k++)
      if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
      {
        mris->vertices[k].obnx = mris->vertices[k].bnx;
        mris->vertices[k].obny = mris->vertices[k].bny;
      }
    for (k=0; k<mris->nvertices; k++)
      if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
      {
        v = &mris->vertices[k];
        n = 1;
        sumx = v->obnx;
        sumy = v->obny;
        for (m=0; m<v->vnum; m++)
          if ((!mris->vertices[v->v[m]].ripflag)&&
              mris->vertices[v->v[m]].border)
          {
            sumx += mris->vertices[v->v[m]].obnx;
            sumy += mris->vertices[v->v[m]].obny;
            n++;
          }
        v->bnx = (n>0)?sumx/n:0;
        v->bny = (n>0)?sumy/n:0;
      }
  }
  for (k=0; k<mris->nvertices; k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
      if (r>0)
      {
        mris->vertices[k].bnx /= r;
        mris->vertices[k].bny /= r;
      }
    }
#endif
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum, n;
  float x, y, z, num;
  int nmarked;

  MRIScopyMarkedToMarked3(mris);
  for (i = 0; i < navgs; i++) {
    if (Gdiag_no >= 0 && (mris->vertices[Gdiag_no].marked == 0)) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];
      VERTEX                * const v  = &mris->vertices         [Gdiag_no];
      printf("v %d @ (%2.1f, %2.1f, %2.1f) will be moved during soap bubble:\n", Gdiag_no, v->x, v->y, v->z);
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        printf("\tnbr %d=%d, marked = %d, (%2.1f, %2.1f, %2.1f)\n", n, vt->v[n], vn->marked, vn->x, vn->y, vn->z);
      }
      DiagBreak();
    }
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      x = y = z = 0;
      num = 0;
      /*      if (v->marked == 2)*/
      {
        x = v->x;
        y = v->y;
        z = v->z;
        num++; /* account for central vertex */
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag) /* no valid data */
        continue;
        num++;
        x += vn->x;
        y += vn->y;
        z += vn->z;
      }
      if (num > 0) {
        v->tdx = x / num;
        v->tdy = y / num;
        v->tdz = z / num;
        if (!v->marked) {
          nmarked++;
        }
        v->marked = 3; /* need modification */
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      if (vno == Gdiag_no) {
        printf("iter %d: moving v %d from (%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
               i,
               vno,
               v->tdx,
               v->tdy,
               v->tdz,
               v->x,
               v->y,
               v->z);
      }
      if (v->marked) {
        v->x = v->tdx;
        v->y = v->tdy;
        v->z = v->tdz;
      }
      if (v->marked == 3) /* needs modification */
      {
        v->marked = 2; /* modified, but not fixed */
      }
    }
    if (nmarked && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
      printf("%d: %d vertices marked\n", i, nmarked);
    }
  }

  MRIScopyMarked3ToMarked(mris);
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISweightedSoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum, n;
  float x, y, z;
  int nmarked;
  double norm, total_norm;

  for (i = 0; i < navgs; i++) {
    if (Gdiag_no >= 0 && (mris->vertices[Gdiag_no].marked == 0)) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];
      VERTEX          const * const v  = &mris->vertices         [Gdiag_no];
      printf("v %d @ (%2.1f, %2.1f, %2.1f) will be moved during soap bubble:\n", Gdiag_no, v->x, v->y, v->z);
      for (n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        printf("\tnbr %d=%d, marked = %d, (%2.1f, %2.1f, %2.1f)\n", n, vt->v[n], vn->marked, vn->x, vn->y, vn->z);
      }
      DiagBreak();
    }
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) {
        v->tdx = v->x;
        v->tdy = v->y;
        v->tdz = v->z;
        continue;
      }
      if (vno == Gdiag_no) {
        DiagBreak();
      }

      total_norm = 1.0;  // central vertex gets unit weight
      x = v->x;
      y = v->y;
      z = v->z;
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        norm = 1.0;
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag)              /* no valid data */
        {
          norm *= 20;
        }
        x += norm * vn->x;
        y += norm * vn->y;
        z += norm * vn->z;
        total_norm += norm;
      }
      if (total_norm > 0) {
        v->tdx = x / total_norm;
        v->tdy = y / total_norm;
        v->tdz = z / total_norm;
        if (!v->marked) {
          nmarked++;
        }
        v->marked = 3; /* need modification */
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      if (vno == Gdiag_no) {
        printf("iter %d: moving v %d from (%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
               i,
               vno,
               v->tdx,
               v->tdy,
               v->tdz,
               v->x,
               v->y,
               v->z);
      }
      if (v->marked) {
        v->x = v->tdx;
        v->y = v->tdy;
        v->z = v->tdz;
      }
      if (v->marked == 3) /* needs modification */
      {
        v->marked = 2; /* modified, but not fixed */
      }
    }
    if (nmarked && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
      printf("%d: %d vertices marked\n", i, nmarked);
    }
    if (DIAG_VERBOSE_ON) MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES);
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISsoapBubbleOrigVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float x, y, z, num;
  int nmarked, num_none_marked = 0;

  /*
    v->marked:

    0 -  never processed
    1 -  fixed value
    2 -  has had value computed via previous soap bubble.
    3 -  has had value computed on this soap bubble, but not yet usable.
  */
  for (i = 0; i < navgs; i++) {
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      x = y = z = 0;
      num = 0;
      if (v->marked == 2) /* computed on previous iteration, use value */
      {
        x = v->origx;
        y = v->origy;
        z = v->origz;
        num++; /* account for central vertex */
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++];                     /* neighboring vertex pointer */
        if (vn->ripflag || !vn->marked || vn->marked > 2) /* no valid data */
        {
          continue;
        }
        num++;
        x += vn->origx;
        y += vn->origy;
        z += vn->origz;
      }
      if (num > 0) {
        v->tdx = x / num;
        v->tdy = y / num;
        v->tdz = z / num;
        if (!v->marked) {
          nmarked++;
        }
        v->marked = 3; /* need modification */
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      if (v->marked) /* update value */
      {
        MRISsetOriginalXYZ(mris, vno, v->tdx, v->tdy, v->tdz);
      }
      if (v->marked == 3) /* needs modification */
      {
        v->marked = 2; /* modified, but not fixed */
      }
    }
    if (Gdiag & DIAG_SHOW) {
      printf("%d: %d vertices marked\n", i, nmarked);
    }
    if (!nmarked && ++num_none_marked > 5) {
      break;
    }
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag || v->marked == 1) {
      continue;
    }
    v->marked = 0;
  }
  return (NO_ERROR);
}

int MRISsoapBubbleTargetVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float x, y, z, num;
  int nmarked, num_none_marked = 0;

  /*
    v->marked:

    0 -  never processed
    1 -  fixed value
    2 -  has had value computed via previous soap bubble.
    3 -  has had value computed on this soap bubble, but not yet usable.
  */
  for (i = 0; i < navgs; i++) {
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no) DiagBreak();
      if (v->ripflag || v->marked == 1) continue;

      x = y = z = 0;
      num = 0;
      if (v->marked == 2) /* computed on previous iteration, use value */
      {
        x = v->targx;
        y = v->targy;
        z = v->targz;
        num++; /* account for central vertex */
      }
      int const * pnb  = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++];                     /* neighboring vertex pointer */
        if (vn->ripflag || !vn->marked || vn->marked > 2) /* no valid data */
        {
          continue;
        }
        num++;
        x += vn->targx;
        y += vn->targy;
        z += vn->targz;
      }
      if (vno == Gdiag_no) DiagBreak();
      if (num > 0) {
        if (vno == Gdiag_no) DiagBreak();
        v->tdx = x / num;
        v->tdy = y / num;
        v->tdz = z / num;
        if (!v->marked) {
          nmarked++;
        }
        v->marked = 3; /* need modification */
      }
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag || v->marked == 1) {
        continue;
      }
      if (v->marked) /* update value */
      {
        v->targx = v->tdx;
        v->targy = v->tdy;
        v->targz = v->tdz;
      }
      if (v->marked == 3) /* needs modification */
      {
        v->marked = 2; /* modified, but not fixed */
      }
    }
    if (Gdiag & DIAG_SHOW) {
      printf("%d: %d vertices marked\n", i, nmarked);
    }
    if (!nmarked && ++num_none_marked > 5) {
      break;
    }
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag || v->marked == 1) {
      continue;
    }
    v->marked = 0;
  }
  return (NO_ERROR);
}


/* extract the non-marked vertices from a surface */
MRIS *MRISextractMarkedVertices(MRIS *mris)
{
  int *face_trans, *vertex_trans;
  FACE *f, *fdst;
  int vno, i, n, fno;

  face_trans = (int *)calloc(mris->nfaces, sizeof(int));
  vertex_trans = (int *)calloc(mris->nvertices, sizeof(int));
  memset(vertex_trans, -1, mris->nvertices * sizeof(int));
  memset(face_trans, -1, mris->nfaces * sizeof(int));
  // create a new surface
  MRIS *mris_corrected = MRISalloc(mris->nvertices, mris->nfaces);
  // keep the extra info into the new one
  mris_corrected->useRealRAS = mris->useRealRAS;
  copyVolGeom(&mris->vg, &mris_corrected->vg);

  mris_corrected->type = MRIS_TRIANGULAR_SURFACE;
  mris_corrected->status         = mris->status;            // this should have been done
  mris_corrected->origxyz_status = mris->origxyz_status;    // this is new
  
  int newNVertices;
  for (newNVertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    /* ignore the marked defect vertices but not the bordering ones */
    if (v->marked) {
      continue;
    }
    VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[newNVertices];
    VERTEX          * const vdst  = &mris_corrected->vertices         [newNVertices];
    
    /* original vertices */
    vdst->x = v->x;
    vdst->y = v->y;
    vdst->z = v->z;
    /* smoothed vertices */
    
    MRISsetOriginalXYZ(mris_corrected, newNVertices,
      v->origx, v->origy, v->origz);
    
    vdst->tx = v->tx;
    vdst->ty = v->ty;
    vdst->tz = v->tz;
    vdst->nx = v->nx;
    vdst->ny = v->ny;
    vdst->nz = v->nz;
    /* canonical vertices */
    vdst->cx = v->cx;
    vdst->cy = v->cy;
    vdst->cz = v->cz;
    vdstt->num = vt->num;
    vdst->val = v->val;
    vdst->val2 = v->val2;
    vdst->valbak = v->valbak;
    vdst->val2bak = v->val2bak;
    vdst->imag_val = v->imag_val;
    vdst->curv = v->curv;
    vdst->curvbak = v->curvbak;
    vdst->stat = v->stat;
    vdst->mean = v->mean;
    vdst->mean_imag = v->mean_imag;
    vdst->std_error = v->std_error;
    vdst->H = v->H;
    vdst->K = v->K;
    vdst->k1 = v->k1;
    vdst->k2 = v->k2;
    vdst->border = 0;
    vertex_trans[vno] = newNVertices++;
  }
  MRIStruncateNVertices(mris_corrected, newNVertices);

  int newNfaces = 0;
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    /* don't update triangle with marked vertices */
    if (triangleMarked(mris, fno)) {
      continue;
    }
    /* initialize face */
    fdst = &mris_corrected->faces[newNfaces];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fdst->v[n] = vertex_trans[f->v[n]];
    }
    face_trans[fno] = newNfaces++;
  }
  MRIStruncateNFaces(mris_corrected, newNfaces); 

  /* now allocate face and neighbor stuff in mris_corrected */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (v->marked) {
      continue;
    }
    if (vertex_trans[vno] < 0 || vertex_trans[vno] >= mris_corrected->nvertices) {
      continue;
    }
    
    int const vno_dst = vertex_trans[vno];
    VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[vno_dst];

    /* count # of good triangles attached to this vertex */
    for (vdstt->num = n = 0; n < vt->num; n++)
      if (triangleMarked(mris, vt->f[n]) == 0) {
        vdstt->num++;
      }
    vdstt->f = (int *)calloc(vdstt->num, sizeof(int));
    vdstt->n = (uchar *)calloc(vdstt->num, sizeof(uchar));
    for (i = n = 0; n < vt->num; n++) {
      if (triangleMarked(mris, vt->f[n])) {
        continue;
      }
      vdstt->n[i] = vt->n[n];
      vdstt->f[i] = face_trans[vt->f[n]];
      i++;
    }

    /* count # of valid neighbors */
    clearVnum(mris_corrected,vno_dst);
    for (n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].marked == 0) {
        addVnum(mris_corrected,vno_dst,1);
      }
    vdstt->nsizeMax = 1;
    vdstt->v = (int *)calloc(vdstt->vnum, sizeof(int));
    for (i = n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].marked == 0) {
        vdstt->v[i++] = vertex_trans[vt->v[n]];
      }
      
    MRIS_setNsizeCur(mris_corrected,vno_dst,1);
  }

  mrisCheckVertexFaceTopology(mris_corrected);

  return mris_corrected;
}

/* extract the main connected component from a triangulation
   use the v->marked to extract the surface */
MRIS *MRISextractMainComponent(MRI_SURFACE *mris, int do_not_extract, int verbose, int *ncpts)
{
  MRIS *mris_out;

  int n, vn0, vn1, vn2, p, count, max_c, max_nbr, found, ncpt;
  int nv, ne, nf, nX, tX, tv, te, tf;
  float cx, cy, cz;

  /* use marked field for counting purposes */
  MRISclearMarks(mris);

  tX = tv = te = tf = 0;
  max_c = 0;
  max_nbr = 0;
  if (verbose) {
    fprintf(WHICH_OUTPUT, "\ncounting number of connected components...");
  }
  for (count = 0, vn0 = 0; vn0 < mris->nvertices; vn0++) {
    VERTEX                * const v  = &mris->vertices         [vn0];
    if (v->marked) {
      continue;
    }
    count++;
    v->marked = count;
    found = 1;
    ncpt = 1;
    while (found) {
      found = 0;
      for (vn1 = 0; vn1 < mris->nvertices; vn1++) {
        VERTEX_TOPOLOGY const * const vp1t = &mris->vertices_topology[vn1];
        VERTEX                * const vp1  = &mris->vertices         [vn1];
        if (vp1->marked) {
          continue;
        }
        /* check neighbors */
        for (p = 0; p < vp1t->vnum; p++) {
          vn2 = vp1t->v[p];
          VERTEX * const vp2 = &mris->vertices[vn2];
          if (vp2->marked == count) {
            vp1->marked = count;
            found = 1;
            ncpt++;
            break;
          }
        }
      }
    }
    if (max_nbr < ncpt) {
      max_c = count;
      max_nbr = ncpt;
    }
    /* computing Euler Number of component */
    ne = 0;
    nv = 0;
    cx = cy = cz = 0.0f;
    for (vn1 = 0; vn1 < mris->nvertices; vn1++) {
      VERTEX_TOPOLOGY const * const vp1t = &mris->vertices_topology[vn1];
      VERTEX                * const vp1  = &mris->vertices         [vn1];
      if (vp1->marked != count) {
        continue;
      }
      ne += vp1t->vnum;
      nv++;
      cx += vp1->x;
      cy += vp1->y;
      cz += vp1->z;
    }
    cx /= nv;  // center of gravity
    cy /= nv;
    cz /= nv;
    ne /= 2;  // half the number of edges
    /* now counting faces */
    nf = 0;
    for (vn1 = 0; vn1 < mris->nfaces; vn1++) {
      if (mris->vertices[mris->faces[vn1].v[0]].marked == count) {
        nf++;
      }
    }
    nX = nv - ne + nf;
    tX += nX;
    tv += nv;
    te += ne;
    tf += nf;
    if (verbose)
      fprintf(stderr,
              "\n   %d voxel in cpt #%d: X=%d [v=%d,e=%d,f=%d] located at (%f, %f, %f)",
              ncpt,
              count,
              nX,
              nv,
              ne,
              nf,
              cx,
              cy,
              cz);
  }
  if (verbose) {
    fprintf(stderr, "\nFor the whole surface: X=%d [v=%d,e=%d,f=%d]", tX, tv, te, tf);
  }

  if (count <= 1) {
    if (verbose) {
      fprintf(WHICH_OUTPUT, "\nOne single component has been found");
    }
    if (!do_not_extract && verbose) {
      fprintf(WHICH_OUTPUT, "\nnothing to do");
    }
  }
  else {
    if (verbose) {
      fprintf(WHICH_OUTPUT, "\n%d components have been found", count);
    }
    if (!do_not_extract && verbose) fprintf(WHICH_OUTPUT, "\nkeeping component #%d with %d vertices", max_c, max_nbr);
  }

  if (ncpts) {
    *ncpts = count;
  }

  /* nothing to do */
  if (do_not_extract) {
    return mris;
  }

  /* set vertex mark to 1 if removed vertex */
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX * const v = &mris->vertices[n];
    if (v->marked == max_c) {
      v->marked = 0;
    }
    else {
      v->marked = 1;
    }
  }

  mris_out = MRISextractMarkedVertices(mris);
  MRISclearMarks(mris);
  return mris_out;
}

int MRIScombine(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total, MRIS_HASH_TABLE *mht, int which)
{
  if (which == VERTEX_COORDS) {
    cheapAssert(mris_src->origxyz_status == mris_total->origxyz_status);
  }
  
  int vno;
  VERTEX *v, *vdst;
  MHT *mht_src = NULL;
  double max_len, mean;

  MRISclearMarks(mris_total);
  for (vno = 0; vno < mris_total->nvertices; vno++) {
    vdst = &mris_total->vertices[vno];
    vdst->d = 0;
  }

  for (vno = 0; vno < mris_src->nvertices; vno++) {
    v = &mris_src->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    float min_dist;
    vdst = MHTfindClosestVertex2(mht, mris_total, mris_src, v, &min_dist);
    if (!vdst) {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map vno %d", vno);
      continue;
    }
    if (vdst - mris_total->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    switch (which) {
      case VERTEX_COORDS: {
        int const vdst_vno = vdst - mris_total->vertices;   // GROSS HACK since MHTfindClosestVertex doesn't return the vno
        MRISsetOriginalXYZ(mris_total, vdst_vno,
          vdst->origx + v->origx,
          vdst->origy + v->origy,
          vdst->origz + v->origz);
        } break;
      case VERTEX_AREA:
        vdst->d += v->origarea;
        break;
      case VERTEX_CURV:
        vdst->d += v->curv;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->d += v->val;
        break;
      case VERTEX_ANNOTATION:
        vdst->annotation = v->annotation;
        break;
    }
  }

  /* normalize by # of times a vertex is mapped to */
  for (vno = 0; vno < mris_total->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vdst = &mris_total->vertices[vno];
    if (vdst->ripflag || !vdst->marked) {
      continue;
    }
    mean = vdst->d / (float)vdst->marked;
    switch (which) {
      case VERTEX_COORDS:
        MRISsetOriginalXYZ(mris_total, vno,
          vdst->origx / (float)vdst->marked,
          vdst->origy / (float)vdst->marked,
          vdst->origz / (float)vdst->marked);
        break;
      case VERTEX_AREA: /* don't normalize by # of vertices mapped!! */
        vdst->origarea += vdst->d;
        vdst->val2 += vdst->d * vdst->d;
        break;
      case VERTEX_CURV:
        vdst->curv += mean;
        vdst->val2 += mean * mean;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->val += mean;
        vdst->val2 += mean * mean;
        break;
    }
  }

  /* sample from dst to source to fill holes */
  for (vno = 0; vno < mris_total->nvertices; vno++) {
    vdst = &mris_total->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (vdst->marked || vdst->ripflag) {
      continue;
    }
    if (!mht_src) {
      double mean, sigma;

      mean = MRIScomputeVertexSpacingStats(mris_src, &sigma, NULL, &max_len, NULL, NULL, CURRENT_VERTICES);
      if (max_len > mean + 3 * sigma) {
        max_len = mean + 3 * sigma;
      }
      mht_src = MHTcreateVertexTable_Resolution(mris_src, CURRENT_VERTICES, 2 * max_len);
    }
    float min_dist;
    v = MHTfindClosestVertex2(mht_src, mris_src, mris_total, vdst, &min_dist);
    if (!v) {
      ErrorPrintf(ERROR_BADPARM, "MRIScombine: cannot map dst vno %d", vno);
      continue;
    }
    if (v - mris_src->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    switch (which) {
      case VERTEX_COORDS:
        MRISsetOriginalXYZ(mris_total, vno, v->origx, v->origy, v->origz);
        break;
      case VERTEX_ANNOTATION:
        vdst->annotation = v->annotation;
        break;
      case VERTEX_AREA:
        vdst->origarea += v->origarea;
        vdst->val2 += (v->origarea * v->origarea);
        break;
      case VERTEX_CURV:
        vdst->curv += v->curv;
        vdst->val2 += (v->curv * v->curv);
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->val += v->val;
        vdst->val2 += (v->val * v->val);
        break;
    }
  }
  if (mht_src) {
    MHTfree(&mht_src);
  }

  return (NO_ERROR);
}

#if 1
int MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, MRIS_HASH_TABLE *mht, int which)
{
  switch (which) {
    case VERTEX_COORDS:
      cheapAssert(mris_dst->origxyz_status == mris_src->origxyz_status);
      break;
  }
  
  int vno;
  MHT *mht_src = NULL;
  double max_len;

  MRISclearWhichAndVal2(mris_dst, which);
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    VERTEX * const vdst = &mris_dst->vertices[vno];
    vdst->d = 0;
    vdst->val2 = 0;
  }

  MRIScomputeVertexSpacingStats(mris_src, NULL, NULL, &max_len, NULL, NULL, CURRENT_VERTICES);
  mht_src = MHTcreateVertexTable_Resolution(mris_src, CURRENT_VERTICES, 2 * max_len);

  /* sample from dst to source to fill holes */
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    VERTEX * const vdst = &mris_dst->vertices[vno];
    if (vdst->ripflag) {
      continue;
    }
    float min_dist;
    VERTEX const * const v = MHTfindClosestVertex2(mht_src, mris_src, mris_dst, vdst, &min_dist);
    if (!v) {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map dst vno %d", vno);
      continue;
    }
    if (v - mris_src->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    vdst->val2 = v->val2;
    switch (which) {
      case VERTEX_COORDS:
        MRISsetOriginalXYZ(mris_dst, vno, v->origx, v->origy, v->origz);
        break;
      case VERTEX_ANNOTATION:
        vdst->annotation = v->annotation;
        break;
      case VERTEX_AREA:
        vdst->origarea = v->origarea;
        break;
      case VERTEX_CURV:
        vdst->curv = v->curv;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->val = v->val;
        break;
    }
  }
  if (mht_src) {
    MHTfree(&mht_src);
  }

  return (NO_ERROR);
}
#else
int MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, MRIS_HASH_TABLE *mht, int which)
{
  int vno;
  VERTEX *v, *vdst;
  MHT *mht_src = NULL;
  double max_len, mean;

  MRISclearMarks(mris_dst);
  MRISclearWhichAndVal2(mris_dst, which);
  MRISclearMarks(mris_src);
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    vdst = &mris_dst->vertices[vno];
    vdst->d = 0;
    vdst->val2 = 0;
  }

  /*
  First determine how much of a fan in there is in the mapping.
  */

  /*
  go through each vertex in the source and see what destination
  vertex it maps to.
  */
  for (vno = 0; vno < mris_src->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    v = &mris_src->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    vdst = MHTfindClosestVertex(mht, mris_dst, v);
    if (!vdst) {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map vno %d", vno);
      continue;
    }
    if (vdst - mris_dst->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    v->marked++;
  }

  /* sample from dst to source  */
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vdst = &mris_dst->vertices[vno];
    if (vdst->marked || vdst->ripflag) {
      continue;
    }
    if (!mht_src) {
      MRIScomputeVertexSpacingStats(mris_src, NULL, NULL, &max_len, NULL, NULL, CURRENT_VERTICES);
      mht_src = MHTcreateVertexTable_Resolution(mris_src, CURRENT_VERTICES, 2 * max_len);
    }
    v = MHTfindClosestVertex(mht_src, mris_src, vdst);
    if (!v) {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map dst vno %d", vno);
      continue;
    }
    if (v - mris_src->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    v->marked++;
  }

  MRISclearMarks(mris_dst);
  /*
  go through each vertex in the source and sample it onto
  the destination surface.
  */
  for (vno = 0; vno < mris_src->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    v = &mris_src->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    vdst = MHTfindClosestVertex(mht, mris_dst, v);
    if (!vdst) {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map vno %d", vno);
      continue;
    }
    if (vdst - mris_dst->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    vdst->val2 += v->val2 / (float)v->marked; /* variances */
    switch (which) {
      case VERTEX_AREA:
        vdst->d += v->origarea / (float)v->marked;
        break;
      case VERTEX_CURV:
        vdst->d += v->curv;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->d += v->val;
        break;
    }
  }

  /* normalize by # of times a vertex is mapped to */
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vdst = &mris_dst->vertices[vno];
    if (vdst->ripflag || !vdst->marked) {
      continue;
    }
    vdst->val2 /= (float)vdst->marked;
    mean = vdst->d / (float)vdst->marked;
    switch (which) {
      case VERTEX_AREA:
        vdst->origarea = mean;
        break;
      case VERTEX_CURV:
        vdst->curv = mean;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->val = mean;
        break;
    }
  }

  /* sample from dst to source to fill holes */
  for (vno = 0; vno < mris_dst->nvertices; vno++) {
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vdst = &mris_dst->vertices[vno];
    if (vdst->marked || vdst->ripflag) {
      continue;
    }
    v = MHTfindClosestVertex(mht_src, mris_src, vdst);
    if (!v) {
      ErrorPrintf(ERROR_BADPARM, "MRISsphericalCopy: cannot map dst vno %d", vno);
      continue;
    }
    if (v - mris_src->vertices == Gdiag_no) {
      DiagBreak();
    }
    vdst->marked++;
    vdst->val2 = v->val2 / (float)v->marked;
    switch (which) {
      case VERTEX_AREA:
        vdst->origarea = v->origarea / (float)v->marked;
        break;
      case VERTEX_CURV:
        vdst->curv = v->curv;
        break;
      case VERTEX_LOGODDS:
      case VERTEX_VALS:
        vdst->val = v->val;
        break;
    }
  }
  if (mht_src) {
    MHTfree(&mht_src);
  }

  return (NO_ERROR);
}
#endif

MRIS *MRISremoveRippedSurfaceElements(MRIS *mris)
{
  int *vertex_trans, *face_trans;
  int vno, fno, i, n, nrippedfaces, nrippedvertices, kept_vertices, kept_faces;

  FACE *f, *fdst;

  fprintf(WHICH_OUTPUT, "building final representation...\n");

  vertex_trans = (int *)malloc(mris->nvertices * sizeof(int));
  face_trans = (int *)malloc(mris->nfaces * sizeof(int));

  memset(vertex_trans, -1, mris->nvertices * sizeof(int));
  memset(face_trans, -1, mris->nfaces * sizeof(int));

  // cout the number of faces and vertices
  for (kept_vertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    kept_vertices++;
  }
  for (kept_faces = fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    kept_faces++;
  }

  // create a new surface
  MRIS *mris_corrected ; //= MRISalloc(kept_vertices, kept_faces);
  mris_corrected = MRISoverAlloc(nint(1.5*kept_vertices), nint(1.5*kept_faces), kept_vertices, kept_faces);
  // keep the extra info into the new one
  mris_corrected->useRealRAS = mris->useRealRAS;
  copyVolGeom(&mris->vg, &mris_corrected->vg);

  mris_corrected->type           = MRIS_TRIANGULAR_SURFACE;
  mris_corrected->status         = mris->status;                // this should have been set, but wasn't
  mris_corrected->origxyz_status = mris->origxyz_status;        // this is new

  int newNVertices = 0;
  for (nrippedvertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX          const * const v  = &mris->vertices         [vno];
    /* ignore the ripped vertices */
    if (v->ripflag) {
      nrippedvertices++;
      continue;
    }
    /* save vertex information */
    // TEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[newNVertices];
    VERTEX          * const vdst  = &mris_corrected->vertices         [newNVertices];

    vdst->x = v->x;
    vdst->y = v->y;
    vdst->z = v->z;
    
    MRISsetOriginalXYZ(mris_corrected, newNVertices, v->origx, v->origy, v->origz);
    
    vdst->tx = v->tx;
    vdst->ty = v->ty;
    vdst->tz = v->tz;
    vdst->nx = v->nx;
    vdst->ny = v->ny;
    vdst->nz = v->nz;
    vdst->cx = v->cx;
    vdst->cy = v->cy;
    vdst->cz = v->cz;
    vdst->val = v->val;
    vdst->val2 = v->val2;
    vdst->valbak = v->valbak;
    vdst->val2bak = v->val2bak;
    vdst->imag_val = v->imag_val;
    vdst->curv = v->curv;
    vdst->curvbak = v->curvbak;
    vdst->stat = v->stat;
    vdst->mean = v->mean;
    vdst->mean_imag = v->mean_imag;
    vdst->std_error = v->std_error;
    vdst->H = v->H;
    vdst->K = v->K;
    vdst->k1 = v->k1;
    vdst->k2 = v->k2;
    vertex_trans[vno] = newNVertices++;
  }
  MRIStruncateNVertices(mris_corrected, newNVertices);

  int newNFaces = 0;
  for (nrippedfaces = newNFaces = fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    /* don't update triangle with marked vertices */
    if (f->ripflag) {
      nrippedfaces++;
      continue;
    }
    /* save face information */
    fdst = &mris_corrected->faces[newNFaces];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fdst->v[n] = vertex_trans[f->v[n]];
      if (mris->vertices[f->v[n]].ripflag) {
        fprintf(stderr,
                "Error with face %d (%d): vertex %d (%d) is ripped\n",
                fno,
                newNFaces,
                f->v[n],
                fdst->v[n]);
      }
    }
    face_trans[fno] = newNFaces++;
  }
  MRIStruncateNFaces(mris_corrected, newNFaces);

  /* now allocate face and neighbor stuff in mris_corrected */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    int const vno_dst = vertex_trans[vno];
    VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[vno_dst];

    /* count # of good triangles attached to this vertex */
    for (vdstt->num = n = 0; n < vt->num; n++)
      if (mris->faces[vt->f[n]].ripflag == 0) {
        vdstt->num++;
      }
    vdstt->f = (int *)calloc(vdstt->num, sizeof(int));
    vdstt->n = (uchar *)calloc(vdstt->num, sizeof(uchar));
    for (i = n = 0; n < vt->num; n++) {
      if (mris->faces[vt->f[n]].ripflag) {
        continue;
      }
      vdstt->n[i] = vt->n[n];
      vdstt->f[i] = face_trans[vt->f[n]];
      i++;
    }
    /* count # of valid neighbors */
    clearVnum(mris_corrected,vno_dst);
    for (n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].ripflag == 0) {
        addVnum(mris_corrected,vno_dst,1);
      }

    vdstt->nsizeMax = 1;
    vdstt->v = (int *)calloc(vdstt->vnum, sizeof(int));
    for (i = n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].ripflag == 0) {
        vdstt->v[i++] = vertex_trans[vt->v[n]];
      }
      
    MRIS_setNsizeCur(mris_corrected, vno_dst, 1);
  }

  free(vertex_trans);
  free(face_trans);

  fprintf(stderr, "%d vertices and %d faces have been removed from triangulation\n", nrippedvertices, nrippedfaces);

  mrisCheckVertexFaceTopology(mris_corrected);
  
  return mris_corrected;
}

static int mrisGraphLaplacian(MRI_SURFACE *mris, int vno, double *pdx, double *pdy, double *pdz, int which)
{
  double lx, ly, lz, w, norm, dx, dy, dz, wtotal;
  int n;

  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  if (vt->vnum == 0) {
    return (NO_ERROR);
  }
  w = 1.0 / (double)vt->vnum;
  for (wtotal = 0.0, lx = ly = lz = 0.0, n = 0; n < vt->vnum; n++) {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];
    dx = (v->x - vn->x);
    dy = (v->y - vn->y);
    dz = (v->z - vn->z);
    switch (which) {
      default:
      case TAUBIN_UNIFORM_WEIGHTS:
        // w = 1/v->vnum above
        break;
      case TAUBIN_INVERSE_WEIGHTS:
        norm = sqrt(dx * dx + dy * dy + dz * dz);

        if (!FZERO(norm)) {
          w = 1.0 / norm;
        }
        else {
          w = 1;
        }
        break;
      case TAUBIN_EDGE_WEIGHTS:
        norm = sqrt(dx * dx + dy * dy + dz * dz);
        w = norm;
        break;
    }

    wtotal += w;
    lx += w * dx;
    ly += w * dy;
    lz += w * dz;
  }
  *pdx = lx / wtotal;
  *pdy = ly / wtotal;
  *pdz = lz / wtotal;

  return (NO_ERROR);
}

/* bigger values of lambda will give more smoothing. Lambda and Mu
   bigger K_bp will also yield more smoothing
#define Lambda .3
#define Mu     (1.0)/((K_bp)-1.0/Lambda)
*/
int MRIStaubinSmooth(MRI_SURFACE *mris, int niters, double lambda, double mu, int which)
{
  int n, vno;
  double dx, dy, dz;
  VERTEX *v;

  if (lambda < 0 || lambda > 1 || mu >= -lambda) {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MRIStaubinSmooth: mu < -lambda < 0 violated"));
  }

  dx = dy = dz = 0;  // to get rid of compiler warning
  for (n = 0; n < niters; n++) {
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (v->ripflag) {
        v->tx = v->x;
        v->ty = v->y;
        v->tz = v->z;
        continue;
      }
      mrisGraphLaplacian(mris, vno, &dx, &dy, &dz, which);
      if (EVEN(n)) {
        v->tx = v->x + lambda * dx;
        v->ty = v->y + lambda * dy;
        v->tz = v->z + lambda * dz;
      }
      else {
        v->tx = v->x + mu * dx;
        v->ty = v->y + mu * dy;
        v->tz = v->z + mu * dz;
      }
    }
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
  }

  return (NO_ERROR);
}


/*-------------------------------------------------------------------
  MRISgaussianSmooth() - perform gaussian smoothing on a spherical
  surface. The gaussian is defined by stddev GStd and is truncated
  at TruncFactor stddevs. Note: this will change the val2bak of all
  the vertices. See also MRISspatialFilter() and MRISgaussianWeights().
  -------------------------------------------------------------------*/
MRI *MRISgaussianSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ, double TruncFactor)
{
  int vtxno1, vtxno2;
  float val;
  MRI *SrcTmp, *GSum, *GSum2, *nXNbrsMRI;
  VERTEX *vtx1;
  double Radius, Radius2, dmax, GVar2, f, d, costheta, theta, g, dotprod;
  int n, err, nXNbrs, *XNbrVtxNo, frame;
  double *XNbrDotProd, DotProdThresh;
  double InterVertexDistAvg, InterVertexDistStdDev;
  double VertexRadiusAvg, VertexRadiusStdDev;

  if (Surf->nvertices != Src->width) {
    printf("ERROR: MRISgaussianSmooth: Surf/Src dimension mismatch\n");
    return (NULL);
  }

  if (Targ == NULL) {
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, Src->nframes);
    if (Targ == NULL) {
      printf("ERROR: MRISgaussianSmooth: could not alloc\n");
      return (NULL);
    }
  }
  else {
    if (Src->width != Targ->width || Src->height != Targ->height || Src->depth != Targ->depth ||
        Src->nframes != Targ->nframes) {
      printf("ERROR: MRISgaussianSmooth: output dimension mismatch\n");
      return (NULL);
    }
    if (Targ->type != MRI_FLOAT) {
      printf("ERROR: MRISgaussianSmooth: structure passed is not MRI_FLOAT\n");
      return (NULL);
    }
  }

  /* Make a copy in case it's done in place */
  SrcTmp = MRIcopy(Src, NULL);

  /* This is for normalizing */
  GSum = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if (GSum == NULL) {
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum\n");
    return (NULL);
  }

  GSum2 = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);
  if (GSum2 == NULL) {
    printf("ERROR: MRISgaussianSmooth: could not alloc GSum2\n");
    return (NULL);
  }

  MRIScomputeMetricProperties(Surf);
  nXNbrsMRI = MRIallocSequence(Src->width, Src->height, Src->depth, MRI_FLOAT, 1);

  vtx1 = &Surf->vertices[0];
  Radius2 = (vtx1->x * vtx1->x) + (vtx1->y * vtx1->y) + (vtx1->z * vtx1->z);
  Radius = sqrt(Radius2);
  dmax = TruncFactor * GStd;  // truncate after TruncFactor stddevs
  GVar2 = 2 * (GStd * GStd);
  f = pow(1 / (sqrt(2 * M_PI) * GStd), 2.0);  // squared for 2D
  DotProdThresh = Radius2 * cos(dmax / Radius) * (1.0001);

  printf(
      "Radius = %g, gstd = %g, dmax = %g, GVar2 = %g, f = %g, dpt = %g\n", Radius, GStd, dmax, GVar2, f, DotProdThresh);

  InterVertexDistAvg = Surf->avg_vertex_dist;
  InterVertexDistStdDev = Surf->std_vertex_dist;
  VertexRadiusAvg = MRISavgVetexRadius(Surf, &VertexRadiusStdDev);

  printf("Total Area = %g \n", Surf->total_area);
  printf("Dist   = %g +/- %g\n", InterVertexDistAvg, InterVertexDistStdDev);
  printf("Radius = %g +/- %g\n", VertexRadiusAvg, VertexRadiusStdDev);

  /* Initialize */
  for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
    MRIFseq_vox(GSum, vtxno1, 0, 0, 0) = 0;
    MRIFseq_vox(GSum2, vtxno1, 0, 0, 0) = 0;
    for (frame = 0; frame < Targ->nframes; frame++) {
      MRIFseq_vox(Targ, vtxno1, 0, 0, frame) = 0;
    }
    Surf->vertices[vtxno1].val2bak = -1;
  }

  /* These are needed by MRISextendedNeighbors()*/
  XNbrVtxNo = (int *)calloc(Surf->nvertices, sizeof(int));
  XNbrDotProd = (double *)calloc(Surf->nvertices, sizeof(double));

  printf("nvertices = %d\n", Surf->nvertices);
  for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
    nXNbrs = 0;
    err =
        MRISextendedNeighbors(Surf, vtxno1, vtxno1, DotProdThresh, XNbrVtxNo, XNbrDotProd, &nXNbrs, Surf->nvertices, 1);
    MRIFseq_vox(nXNbrsMRI, vtxno1, 0, 0, 0) = nXNbrs;

    if (vtxno1 % 10000 == 0 && Gdiag_no > 0) {
      printf("vtxno1 = %d, nXNbrs = %d\n", vtxno1, nXNbrs);
      fflush(stdout);
    }

    for (n = 0; n < nXNbrs; n++) {
      vtxno2 = XNbrVtxNo[n];
      dotprod = XNbrDotProd[n];
      costheta = dotprod / Radius2;

      // cos theta might be slightly > 1 due to precision
      if (costheta > +1.0) {
        costheta = +1.0;
      }
      if (costheta < -1.0) {
        costheta = -1.0;
      }

      // Compute the angle between the vertices
      theta = acos(costheta);

      /* Compute the distance bet vertices along the surface of the sphere */
      d = Radius * theta;

      /* Compute weighting factor for this distance */
      g = f * exp(-(d * d) / (GVar2)); /* f not really nec */
      // ga = g * Surf->vertices[vtxno2].area;

      if (vtxno1 == 10000 && 1) {
        printf("%d %d %g %g %g %g %g\n", vtxno1, vtxno2, dotprod, costheta, theta, d, g);
        fflush(stdout);
      }

      // MRIFseq_vox(GSum,vtxno1,0,0,0)  += g;
      // MRIFseq_vox(GSum2,vtxno1,0,0,0) += (g*g);

      for (frame = 0; frame < Targ->nframes; frame++) {
        val = g * MRIFseq_vox(SrcTmp, vtxno2, 0, 0, frame);
        MRIFseq_vox(Targ, vtxno1, 0, 0, frame) += val;
      }

    } /* end loop over vertex2 */

  } /* end loop over vertex1 */

  /* Normalize */
  if (0) {
    for (vtxno1 = 0; vtxno1 < Surf->nvertices; vtxno1++) {
      vtx1 = &Surf->vertices[vtxno1];
      g = MRIFseq_vox(GSum, vtxno1, 0, 0, 0);
      MRIFseq_vox(GSum2, vtxno1, 0, 0, 0) /= (g * g);

      for (frame = 0; frame < Targ->nframes; frame++) {
        val = MRIFseq_vox(Targ, vtxno1, 0, 0, frame);
        MRIFseq_vox(Targ, vtxno1, 0, 0, frame) = val / g;
      }
    }
  }

  // MRIwrite(GSum,"gsum.mgh");
  // MRIwrite(GSum2,"gsum2.mgh");
  // MRIwrite(nXNbrsMRI,"nxnbrs.mgh");

  MRIfree(&SrcTmp);
  MRIfree(&GSum);
  MRIfree(&GSum2);
  MRIfree(&nXNbrsMRI);

  free(XNbrVtxNo);
  free(XNbrDotProd);

  return (Targ);
}

int MRISupsampleIco(MRI_SURFACE *mris, MRI_SURFACE *mris_new)
{
  MRISclearMarks(mris_new);

  mris_new->origxyz_status = mris->origxyz_status;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const vold = &mris->vertices[vno];
    VERTEX * const vnew = &mris_new->vertices[vno];

    vnew->x = vold->x;
    vnew->y = vold->y;
    vnew->z = vold->z;
    
    MRISsetOriginalXYZ(mris_new, vno, vold->origx, vold->origy, vold->origz);
    
    vnew->marked = 1;
  }

  MRISsoapBubbleVertexPositions(mris_new, 100);
  MRISsoapBubbleOrigVertexPositions(mris_new, 100);
  copyVolGeom(&mris->vg, &mris_new->vg);
  return (NO_ERROR);
}

#define MAX_INT_REMOVAL_NEIGHBORS 5
int MRISremoveIntersections(MRI_SURFACE *mris, int FillHoles)
{
  int n, num, writeit = 0, old_num, nbrs, min_int, no_progress = 0;

  n = 0;

  num = mrisMarkIntersections(mris,FillHoles);
  if (num == 0) return (NO_ERROR);
  printf("removing intersecting faces\n");
  min_int = old_num = mris->nvertices;
  MRISsaveVertexPositions(mris, TMP2_VERTICES);
  nbrs = 0;
  MRISclearMarks(mris);
  num = mrisMarkIntersections(mris,FillHoles);
  while (num > 0) {
    if ((num > old_num) || ((num == old_num) && (no_progress >= nbrs)))  // couldn't remove any
    {
      no_progress++;
      printf("step %d with no progress (num=%d, old_num=%d)\n", no_progress, num, old_num);

      // couldn't make any more progress with current size of neighborhood, expand, reset or quit
      if (nbrs >= MAX_INT_REMOVAL_NEIGHBORS)  // don't let neighborhood get too big
      {
        if (no_progress >= nbrs)  // couldnt' find any better configurations - give up
          break;
        printf("max nbrs reached, resetting neighborhood size\n");
        nbrs = 0;
      }
      else {
        // disable this code as it was causing the surface to do wacky things. Just try smoothing
        // vertices that are actually intersecting
        //    nbrs++ ;
        //    printf("expanding nbhd size to %d\n", nbrs);
        if (no_progress > 15) break;
      }
    }
    else  // num didn't get bigger
    {
      if (num < old_num)  // num actually decreased
        no_progress = 0;
      else
        no_progress++;
      if (num < min_int) {
        min_int = num;
        MRISsaveVertexPositions(mris, TMP2_VERTICES);
      }
    }

    MRISdilateMarked(mris, nbrs);
    old_num = num;

    printf("%03d: %d intersecting\n", n, num);
    MRISnotMarked(mris);  // turn off->on and on->off so soap bubble is correct (marked are fixed)
    MRISsoapBubbleVertexPositions(mris, 100);
    if (writeit) {
      char fname[STRLEN];
      sprintf(fname, "%s.avg%03d", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", n);
      MRISwrite(mris, fname);
    }
    if (n++ > 100)  // don't let it go forever
      break;

    MRISclearMarks(mris);
    num = mrisMarkIntersections(mris,FillHoles);
  }

  if (num > min_int) {
    MRISrestoreVertexPositions(mris, TMP2_VERTICES);  // the least number of negative vertices
    num = mrisMarkIntersections(mris,0);
  }
  printf("terminating search with %d intersecting\n", num);
  return (NO_ERROR);
}


int MRISaverageMarkedVertexPositions(MRI_SURFACE *mris, int navgs)
{
  int i, vno, vnb, vnum;
  float x, y, z, num;
  int nmarked;

  for (i = 0; i < navgs; i++) {
    for (nmarked = vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      x = y = z = 0;
      num = 0;
      x = v->x;
      y = v->y;
      z = v->z;
      num++; /* account for central vertex */
      int const *pnb = vt->v;
      vnum = vt->vnum;
      for (vnb = 0; vnb < vnum; vnb++) {
        VERTEX const * const vn = &mris->vertices[*pnb++]; /* neighboring vertex pointer */
        if (vn->ripflag)              /* no valid data */
        {
          continue;
        }
        num++;
        x += vn->x;
        y += vn->y;
        z += vn->z;
      }
      v->tdx = x / num;
      v->tdy = y / num;
      v->tdz = z / num;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v  = &mris->vertices[vno];
      if (v->ripflag || v->marked == 0) {
        continue;
      }
      v->x = v->tdx;
      v->y = v->tdy;
      v->z = v->tdz;
    }
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISripMedialWall(MRI_SURFACE *mris)
{
  int vno, med_index, unknown_index;
  VERTEX *v;
  int structure;

  printf("ripping medial wall...\n");
  CTABfindName(mris->ct, "Unknown", &unknown_index);
  CTABfindName(mris->ct, "Medial_wall", &med_index);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    CTABfindAnnotation(mris->ct, v->annotation, &structure);
    if (structure == unknown_index || structure == med_index) {
      v->ripflag = 1;
    }
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISzeroMedialWallCurvature(MRI_SURFACE *mris)
{
  int vno, med_index, unknown_index;
  VERTEX *v;
  int structure;

  printf("erasing medial wall curvatures...\n");
  CTABfindName(mris->ct, "Unknown", &unknown_index);
  CTABfindName(mris->ct, "Medial_wall", &med_index);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    CTABfindAnnotation(mris->ct, v->annotation, &structure);
    if (structure == unknown_index || structure == med_index) {
      v->curv = 0;
    }
  }
  MRISsetRipInFacesWithRippedVertices(mris);
  return (NO_ERROR);
}


void UpdateMRIS(MRI_SURFACE *mris, const char *fname)
{
  mrisCompleteTopology(mris);
  mrisComputeVertexDistances(mris);
  if (fname) {
    mrisReadTransform(mris, fname);
  }
  mris->radius = MRISaverageRadius(mris);
  MRIScomputeMetricProperties(mris);
  MRISstoreCurrentPositions(mris);
}


int MRISwriteCoordsToIco(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, int which_vertices)
{
  int vno, fno, debug;
  VERTEX *v;
  FACE *face;
  double fdist;
  float x, y, z;

  MHT* mht = MHTcreateFaceTable_Resolution(mris, CANONICAL_VERTICES, 1.0);

  for (vno = 0; vno < mris_ico->nvertices; vno++) {
    v = &mris_ico->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }

    MHTfindClosestFaceGeneric(mht, mris, v->cx, v->cy, v->cz, 1000, -1, 1, &face, &fno, &fdist);
    if (face == NULL) {
      DiagBreak();
      continue;
    }
    debug = (face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no);
    MRISsampleFaceCoordsCanonical(mht, mris, v->cx, v->cy, v->cz, which_vertices, &x, &y, &z);
    if (debug) {
      printf("ico vertex %d, coords = (%2.2f, %2.2f, %2.2f)\n", vno, x, y, z);
    }
    switch (which_vertices) {
      case PIAL_VERTICES:
        v->pialx = x;
        v->pialy = y;
        v->pialz = z;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "%s: unsupported vertex set %d", __MYFUNCTION__, which_vertices);
        break;
    }
  }
  MHTfree(&mht);
  return (NO_ERROR);
}


int MRISrepositionSurface(
    MRI_SURFACE *mris, MRI *mri, int *target_vnos, float *target_vals, int nv, int nsize, double sigma, int flags)
{
  int vno, n;
  INTEGRATION_PARMS parms;
  VERTEX *v;

  printf("flags = %x, size = %ld\n", flags, (long)sizeof(flags));

  parms.fill_interior = 0;
  parms.projection = NO_PROJECTION;
  parms.tol = 1e-4;
  parms.dt = 0.5f;
  parms.base_dt = .1;
  parms.l_spring = 1.0f;
  parms.l_curv = 1.0;
  parms.l_intensity = 0.2;
  parms.l_spring = 0.0f;
  parms.l_curv = 1.0;
  parms.l_intensity = 0.2;
  parms.l_tspring = 1.0f;
  parms.l_nspring = 0.5f;

  parms.flags |= flags;
  if (flags & IPFLAG_FORCE_GRADIENT_IN) {
    parms.grad_dir = -1;
    printf("forcing gradient in\n");
  }
  else if (flags & IPFLAG_FORCE_GRADIENT_IN) {
    parms.grad_dir = 1;
    printf("forcing gradient out\n");
  }
  parms.niterations = 0;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM;
  parms.momentum = 0.0 /*0.8*/;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.l_surf_repulse = 0.0;
  parms.l_repulse = 5;
  parms.niterations = 100;
  sprintf(parms.base_name, "nudge");

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].ripflag = 1;
  }
  for (n = 0; n < nv; n++) {
    mris->vertices[target_vnos[n]].ripflag = 0;
  }
  MRISerodeRipped(mris, nsize);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->val = target_vals[0];
    v->val2 = sigma;
    v->marked = 1;
  }

  MRISpositionSurface(mris, mri, mri, &parms);
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].ripflag = 0;
  }
  return (NO_ERROR);
}


int MRISrepositionSurfaceToCoordinate(
    MRI_SURFACE *mris, MRI *mri, int target_vno, float tx, float ty, float tz, int nsize, double sigma, int flags)
{
  int vno;
  INTEGRATION_PARMS parms;
  VERTEX *v;
  double xv, yv, zv, val;

#define SCALE .01

  printf("MRISrepositionSurfaceToCoordinate(%d, %f, %f, %f, %d, %f, %x)\n",
	 target_vno, tx, ty, tz, nsize, sigma, flags);
  parms.fill_interior = 0;
  parms.projection = NO_PROJECTION;
  parms.tol = 1e-4;
  parms.n_averages = sigma * sigma * M_PI / 2;
  parms.dt = 0.5f;
  parms.base_dt = .1;
  parms.l_spring = SCALE * 1.0f;
  parms.l_curv = SCALE * 1.0;
  parms.l_location = 10;
  parms.flags = flags;
  parms.l_spring = 0.0f;
  parms.l_intensity = 0 * SCALE * 0.2;
  parms.l_tspring = SCALE * .1f;
  parms.l_nspring = SCALE * 0.5f;

  parms.niterations = 0;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM;
  parms.momentum = 0.0 /*0.8*/;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
  parms.l_surf_repulse = 0.0;
  //  parms.l_repulse = SCALE*5 ;
  parms.niterations = 100;
  sprintf(parms.base_name, "nudge");

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].ripflag = 1;
  }
  mris->vertices[target_vno].ripflag = 0;
#if 0
  for (n = 0 ; n < nv ; n++)
  {
    mris->vertices[target_vnos[n]].ripflag = 0 ;
  }
#endif
  MRISerodeRipped(mris, nsize);
  v = &mris->vertices[target_vno];
  v->targx = tx;
  v->targy = ty;
  v->targz = tz;
  v->val2 = sigma;
  MRISsurfaceRASToVoxelCached(mris, mri, tx, ty, tz, &xv, &yv, &zv);
  MRIsampleVolume(mri, xv, yv, zv, &val);
  printf("volume @ (%2.1f, %2.1f, %2.1f) = %2.1f\n", xv, yv, zv, val);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *v, *vn;

    vn = &mris->vertices[vno];
    if (vn->ripflag == 0) {
      vn->val = val;
      if (vno != target_vno)  // give it a target which maintains its relative position
      {
        v = &mris->vertices[target_vno];
        vn->targx = tx + (vn->x - v->x);
        vn->targy = ty + (vn->y - v->y);
        vn->targz = tz + (vn->z - v->z);
        printf("setting target for vertex %d to (%2.1f %2.1f %2.1f)\n", vno, vn->targx, vn->targy, vn->targz);
      }
    }
  }

  /*
  Gdiag |= DIAG_SHOW ;
  Gdiag_no = target_vno ;
  */

  MRISpositionSurface(mris, mri, mri, &parms);
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].ripflag = 0;
  }
  return (NO_ERROR);
}

static void MRIScopyOneVertex(MRI_SURFACE *out_mris, int const out_vno, MRI_SURFACE *inp_mris, int const inp_vno) {
    VERTEX_TOPOLOGY       * const out_vt = &out_mris->vertices_topology[out_vno];
    VERTEX_TOPOLOGY const * const inp_vt = &inp_mris->vertices_topology[inp_vno];
    VERTEX                * const out_v  = &out_mris->vertices         [out_vno];
    VERTEX          const * const inp_v  = &inp_mris->vertices         [inp_vno];
    memmove(out_v,  inp_v,  sizeof(VERTEX));
#ifdef SEPARATE_VERTEX_TOPOLOGY
    memmove(out_vt, inp_vt, sizeof(VERTEX_TOPOLOGY));
#endif
    out_vt->v = (int *)calloc(out_vt->vtotal, sizeof(int));
    if (out_vt->v == NULL) ErrorExit(ERROR_NOMEMORY, "MRISconcat: could not allocate %dth vertex array", out_vno);
    int n;
    for (n = 0; n < out_vt->vtotal; n++) out_vt->v[n] = inp_vt->v[n];
    out_vt->f = (int *)calloc(out_vt->num, sizeof(int));
    if (out_vt->f == NULL) ErrorExit(ERROR_NOMEMORY, "MRISconcat: could not allocate %dth face array", out_vno);
    for (n = 0; n < out_vt->vtotal; n++) out_vt->f[n] = inp_vt->f[n];
}

MRI_SURFACE *MRISconcat(MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris)
{
  int vno, n;
  int fno;
  FACE *f, *fo;

  if (mris == NULL) mris = MRISalloc(mris1->nvertices + mris2->nvertices, mris1->nfaces + mris2->nfaces);

  for (vno = 0; vno < mris1->nvertices; vno++) {
    MRIScopyOneVertex(mris, vno, mris1, vno);
  }

  for (fno = 0; fno < mris1->nfaces; fno++) {
    f  = &mris ->faces[fno];
    fo = &mris1->faces[fno];
    memmove(f, fo, sizeof(FACE));
  }
  for (vno = mris1->nvertices; vno < mris->nvertices; vno++) {
    MRIScopyOneVertex(mris, vno, mris2, vno - mris1->nvertices);
  }
  for (fno = mris1->nfaces; fno < mris->nfaces; fno++) {
    f  = &mris ->faces[fno];
    fo = &mris2->faces[fno - mris1->nfaces];
    for (n = 0; n < VERTICES_PER_FACE; n++) f->v[n] = fo->v[n] + mris1->nvertices;
  }

  if (mris1->hemisphere != mris2->hemisphere)
    mris->hemisphere = BOTH_HEMISPHERES;
  else
    mris->hemisphere = mris1->hemisphere;
  mris->type = mris1->type;
  mris->nsize = mris1->nsize;
  MRISsetNeighborhoodSizeAndDist(mris, mris->nsize);

  //memmove(&mris->vg, &mris1->vg, sizeof(mris1->vg));
  mris->vg = mris1->vg;
  MRIScomputeMetricProperties(mris);
  MRIScomputeSecondFundamentalForm(mris);
  strcpy(mris->fname, mris1->fname);

  mrisCheckVertexFaceTopology(mris);
  
  return (mris);
}

void MRISmapOntoSphere(MRIS *mris)
{
  MRISprojectOntoSphere(mris, mris, 100.0f);

  MRISsmoothOnSphere(mris, 10000);
}

/*!
  \fn int MRISfixAverageSurf7(MRIS *surf7)
  \brief This fixeds a problem with ico7 average surfaces as created
  by mris_make_average_surface where two vertices (0 and 40969) have
  the same xyz coordinates. The fix is to move vertex 40969 to a point
  half way between itself and its nearest neighbor. mris_make_average_surface
  has since been fixed.
*/
int MRISfixAverageSurf7(MRIS *surf7)
{

  if(surf7->nvertices != 163842){
    printf("ERROR: MRISfixAverageSurf7(): must be ico7\n");
    return(1);
  }

  VERTEX          const * const v0  = &(surf7->vertices         [    0]);
  VERTEX_TOPOLOGY const * const v1t = &(surf7->vertices_topology[40969]);
  VERTEX                * const v1  = &(surf7->vertices         [40969]);

  // Make sure that the xyz are the same at these two vertices
  {
    double const dx = (v0->x - v1->x);
    double const dy = (v0->y - v1->y);
    double const dz = (v0->z - v1->z);
    double const d = sqrt(dx*dx + dy*dy + dz*dz);
    if(d > .001){
      printf("INFO: MRISfixAverageSurf7(): this surface appears to have been fixed already\n");
      printf(" v0 (%g,%g,%g) v40969 (%g,%g,%g), d = %g\n",v0->x,v0->y,v0->z,v1->x,v1->y,v1->z,d);
      return(1);
    }
  }

  // Go through the neighbors of 40969 and find the closest vertex. Note: exclude any vertex
  // with a distance of 0 (because this is what the problem is)
  double dmin = 10e10;
  double xmin = dmin, ymin = dmin, zmin = dmin;
  int n;
  for(n=0; n < v1t->num; n++){
    VERTEX const * const vn = &(surf7->vertices[v1t->v[n]]);	
    double const dx = (vn->x - v1->x);
    double const dy = (vn->y - v1->y);
    double const dz = (vn->z - v1->z);
    double const d = sqrt(dx*dx + dy*dy + dz*dz);
    if(d<dmin && d > .001){
      dmin = d;
      xmin = vn->x;
      ymin = vn->y;
      zmin = vn->z;
    }
  }

  // Move vertex 40969 half way between itself and its nearest neighbor
  v1->x = (v1->x + xmin)/2;
  v1->y = (v1->y + ymin)/2;
  v1->z = (v1->z + zmin)/2;

  return(0);
}


/*!
  \fn MRIS *MRISsortVertices(MRIS *mris0)
  \brief Creates a new surface with the vertices sorted by y,x,z; the
  faces are also sorted. The purpose of this function is to run after
  mris_decimate to make the output deterministic. mris_decimate will
  always produce the same surface xyz, but the vertices and faces may
  be sorted differently, which is pretty annoying. The input surface
  is cloned, but one should not trust anything in the vertex strcture
  except the xyz and neighboring vertices. The surface is actually
  written to disk and read back in to clean out all the elements.
  There may still be some non-deterministic behavior. For example,
  when the orig.nofix is input, mris_decimate can produce different
  surfaces (same number of vertices, but not in the same places). Not
  sure why, but probably because the lengths of the edges are all
  either 1 or sqrt(2) thus creating some abiguity which is handled
  differently on different runs.
 */
MRIS *MRISsortVertices(MRIS *mris0)
{
  VERTEX_SORT * vtxsort = (VERTEX_SORT *) calloc(mris0->nvertices,sizeof(VERTEX_SORT));
  int nthvtx;
  for(nthvtx = 0; nthvtx < mris0->nvertices; nthvtx++){
    vtxsort[nthvtx].vtxno = nthvtx;
    vtxsort[nthvtx].x = mris0->vertices[nthvtx].x;
    vtxsort[nthvtx].y = mris0->vertices[nthvtx].y;
    vtxsort[nthvtx].z = mris0->vertices[nthvtx].z;
  }
  qsort(vtxsort, mris0->nvertices, sizeof(VERTEX_SORT), CompareVertexCoords);

  // Create a LUT to quickly map from old to new
  int* const lut = (int *) calloc(mris0->nvertices,sizeof(int));
  for(nthvtx = 0; nthvtx < mris0->nvertices; nthvtx++){
    lut[vtxsort[nthvtx].vtxno] = nthvtx;
  }

  // Create a new surface. Cloning is tricky here because it will
  // copy all the elements from the vertices and faces. But when the
  // vertex identity is reassigned, most of the elements will be out
  // of synch. So below, the surface is written out and read back in
  // to clean it. When writing, only xyz and face vertices are kept.
  MRIS* mris = MRISclone(mris0);

  // Change the vertex xyz and neighbors
  for(nthvtx = 0; nthvtx < mris0->nvertices; nthvtx++){
    VERTEX_TOPOLOGY const * const vtxtold = &(mris0->vertices_topology[vtxsort[nthvtx].vtxno]);
    VERTEX          const * const vtxold  = &(mris0->vertices         [vtxsort[nthvtx].vtxno]);
    VERTEX_TOPOLOGY       * const vtxtnew = &(mris ->vertices_topology        [nthvtx]);
    VERTEX                * const vtxnew  = &(mris ->vertices                 [nthvtx]);
    vtxnew->x = vtxold->x;
    vtxnew->y = vtxold->y;
    vtxnew->z = vtxold->z;
    modVnum(mris,nthvtx,vtxtold->vnum,true); // number of neighboring vertices
    // Now copy the neighbors
    if(vtxtnew->v) free(vtxtnew->v);
    vtxtnew->v = (int *)calloc(vtxtnew->vnum, sizeof(int));
    int n;
    for(n=0; n < vtxtold->vnum; n++){
      int vno_old = vtxtold->v[n];
      vtxtnew->v[n] = lut[vno_old];
    }
    // The order of the neighbors appears to be the deterministic
    // after decimation but this is added to make sure.
    qsort(vtxtnew->v, vtxtnew->vnum, sizeof(int), compare_ints);
  }
  freeAndNULL(vtxsort);

  // Change the face vertex numbers to match the new vertex
  // order. Also change the order of the vertices at each face. The
  // face order is changed in such a way that the minimum vertex
  // number is first but the vertex order is same. Eg, if the original
  // vertex order was 2-1-3 the new order would be 1-3-2 (NOT
  // 1-2-3). This is necessary to keep the direction of the normal
  // correct.
  int nthface;
  for(nthface = 0; nthface < mris0->nfaces; nthface++){
    FACE const * const f = &(mris0->faces[nthface]); // old and new faces are the same
    int vno_new_min = 0;
    int n_vno_new_min = 0;
    int n;
    int vnolist[3];
    for(n=0; n < 3; n++){
      int vno_old = f->v[n];
      int vno_new = lut[vno_old];
      vnolist[n] = vno_new;
      if(vno_new_min > vno_new){
	vno_new_min = vno_new;
	n_vno_new_min = n;
      }
    }
    // Now change to order so that the min vertex number is first.
    for(n=0; n < 3; n++){
      int m = (n_vno_new_min + n)%3;
      mris->faces[nthface].v[n] = vnolist[m];
    }    
  }
  free(lut);

  // Now sort the faces to be in deterministic order
  FACE_SORT * facesort = (FACE_SORT *) calloc(mris->nfaces,sizeof(FACE_SORT));
  for(nthface = 0; nthface < mris->nfaces; nthface++){
    facesort[nthface].faceno = nthface;
    facesort[nthface].v0 = mris->faces[nthface].v[0];
    facesort[nthface].v1 = mris->faces[nthface].v[1];
    facesort[nthface].v2 = mris->faces[nthface].v[2];
  }
  qsort(facesort, mris->nfaces, sizeof(FACE_SORT), CompareFaceVertices);
  for(nthface = 0; nthface < mris->nfaces; nthface++){
    FACE* const f = &(mris->faces[nthface]);
    f->v[0] = facesort[nthface].v0;
    f->v[1] = facesort[nthface].v1;
    f->v[2] = facesort[nthface].v2;
  }
  free(facesort);

  mrisCheckVertexFaceTopology(mris);

  // Create a temporary file to write the output and then read it back in
  // to get the decimated surface. The reason I do this is because there
  // was no clear API call in mrisurf.c that would allow me to muck with
  // the vertex/face array and properly calculate everything else in the
  // structure.  I felt this was the safest way to make sure everything
  // in the surface got recalculated properly.
  std::string tmpfile = makeTempFile();
  MRISwrite(mris, tmpfile.c_str());
  MRISfree(&mris);
  mris = MRISread(tmpfile.c_str());
  remove(tmpfile.c_str());

  return(mris);
}



int MRISrigidBodyAlignLocal(MRI_SURFACE *mris, INTEGRATION_PARMS *old_parms)
{
  int steps;
  INTEGRATION_PARMS parms;

  /* dx,dy,dz interpreted as rotations in applyGradient when status is rigid */
  auto const old_status = mris->status; /* okay, okay, this is a hack too... */
  mris->status = MRIS_RIGID_BODY;
  parms.integration_type = INTEGRATE_LM_SEARCH;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE;

  parms.mrisp_template = old_parms->mrisp_template;
  INTEGRATION_PARMS_copyFp(&parms, old_parms);
  parms.niterations = 25;
  parms.frame_no = old_parms->frame_no;
  parms.mrisp = old_parms->mrisp;
  parms.tol = old_parms->tol;
  parms.l_pcorr = 1.0f;
  parms.dt = old_parms->dt;
  /*  parms.integration_type = old_parms->integration_type ;*/
  parms.momentum = old_parms->momentum;
  parms.write_iterations = old_parms->write_iterations;
  parms.start_t = old_parms->start_t;
  strcpy(parms.base_name, old_parms->base_name);

  steps = MRISintegrate(mris, &parms, 0);
  old_parms->start_t += steps;
  mris->status = old_status;
  if (Gdiag & DIAG_WRITE)
    fprintf(old_parms->fp, "rigid alignment complete, sse = %2.3f\n", MRIScomputeSSE(mris, old_parms));
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    float area_rms, angle_rms, curv_rms, dist_rms, corr_rms, rms;

    rms = mrisComputeError(mris, &parms, &area_rms, &angle_rms, &curv_rms, &dist_rms, &corr_rms);
    fprintf(stdout, "rms = %2.3f, corr_rms = %2.3f ", rms, corr_rms);
    rms = mrisComputeError(mris, old_parms, &area_rms, &angle_rms, &curv_rms, &dist_rms, &corr_rms);
    fprintf(stdout, "(%2.3f, %2.3f)\n", rms, corr_rms);
  }
  return (NO_ERROR);
}

int MRISrigidBodyAlignVectorLocal(MRI_SURFACE *mris, INTEGRATION_PARMS *old_parms)
{
  int n, steps;
  INTEGRATION_PARMS parms;

  /* dx,dy,dz interpreted as rotations in applyGradient when status is rigid */
  auto const old_status = mris->status; /* okay, okay, this is a hack too... */
  mris->status = MRIS_RIGID_BODY;
  parms.integration_type = INTEGRATE_LM_SEARCH;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE;

  parms.mrisp_template = old_parms->mrisp_template;
  INTEGRATION_PARMS_copyFp(&parms, old_parms);
  parms.niterations = 25;
  parms.frame_no = old_parms->frame_no;
  parms.mrisp = old_parms->mrisp;
  parms.tol = old_parms->tol;

  parms.nfields = old_parms->nfields;
  parms.flags &= IP_USE_MULTIFRAMES;
  for (n = 0; n < MNOFIV; n++) {
    memmove(&parms.fields[n], &old_parms->fields[n], sizeof(FIELD_LABEL));
  }

  parms.l_pcorr = parms.l_corr = 0.0f;
  parms.dt = old_parms->dt;
  /*  parms.integration_type = old_parms->integration_type ;*/
  parms.momentum = old_parms->momentum;
  parms.write_iterations = old_parms->write_iterations;
  parms.start_t = old_parms->start_t;
  strcpy(parms.base_name, old_parms->base_name);

  steps = MRISintegrate(mris, &parms, 0);
  old_parms->start_t += steps;
  mris->status = old_status;
  if (Gdiag & DIAG_WRITE)
    fprintf(old_parms->fp, "rigid alignment complete, sse = %2.3f\n", MRIScomputeSSE(mris, old_parms));
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    float area_rms, angle_rms, curv_rms, dist_rms, corr_rms, rms;

    rms = mrisComputeError(mris, &parms, &area_rms, &angle_rms, &curv_rms, &dist_rms, &corr_rms);
    fprintf(stdout, "rms = %2.3f, corr_rms = %2.3f ", rms, corr_rms);
    rms = mrisComputeError(mris, old_parms, &area_rms, &angle_rms, &curv_rms, &dist_rms, &corr_rms);
    fprintf(stdout, "(%2.3f, %2.3f)\n", rms, corr_rms);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define STARTING_ANGLE RADIANS(16.0f)
#define ENDING_ANGLE RADIANS(4.0f)
#define NANGLES 8

int MRISrigidBodyAlignGlobal(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, float min_degrees, float max_degrees, int nangles)
{
  bool const tracing             = (Gdiag & DIAG_SHOW);
  bool const tracingWithSnapshot = tracing && (Gdiag & DIAG_WRITE);
    
  float const min_radians = RADIANS(min_degrees);
  float const max_radians = RADIANS(max_degrees);

  printf("Starting MRISrigidBodyAlignGlobal()\n");
  Timer mytimer;

  int const old_norm = parms->abs_norm;
  parms->abs_norm = 1;

  mrisOrientSurface(mris);

  auto const old_status = mris->status;
  mris->status = MRIS_RIGID_BODY;

  if (!parms->start_t) {
    mrisLogStatus(mris, parms, stdout, 0.0f, -1);
    if (tracingWithSnapshot) {
      if (parms->fp)
	mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
      if (parms->write_iterations > 0) {
        mrisWriteSnapshot(mris, parms, 0);
      }
    }
  }

  static bool 
    once,
    use_old,
    use_new;

  if (!once) { once = true;
    use_old = !!getenv("FREESURFER_MRISrigidBodyAlignGlobal_useOld");
    use_new = !!getenv("FREESURFER_MRISrigidBodyAlignGlobal_useNew") || !use_old ;
  }

  double new_mina = 666.0, new_minb = 666.0, new_ming = 666.0, new_sse = 666.0;

  //double const ext_sse = gMRISexternalSSE ? (*gMRISexternalSSE)(mris, parms) : 0.0;

  if (use_new) {
  
    printf("Starting new MRISrigidBodyAlignGlobal_findMinSSE()\n");
    Timer new_timer;
    
    // This does not modify either mris or params until after the old code has executed
    //
    double ext_sse = 0.0;  
    //double ext_sse = (gMRISexternalSSE) ? (*gMRISexternalSSE)(mris, parms) : 0.0;
    MRISrigidBodyAlignGlobal_findMinSSE(
        &new_mina, &new_minb, &new_ming, &new_sse,
        mris,
        parms,
        min_radians,
        max_radians,
        ext_sse,
        nangles); 

    parms->start_t += 1.0f;
    parms->t       += 1.0f;

    int msec = new_timer.milliseconds();
    printf("  new MRISrigidBodyAlignGlobal_findMinSSE"
      " min @ (%2.2f, %2.2f, %2.2f) sse = %2.1f, elapsed since starting=%6.4f min\n",
      (float)DEGREES(new_mina), (float)DEGREES(new_minb), (float)DEGREES(new_ming), new_sse,
      msec / (1000 * 60.0));
  }
  
  double old_mina = 666.0, old_minb = 666.0, old_ming = 666.0, old_sse = 666.0;
  if (use_old) {

    printf("Starting old MRISrigidBodyAlignGlobal_findMinSSE()\n");
    Timer old_timer;
 
    // Note: 
    //      This used to do a series of smaller and smaller rotations to mris_sphere, rotating the sphere to the minimum each time.
    //      This accumulated errors, and did an unnecessary rotation.
    //      Now it converges the three angles separately...
    // 
    MRISsaveVertexPositions(mris, TMP_VERTICES);

    double center_a = 0.0, center_b = 0.0, center_g = 0.0;

    double center_sse = -1.0;  // not known
    bool   center_sse_known = false;
    
    double radians;
    for (radians = max_radians; radians >= min_radians; radians /= 2.0f) {

      bool trace = tracing;
    
      if (tracing) {
        fprintf(stdout, 
          "scanning %2.2f degree nbhd of (%2.2f, %2.2f, %2.2f), min sse = %2.2f\n", 
          (float)DEGREES(radians),
          (float)DEGREES(center_a),  (float)DEGREES(center_b), (float)DEGREES(center_g),
          center_sse_known ? (float)(center_sse) : -666.0);
      }

      double const delta = 2 * radians / (float)nangles;
      double mina = center_a, minb = center_b, ming = center_g, min_sse = center_sse;
      bool min_known = center_sse_known;
      
      double alpha, beta, gamma;
#if 0
      // The oldest code went in this order
      // but the order should not matter
      // and the cache behavior and the new code both think that alpha should the innermost
      //
      for (alpha = center_a - radians; alpha <= center_a + radians; alpha += delta) {
        for (beta = center_b - radians; beta <= center_b + radians; beta += delta) {
          for (gamma = center_g - radians; gamma <= center_g + radians; gamma += delta) {
#else
      for (gamma = center_g - radians; gamma <= center_g + radians; gamma += delta) {
        for (beta = center_b - radians; beta <= center_b + radians; beta += delta) {
          for (alpha = center_a - radians; alpha <= center_a + radians; alpha += delta) {
#endif     
            VERTEX const * vertex0 = &mris->vertices[0];
            float vertex0_x = vertex0->x, vertex0_y = vertex0->y, vertex0_z = vertex0->z;
            MRISrotate(mris, mris, center_a + alpha, center_b + beta, center_g + gamma);
            
            if (trace) {
              fprintf(stdout, "%s:%d rotated (%g,%g,%g) by (a:%g, b:%g, g:%g) to (%g,%g,%g)\n", __FILE__, __LINE__, 
                vertex0_x,vertex0_y,vertex0_z, 
                center_a + alpha, center_b + beta, center_g + gamma,
                vertex0->x,vertex0->y,vertex0->z); 
            }
            
            double sse = mrisComputeCorrelationError(mris, parms, 1);
            if (gMRISexternalSSE)
            {
              double ext_sse = (*gMRISexternalSSE)(mris, parms) ;
              sse += ext_sse ;
            }
            if (trace) fprintf(stdout, "%s:%d sse:%g\n", __FILE__, __LINE__, sse);
            
            MRISrestoreVertexPositions(mris, TMP_VERTICES);
            
            if (!min_known || sse < min_sse) {
              min_known = true;
              mina    = center_a + alpha;
              minb    = center_b + beta;
              ming    = center_g + gamma;
              min_sse = sse;
            }
            
            trace = false;

            if (false && tracing) {
              fprintf(stdout, "\r  gamma "
                "min @ (%2.2f, %2.2f, %2.2f) sse:%2.1f   try @ (%+2.2f, %+2.2f, %+2.2f) sse:%2.2f",
                (float)DEGREES(mina),  (float)DEGREES(minb), (float)DEGREES(ming),   (float)(min_sse),
                (float)DEGREES(alpha), (float)DEGREES(beta), (float)DEGREES(gamma),  (float)(    sse));
              fflush(stdout);
            }
          }  // alpha
          if (false && tracing) {
            fprintf(stdout, "\r  beta "
              "min @ (%2.2f, %2.2f, %2.2f) sse:%2.1f   try @ (%+2.2f, %+2.2f, *)",
              (float)DEGREES(mina),  (float)DEGREES(minb), (float)DEGREES(ming),   (float)(min_sse),
              (float)DEGREES(alpha), (float)DEGREES(beta));
          }
        }    // beta
        if (false && tracing) {
          fprintf(stdout, "\n");
        }
      }      // gamma
    
      int msec = old_timer.milliseconds();
      printf("  d=%4.2f min @ (%2.2f, %2.2f, %2.2f) sse = %2.1f, elapsed since starting=%6.4f min\n",
        (float)DEGREES(radians), (float)DEGREES(mina), (float)DEGREES(minb), (float)DEGREES(ming),
        (float)(min_sse),
        msec / (1000 * 60.0));
      fflush(stdout);

      center_a = mina;
      center_b = minb;
      center_g = ming;
      center_sse = min_sse;
      center_sse_known = min_known;
      
      parms->start_t += 1.0f;
      parms->t       += 1.0f;

      if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
        mrisWriteSnapshot(mris, parms, parms->start_t);
      }
      if (Gdiag & DIAG_WRITE) {
        mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
      }
      if (Gdiag & DIAG_SHOW) {
        mrisLogStatus(mris, parms, stdout, 0.0f, -1);
      }
    }   // radians

    old_mina = center_a; old_minb = center_b, old_ming = center_g, old_sse = center_sse;

    int msec = old_timer.milliseconds();
    printf("  old MRISrigidBodyAlignGlobal_findMinSSE"
      " min @ (%2.2f, %2.2f, %2.2f) sse = %2.1f, elapsed since starting=%6.4f min\n",
      (float)DEGREES(old_mina), (float)DEGREES(old_minb), (float)DEGREES(old_ming), old_sse,
      msec / (1000 * 60.0));

  }  // use_old

  // Rotate to the best position
  //
  if (use_old) { 
    new_mina = old_mina; new_minb = old_minb; new_ming = old_ming; new_sse =  old_sse; 
    
    // TODO compare with new once the new is implemented
  }

  if (tracing) {
    fprintf(stdout,
      "rotating brain by (%2.2f, %2.2f, %2.2f), final sse: %2.2f\n",
      (float)DEGREES(new_mina),
      (float)DEGREES(new_minb),
      (float)DEGREES(new_ming),
      (float)(new_sse));
  }
  if (tracingWithSnapshot) {
    fprintf(parms->fp,
      "rotating brain by (%2.2f, %2.2f, %2.2f), final sse: %2.2f\n",
      (float)DEGREES(new_mina),
      (float)DEGREES(new_minb),
      (float)DEGREES(new_ming),
      (float)(new_sse));
  }
  
  MRISrotate(mris, mris, new_mina, new_minb, new_ming);

  mris->status    = old_status;
  parms->abs_norm = old_norm;

  int msec = mytimer.milliseconds();
  printf("MRISrigidBodyAlignGlobal() done %6.2f min\n", msec / (1000 * 60.0));

  return (NO_ERROR);
}

int MRISrigidBodyAlignVectorGlobal(
    MRI_SURFACE *mris, INTEGRATION_PARMS *parms, float min_degrees, float max_degrees, int nangles)
{
  double alpha, beta, gamma, degrees, delta, mina, minb, ming, sse, min_sse;
  auto const old_status = mris->status;

  min_degrees = RADIANS(min_degrees);
  max_degrees = RADIANS(max_degrees);
  mrisOrientSurface(mris);
  mris->status = MRIS_RIGID_BODY;
  if (!parms->start_t) {
    mrisLogStatus(mris, parms, stdout, 0.0f, -1);
    if (Gdiag & DIAG_WRITE) {
      mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
      if (parms->write_iterations > 0) {
        mrisWriteSnapshot(mris, parms, 0);
      }
    }
  }
  for (degrees = max_degrees; degrees >= min_degrees; degrees /= 2.0f) {
    mina = minb = ming = 0.0;
    min_sse = mrisComputeVectorCorrelationError(mris, parms, 1); /* was 0 !!!! */
    /*                if (gMRISexternalSSE) */
    /*                        min_sse += (*gMRISexternalSSE)(mris, parms) ; */
    delta = 2 * degrees / (float)nangles;
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "scanning %2.2f degree nbhd, min sse = %2.2f\n", (float)DEGREES(degrees), (float)min_sse);
    for (alpha = -degrees; alpha <= degrees; alpha += delta) {
      for (beta = -degrees; beta <= degrees; beta += delta) {
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout,
                  "\r(%+2.2f, %+2.2f, %+2.2f), "
                  "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                  (float)DEGREES(alpha),
                  (float)DEGREES(beta),
                  (float)DEGREES(-degrees),
                  (float)DEGREES(mina),
                  (float)DEGREES(minb),
                  (float)DEGREES(ming),
                  (float)min_sse);

        for (gamma = -degrees; gamma <= degrees; gamma += delta) {
          MRISsaveVertexPositions(mris, TMP_VERTICES);
          MRISrotate(mris, mris, alpha, beta, gamma);
          sse = mrisComputeVectorCorrelationError(mris, parms, 1); /* was 0 !!!! */
          /* if (gMRISexternalSSE) */
          /*   sse += (*gMRISexternalSSE)(mris, parms) ; */
          MRISrestoreVertexPositions(mris, TMP_VERTICES);
          if (sse < min_sse) {
            mina = alpha;
            minb = beta;
            ming = gamma;
            min_sse = sse;
          }
#if 0
          if (Gdiag & DIAG_SHOW)
            fprintf(stdout, "\r(%+2.2f, %+2.2f, %+2.2f), "
                    "min @ (%2.2f, %2.2f, %2.2f) = %2.1f   ",
                    (float)DEGREES(alpha), (float)DEGREES(beta), (float)
                    DEGREES(gamma), (float)DEGREES(mina),
                    (float)DEGREES(minb), (float)DEGREES(ming),(float)min_sse);
#endif
        }
      }
    }
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "\n");
    }
    if (!FZERO(mina) || !FZERO(minb) || !FZERO(ming)) {
      MRISrotate(mris, mris, mina, minb, ming);
      sse = mrisComputeVectorCorrelationError(mris, parms, 1); /* was 0 !!!! */
      /* if (gMRISexternalSSE) */
      /*   sse += (*gMRISexternalSSE)(mris, parms) ; */
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout,
                "min sse = %2.2f at (%2.2f, %2.2f, %2.2f)\n",
                sse,
                (float)DEGREES(mina),
                (float)DEGREES(minb),
                (float)DEGREES(ming));
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "rotating brain by (%2.2f, %2.2f, %2.2f), sse: %2.2f\n",
                (float)DEGREES(mina),
                (float)DEGREES(minb),
                (float)DEGREES(ming),
                (float)sse);
      parms->start_t += 1.0f;
      parms->t += 1.0f;
      if (Gdiag & DIAG_WRITE && parms->write_iterations > 0) {
        mrisWriteSnapshot(mris, parms, parms->start_t);
      }
      if (Gdiag & DIAG_WRITE) {
        mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
      }
      if (Gdiag & DIAG_SHOW) {
        mrisLogStatus(mris, parms, stdout, 0.0f, -1);
      }
    }
  }

  mris->status = old_status;
  return (NO_ERROR);
}


#if 0
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Given an intersection point int_pt in which a vertex intersects
  a face, compute the analagous location that the vertex should be
  placed in the faces original coordnates.

  note that fno refers to a face in mris and
  v is a vertex in mris_ico (not given), NOT mris.
  ------------------------------------------------------*/
static int mrisPlaceVertexInOrigFace(MRIS * const mris_vno, int const vno, MRIS * const mris, int const fno)
{
  VERTEX * const v = &mris_vno->vertices[vno];
  
  double U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], l1[3], l2[3], l2_len, l1_len, l_len, P[3], theta1, theta2, dot,
      theta_ratio, len_scale, e1[3], e2[3], etmp[3], x, y;
  int ret;

  e2[0] = 0.0;
  e2[1] = 0.0;
  e2[2] = 0.0;

  /* first compute point where normal to vertex intersects face */
  dir[0] = v->nx;
  dir[1] = v->ny;
  dir[2] = v->nz;
  pt[0] = v->x;
  pt[1] = v->y;
  pt[2] = v->z;
  load_triangle_vertices(mris, fno, U0, U1, U2, CURRENT_VERTICES);
  ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt);
  if (ret == 0) /* try in negative of normal direction */
  {
    dir[0] = -v->nx;
    dir[1] = -v->ny;
    dir[2] = -v->nz;
    pt[0] = v->x;
    pt[1] = v->y;
    pt[2] = v->z;
    load_triangle_vertices(mris, fno, U0, U1, U2, CURRENT_VERTICES);
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt);
    if (ret == 0) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisPlaceVertexInOrigFace: v does not intersect face!"));
  }

  /* now normalize the edges (l1 and l2) of the current triangle */
  SUB(l1, U1, U0);
  SUB(l2, U2, U0);
  SUB(P, int_pt, U0);
  l1_len = VLEN(l1);
  SCALAR_MUL(l1, 1.0 / l1_len, l1);
  l2_len = VLEN(l2);
  SCALAR_MUL(l2, 1.0 / l2_len, l2);
  l_len = VLEN(P);
  SCALAR_MUL(P, 1.0 / l_len, P);

  /*
    compute the angle between the two legs, and between P and l1.
    The ratio of these two angles will be used to place the vertex
    in the original triangle.
  */
  dot = DOT(l1, P);
  theta1 = acos(dot);
  if (theta1 < 0) {
    theta1 += 2 * PI;
  }
  dot = DOT(l1, l2);
  theta2 = acos(dot);
  if (theta2 < 0) {
    theta2 += 2 * PI;
  }
  if (!DZERO(theta2)) {
    theta_ratio = theta1 / theta2;
  }
  else /* degenerate triangle */
  {
    theta_ratio = 0;
  }

  /*
    express the ratio of the length of the line segment P-U0 as
    a scaled linear combination of the l1 and l2 (the legs of the
    triangle), where the relative weighting is based on the theta
    ratio. This will allow us to use the ratio and the original
    lengths of the legs to place the point in the corresponding location
    in the original triangle.
  */
  len_scale = l_len / (((1 - theta_ratio) * l1_len) + theta_ratio * l2_len);
  load_orig_triangle_vertices(mris, fno, U0, U1, U2);
  SUB(l1, U1, U0);
  SUB(l2, U2, U0);
  l1_len = VLEN(l1);
  SCALAR_MUL(l1, 1.0 / l1_len, l1);
  l2_len = VLEN(l2);
  SCALAR_MUL(l2, 1.0 / l2_len, l2);
  l_len = VLEN(P);
  SCALAR_MUL(P, 1.0 / l_len, P);

  /* compute angle between original legs */
  dot = DOT(l1, l2);
  theta2 = acos(dot);
  if (theta2 < 0) {
    theta2 += 2 * PI;
  }

  theta1 = theta_ratio * theta2; /* analogous angle in orig triangle */

  /* construct basis vector for plane defined by original triangle */
  SCALAR_MUL(l1, 1.0, e1); /* 1st basis vector is just l1 */
  CROSS(l1, l2, etmp);     /* vector orthogonal to l1 and l2 */
  CROSS(etmp, l1, e2);
  SCALAR_MUL(e2, (1.0 / VLEN(e2)), e2);

  /*
    express length of line segment in original triangle as a linear
    combination of the original leg lengths, using the same weighting.
  */
  l_len = len_scale * (((1 - theta_ratio) * l1_len) + theta_ratio * l2_len);

  /* rotate e1 by theta1 and scale it by l_len */
  x = l_len * cos(theta1);
  y = l_len * sin(theta1);

  /* express it in the global coordinate system */
  SCALAR_MUL(e1, x, e1);
  SCALAR_MUL(e2, y, e2);
  ADD(e1, e2, P);
  ADD(P, U0, P);

  cheapAssert(mris->status == mris->origxyz_status);
  MRISsetOriginalXYZ(mris_vno, vno, P[0], P[1], P[2]);

  return (NO_ERROR);
}


int MRISinverseSphericalMap(MRIS *mris, MRIS *mris_ico)
{
  double r;
  int fno, vno, num_ambiguous = 0, nfound, i, max_count;
  short *vcount = (short *)calloc(mris->nfaces, sizeof(short));

  MRISfreeDistsButNotOrig(mris);

  /* make sure they are they same size */
  r = MRISaverageRadius(mris_ico);
  MRISscaleBrain(mris_ico, mris_ico, 100.0 / r);
  r = MRISaverageRadius(mris);
  MRISscaleBrain(mris, mris, 100.0 / r);
  MRISstoreMetricProperties(mris);

  /*
    orig       positions are on cortical surface
    current    positions are on sphere.
  */
  MHT *mht = MHTcreateFaceTable(mris);

  /*
    for each vertex on the icosahedral surface, find what face it lies
    in on the spherical representation of the cortex. If it only lies
    within a single face, position it at the centroid of the original
    position of the face, and mark it as positioned.
  */
  MRISclearMarks(mris_ico);
  MRISclearMarks(mris);
  for (vno = 0; vno < mris_ico->nvertices; vno++) {
    VERTEX * const v = &mris_ico->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    fno = mrisFindUnambiguousFace(mris, mht, v, &nfound);
    if (fno >= 0) {
      vcount[fno]++;
      mrisPlaceVertexInOrigFace(mris_ico, vno, mris, fno);
      if (vno == Gdiag_no) {
        fprintf(stdout, "vertex %d maps to face %d at (%2.1f, %2.1f, %2.1f)\n", vno, fno, v->origx, v->origy, v->origz);
        mrisDumpFace(mris, fno, stderr);
      }
      v->marked = 1;
    }
    else {
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "v %d maps to %d faces\n", vno, nfound);
      }
      num_ambiguous++;
    }
  }
  fprintf(stdout, "%d non-invertible locations found - resolving ambiguity\n", num_ambiguous);

  MRISsoapBubbleOrigVertexPositions(mris_ico, 100);
  for (vno = 0; vno < mris_ico->nvertices; vno++) {
    VERTEX * const v = &mris_ico->vertices[vno];
    if (v->marked || v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    fno = mrisChooseFace(mris, mht, v);
    if (fno < 0)
      ErrorPrintf(ERROR_BADPARM, "unable to find face for ico vertex %d!!!\n", vno);
    else {
      mrisPlaceVertexInOrigFace(mris, vno, mris, fno);
      vcount[fno]++;
    }
  }
  for (max_count = i = 0; i < 50; i++) {
    for (nfound = fno = 0; fno < mris->nfaces; fno++) {
      if (vcount[fno] == i) {
        nfound++;
      }
    }
    if (nfound) {
      if (i > max_count) {
        max_count = i;
      }
      fprintf(stdout, "%d mappings to a single face %d times.\n", i, nfound);
    }
  }

  fprintf(stdout, "faces mapped to %d times: \n", max_count);
  for (fno = 0; fno < mris->nfaces; fno++) {
    if (vcount[fno] == max_count)
      fprintf(stdout, "\t%d (%d, %d, %d)\n", fno, mris->faces[fno].v[0], mris->faces[fno].v[1], mris->faces[fno].v[2]);
  }
  MHTfree(&mht);
  free(vcount);

  return (NO_ERROR);
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Add the triangle vno0 --> vno1 --> vno2 to the tessellation
  ------------------------------------------------------*/
int mrisAddFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2)
{
  int n, fno, ilist[1000], n0, n1;
  FACE *f;
  uchar ulist[1000];
  float norm[3], dot, cx, cy, cz;

  if (vno0 < 0 || vno1 < 0 || vno2 < 0)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisAddFace(%d,%d,%d)!\n", vno0, vno1, vno2));
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "adding face (%d, %d, %d)\n", vno0, vno1, vno2);
  }
  if (vno0 == 6160 && vno1 == 6189 && vno2 == 6176) {
    DiagBreak();
  }

  fno = mris->nfaces;
  MRISgrowNFaces(mris, fno+1);
  if (fno == Gdiag_no) {
    DiagBreak();
  }

  cheapAssert(mrisCanAttachFaceToVertices(mris, vno0, vno1, vno2));
  f = &mris->faces[fno];
  f->v[0] = vno0;
  f->v[1] = vno1;
  f->v[2] = vno2;

  for (n = 0; n < VERTICES_PER_FACE; n++) {
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[f->v[n]];
    
    if (v->num >= 255) {
      continue;
    }
    memmove(ilist, v->f, v->num * sizeof(int));
    ilist[v->num++] = fno;
    if (v->f) {
      free(v->f);
    }
    v->f = (int *)calloc(v->num, sizeof(int));
    if (!v->f) ErrorExit(ERROR_NOMEMORY, "mrisAddFace: could not allocate face list");
    memmove(v->f, ilist, v->num * sizeof(int));

    memmove(ulist, v->n, (v->num - 1) * sizeof(uchar));
    ulist[v->num - 1] = n;
    if (v->n) {
      free(v->n);
    }
    v->n = (uchar *)calloc(v->num, sizeof(uchar));
    if (!v->n) {
      ErrorExit(ERROR_NOMEMORY, "mrisAddFace: could not allocate n list");
    }
    memmove(v->n, ulist, v->num * sizeof(uchar));
  }

  mrisCheckVertexFaceTopology(mris);

  mrisCalculateCanonicalFaceCentroid(mris, fno, &cx, &cy, &cz);
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    mrisNormalFace(mris, fno, n, norm); /* compute face normal */
    dot = norm[0] * cx + norm[1] * cy + norm[2] * cz;

    if (dot < 0) /* they disagree - change order of vertices in face */
    {
      n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
      n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
      vno0 = f->v[n0];
      vno1 = f->v[n1];
      f->v[n0] = vno1;
      f->v[n1] = vno0;
      mrisSetVertexFaceIndex(mris, vno0, fno);
      mrisSetVertexFaceIndex(mris, vno1, fno);
    }
  }

  return (NO_ERROR);
}

/*!
  \fn int MRISshrinkFaces(MRIS *surf, double zthresh, int nmax)
  \brief Iteratively shrinks triangles by moving the three vertices
  toward the centroid to make the size of the triangle equal to the
  mean triangle size (prior to any shrinkage). Any triangle whose area
  is > mean+zthresh*std is shrunk. The process is iterative and
  finishes when all traingle areas are below threshold. A maximum
  number of iterations can be set with nmax (set to <=0 to have no
  limit).  Not clear whether ripped vertices or faces are
  excluded. This process can result in an intersection. Returns the
  number of iterations.
 */
int MRISshrinkFaces(MRIS *surf, double zthresh, int nmax)
{
  double areamean, areastd, areamin, areamax, areathresh;
  double *stats;
  FACE *f;

  MRIScomputeMetricProperties(surf);

  // Does not seem to work
  //areamean =  MRIScomputeFaceAreaStats(surf, &areastd, &areamin, &areamax);

  stats = MRIStriangleAreaStats(surf, NULL, NULL);
  printf("Shrink Pre:  %d %g %g %g %g\n",(int)stats[0],stats[1],stats[2],stats[3],stats[4]);
  areamean = stats[1];
  areastd  = stats[2];
  areamin  = stats[3];
  areamax  = stats[4];

  areathresh = areamean + zthresh*areastd;
  printf("MRISshrinkFaces(): %g %g %g %g %g\n",areamean,areastd,zthresh,areathresh,areamax);

  if(areamax < areathresh) return(0);

  // Continue to reduce face sizes until the largest goes below threshold
  int faceno, facenomax, iter=0;
  double localmax;
  while(areamax > areathresh){
    if(nmax > 0 && iter > nmax) break;
    iter++;
    // Go through all the faces and find the one with largest area
    facenomax = 0;
    areamax = 0;
    for(faceno=0; faceno < surf->nfaces; faceno++){
      f = &(surf->faces[faceno]);
      if(f->area > areamax){
	facenomax = faceno;
	areamax = f->area;
      }
    }
    // Shrink the biggest face
    f = &(surf->faces[facenomax]);
    localmax = MRISshrinkFace(surf, facenomax, areamean/f->area);
    // Check whether a neighboring triangle was made bigger than the thresh
    if(localmax > areamax) areamax = localmax;
  }

  printf("iter = %d\n",iter);
  MRIScomputeMetricProperties(surf);

  stats = MRIStriangleAreaStats(surf, NULL, stats);
  printf("Shrink Post: %d %g %g %g %g\n",(int)stats[0],stats[1],stats[2],stats[3],stats[4]);
  free(stats);

  return(iter);
}

/*!
\fn int MRISshrinkFaceCorner(MRIS *surf, int faceno, int nthv, double dist, double *vx, double *vy, double *vz)
\brief Shinks the given corner of a given face by moving the corner towards the centoid by dist where
dist is a value between 0 (don't move it at all) and 1 (move all the way to the centroid). This function
just computes the new vertex location and does not change the actual vertex xyz. The area of the triangle
will shrink by (1-dist)^2 assuming that all vertices of the face are moved in the same way.
*/
int MRISshrinkFaceCorner(MRIS *surf, int faceno, int nthv, double dist, double *vx, double *vy, double *vz)
{
  FACE *f;
  f = &(surf->faces[faceno]);
  VERTEX *v0, *v1, *v2;
  double cx,cy,cz, wx,wy,wz, ux,uy,uz, m;
  int n0, n1, n2;

  n0 = nthv;
  n1 = n0 + 1;
  if(n1 >= 3) n1 = 0;
  n2 = n1 + 1;
  if(n2 >= 3) n2 = 0;

  v0 = &(surf->vertices[f->v[n0]]);// this vertex
  v1 = &(surf->vertices[f->v[n1]]);
  v2 = &(surf->vertices[f->v[n2]]);

  // Compute the centroid.
  // This only needs to be done once per face, so a little inefficient
  cx = (v0->x + v1->x + v2->x)/3.0;
  cy = (v0->y + v1->y + v2->y)/3.0;
  cz = (v0->z + v1->z + v2->z)/3.0;

  // Vector from this vertex to centroid
  wx = cx - v0->x;
  wy = cy - v0->y;
  wz = cz - v0->z;
  // Distance from this vertex to centroid
  m = sqrt(wx*wx + wy*wy + wz*wz);
  // Unit vector from this vertex to centroid
  ux = wx/m;
  uy = wy/m;
  uz = wz/m;

  // Compute new location for this vertex that is dist fraction from
  // this vertex to the centroid. Does not change face structure
  *vx = (v0->x + dist*m*ux);
  *vy = (v0->y + dist*m*uy);
  *vz = (v0->z + dist*m*uz);

  return(0);
}

/*!
  \fn double MRISshrinkFace(MRIS *surf, int faceno, double newareafraction)
  \brief Shrinks the given face to be newareafraction times the original area.
  The f->area of the face and all affected neighboring faces are updated
  so that MRIScomputeMetricProperties() does not need to be run to update
  the face areas (but will be for other metric properties). It returns
  the maximum area of all the affected faces after shrinking the given
  face (other faces will have gotten bigger).
 */
double MRISshrinkFace(MRIS *surf, int faceno, double newareafraction)
{
  FACE *f;
  f = &(surf->faces[faceno]);
  double vx0, vy0, vz0, vx1, vy1, vz1, vx2, vy2, vz2;
  VERTEX *v0, *v1, *v2;
  double dist,maxarea;

  // Compute the distance factor needed to realize the area
  dist = 1-sqrt(newareafraction);

  // Compute the coords for each vertex to shrink
  MRISshrinkFaceCorner(surf, faceno, 0, dist, &vx0, &vy0, &vz0);
  MRISshrinkFaceCorner(surf, faceno, 1, dist, &vx1, &vy1, &vz1);
  MRISshrinkFaceCorner(surf, faceno, 2, dist, &vx2, &vy2, &vz2);

  // Now update the vertex positions
  v0 = &(surf->vertices[f->v[0]]);
  v1 = &(surf->vertices[f->v[1]]);
  v2 = &(surf->vertices[f->v[2]]);

  v0->x = vx0;
  v0->y = vy0;
  v0->z = vz0;
  
  v1->x = vx1;
  v1->y = vy1;
  v1->z = vz1;
  
  v2->x = vx2;
  v2->y = vy2;
  v2->z = vz2;

  // Now update the affected neighboring faces
  // First get a list of neighboring faces
  int n, nthface, facenolist[1000], nfacenolist=0;
  VERTEX_TOPOLOGY *vt;
  for(n=0; n < 3; n++){
    vt = &(surf->vertices_topology[f->v[n]]);
    for(nthface = 0; nthface < vt->vnum; nthface++){
      facenolist[nfacenolist] = vt->f[nthface];
      nfacenolist++;
    }
  }
  // Make sure they are unique
  int nunique, *ulist;
  ulist = unqiue_int_list(facenolist, nfacenolist, &nunique);
  // Now change the area of each face
  maxarea = 0;
  for(nthface = 0; nthface < nunique; nthface++){
    faceno = ulist[nthface];
    f = &(surf->faces[faceno]);
    f->area = fabs(mrisComputeArea(surf, faceno, 0))/2.0; //0 does not matter
    // Not sure why it needs to be div by 2, except that it works
    if(maxarea < f->area) maxarea = f->area;
  }
  free(ulist);

  return(maxarea);
}

