#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
#define COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND
/*
 * @file utilities operating on Original
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_io.h"


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


/* project onto the sphere of radius DEFAULT_RADIUS */
static void sphericalProjection(float xs, float ys, float zs, float *xd, float *yd, float *zd)
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


static int mrisSphericalProjection(MRIS *mris)
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

    sphericalProjection(v->x, v->y, v->z, &v->x, &v->y, &v->z);
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

      sphericalProjection(v->x,v->y,v->z,&v->x,&v->y,&v->z);

      if(n == 88 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 89 )
      fprintf(stderr,"af 2 sp: vertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
      if(n == 209 )
      fprintf(stderr,"af 2 sp: nvertex %d (%f,%f,%f)\n",n,v->x,v->y,v->z);
    */
  }
  return NO_ERROR;
}

static int mris_project_point_into_face(
    MRI_SURFACE *mris, FACE *face, int which, double x, double y, double z, double *px, double *py, double *pz)
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

  return (NO_ERROR);
}


static int mrisProjectOntoSurface(MRI_SURFACE *mris, int which_vertices)
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

  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisProjectSurface(MRI_SURFACE *mris)
{
  /*  MRISupdateSurface(mris) ;*/
  switch (mris->status) {
    case MRIS_PLANE:
      MRISflattenPatch(mris);
      break;
    case MRIS_SPHERICAL_PATCH:
      mrisSphericalProjection(mris);
      break;
    case MRIS_PARAMETERIZED_SPHERE:
      MRISprojectOntoSphere(mris, mris, mris->radius);
      break;
    case MRIS_SPHERE:
      MRISprojectOntoSphere(mris, mris, mris->radius);
      break;
    case MRIS_ELLIPSOID:
      MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
      break;
    //        case PROJECT_PLANE:
    /*    mrisOrientPlane(mris) ;*/
    //                break ;
    case MRIS_RIGID_BODY:
      /*    MRISprojectOntoSphere(mris, mris, mris->radius) ;*/
      mris->status = MRIS_RIGID_BODY;
      break;
    case MRIS_PIAL_SURFACE:
      mrisProjectOntoSurface(mris, PIAL_VERTICES);
      break;
    default:
      break;
  }
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
      if (!isfinite(v->tx)) {
        DiagBreak();
      }
    }

    for (n = 0; n < mris->nvertices; n++) {
      VERTEX * const v  = &mris->vertices[n];
      sphericalProjection(v->tx, v->ty, v->tz, &v->x, &v->y, &v->z);
      if (!isfinite(v->x)) {
        DiagBreak();
      }
    }
  }

  return NO_ERROR;
}


static int mrisAssignFaces(MRI_SURFACE *mris, MHT *mht, int which_vertices)
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
  
  return (NO_ERROR);
}


int mrisComputePlaneTerm(MRI_SURFACE *mris, double l_plane, double l_spacing)
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int mrisComputeSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, dist_scale;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

#if METRIC_SCALE
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
#else
  dist_scale = 1.0;
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->border && !v->neg) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    sx = sy = sz = 0.0;
    n = 0;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n > 0) {
      sx = dist_scale * sx / n;
      sy = dist_scale * sy / n;
      sz = dist_scale * sz / n;
    }

    sx *= l_spring;
    sy *= l_spring;
    sz *= l_spring;
    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d spring term:         (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  Note: this function assumes that the mris surface
  has the original (i.e. after
  global rotational alignment) spherical coordinates in the TMP2_VERTICES
  ------------------------------------------------------*/
int mrisComputeLaplacianTerm(MRI_SURFACE *mris, double l_lap)
{
  int vno, n, m;
  float x, y, z, vx, vy, vz, vnx, vny, vnz, dx, dy, dz;

  if (FZERO(l_lap)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->border && !v->neg) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    n = 0;
    vx = v->x - v->tx2;
    vy = v->y - v->ty2;
    vz = v->z - v->tz2;
    dx = dy = dz = 0.0f;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        vnx = vn->x - vn->tx2;
        vny = vn->y - vn->ty2;
        vnz = vn->z - vn->tz2;
        dx += (vnx - vx);
        dy += (vny - vy);
        dz += (vnz - vz);
        if ((x == Gx && y == Gy && z == Gz) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)
          printf(
              "\tvertex %d: V=(%2.2f,%2.2f,%2.2f), "
              "DX=(%2.2f,%2.2f,%2.2f)\n",
              vno,
              vnx,
              vny,
              vnz,
              vnx - vx,
              vny - vy,
              vnz - vz);
        n++;
      }
    }
    if (n > 0) {
      dx = dx * l_lap / n;
      dy = dy * l_lap / n;
      dz = dz * l_lap / n;
    }

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no) {
      printf("l_lap: v %d: DX=(%2.2f,%2.2f,%2.2f)\n", vno, dx, dy, dz);
    }
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Compute a spring term, and normalize it by removing the
  average normal component which typically forces the surface
  to shrink.
  ------------------------------------------------------*/
int mrisComputeNormalizedSpringTerm(MRI_SURFACE *const mris, double const l_spring)
{
  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  float dist_scale_init;
#if METRIC_SCALE
  if (mris->patch) {
    dist_scale_init = 1.0;
  }
  else {
    dist_scale_init = sqrt(mris->orig_area / mris->total_area);
  }
#else
  dist_scale_init = 1.0;
#endif
  const float dist_scale = dist_scale_init;

  const double num = (double)MRISvalidVertices(mris);

  double dot_total = 0.0;
  int vno;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : dot_total)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    float sx = 0.0, sy = 0.0, sz = 0.0;
    int n = 0;
    int m;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX *vn = &mris->vertices[vt->v[m]];
      if (vn->ripflag) continue;

      // Almost all the time in this loop is spent in the above conditions
      // The following is NOT where the time goes!
      //
      sx += vn->x - v->x;
      sy += vn->y - v->y;
      sz += vn->z - v->z;

      n++;
    }
    if (n == 0) continue;

    float multiplier = dist_scale / n;
    sx *= multiplier;
    sy *= multiplier;
    sz *= multiplier;

    dot_total += l_spring * (v->nx * sx + v->ny * sy + v->nz * sz);
    v->dx += l_spring * sx;
    v->dy += l_spring * sy;
    v->dz += l_spring * sz;

    if (vno == Gdiag_no)
      fprintf(
          stdout, "v %d spring norm term: (%2.3f, %2.3f, %2.3f)\n", vno, l_spring * sx, l_spring * sy, l_spring * sz);

    ROMP_PFLB_end
  }
  ROMP_PF_end

  float const dot_avg = dot_total / num;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX *v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    v->dx -= dot_avg * v->nx;
    v->dy -= dot_avg * v->ny;
    v->dz -= dot_avg * v->nz;
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int mrisComputeConvexityTerm(MRI_SURFACE *mris, double l_convex)
{
  int vno;

  if (FZERO(l_convex)) {
    return (NO_ERROR);
  }

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    int n, m;
    float sx, sy, sz, nx, ny, nz, nc, x, y, z;

    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    nx = vertex->nx;
    ny = vertex->ny;
    nz = vertex->nz;
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
    nc = sx * nx + sy * ny + sz * nz; /* projection onto normal */
    if (nc < 0) {
      nc = 0;
    }
    sx = nc * nx; /* move in normal direction */
    sy = nc * ny;
    sz = nc * nz;

    vertex->dx += l_convex * sx;
    vertex->dy += l_convex * sy;
    vertex->dz += l_convex * sz;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d convexity term: (%2.3f, %2.3f, %2.3f)\n", vno, l_convex * sx, l_convex * sy, l_convex * sz);
      
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, nx, ny, nz, nc, x, y, z;

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

    nx = vertex->nx;
    ny = vertex->ny;
    nz = vertex->nz;
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
    nc = sx * nx + sy * ny + sz * nz; /* projection onto normal */
    sx = l_spring * nc * nx;          /* move in normal direction */
    sy = l_spring * nc * ny;
    sz = l_spring * nc * nz;

    vertex->dx += sx;
    vertex->dy += sy;
    vertex->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, nc;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->border && !v->neg) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    sx = sy = sz = 0.0;
    n = 0;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n > 0) {
      sx = sx / n;
      sy = sy / n;
      sz = sz / n;
    }

    nc = sx * v->nx + sy * v->ny + sz * v->nz; /* projection onto normal */
    sx = l_spring * (sx - nc * v->nx);         /* remove  normal component
                                                            and then scale */
    sy = l_spring * (sy - nc * v->ny);
    sz = l_spring * (sz - nc * v->nz);

    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) printf("v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  File Format is:

  name x y z
  ------------------------------------------------------*/
int mrisComputeNonlinearTangentialSpringTerm(MRI_SURFACE *mris, double l_spring, double min_dist)
{
  int vno, m, n;
  float sx, sy, sz, x, y, z, dx, dy, dz;
  double d, scale;

  if (FZERO(l_spring)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    sx = sy = sz = 0.0;
    for (dx = dy = dz = 0.0, n = m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        sx = x - vn->x;
        sy = y - vn->y;
        sz = z - vn->z;  // move away from nbr
        d = sqrt(sx * sx + sy * sy + sz * sz);
        if (d < min_dist) {
          scale = (min_dist - d) / min_dist;
          d = scale * (v->e1x * sx + v->e1y * sy + v->e1z * sz);
          dx += v->e1x * d;
          dy += v->e1y * d;
          dz += v->e1z * d;
          d = scale * (v->e2x * sx + v->e2y * sy + v->e2z * sz);
          dx += v->e2x * d;
          dy += v->e2y * d;
          dz += v->e2z * d;
          if (vno == Gdiag_no) {
            DiagBreak();
          }
          n++;
        }
      }
    }
    dx *= l_spring;
    dy *= l_spring;
    dz *= l_spring;
    if (vno == Gdiag_no && n > 0)
      printf("v %d nonlinear spring tangent term: (%2.3f, %2.3f, %2.3f)\n", vno, dx, dy, dz);
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
  }

  return (NO_ERROR);
}

int mrisComputeLinkTerm(MRI_SURFACE *mris, double l_link, int pial)
{
  int vno;
  VERTEX *v;
  float dx, dy, dz, lx, ly, lz, len;

  if (FZERO(l_link)) {
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

    lx = v->pialx - v->origx;
    ly = v->pialy - v->origy;
    lz = v->pialz - v->origz;
    len = sqrt(lx * lx + ly * ly + lz * lz);
    if (len < .25) /* can't accurately estimate vector
                                    connecting white and pial */
    {
      continue;
    }
    lx /= len;
    ly /= len;
    lz /= len;

    dx = l_link * (v->nx - lx);
    dy = l_link * (v->ny - ly);
    dz = l_link * (v->nz - lz);

    if (pial == 0) {
      dx *= -1;
      dy *= -1;
      dz *= -1;
    }

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d link %s term: (%2.3f, %2.3f, %2.3f), "
              "Nl=(%2.1f, %2.1f, %2.1f), Ns=(%2.1f, %2.1f, %2.1f), "
              "dot=%2.3f\n",
              vno,
              pial ? "pial" : "white",
              dx,
              dy,
              dz,
              lx,
              ly,
              lz,
              v->nx,
              v->ny,
              v->nz,
              lx * v->nx + ly * v->ny + lz * v->nz);
  }

  return (NO_ERROR);
}


//  these versions use a full 2D fit y = a x^2 + b y^2 + c x + d y + e, then use e as the error term

/*-----------------------------------------------------
  Description
  Fit a 2-d quadratic to the surface locally and compute the SSE as
  the square of the constant term (the distance the quadratic fit surface
  is from going through the central vertex).  Move the
  vertex in the normal direction to improve the fit.
  ------------------------------------------------------*/
int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)  // BEVIN mris_make_surfaces 4
{
  if (FZERO(l_curv)) {
    return (NO_ERROR);
  }

  mrisComputeTangentPlanes(mris);

  typedef struct Reused {
    VECTOR * v_n  ;
    VECTOR * v_P  ;
    VECTOR * v_e1 ;
    VECTOR * v_e2 ;
    VECTOR * v_nbr;
  } Reused;
  
#ifdef HAVE_OPENMP
  int const maxThreads = omp_get_max_threads();
#else
  int const maxThreads = 1;
#endif
  Reused* reusedByThread = (Reused*)calloc(maxThreads, sizeof(Reused));

  int vno;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#ifdef HAVE_OPENMP
  int const tid        = omp_get_thread_num();
#else
  int const tid        = 0;
#endif

    Reused* reused = reusedByThread + tid;
    #define REUSE(NAME,DIM) \
        VECTOR* NAME = reused->NAME; if (!NAME) NAME = reused->NAME = VectorAlloc(DIM, MATRIX_REAL);
    REUSE(v_n   ,3)
    REUSE(v_P   ,5)
    REUSE(v_e1  ,3)
    REUSE(v_e2  ,3)
    REUSE(v_nbr ,3)
    #undef REUSE

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    FILE *fp = NULL;
    if (vno == Gdiag_no) fp = fopen("qcurv.dat", "w");

    VECTOR* v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL);    /* heights above TpS */
    VECTOR_LOAD(v_n,  v->nx,  v->ny,  v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);

    MATRIX* m_X = MatrixAlloc(vt->vtotal, 5, MATRIX_REAL); /* 2-d quadratic fit */

    int n;
    for (n = 0; n < vt->vtotal; n++) /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      float ui = V3_DOT(v_e1, v_nbr);
      float vi = V3_DOT(v_e2, v_nbr);

      *MATRIX_RELT(m_X, n + 1, 1) = ui * ui;
      *MATRIX_RELT(m_X, n + 1, 2) = vi * vi;
      *MATRIX_RELT(m_X, n + 1, 3) = ui;
      *MATRIX_RELT(m_X, n + 1, 4) = vi;
      *MATRIX_RELT(m_X, n + 1, 5) = 1;
      if (vno == Gdiag_no) fprintf(fp, "%d %f %f %f\n", vt->v[n], ui, vi, VECTOR_ELT(v_Y, n + 1));
    }

    if (vno == Gdiag_no) fclose(fp);

    MATRIX *m_X_inv = MatrixPseudoInverse(m_X, NULL);
    if (!m_X_inv) {
      MatrixFree(&m_X);
      VectorFree(&v_Y);
      continue;
    }
    
    v_P = MatrixMultiply(m_X_inv, v_Y, v_P);
    //oat a = VECTOR_ELT(v_P, 1);
    float e = VECTOR_ELT(v_P, 5);
    e *= l_curv;

    v->dx += e * v->nx;
    v->dy += e * v->ny;
    v->dz += e * v->nz;

    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "e=%2.1f\n",
              vno,
              e * v->nx,
              e * v->ny,
              e * v->nz,
              e);
              
    VectorFree(&v_Y);
    MatrixFree(&m_X);
    MatrixFree(&m_X_inv);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  { int tid;
    for (tid = 0; tid < maxThreads; tid++) {
      Reused* reused = reusedByThread + tid;
      if (reused->v_n)   VectorFree(&reused->v_n);
      if (reused->v_e1)  VectorFree(&reused->v_e1);
      if (reused->v_e2)  VectorFree(&reused->v_e2);
      if (reused->v_nbr) VectorFree(&reused->v_nbr);
      if (reused->v_P)   VectorFree(&reused->v_P);
  } }
  
  return (NO_ERROR);
}


// these versions use a 1d fit y = a*r^2 + b

/*-----------------------------------------------------
  Description
  Fit a 1-d quadratic to the surface locally and move the
  vertex in the normal direction to improve the fit.
  ------------------------------------------------------*/
int mrisComputeSurfaceNormalIntersectionTerm(MRI_SURFACE *mris, MHT *mht, double l_norm, double max_dist)
{
  int vno;
  double step;

  if (FZERO(l_norm)) {
    return (NO_ERROR);
  }

  step = mht->vres / 2;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX *v;
    double d, dist, dx, dy, dz, scale;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    if (v->ripflag) continue;
    dist = max_dist;
    for (dist = 0; dist < 2 * max_dist; dist += step) {
      d = dist;
      if (mrisDirectionTriangleIntersection(mris, v->x, v->y, v->z, v->nx, v->ny, v->nz, mht, &d, vno)) {
        if (d < max_dist) {
          if (vno == Gdiag_no) printf("v %d surface self intersection at distance %2.3f\n", vno, d);
          scale = (max_dist - d) / max_dist;
#define MAX_SCALE 3
          scale = MIN(MAX_SCALE, l_norm * exp(1.0 / (20 * (1 - scale))));
#if 0
      dx = (d-max_dist)/max_dist*l_norm*v->nx ;
      dy = (d-max_dist)/max_dist*l_norm*v->ny ;
      dz = (d-max_dist)/max_dist*l_norm*v->nz ;
#else
          dx = -scale * v->nx;
          dy = -scale * v->ny;
          dz = -scale * v->nz;
#endif
          if (vno == Gdiag_no) printf("\t intersection scale %2.2f, vector (%2.2f, %2.2f, %2.2f)\n", scale, dx, dy, dz);
          v->dx += dx;
          v->dy += dy;
          v->dz += dz;
          break;
        }
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  return (NO_ERROR);
}


double mrisComputeDuraError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double dura_thresh = parms->dura_thresh, sse;
  MRI *mri_dura = parms->mri_dura;
  int vno;
  VERTEX *v;
  float x, y, z;
  double val0, xw, yw, zw, delV;

  if (FZERO(parms->l_dura)) {
    return (0.0);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    MRISvertexToVoxel(mris, v, mri_dura, &xw, &yw, &zw);
    MRIsampleVolume(mri_dura, xw, yw, zw, &val0);
    if (val0 < dura_thresh) {
      continue;  // no effect
    }

    delV = dura_thresh - val0;
    sse += delV * delV;
  }

  return (sse);
}


int mrisComputeDuraTerm(MRI_SURFACE *mris, double l_dura, MRI *mri_dura, double dura_thresh)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz, dx, dy, dz;
  double val0, xw, yw, zw, del, delI, delV;

  if (FZERO(l_dura)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    MRISvertexToVoxel(mris, v, mri_dura, &xw, &yw, &zw);
    MRIsampleVolume(mri_dura, xw, yw, zw, &val0);
    if (val0 < dura_thresh) {
      continue;  // no effect
    }

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

    delV = dura_thresh - val0;
    delI = 1;

    del = l_dura * delV * delI;

    dx = nx * del;
    dy = ny * del;
    dz = nz * del;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;

    if (vno == Gdiag_no) {
      double xwi, ywi, zwi, xwo, ywo, zwo, val_inside, val_outside;

      x = v->x;
      y = v->y;
      z = v->z;

      /* sample outward from surface */
      xw = x + nx;
      yw = y + ny;
      zw = z + nz;
      MRISsurfaceRASToVoxel(mris, mri_dura, xw, yw, zw, &xwo, &ywo, &zwo);
      MRIsampleVolume(mri_dura, xw, yw, zw, &val_outside);

      /* sample inward from surface */
      xw = x - nx;
      yw = y - ny;
      zw = z - nz;
      MRISsurfaceRASToVoxel(mris, mri_dura, xw, yw, zw, &xwi, &ywi, &zwi);
      MRIsampleVolume(mri_dura, xw, yw, zw, &val_inside);

      MRISsurfaceRASToVoxel(mris, mri_dura, x, y, z, &xw, &yw, &zw);
      fprintf(stdout,
              "D(%2.1f,%2.1f,%2.1f)=%2.1f, Do(%2.1f,%2.1f,%2.1f)=%2.1f, "
              "Di(%2.1f,%2.1f,%2.1f)=%2.1f\n",
              xw,
              yw,
              zw,
              val0,
              xwo,
              ywo,
              zwo,
              val_outside,
              xwi,
              ywi,
              zwi,
              val_inside);
      fprintf(stdout,
              "v %d dura term:      (%2.3f, %2.3f, %2.3f), "
              "delV=%2.1f, delI=%2.0f\n",
              vno,
              dx,
              dy,
              dz,
              delV,
              delI);
    }
  }

  return (NO_ERROR);
}


#define NORMAL_MOVEMENT 0.1
#define NSAMPLES 15
#define SAMPLE_DISTANCE 0.1

/*!
  \fn static int mrisComputeIntensityTerm()
  \brief Computes the step needed to minimize the intensity term. The
  intensity term is a target intensity value as indicated by v->val
  (eg, see MRIScomputeBorderValues_new()).  The error is the
  difference between the actual intensity at the vertex and
  v->val. v->{dx,dy,dz} are incremented. v->sigma is used when
  computing the direction of the gradient of the intensity along the
  normal. There is a limit on the maximum steps size controlled by a
  hidden parameter.
 */
int mrisComputeIntensityTerm(MRI_SURFACE *mris,
                                    double l_intensity,
                                    MRI *mri_brain,
                                    MRI *mri_smooth,
                                    double sigma_global,
                                    INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz, dx, dy, dz;
  double val0, xw, yw, zw, del, val_outside, val_inside, delI, delV, k, ktotal_outside, xvi, yvi, zvi, interior,
      ktotal_inside;
  double sigma;
  MRI *mri_interior;

  if (FZERO(l_intensity)) {
    return (NO_ERROR);
  }

  if (parms->grad_dir == 0 && parms->fill_interior)  // create binary mask of interior of surface
  {
    mri_interior = MRISfillInterior(mris, 0.5, NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRIwrite(mri_interior, "int.mgz");
    }
  }
  else {
    mri_interior = NULL;
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    // Sample the volume at the vertex
    MRISvertexToVoxel(mris, v, mri_brain, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val0);

    sigma = v->val2; // smoothing level for this vertex 
    if (FZERO(sigma)) sigma = sigma_global;
    if (FZERO(sigma)) sigma = 0.25;

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

    /* compute intensity gradient along the normal. Only used to get the right sign */
    if (parms->grad_dir == 0) {
      double dist, val, step_size;
      int n;

      // Hidden parameter used to compute the step size
      step_size = MIN(sigma / 2, MIN(mri_brain->xsize, MIN(mri_brain->ysize, mri_brain->zsize)) * 0.5);
      ktotal_inside = ktotal_outside = 0.0;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size; dist <= 2 * sigma; dist += step_size, n++) {
        k = exp(-dist * dist / (2 * sigma * sigma));
        xw = x + dist * nx;
        yw = y + dist * ny;
        zw = z + dist * nz;
        if (mri_interior) {
          MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
          MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
        }

        if (mri_interior == NULL || interior < .9) {
          ktotal_outside += k;
          MRISsurfaceRASToVoxelCached(mris, mri_brain, xw, yw, zw, &xw, &yw, &zw);
          MRIsampleVolume(mri_brain, xw, yw, zw, &val);
          val_outside += k * val;
        }
        else {
          DiagBreak();
        }

        xw = x - dist * nx;
        yw = y - dist * ny;
        zw = z - dist * nz;
        if (mri_interior) {
          MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
          MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
        }

        if (mri_interior == NULL || interior > 0) {
          MRISsurfaceRASToVoxelCached(mris, mri_brain, xw, yw, zw, &xw, &yw, &zw);
          MRIsampleVolume(mri_brain, xw, yw, zw, &val);
          val_inside += k * val;
          ktotal_inside += k;
        }
        else {
          DiagBreak();
        }
      }
      if (ktotal_inside > 0) {
        val_inside /= (double)ktotal_inside;
      }
      if (ktotal_outside > 0) {
        val_outside /= (double)ktotal_outside;
      }
    }
    else  // don't compute gradient - assume
    {
      val_outside = parms->grad_dir;
      val_inside = -parms->grad_dir;
    }

    // Difference between target intensity and actual intensity
    delV = v->val - val0; 
    // Dont allow the difference to be greater than 5 or less than -5
    // Hidden parameter 5
    if (delV > 5)
      delV = 5;
    else if (delV < -5) 
      delV = -5;

    // Gradient of the intensity at this location wrt a change along the normal
    delI = (val_outside - val_inside) / 2.0;
    // Change delI into +1 or -1
    if (!FZERO(delI))  delI /= fabs(delI);
    else               delI = -1; /* intensities tend to increase inwards */

    // Weight intensity error by cost weighting
    del = l_intensity * delV * delI;

    // Set to push vertex in the normal direction by this amount
    dx = nx * del;
    dy = ny * del;
    dz = nz * del;
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;

    if (vno == Gdiag_no) {
      double xwi, ywi, zwi, xwo, ywo, zwo;

      x = v->x;
      y = v->y;
      z = v->z;

      /* sample outward from surface */
      xw = x + mri_smooth->xsize * nx;
      yw = y + mri_smooth->ysize * ny;
      zw = z + mri_smooth->zsize * nz;
      MRISsurfaceRASToVoxelCached(mris, mri_smooth, xw, yw, zw, &xwo, &ywo, &zwo);
      /* sample inward from surface */
      xw = x - mri_smooth->xsize * nx;
      yw = y - mri_smooth->ysize * ny;
      zw = z - mri_smooth->zsize * nz;
      MRISsurfaceRASToVoxelCached(mris, mri_smooth, xw, yw, zw, &xwi, &ywi, &zwi);
      MRISsurfaceRASToVoxelCached(mris, mri_smooth, x, y, z, &xw, &yw, &zw);
      fprintf(stdout,
              "I(%2.1f,%2.1f,%2.1f)=%2.1f, Io(%2.1f,%2.1f,%2.1f)=%2.1f, "
              "Ii(%2.1f,%2.1f,%2.1f)=%2.1f\n",
              xw,
              yw,
              zw,
              val0,
              xwo,
              ywo,
              zwo,
              val_outside,
              xwi,
              ywi,
              zwi,
              val_inside);
      if (val_inside < -20 && val_outside > -2) DiagBreak();
      fprintf(stdout,
              "v %d intensity term:      (%2.3f, %2.3f, %2.3f), "
              "delV=%2.1f, delI=%2.0f, sigma=%2.1f, target=%2.1f\n",
              vno,
              dx,
              dy,
              dz,
              delV,
              delI,
              sigma,
              v->val);
    } // end diag
  } // loop over vertices

  if (mri_interior) {
    MRIfree(&mri_interior);
  }
  return (NO_ERROR);
}

int mrisComputeTargetLocationTerm(MRI_SURFACE *mris, double l_location, INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  double dx, dy, dz, norm;

  if (FZERO(l_location)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();

    dx = v->targx - v->x;
    dy = v->targy - v->y;
    dz = v->targz - v->z;

    norm = sqrt(dx * dx + dy * dy + dz * dz);
#define LOCATION_MOVE_LEN 0.25
    if (norm > LOCATION_MOVE_LEN)  // so things move at the same speed
    {
      dx /= norm;
      dy /= norm;
      dz /= norm;
      dx *= LOCATION_MOVE_LEN;
      dy *= LOCATION_MOVE_LEN;
      dz *= LOCATION_MOVE_LEN;
    }

    if (vno == Gdiag_no) {
      fprintf(stdout,
              "l_location: targ (%2.1f, %2.1f, %2.f), "
              "current (%2.1f, %2.1f, %2.1f), "
              "del (%2.1f, %2.1f, %2.1f), norm=%2.1f, dot=%2.3f\n",
              v->targx,
              v->targy,
              v->targz,
              v->x,
              v->y,
              v->z,
              l_location * dx,
              l_location * dy,
              l_location * dz,
              norm,
              dx * v->nx + dy * v->ny + dz * v->nz);
    }
    if (!devFinite(dx) || !devFinite(dy) || !devFinite(dz)) {
      DiagBreak();
    }

    v->dx += l_location * dx;
    v->dy += l_location * dy;
    v->dz += l_location * dz;
  }

  return (NO_ERROR);
}

int mrisComputeIntensityTerm_mef(MRI_SURFACE *mris,
                                        double l_intensity,
                                        MRI *mri_30,
                                        MRI *mri_5,
                                        double sigma_global,
                                        float weight30,
                                        float weight5,
                                        INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz, dx, dy, dz;
  double val0, xw, yw, zw, del, val_outside, val_inside, delI, delV, k, ktotal_outside, ktotal_inside, interior, xvi,
      yvi, zvi;
  double sigma;
  MRI *mri_interior;

  if (FZERO(l_intensity)) {
    return (NO_ERROR);
  }

  if (parms->grad_dir == 0 && parms->fill_interior)  // create binary mask of interior of surface
  {
    mri_interior = MRISfillInterior(mris, 0.5, NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRIwrite(mri_interior, "int.mgz");
    }
  }
  else {
    mri_interior = NULL;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

    sigma = v->val2;
    if (FZERO(sigma)) {
      sigma = sigma_global;
    }
    if (!FZERO(weight30)) {
      MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw);
      MRIsampleVolume(mri_30, xw, yw, zw, &val0);

      /* compute intensity gradient using smoothed volume */

      if (parms->grad_dir == 0) {
        double dist, val, step_size;
        int n;

        step_size = MIN(sigma / 2, MIN(mri_30->xsize, MIN(mri_30->ysize, mri_30->zsize)) * 0.5);
        ktotal_outside = ktotal_inside = 0.0;
        for (n = 0, val_outside = val_inside = 0.0, dist = step_size; dist <= 2 * sigma; dist += step_size, n++) {
          k = exp(-dist * dist / (2 * sigma * sigma));
          xw = x + dist * nx;
          yw = y + dist * ny;
          zw = z + dist * nz;
          if (mri_interior) {
            MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
            MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
          }

          if (mri_interior == NULL || interior < .9) {
            MRIsampleVolume(mri_30, xw, yw, zw, &val);
            val_outside += k * val;
            ktotal_outside += k;
          }
          else {
            DiagBreak();
          }

          xw = x - dist * nx;
          yw = y - dist * ny;
          zw = z - dist * nz;
          if (mri_interior) {
            MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
            MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
          }
          if (mri_interior == NULL || interior > 0) {
            MRIsampleVolume(mri_30, xw, yw, zw, &val);
            val_inside += k * val;
            ktotal_inside += k;
          }
          else {
            DiagBreak();
          }
        }
        if (ktotal_inside > 0) {
          val_inside /= (double)ktotal_inside;
        }
        if (ktotal_outside > 0) {
          val_outside /= (double)ktotal_outside;
        }
      }
      else  // don't compute gradient - assume
      {
        val_outside = parms->grad_dir;
        val_inside = -parms->grad_dir;
      }

      delV = v->val - val0;
      delI = (val_outside - val_inside) / 2.0;

      if (!FZERO(delI)) {
        delI /= fabs(delI);
      }
      else {
        delI = -1;  // intensity tends to increase inwards for flash30
      }

      if (delV > 5) {
        delV = 5;
      }
      else if (delV < -5) {
        delV = -5;
      }

      del = l_intensity * delV * delI;

      dx = nx * del * weight30;
      dy = ny * del * weight30;
      dz = nz * del * weight30;

      if (dx * nx + dy * ny + dz * nz < 0) {
        DiagBreak();
      }

      if (vno == Gdiag_no) {
        double xwi, ywi, zwi, xwo, ywo, zwo;

        x = v->x;
        y = v->y;
        z = v->z;

        /* sample outward from surface */
        xw = x + nx;
        yw = y + ny;
        zw = z + nz;
        MRISsurfaceRASToVoxelCached(mris, mri_30, xw, yw, zw, &xwo, &ywo, &zwo);
        /* sample inward from surface */
        xw = x - nx;
        yw = y - ny;
        zw = z - nz;
        MRISsurfaceRASToVoxelCached(mris, mri_30, xw, yw, zw, &xwi, &ywi, &zwi);
        fprintf(stdout,
                "I30(%2.1f,%2.1f,%2.1f)=%2.1f, Io(%2.1f,%2.1f,%2.1f)=%2.1f, "
                "Ii(%2.1f,%2.1f,%2.1f)=%2.1f\n",
                xw,
                yw,
                zw,
                val0,
                xwo,
                ywo,
                zwo,
                val_outside,
                xwi,
                ywi,
                zwi,
                val_inside);
        fprintf(stdout,
                "v %d I30 intensity term:  (%2.3f, %2.3f, %2.3f), "
                "targ=%2.1f, delV=%2.1f, delI=%2.0f, dot=%2.1f\n",
                vno,
                dx,
                dy,
                dz,
                v->val,
                delV,
                delI,
                dx * nx + dy * ny + dz * nz);
      }
      v->dx += dx;
      v->dy += dy;
      v->dz += dz;
    }

    // now compute flash5
    /* compute intensity gradient using smoothed volume */
    if (!FZERO(weight5)) {
      MRISsurfaceRASToVoxelCached(mris, mri_5, x, y, z, &xw, &yw, &zw);
      MRIsampleVolume(mri_5, xw, yw, zw, &val0);
      if (parms->grad_dir == 0) {
        double dist, val, step_size;
        int n;

        step_size = MIN(sigma / 2, MIN(mri_5->xsize, MIN(mri_5->ysize, mri_5->zsize)) * 0.5);
        ktotal_inside = ktotal_outside = 0.0;
        for (n = 0, val_outside = val_inside = 0.0, dist = step_size; dist <= 2 * sigma; dist += step_size, n++) {
          k = exp(-dist * dist / (2 * sigma * sigma));

          xw = x + dist * nx;
          yw = y + dist * ny;
          zw = z + dist * nz;
          if (mri_interior) {
            MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
            MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
          }
          if (mri_interior == NULL || interior < .9) {
            MRISsurfaceRASToVoxelCached(mris, mri_5, xw, yw, zw, &xw, &yw, &zw);
            MRIsampleVolume(mri_5, xw, yw, zw, &val);
            val_outside += k * val;
            ktotal_outside += k;
          }

          xw = x - dist * nx;
          yw = y - dist * ny;
          zw = z - dist * nz;
          if (mri_interior) {
            MRISsurfaceRASToVoxelCached(mris, mri_interior, xw, yw, zw, &xvi, &yvi, &zvi);
            MRIsampleVolume(mri_interior, xvi, yvi, zvi, &interior);
          }
          if (mri_interior == NULL || interior > 0) {
            MRISsurfaceRASToVoxelCached(mris, mri_5, xw, yw, zw, &xw, &yw, &zw);
            MRIsampleVolume(mri_5, xw, yw, zw, &val);
            val_inside += k * val;
            ktotal_inside += k;
          }
        }
        if (ktotal_inside > 0) {
          val_inside /= (double)ktotal_inside;
        }
        if (ktotal_outside > 0) {
          val_outside /= (double)ktotal_outside;
        }
      }
      else  // don't compute gradient - assume
      {
        val_outside = parms->grad_dir;
        val_inside = -parms->grad_dir;
      }

      delV = v->valbak - val0;
      delI = (val_outside - val_inside) / 2.0;

      if (!FZERO(delI)) {
        delI /= fabs(delI);
      }
      else {
        delI = 0;
      }

      if (delV > 5) {
        delV = 5;
      }
      else if (delV < -5) {
        delV = -5;
      }

      del = l_intensity * delV * delI;

      dx = nx * del * weight5;
      dy = ny * del * weight5;
      dz = nz * del * weight5;

      if (vno == Gdiag_no && !FZERO(weight5)) {
        double xwi, ywi, zwi, xwo, ywo, zwo;

        x = v->x;
        y = v->y;
        z = v->z;

        /* sample outward from surface */
        xw = x + nx;
        yw = y + ny;
        zw = z + nz;
        MRISsurfaceRASToVoxelCached(mris, mri_30, xw, yw, zw, &xwo, &ywo, &zwo);
        /* sample inward from surface */
        xw = x - nx;
        yw = y - ny;
        zw = z - nz;
        MRISsurfaceRASToVoxelCached(mris, mri_30, xw, yw, zw, &xwi, &ywi, &zwi);
        MRISsurfaceRASToVoxelCached(mris, mri_30, x, y, z, &xw, &yw, &zw);
        fprintf(stdout,
                "I5(%2.1f,%2.1f,%2.1f)=%2.1f, Io(%2.1f,%2.1f,%2.1f)=%2.1f, "
                "Ii(%2.1f,%2.1f,%2.1f)=%2.1f\n",
                xw,
                yw,
                zw,
                val0,
                xwo,
                ywo,
                zwo,
                val_outside,
                xwi,
                ywi,
                zwi,
                val_inside);
        fprintf(stdout,
                "v %d I5 intensity term:   (%2.3f, %2.3f, %2.3f), "
                "delV=%2.1f, delI=%2.0f, dot=%2.1f\n",
                vno,
                dx,
                dy,
                dz,
                delV,
                delI,
                dx * nx + dy * ny + dz * nz);
      }
      v->dx += dx;
      v->dy += dy;
      v->dz += dz;  // add flash5 component
    }
  }

  if (mri_interior) {
    MRIfree(&mri_interior);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Compute the effects of the gradient of the distance term
  ------------------------------------------------------*/
int mrisComputeDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  if (!(mris->dist_alloced_flags & 1)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceTerm","should have computed distances already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeVertexDistances(mris);
    }
  }
  if (!(mris->dist_alloced_flags & 2)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceTerm","should have computed dist_origs already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeOriginalVertexDistances(mris);
    }
  }

  float l_dist, scale, norm;
  int vno, tno;
  int diag_vno1, diag_vno2;
  char *cp;
  VECTOR *v_y[_MAX_FS_THREADS], *v_delta[_MAX_FS_THREADS], *v_n[_MAX_FS_THREADS];

  if ((cp = getenv("VDIAG1")) != NULL) {
    diag_vno1 = atoi(cp);
  }
  else {
    diag_vno1 = -1;
  }
  if ((cp = getenv("VDIAG2")) != NULL) {
    diag_vno2 = atoi(cp);
  }
  else {
    diag_vno2 = -1;
  }

  if (!FZERO(parms->l_nldist)) {
    mrisComputeNonlinearDistanceTerm(mris, parms);
  }

  l_dist = parms->l_dist;
  if (DZERO(l_dist)) {
    return (NO_ERROR);
  }

  norm = 1.0f / mris->avg_nbrs;

#if METRIC_SCALE
  if (mris->patch) {
    scale = 1.0f;
  }
  else if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_SPHERE) {
    scale = sqrt(mris->orig_area / mris->total_area);
  }
  else
    scale = mris->neg_area < mris->total_area ? sqrt(mris->orig_area / (mris->total_area - mris->neg_area))
                                              : sqrt(mris->orig_area / mris->total_area);
#else
  scale = 1.0f;
#endif

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "distance scale = %2.3f\n", scale);
  }
  for (tno = 0; tno < _MAX_FS_THREADS; tno++) {
    v_n[tno] = VectorAlloc(3, MATRIX_REAL);
    v_y[tno] = VectorAlloc(3, MATRIX_REAL);
    v_delta[tno] = VectorAlloc(3, MATRIX_REAL);
  }
// need to make v_n etc. into arrays and use tids
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(assume_reproducible) schedule(static, 1)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    int vnum, n;
    float d0, dt, delta, nc, vsmooth = 1.0;

#ifdef HAVE_OPENMP
    // thread ID
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    vnum = vt->vtotal;
    if (v->ripflag || vnum <= 0) continue;

    if (v->border) DiagBreak();

    V3_CLEAR(v_delta[tid]);
    VECTOR_LOAD(v_n[tid], v->nx, v->ny, v->nz);

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "computing distance term for v %d @ (%2.2f, %2.2f, %2.2f)\n", vno, v->x, v->y, v->z);

    for (n = 0; n < vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) continue;

      float const dist_orig_n = !v->dist_orig ? 0.0 : v->dist_orig[n];

      d0 = dist_orig_n / scale;
      dt = v->dist[n];
      delta = dt - d0;
      VECTOR_LOAD(v_y[tid], vn->x - v->x, vn->y - v->y, vn->z - v->z);
      if ((V3_LEN_IS_ZERO(v_y[tid]))) continue;

      V3_NORMALIZE(v_y[tid], v_y[tid]); /* make it a unit vector */
      V3_SCALAR_MUL(v_y[tid], delta, v_y[tid]);
      V3_ADD(v_y[tid], v_delta[tid], v_delta[tid]);
      if (vno == Gdiag_no && Gdiag & DIAG_SHOW & 0)  // deverbosified by dng
        fprintf(stdout,
                "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
                "d0 %2.2f, dt %2.2f, "
                "delta %2.3f\n\tdx=%2.3f, %2.3f, %2.3f)\n",
                n,
                vt->v[n],
                vn->x,
                vn->y,
                vn->z,
                d0,
                dt,
                delta,
                V3_X(v_y[tid]),
                V3_Y(v_y[tid]),
                V3_Z(v_y[tid]));
      if ((vno == diag_vno1 && vt->v[n] == diag_vno2) || (vno == diag_vno2 && vt->v[n] == diag_vno1))
        printf(
            "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
            "d0 %2.2f, dt %2.2f, "
            "delta %2.3f\n\ty=%2.3f, %2.3f, %2.3f)\n",
            n,
            vt->v[n],
            vn->x,
            vn->y,
            vn->z,
            d0,
            dt,
            delta,
            V3_X(v_y[tid]),
            V3_Y(v_y[tid]),
            V3_Z(v_y[tid]));
    }

    V3_SCALAR_MUL(v_delta[tid], norm, v_delta[tid]);

    if ((vno == Gdiag_no || vno == diag_vno1 || vno == diag_vno1) && Gdiag & DIAG_SHOW)
      fprintf(
          stdout, "total delta=(%2.3f, %2.3f, %2.3f)\n", V3_X(v_delta[tid]), V3_Y(v_delta[tid]), V3_Z(v_delta[tid]));
    /* take out normal component */
    nc = V3_DOT(v_n[tid], v_delta[tid]);
    V3_SCALAR_MUL(v_n[tid], -nc, v_n[tid]);
    V3_ADD(v_delta[tid], v_n[tid], v_delta[tid]);

    if (parms->vsmoothness) vsmooth = (1.0 - parms->vsmoothness[vno]);

    v->dx += l_dist * V3_X(v_delta[tid]);
    v->dy += l_dist * V3_Y(v_delta[tid]);
    v->dz += l_dist * V3_Z(v_delta[tid]);
    if (vno == Gdiag_no || vno == diag_vno1 || vno == diag_vno1)
      fprintf(stdout, "v %d, distance term: (%2.3f, %2.3f, %2.3f)\n", vno, v->dx, v->dy, v->dz);
    if (Gdiag_no == vno) {
      FILE *fp;
      char fname[STRLEN];
      
      int i;
      static int iter = 0;

      sprintf(fname, "v%d_dist_%04d.log", Gdiag_no, iter++);
      fp = fopen(fname, "w");
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];
      VERTEX          const * const v  = &mris->vertices         [Gdiag_no];
      for (i = 0; i < vt->vtotal; i++) {
        float const dist_orig_i = !v->dist_orig ? 0.0 : v->dist_orig[i];
        fprintf(fp,
                "%03d: %05d, %f   %f   %f\n",
                i,
                vt->v[i],
                dist_orig_i,
                v->dist[i],
                v->dist[i] - dist_orig_i / scale);
      }
      fclose(fp);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  for (tno = 0; tno < _MAX_FS_THREADS; tno++) {
    VectorFree(&v_n[tno]);
    VectorFree(&v_y[tno]);
    VectorFree(&v_delta[tno]);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Compute the effects of the gradient of the distance term
  ------------------------------------------------------*/
int mrisComputeNonlinearDistanceTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VECTOR *v_y, *v_delta, *v_n;
  float l_dist, d0, dt, delta, nc, scale, norm, ratio;
  int vno, n, vnum;

  l_dist = parms->l_nldist;
  if (FZERO(l_dist)) {
    return (NO_ERROR);
  }

  v_n = VectorAlloc(3, MATRIX_REAL);
  v_y = VectorAlloc(3, MATRIX_REAL);
  v_delta = VectorAlloc(3, MATRIX_REAL);
  norm = 1.0f / mris->avg_nbrs;

#if METRIC_SCALE
  if (mris->patch) {
    scale = 1.0f;
  }
  else if (mris->status == MRIS_PARAMETERIZED_SPHERE) {
    scale = sqrt(mris->orig_area / mris->total_area);
  }
  else
    scale = mris->neg_area < mris->total_area ? sqrt(mris->orig_area / (mris->total_area - mris->neg_area))
                                              : sqrt(mris->orig_area / mris->total_area);
#else
  scale = 1.0f;
#endif
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "distance scale = %2.3f\n", scale);
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    vnum = vt->vtotal;
    if (v->ripflag || vnum <= 0) {
      continue;
    }

    if (v->border) {
      DiagBreak();
    }
    V3_CLEAR(v_delta);
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz);

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "computing distance term for v %d @ (%2.2f, %2.2f, %2.2f)\n", vno, v->x, v->y, v->z);

    for (n = 0; n < vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      d0 = v->dist_orig[n];
      dt = scale * v->dist[n];
      delta = dt - d0;
      VECTOR_LOAD(v_y, vn->x - v->x, vn->y - v->y, vn->z - v->z);
      if ((V3_LEN_IS_ZERO(v_y))) {
        continue;
      }
      V3_NORMALIZE(v_y, v_y); /* make it a unit vector */
      if (!FZERO(d0)) {
        ratio = dt / d0;
        delta *= 1 / (1 + exp(-1 * ratio));
      }
      V3_SCALAR_MUL(v_y, delta, v_y);
      V3_ADD(v_y, v_delta, v_delta);
      if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
        fprintf(stdout,
                "nbr %d (%6.6d) @ (%2.2f, %2.2f, %2.2f), "
                "d0 %2.2f, dt %2.2f, "
                "delta %2.3f\n\ty=%2.3f, %2.3f, %2.3f)\n",
                n,
                vt->v[n],
                vn->x,
                vn->y,
                vn->z,
                d0,
                dt,
                delta,
                V3_X(v_y),
                V3_Y(v_y),
                V3_Z(v_y));
    }

    V3_SCALAR_MUL(v_delta, norm, v_delta);

    if (vno == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, "total delta=(%2.3f, %2.3f, %2.3f)\n", V3_X(v_delta), V3_Y(v_delta), V3_Z(v_delta));
    /* take out normal component */
    nc = V3_DOT(v_n, v_delta);
    V3_SCALAR_MUL(v_n, -nc, v_n);
    V3_ADD(v_delta, v_n, v_delta);

    v->dx += l_dist * V3_X(v_delta);
    v->dy += l_dist * V3_Y(v_delta);
    v->dz += l_dist * V3_Z(v_delta);
    if (vno == Gdiag_no) fprintf(stdout, "v %d, distance term: (%2.3f, %2.3f, %2.3f)\n", vno, v->dx, v->dy, v->dz);
  }

  VectorFree(&v_n);
  VectorFree(&v_y);
  VectorFree(&v_delta);

  return (NO_ERROR);
}


int mrisComputeIntensityGradientTerm(MRI_SURFACE *mris, double l_grad, MRI *mri_brain, MRI *mri_smooth)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz;
  double val0, mag0, xw, yw, zw, del, mag_outside, mag_inside, delI, delV, dx, dy, dz, val_outside, val_inside,
      val_dist, dn, xw1, yw1, zw1;

  if (FZERO(l_grad)) {
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

    /*
      sample intensity value and derivative in normal direction
      at current point.
    */
    x = v->x + v->nx;
    y = v->y + v->ny;
    z = v->z + v->nz;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw1, &yw1, &zw1);
    x = v->x;
    y = v->y;
    z = v->z;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    nx = xw1 - xw;
    ny = yw1 - yw;
    nz = zw1 - zw;
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val0);
    mag0 = sqrt(dx * dx + dy * dy + dz * dz);
    MRIsampleVolumeDerivative(mri_smooth, xw, yw, zw, nx, ny, nz, &dn);

    /* compute intensity gradient using smoothed volume */

    /* sample outward from surface */
    xw = x + nx;
    yw = y + ny;
    zw = z + nz;
    MRISsurfaceRASToVoxelCached(mris, mri_smooth, xw, yw, zw, &xw, &yw, &zw);
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag_outside = sqrt(dx * dx + dy * dy + dz * dz);
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_outside);

    /* sample inward from surface */
    xw = x - nx;
    yw = y - ny;
    zw = z - nz;
    MRISsurfaceRASToVoxelCached(mris, mri_smooth, xw, yw, zw, &xw, &yw, &zw);
    MRIsampleVolumeGradient(mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag_inside = sqrt(dx * dx + dy * dy + dz * dz);
    MRIsampleVolume(mri_smooth, xw, yw, zw, &val_inside);

    if (mag_outside > mag_inside) /* gradients suggest edge is outwards */
    {
      val_dist = val_outside - v->val;
    }
    else /* gradients suggest edge is inwards */
    {
      val_dist = val_inside - v->val;
    }

#if 0
    /* diminish the effect of gradients that are in locations whose
       intensity values are far from the target.
    */
    val_dist = 1 / (1 + val_dist*val_dist) ;
#else
    /* check to make sure that gradient is changing in right direction */
    val_dist = 1;
    /* is this right??? Used to be val0 > v->val, what about dn < 0?? */
    if (val0 > v->val) /* current point is brighter than target */
    {
      /* dn > 0 implies should move in, but gradient mag points out */
      if (((mag_inside < mag_outside) && dn > 0) || ((mag_inside > mag_outside) && dn < 0)) {
        val_dist = 0;
      }
    }
    else /* current point is dimmer than target */
    {
      /* dn > 0 implies should move out, but gradient mag points in */
      if (((mag_inside > mag_outside) && dn > 0) || ((mag_inside < mag_outside) && dn < 0)) {
        val_dist = 0;
      }
    }
#endif

    delV = 1.0f /*v->mean - mag0*/;
    delI = (mag_outside - mag_inside) / 2.0;
#if 1
    if (!FZERO(delI)) {
      delI /= fabs(delI);
    }
#endif
    del = val_dist * l_grad * delV * delI;
    dx = v->nx * del;
    dy = v->ny * del;
    dz = v->nz * del;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d intensity gradient term: (%2.3f, %2.3f, %2.3f) "
              "(mag = %2.1f, [%2.1f,%2.1f])\n",
              vno,
              v->dx,
              v->dy,
              v->dz,
              mag0,
              mag_inside,
              mag_outside);
  }

  return (NO_ERROR);
}


double MRIScomputeCorrelationError(MRI_SURFACE *mris, MRI_SP *mrisp_template, int fno)
{
  INTEGRATION_PARMS parms;
  float error;

  if (!mrisp_template) {
    return (0.0);
  }

  memset(&parms, 0, sizeof(parms));
  parms.mrisp_template = mrisp_template;
  parms.l_corr = 1.0f;
  parms.frame_no = fno;
  error = mrisComputeCorrelationError(mris, &parms, 1);
  return (sqrt(error / (double)MRISvalidVertices(mris)));
}

double mrisComputeCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds) {
    return mrisComputeCorrelationErrorTraceable(mris, parms, use_stds, false);
}

double mrisComputeCorrelationErrorTraceable(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace)
{
  float l_corr;

  l_corr = parms->l_corr + parms->l_pcorr; /* only one will be nonzero */
  if (FZERO(l_corr)) {
    return (0.0);
  }

  double sse = 0.0;
  
#ifdef BEVIN_MRISCOMPUTECORRELATIONERROR_REPRODUCIBLE

  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse  ROMP_PARTIALSUM(0)
    
#else
  int vno;

  ROMP_PF_begin         // Important during mris_register
 
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse)
#endif

  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif
    
    bool const vertexTrace = trace && (vno == 0);
    
    VERTEX *v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      ROMP_PF_continue;
    }

    double src, target, delta, std;
    float x, y, z;

    x = v->x;
    y = v->y;
    z = v->z;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0, vertexTrace) ;
#else
    src = v->curv;
#endif
    target = MRISPfunctionValTraceable(parms->mrisp_template, mris, x, y, z, parms->frame_no, vertexTrace);
#define DEFAULT_STD 4.0f
#define DISABLE_STDS 0
#if DISABLE_STDS
    std = 1.0f;
#else
    std = MRISPfunctionValTraceable(parms->mrisp_template, mris, x, y, z, parms->frame_no + 1, vertexTrace);
    std = sqrt(std);
    if (FZERO(std)) {
      std = DEFAULT_STD /*FSMALL*/;
    }
    if (!use_stds) {
      std = 1.0f;
    }
#endif
    delta = (src - target) / std;
    if (!isfinite(target) || !isfinite(delta)) {
      DiagBreak();
    }
    if (parms->geometry_error) {
      parms->geometry_error[vno] = (delta * delta);
    }
    if (parms->abs_norm) {
      sse += fabs(delta);
    }
    else {
      sse += delta * delta;
    }
#ifdef BEVIN_MRISCOMPUTECORRELATIONERROR_REPRODUCIBLE

    #undef sse
  #include "romp_for_end.h"

#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

  return (sse);
}

double mrisComputeVectorCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds)
{
  double src, target, sse, delta, std;
  VERTEX *v;
  int n, vno, fno;
  float x, y, z, l_corr;

  double *vals, *corrs, *sses;
  int nframes, *frames, *ind;

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    nframes++;
  }
  if (!nframes) {
    return 0.0;
  }

  corrs = (double *)malloc(nframes * sizeof(double));
  sses = (double *)malloc(nframes * sizeof(double));
  ind = (int *)malloc(nframes * sizeof(int));

  vals = (double *)malloc(2 * nframes * sizeof(double)); /* include the variances */
  frames = (int *)malloc(2 * nframes * sizeof(int));
  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;
    frames[2 * nframes] = fno;
    frames[2 * nframes + 1] = fno + 1;
    ind[nframes] = n;
    corrs[nframes] = l_corr;
    nframes++;
  }

  memset(sses, 0, nframes * sizeof(double)); /* set the vals to zero */
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    x = v->x;
    y = v->y;
    z = v->z;
    /* get the template values (variance included) */
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x, y, z, frames, 2 * nframes, vals);
#define DEFAULT_STD 4.0f
#if DISABLE_STDS
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = 1.0f;
      delta = (src - target) / std;
      sses[n] += delta * delta;
    }
#else
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = sqrt(vals[2 * n + 1]);
      if (FZERO(std)) {
        std = DEFAULT_STD /*FSMALL*/; /* to be checked */
      }
      if (!use_stds) {
        std = 1.0f;
      }

      if (parms->fields[ind[n]].type) {
        std = MAX(0.01, std);
      }

      delta = (src - target) / std;

      // if(vno==0) fprintf(stderr,"[%d : %f , %f , %f ]",n, src,target,std);

      sses[n] += delta * delta;
    }
#endif
  }
  sse = 0.0f;
  for (n = 0; n < nframes; n++) {
    // fprintf(stderr,"(%d,%f,%f -> %f)\n",n,corrs[n],sses[n],corrs[n]*sses[n]);
    parms->fields[ind[n]].sse = sses[n];
    sse += corrs[n] * sses[n];
  }

  free(corrs);
  free(sses);
  free(ind);
  free(frames);
  free(vals);

  return (sse);

#if 0
  sse=0.0f; /* compiler warnings */
  for (n=0 ; n < parms->nfields ; n++)
  {

    l_corr = parms->fields[n].l_corr+parms->fields[n].l_pcorr ;

    if (FZERO(l_corr))
    {
      continue;
    }

    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }

      x = v->x ;
      y = v->y ;
      z = v->z ;
#if 0
      src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
      src = ((VALS_VP*)v->vp)->vals[n] ;
#endif
      target = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, fno) ;
#define DEFAULT_STD 4.0f
#define DISABLE_STDS 0
#if DISABLE_STDS
      std = 1.0f ;
#else
      std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno+1);
      std = sqrt(std) ;
      if (FZERO(std))
      {
        std = DEFAULT_STD /*FSMALL*/ ;
      }
      if (!use_stds)
      {
        std = 1.0f ;
      }
#endif
      delta = (src - target) / std ;
      if (!finite(target) || !finite(delta))
      {
        DiagBreak() ;
      }
      sse += l_corr * delta * delta ;
    }
  }
  return(sse) ;
#endif
}


int mrisComputeCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double du, dv, up1, um1, vp1, vm1, delta, src, target, mag, max_mag, l_corr;
  VERTEX *v;
  int vno, fno;
  float x, y, z, e1x, e1y, e1z, e2x, e2y, e2z, ux, uy, uz, vx, vy, vz, std, coef, vsmooth = 1.0;
  double d_dist = D_DIST * mris->avg_vertex_dist;

  l_corr = parms->l_corr;
  if (FZERO(l_corr)) {
    return (NO_ERROR);
  }
  fno = parms->frame_no;
  mrisComputeTangentPlanes(mris);
  max_mag = 0.0f;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRISPwrite(parms->mrisp_template, "temp.hipl");
    MRISPwrite(parms->mrisp, "srf.hipl");
  }
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
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
    src = v->curv;
#endif
    target = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, parms->frame_no);
    std = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, parms->frame_no + 1);
    std = sqrt(std);
    if (FZERO(std)) {
      std = DEFAULT_STD /*FSMALL*/;
    }
#if DISABLE_STDS
    std = 1.0f;
#endif
    delta = (target - src);
    if (parms->vsmoothness) {
      vsmooth = 1.0 + parms->vsmoothness[vno];
    }
    coef = vsmooth * delta * l_corr / std;

    /* now compute gradient of template w.r.t. a change in vertex position */

    /*
      sample the curvature functions along the tangent plane axes and
      compute the derivates using them.
    */
    ux = e1x * d_dist;
    uy = e1y * d_dist;
    uz = e1z * d_dist;
    vx = e2x * d_dist;
    vy = e2y * d_dist;
    vz = e2z * d_dist;

#if 0
    /* compute src term */
    up1 = MRISPfunctionVal(parms->mrisp, mris, x+ux, y+uy, z+uz, fno) ;
    um1 = MRISPfunctionVal(parms->mrisp, mris, x-ux, y-uy, z-uz, fno) ;
    vp1 = MRISPfunctionVal(parms->mrisp, mris, x+vx, y+vy, z+vz, fno) ;
    vm1 = MRISPfunctionVal(parms->mrisp, mris, x-vx, y-vy, z-vz, fno) ;
    du = (up1 - um1) / (2 * d_dist) ;
    dv = (vp1 - vm1) / (2 * d_dist) ;
    v->dx += coef * (du*e1x + dv*e2x) ;  /* in negative of grad. direction */
    v->dy += coef * (du*e1y + dv*e2y) ;
    v->dz += coef * (du*e1z + dv*e2z) ;
#endif

    /* compute target term */
    up1 = MRISPfunctionVal(parms->mrisp_template, mris, x + ux, y + uy, z + uz, fno);
    um1 = MRISPfunctionVal(parms->mrisp_template, mris, x - ux, y - uy, z - uz, fno);
    vp1 = MRISPfunctionVal(parms->mrisp_template, mris, x + vx, y + vy, z + vz, fno);
    vm1 = MRISPfunctionVal(parms->mrisp_template, mris, x - vx, y - vy, z - vz, fno);
    du = (up1 - um1) / (2 * d_dist);
    dv = (vp1 - vm1) / (2 * d_dist);
    v->dx -= coef * (du * e1x + dv * e2x);
    v->dy -= coef * (du * e1y + dv * e2y);
    v->dz -= coef * (du * e1z + dv * e2z);

    mag = sqrt(v->dx * v->dx + v->dy * v->dy + v->dz * v->dz);
    if (mag > max_mag) {
      max_mag = mag;
    }
    if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz)) {
      DiagBreak();
      ErrorExit(ERROR_BADPARM, "mrisComputeCorrelationTerm: delta is not finite at vno %d", vno);
    }
    if (vno == Gdiag_no)
      printf("l_corr(%d): dx = (%2.2f, %2.2f, %2.2f), vsmooth = %2.3f\n",
             vno,
             coef * (du * e1x + dv * e2x),
             coef * (du * e1y + dv * e2y),
             coef * (du * e1z + dv * e2z),
             vsmooth);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag);
  }

  return (NO_ERROR);
}


int mrisComputeVectorCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double du, dv, up1, um1, vp1, vm1, delta, src, target, mag, max_mag, l_corr;
  VERTEX *v;
  int vno, fno, n;
  float x, y, z, e1x, e1y, e1z, e2x, e2y, e2z, ux, uy, uz, vx, vy, vz, std;
  double *vals, *upvals, *umvals, *vpvals, *vmvals, *corrs;
  int nframes, *frames, *wv_frames, *ind;
  double d_dist = D_DIST * mris->avg_vertex_dist;

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr;
    if (FZERO(l_corr)) {
      continue;
    }
    nframes++;
  }
  if (!nframes) {
    return NO_ERROR; /* no frames */
  }

  corrs = (double *)malloc(nframes * sizeof(double)); /* correlation coefficients */
  ind = (int *)malloc(nframes * sizeof(int));

  vals = (double *)malloc(2 * nframes * sizeof(double)); /* include the variances */
  frames = (int *)malloc(2 * nframes * sizeof(int));     /* include the variances */
  upvals = (double *)malloc(nframes * sizeof(double));   /* without the variances */
  umvals = (double *)malloc(nframes * sizeof(double));
  vpvals = (double *)malloc(nframes * sizeof(double));
  vmvals = (double *)malloc(nframes * sizeof(double));
  wv_frames = (int *)malloc(nframes * sizeof(int));

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr;
    if (FZERO(l_corr)) {
      continue;
    }
    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;
    frames[2 * nframes] = fno;         /* mean field */
    frames[2 * nframes + 1] = fno + 1; /* variance field */
    wv_frames[nframes] = fno;
    ind[nframes] = n;
    corrs[nframes] = l_corr;
    nframes++;
  }

  mrisComputeTangentPlanes(mris);
  max_mag = 0.0f;

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
    target = MRISPfunctionVectorVals(parms->mrisp_template, mris, x, y, z, frames, 2 * nframes, vals);
    /* now compute gradient of template w.r.t. a change in vertex position */
    /*
      sample the fields along the tangent plane axes and
      compute the derivates using them.
    */
    ux = e1x * d_dist;
    uy = e1y * d_dist;
    uz = e1z * d_dist;
    vx = e2x * d_dist;
    vy = e2y * d_dist;
    vz = e2z * d_dist;
    /* compute target terms */
    up1 = MRISPfunctionVectorVals(parms->mrisp_template, mris, x + ux, y + uy, z + uz, wv_frames, nframes, upvals);
    um1 = MRISPfunctionVectorVals(parms->mrisp_template, mris, x - ux, y - uy, z - uz, wv_frames, nframes, umvals);
    vp1 = MRISPfunctionVectorVals(parms->mrisp_template, mris, x + vx, y + vy, z + vz, wv_frames, nframes, vpvals);
    vm1 = MRISPfunctionVectorVals(parms->mrisp_template, mris, x - vx, y - vy, z - vz, wv_frames, nframes, vmvals);

#if DISABLE_STDS
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = 1.0f;
      delta = corrs[n] * (target - src) / std;

      du = (upvals[n] - umvals[n]) / (2 * d_dist);
      dv = (vpvals[n] - vmvals[n]) / (2 * d_dist);

      v->dx -= delta * (du * e1x + dv * e2x);
      v->dy -= delta * (du * e1y + dv * e2y);
      v->dz -= delta * (du * e1z + dv * e2z);
    }
#else
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = sqrt(vals[2 * n + 1]);
      if (FZERO(std)) {
        std = DEFAULT_STD /*FSMALL*/; /* to be checked */
      }
      delta = corrs[n] * (target - src) / std;

      du = (upvals[n] - umvals[n]) / (2 * d_dist);
      dv = (vpvals[n] - vmvals[n]) / (2 * d_dist);

      v->dx -= delta * (du * e1x + dv * e2x);
      v->dy -= delta * (du * e1y + dv * e2y);
      v->dz -= delta * (du * e1z + dv * e2z);
    }
#endif
    mag = sqrt(v->dx * v->dx + v->dy * v->dy + v->dz * v->dz);
    if (mag > max_mag) {
      max_mag = mag;
    }
    if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz)) {
      DiagBreak();
      ErrorExit(ERROR_BADPARM, "mrisComputeVectorCorrelationTerm: delta is not finite at vno %d", vno);
    }
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag);
  }

  free(corrs);
  free(ind);
  free(vals);
  free(upvals);
  free(umvals);
  free(vpvals);
  free(vmvals);
  free(frames);
  free(wv_frames);

#if 0
  mrisComputeTangentPlanes(mris) ;
  max_mag = 0.0f ;
  for (n=0 ; n < parms->nfields ; n++)
  {

    l_corr = parms->fields[n].l_corr ;

    if (FZERO(l_corr))
    {
      continue;
    }

    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;

    max_mag = 0.0f ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }

      e1x = v->e1x ;
      e1y = v->e1y ;
      e1z = v->e1z ;
      e2x = v->e2x ;
      e2y = v->e2y ;
      e2z = v->e2z ;
      x = v->x ;
      y = v->y ;
      z = v->z ;
#if 0
      src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
      src = ((VALS_VP*)v->vp)->vals[n] ;
#endif
      target = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno);
      std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno+1);
      std = sqrt(std) ;
      if (FZERO(std))
      {
        std = DEFAULT_STD /*FSMALL*/ ;
      }
#if DISABLE_STDS
      std = 1.0f ;
#endif
      delta = (target-src) ;
      coef = delta * l_corr / std ;

      /* now compute gradient of template w.r.t. a change in vertex position */

      /*
        sample the curvature functions along the tangent plane axes and
        compute the derivates using them.
      */
      ux = e1x*d_dist ;
      uy = e1y*d_dist ;
      uz = e1z*d_dist ;
      vx = e2x*d_dist ;
      vy = e2y*d_dist ;
      vz = e2z*d_dist ;

#if 0
      /* compute src term */
      up1 = MRISPfunctionVal(parms->mrisp, mris, x+ux, y+uy, z+uz, fno) ;
      um1 = MRISPfunctionVal(parms->mrisp, mris, x-ux, y-uy, z-uz, fno) ;
      vp1 = MRISPfunctionVal(parms->mrisp, mris, x+vx, y+vy, z+vz, fno) ;
      vm1 = MRISPfunctionVal(parms->mrisp, mris, x-vx, y-vy, z-vz, fno) ;
      du = (up1 - um1) / (2 * d_dist) ;
      dv = (vp1 - vm1) / (2 * d_dist) ;
      v->dx += coef * (du*e1x + dv*e2x) ;  /* in negative of grad. direction */
      v->dy += coef * (du*e1y + dv*e2y) ;
      v->dz += coef * (du*e1z + dv*e2z) ;
#endif

      /* compute target term */
      up1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x+ux, y+uy, z+uz, fno);
      um1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x-ux, y-uy, z-uz, fno);
      vp1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x+vx, y+vy, z+vz, fno);
      vm1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x-vx, y-vy, z-vz, fno);
      du = (up1 - um1) / (2 * d_dist) ;
      dv = (vp1 - vm1) / (2 * d_dist) ;
      v->dx -= coef * (du*e1x + dv*e2x) ;
      v->dy -= coef * (du*e1y + dv*e2y) ;
      v->dz -= coef * (du*e1z + dv*e2z) ;

      mag = sqrt(v->dx*v->dx + v->dy*v->dy + v->dz*v->dz) ;
      if (mag > max_mag)
      {
        max_mag = mag ;
      }
      if (!finite(v->dx) || !finite(v->dy) || !finite(v->dz))
      {
        DiagBreak() ;
        ErrorExit
        (ERROR_BADPARM,
         "mrisComputeVectorCorrelationTerm: delta is not finite at vno %d",
         vno) ;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag) ;
  }
#endif

  return (NO_ERROR);
}


#define D_ANGLE RADIANS(0.25)

int mrisComputePolarCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double ap1, am1, da, bp1, bm1, db, gp1, gm1, dg, delta, src, target, mag, max_mag;
  VERTEX *v;
  int vno, fno;
  float x, y, z, std, coef, dx, dy, dz, nv, r;

  if (FZERO(parms->l_pcorr)) {
    return (NO_ERROR);
  }
  fno = parms->frame_no;
  max_mag = 0.0f;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRISPwrite(parms->mrisp_template, "temp.hipl");
    /*    MRISPwrite(parms->mrisp, "srf.hipl") ;*/
  }
  mris->gamma = mris->beta = mris->alpha = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
    src = v->curv;
#endif
    target = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, parms->frame_no);
    std = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, parms->frame_no + 1);
    std = sqrt(std);
    if (FZERO(std)) {
      std = DEFAULT_STD /*FSMALL*/;
    }
#if DISABLE_STDS
    std = 1.0f;
#endif
    delta = (target - src);
    coef = delta / std;

    /* now compute gradient of template w.r.t. a change in vertex position */

    /*
      compute the gradient by using differential rotations in the 3
      rotational directions using the associated skew-symmetric
      differential rotation matrices
    */
    /* compute alpha term - differential rotation around z axis */
    dx = y * D_ANGLE;
    dy = -x * D_ANGLE;
    dz = 0;
    am1 = MRISPfunctionVal(parms->mrisp_template, mris, x - dx, y - dy, z - dz, 0);
    ap1 = MRISPfunctionVal(parms->mrisp_template, mris, x + dx, y + dy, z + dz, 0);
    da = (ap1 - am1) / (2 * D_ANGLE);

    /* compute beta term - differential rotation around y axis */
    dx = -z * D_ANGLE;
    dy = 0;
    dz = x * D_ANGLE;
    bm1 = MRISPfunctionVal(parms->mrisp_template, mris, x - dx, y - dy, z - dz, 0);
    bp1 = MRISPfunctionVal(parms->mrisp_template, mris, x + dx, y + dy, z + dz, 0);
    db = (bp1 - bm1) / (2 * D_ANGLE);

    /* compute gamma term - differential rotation around x axis */
    dx = 0;
    dy = -z * D_ANGLE;
    dz = y * D_ANGLE;
    gm1 = MRISPfunctionVal(parms->mrisp_template, mris, x - dx, y - dy, z - dz, 0);
    gp1 = MRISPfunctionVal(parms->mrisp_template, mris, x + dx, y + dy, z + dz, 0);
    dg = (gp1 - gm1) / (2 * D_ANGLE);

    mris->gamma -= coef * dg; /* around x-axis */
    mris->beta -= coef * db;  /* around y-axis */
    mris->alpha -= coef * da; /* around z-axis */

    mag = sqrt(v->dx * v->dx + v->dy * v->dy + v->dz * v->dz);
    if (mag > max_mag) {
      max_mag = mag;
    }
    if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz)) {
      DiagBreak();
      ErrorExit(ERROR_BADPARM, "mrisComputePolarCorrelationTerm: delta is not finite at vno %d", vno);
    }
  }

  nv = MRISvalidVertices(mris);
  r = mris->radius;
  mris->alpha /= nv;
  mris->beta /= nv;
  mris->gamma /= nv;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->dx = r * mris->alpha;
    v->dy = r * mris->beta;
    v->dz = r * mris->gamma;
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "max gradient magnitude = %2.5f\n", max_mag);
  }

  return (NO_ERROR);
}

int mrisComputePolarVectorCorrelationTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double da, db, dg, delta, src, target, l_pcorr;
  VERTEX *v;
  int vno, fno, n;
  float x, y, z, std, dx, dy, dz, nv, r;

  double *vals, *dalpha, *dbeta, *dgamma, *apvals, *amvals, *bpvals, *bmvals, *gpvals, *gmvals, *corrs;
  int nframes, *frames, *wv_frames, *ind;

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_pcorr = parms->fields[n].l_pcorr;
    if (FZERO(l_pcorr)) {
      continue;
    }
    nframes++;
  }
  if (!nframes) {
    return NO_ERROR; /* no frames */
  }

  corrs = (double *)malloc(nframes * sizeof(double));  /* correlation coefficients */
  dalpha = (double *)malloc(nframes * sizeof(double)); /* gradient coefficients */
  dbeta = (double *)malloc(nframes * sizeof(double));  /* gradient coefficients */
  dgamma = (double *)malloc(nframes * sizeof(double)); /* gradient coefficients */
  ind = (int *)malloc(nframes * sizeof(int));

  vals = (double *)malloc(2 * nframes * sizeof(double)); /* include the variances */
  frames = (int *)malloc(2 * nframes * sizeof(int));     /* include the variances */
  apvals = (double *)malloc(nframes * sizeof(double));   /* without the variances */
  amvals = (double *)malloc(nframes * sizeof(double));
  bpvals = (double *)malloc(nframes * sizeof(double));
  bmvals = (double *)malloc(nframes * sizeof(double));
  gpvals = (double *)malloc(nframes * sizeof(double));
  gmvals = (double *)malloc(nframes * sizeof(double));
  wv_frames = (int *)malloc(nframes * sizeof(int));

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_pcorr = parms->fields[n].l_pcorr;
    if (FZERO(l_pcorr)) {
      continue;
    }
    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;
    frames[2 * nframes] = fno;         /* mean field */
    frames[2 * nframes + 1] = fno + 1; /* variance field */
    wv_frames[nframes] = fno;
    ind[nframes] = n;
    corrs[nframes] = l_pcorr;
    nframes++;
  }

  memset(dalpha, 0, nframes * sizeof(double));
  memset(dbeta, 0, nframes * sizeof(double));
  memset(dgamma, 0, nframes * sizeof(double));
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    target = MRISPfunctionVectorVals(parms->mrisp_template, mris, x, y, z, frames, 2 * nframes, vals);
    /*
      compute the gradient by using differential rotations in the 3
      rotational directions using the associated skew-symmetric
      differential rotation matrices
    */
    /* compute alpha term - differential rotation around z axis */
    dx = y * D_ANGLE;
    dy = -x * D_ANGLE;
    dz = 0;
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x - dx, y - dy, z - dz, wv_frames, nframes, amvals);
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x + dx, y + dy, z + dz, wv_frames, nframes, apvals);

    /* compute beta term - differential rotation around y axis */
    dx = -z * D_ANGLE;
    dy = 0;
    dz = x * D_ANGLE;
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x - dx, y - dy, z - dz, wv_frames, nframes, bmvals);
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x + dx, y + dy, z + dz, wv_frames, nframes, bpvals);

    /* compute gamma term - differential rotation around x axis */
    dx = 0;
    dy = -z * D_ANGLE;
    dz = y * D_ANGLE;
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x - dx, y - dy, z - dz, wv_frames, nframes, gmvals);
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x + dx, y + dy, z + dz, wv_frames, nframes, gpvals);

#if DISABLE_STDS
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = 1.0f;
      delta = (target - src) / std;

      da = (apvals[n] - amvals[n]) / (2 * D_ANGLE);
      dalpha[n] += delta * da;

      db = (bpvals[n] - bmvals[n]) / (2 * D_ANGLE);
      dbeta[n] += delta * db;

      dg = (gpvals[n] - gmvals[n]) / (2 * D_ANGLE);
      dgamma[n] += delta * dg;
    }
#else
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = sqrt(vals[2 * n + 1]);
      if (FZERO(std)) {
        std = DEFAULT_STD /*FSMALL*/; /* to be checked */
      }
      delta = (target - src) / std;

      da = (apvals[n] - amvals[n]) / (2 * D_ANGLE);
      dalpha[n] += delta * da;

      db = (bpvals[n] - bmvals[n]) / (2 * D_ANGLE);
      dbeta[n] += delta * db;

      dg = (gpvals[n] - gmvals[n]) / (2 * D_ANGLE);
      dgamma[n] += delta * dg;
    }
#endif
  }

  mris->gamma = mris->beta = mris->alpha = 0;
  for (n = 0; n < nframes; n++) {
    mris->gamma -= corrs[n] * dgamma[n]; /* around x-axis */
    mris->beta -= corrs[n] * dbeta[n];   /* around y-axis */
    mris->alpha -= corrs[n] * dalpha[n]; /* around z-axis */
  }

  nv = MRISvalidVertices(mris);
  r = mris->radius;
  mris->alpha /= nv;
  mris->beta /= nv;
  mris->gamma /= nv;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    v->dx = r * mris->alpha;
    v->dy = r * mris->beta;
    v->dz = r * mris->gamma;
  }

  free(corrs);
  free(dalpha);
  free(dbeta);
  free(dgamma);
  free(ind);
  free(vals);
  free(frames);
  free(apvals);
  free(amvals);
  free(bpvals);
  free(bmvals);
  free(gpvals);
  free(gmvals);
  free(wv_frames);

#if 0

  for (n=0 ; n < parms->nfields ; n++)
  {

    l_pcorr = parms->fields[n].l_pcorr ;
    if (FZERO(l_pcorr))
    {
      continue;
    }

    fno = parms->fields[n].frame * IMAGES_PER_SURFACE ;

    mris->gamma = mris->beta = mris->alpha = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }

      x = v->x ;
      y = v->y ;
      z = v->z ;
#if 0
      src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
      src = ((VALS_VP*)v->vp)->vals[n] ;
#endif
      target = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno);
      std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno+1);
      std = sqrt(std) ;
      if (FZERO(std))
      {
        std = DEFAULT_STD /*FSMALL*/ ;
      }
#if DISABLE_STDS
      std = 1.0f ;
#endif
      delta = (target-src) ;
      coef = delta * l_pcorr / std ;

      /* now compute gradient of template w.r.t. a change in vertex position */

      /*
        compute the gradient by using differential rotations in the 3
        rotational directions using the associated skew-symmetric
        differential rotation matrices
      */
      /* compute alpha term - differential rotation around z axis */
      dx = y*D_ANGLE ;
      dy = -x*D_ANGLE ;
      dz = 0 ;
      am1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
      ap1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
      da = (ap1 - am1) / (2*D_ANGLE);

      /* compute beta term - differential rotation around y axis */
      dx = -z*D_ANGLE ;
      dy = 0 ;
      dz = x*D_ANGLE ;
      bm1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
      bp1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
      db = (bp1 - bm1) / (2*D_ANGLE);

      /* compute gamma term - differential rotation around x axis */
      dx = 0 ;
      dy = -z*D_ANGLE ;
      dz = y*D_ANGLE ;
      gm1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x-dx, y-dy, z-dz, 0) ;
      gp1 =
        MRISPfunctionVal(parms->mrisp_template, mris, x+dx, y+dy, z+dz, 0) ;
      dg = (gp1 - gm1) / (2*D_ANGLE);

      mris->gamma -= coef * dg ;   /* around x-axis */
      mris->beta  -= coef * db ;   /* around y-axis */
      mris->alpha -= coef * da ;   /* around z-axis */

    }
  }
  nv = MRISvalidVertices(mris) ;
  r = mris->radius ;
  mris->alpha /= nv ;
  mris->beta /= nv ;
  mris->gamma /= nv ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    v->dx = r*mris->alpha ;
    v->dy = r*mris->beta ;
    v->dz = r*mris->gamma ;
  }
#endif

  return (NO_ERROR);
}

int mrisComputeExpandwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap)
{
  int vno;
  double xw, yw, zw, x, y, z, val, dx, dy, dz;
  VERTEX *v;
  float min_val, max_val, target_val, delta;

  if (FZERO(l_expandwrap)) {
    return (NO_ERROR);
  }

  MRIvalRange(mri_brain, &min_val, &max_val);
  target_val = (min_val + max_val) / 2;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    target_val = v->val;
    x = v->x;
    y = v->y;
    z = v->z;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val);
    delta = (val - target_val);
    dx = -delta * v->nx * l_expandwrap;
    dy = -delta * v->ny * l_expandwrap;
    dz = -delta * v->nz * l_expandwrap;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d expandwrap term: (%2.3f, %2.3f, %2.3f), "
              "target %2.1f, MRI %2.1f, del=%2.1f, "
              "N=(%2.1f, %2.1f, %2.1f)\n",
              vno,
              dx,
              dy,
              dz,
              target_val,
              val,
              delta,
              v->nx,
              v->ny,
              v->nz);
  }
  return (NO_ERROR);
}

double mrisComputeExpandwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap, double target_radius)
{
  int vno;
  double xw, yw, zw, x, y, z, val, dx, dy, dz, sse, error, dist;
  VERTEX *v;
  float min_val, max_val, target_val, delta;

  if (FZERO(l_expandwrap)) {
    return (NO_ERROR);
  }

  mrisComputeSurfaceDimensions(mris);
  MRIvalRange(mri_brain, &min_val, &max_val);
  target_val = (min_val + max_val) / 2;
  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    target_val = v->val;
    x = v->x;
    y = v->y;
    z = v->z;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val);
    delta = (val - target_val);
    if (val < 0.25 * target_val) {
      dx = x - mris->xctr;
      dy = y - mris->yctr;
      dz = z - mris->zctr;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      error = (target_radius - dist);
      sse += error * error;
    }
    else {
      error = 0.0;
    }

    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d expandwrap error: %2.3f, target %2.1f, "
              "MRI %2.1f, del=%2.1f, \n",
              vno,
              error,
              target_val,
              val,
              delta);
  }
  return (sse);
}

int mrisComputeShrinkwrapTerm(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap)
{
  int vno;
  double xw, yw, zw, x, y, z, val, dx, dy, dz;
  VERTEX *v;
  float min_val, max_val, target_val, delta;

  if (FZERO(l_shrinkwrap)) {
    return (NO_ERROR);
  }

  MRIvalRange(mri_brain, &min_val, &max_val);
  target_val = (min_val + max_val) / 2;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    target_val = v->val;
    x = v->x;
    y = v->y;
    z = v->z;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val);

    delta = (val - target_val);

    dx = delta * v->nx * l_shrinkwrap;
    dy = delta * v->ny * l_shrinkwrap;
    dz = delta * v->nz * l_shrinkwrap;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d shrinkwrap term: (%2.3f, %2.3f, %2.3f), "
              "target %2.1f, MRI %2.1f, del=%2.1f, "
              "N=(%2.1f, %2.1f, %2.1f)\n",
              vno,
              dx,
              dy,
              dz,
              target_val,
              val,
              delta,
              v->nx,
              v->ny,
              v->nz);
  }
  return (NO_ERROR);
}

double mrisComputeShrinkwrapError(MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap)
{
#if 0
  static int iter = 100 ;
  int    vno ;
  double   xw, yw, zw, x, y, z, val ;
  VERTEX *v ;
  float  min_val, max_val, target_val, error ;
  double sse ;

  MRIvalRange(mri_brain, &min_val, &max_val) ;
  target_val = (min_val + max_val) / 2 ;
  sse = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    sse += iter ;
  }
  iter-- ;
  return(sse) ;
#else
  return (0.0);
#endif
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Apply a uniform outward expansion force.
  ------------------------------------------------------*/
int mrisComputeSphereTerm(MRI_SURFACE *mris, double l_sphere, float radius, int explode_flag)
{
  int vno;
  VERTEX *v;
  float r, x, y, z, x0, y0, z0;

  if (FZERO(l_sphere)) {
    return (0.0f);
  }

#if 1
  x0 = (mris->xlo + mris->xhi) / 2.0f;
  y0 = (mris->ylo + mris->yhi) / 2.0f;
  z0 = (mris->zlo + mris->zhi) / 2.0f;
#else
  x0 = mris->x0;
  y0 = mris->y0;
  z0 = mris->z0;
#endif

#if 0
  // make sure center is inside surface, otherwise move it
  {
    float dot, dx, dy, dz ;
    static int iter = 0 ;

    vno = MRISfindClosestVertex(mris, x0, y0, z0, &r, CURRENT_VERTICES) ;
    v = &mris->vertices[vno] ;
    dx = x0-v->x ;
    dy = y0-v->y ;
    dz = z0-v->z ;
    dot = v->nx*dx + v->ny*dy + v->nz*dz ;
    if (iter < 20 && dot > 0) // outside surface!
    {
      if (iter > 500)
      {
        DiagBreak() ;
      }
      printf("centroid (%2.1f, %2.1f, %2.1f) outside surface "
             "(dot = %2.1f, v %d = (%2.1f, %2.1f, %2.1f) n = "
             "(%2.1f, %2.1f, %2.1f)\n",
             x0, y0, z0, dot, vno, v->x, v->y, v->z, v->nx, v->ny, v->nz) ;
      x0 = v->x-0.5*v->nx ;
      y0 = v->y-(0.5*v->ny) ;
      z0 = v->z-(0.5*v->nz) ;
      mris->x0 = x0 ;
      mris->y0 = y0 ;
      mris->z0 = z0 ;
      print("moving centroid to (%2.1f, %2.1f, %2.1f)\n", x0, y0, z0) ;
      vno = MRISfindClosestVertex(mris, x0, y0, z0, &r, CURRENT_VERTICES) ;
      v = &mris->vertices[vno] ;
      dx = x0-v->x ;
      dy = y0-v->y ;
      dz = z0-v->z ;
      dot = v->nx*dx + v->ny*dy + v->nz*dz ;
      if (dot > 0)
      {
        printf("still outside!!!!!\n") ;
      }
    }
    iter++ ;
  }
#endif

  if (radius < 0) {
    radius = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      if (vno == Gdiag_no) {
        DiagBreak();
      }

      x = v->x - x0;
      y = v->y - y0;
      z = v->z - z0;
      r = sqrt(x * x + y * y + z * z);
      if (r > radius) {
        radius = r;
      }
    }
    radius++;
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x - x0;
    y = v->y - y0;
    z = v->z - z0;
    r = sqrt(x * x + y * y + z * z);
    if (FZERO(r)) {
      continue;
    }
    x /= r;
    y /= r;
    z /= r; /* normal direction */
    //    x = v->nx ; y = v->ny ; z = v->nz ;
    r = (radius - r) / radius;
    if (explode_flag)
      r = 1 ;

    if (vno == Gdiag_no && (r * l_sphere * x * v->nx + r * l_sphere * y * v->ny + r * l_sphere * z * v->nz < 0)) {
      DiagBreak();
    }
    v->dx += r * l_sphere * x;
    v->dy += r * l_sphere * y;
    v->dz += r * l_sphere * z;
    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d sphere   "
              " term: (%2.3f, %2.3f, %2.3f), r=%2.2f\n",
              vno,
              v->dx,
              v->dy,
              v->dz,
              r);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Apply a uniform outward expansion force.
  ------------------------------------------------------*/
int mrisComputeExpansionTerm(MRI_SURFACE *mris, double l_expand)
{
  int vno;
  VERTEX *v;

  if (FZERO(l_expand)) {
    return (0.0f);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    v->dx += l_expand * v->nx;
    v->dy += l_expand * v->ny;
    v->dz += l_expand * v->nz;
    if (vno == Gdiag_no)
      printf("v %d expansion term: (%2.3f, %2.3f, %2.3f)\n", vno, l_expand * v->nx, l_expand * v->ny, l_expand * v->nz);
  }

  return (NO_ERROR);
}

int mrisComputeBorderTerm(MRI_SURFACE *mris, double l_border)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, dist_scale;

  if (FZERO(l_border)) {
    return (NO_ERROR);
  }

  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }

  MRIScopyMarkedToMarked3(mris);
  MRISclearMarks(mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];

    if (v->ripflag == 0) {
      continue;
    }
    for (m = 0; m < vt->vtotal; m++) {
      VERTEX * const vn = &mris->vertices[vt->v[m]];
      vn->marked = 1;
    }
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->border && !v->neg) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    sx = sy = sz = 0.0;
    n = 0;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (vn->marked)  // move towards ripped vertices
                       // that represent the border of this region
      {
        sx += .5 * (vn->x - x);
        sy += .5 * (vn->y - y);
        sz += .5 * (vn->z - z);
        n++;
      }
    }
    if (n > 0) {
      sx = dist_scale * sx / n;
      sy = dist_scale * sy / n;
      sz = dist_scale * sz / n;
    }

    sx *= l_border;
    sy *= l_border;
    sz *= l_border;
    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d border term:         (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  MRIScopyMarked3ToMarked(mris);
  return (NO_ERROR);
}

int mrisComputeMaxSpringTerm(MRI_SURFACE *mris, double l_max_spring)
{
  int vno, n, m, m_max;
  float dx, dy, dz, x, y, z, dist_scale, dist, max_dist;

  if (FZERO(l_max_spring)) {
    return (NO_ERROR);
  }

  if (mris->patch) {
    dist_scale = l_max_spring;
  }
  else {
    dist_scale = l_max_spring * sqrt(mris->orig_area / mris->total_area);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    n = 0;
    m_max = 0;
    max_dist = 0;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      dx = (vn->x - x);
      dy = (vn->y - y);
      dz = (vn->z - z);
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (dist >= max_dist) {
        m_max = m;
        max_dist = dist;
      }
    }

    VERTEX const * const vn = &mris->vertices[vt->v[m_max]];
    dx = (vn->x - x);
    dy = (vn->y - y);
    dz = (vn->z - z);
    dx *= dist_scale;
    dy *= dist_scale;
    dz *= dist_scale;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no) {
      fprintf(stdout,
              "v %d spring max term:     "
              "(%2.3f, %2.3f, %2.3f)\n",
              vno,
              dx,
              dy,
              dz);
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
double mrisComputeThicknessParallelEnergy(MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms)
{
  int vno, max_vno;
  double sse_tparallel, max_inc;
  static int cno = 0;
  static double last_sse[MAXVERTICES];

  if (FZERO(l_thick_parallel)) {
    return (0.0);
  }
  if (cno == 0) {
    memset(last_sse, 0, sizeof(last_sse));
  }
  cno++;

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time

  max_inc = 0;
  max_vno = 0;
  sse_tparallel = 0.0;
  ROMP_PF_begin
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental) reduction(+:sse_tparallel)
  // endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    
    double sse;

    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();

    sse = mrisSampleParallelEnergy(mris, vno, parms, v->x, v->y, v->z);
    if ((vno < MAXVERTICES) && (sse > last_sse[vno] && cno > 1 && vno == Gdiag_no)) DiagBreak();

    if ((vno < MAXVERTICES) && (sse > last_sse[vno] && cno > 1)) {
      if (sse - last_sse[vno] > max_inc) {
        max_inc = sse - last_sse[vno];
        max_vno = vno;
      }
      DiagBreak();
    }

    if (vno < MAXVERTICES) last_sse[vno] = sse;
    sse_tparallel += sse;
    if (vno == Gdiag_no) {
      printf("E_parallel: vno = %d, E = %f\n", vno, sse);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  sse_tparallel /= 2;
  return (sse_tparallel);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeAshburnerTriangleEnergy(MRI_SURFACE *mris,
                                                 double l_ashburner_triangle,
                                                 INTEGRATION_PARMS *parms)
{
  int vno;
  double sse_ashburner;

  if (FZERO(l_ashburner_triangle)) return (0.0);

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time
  sse_ashburner = 0.0;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse_ashburner)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    double sse;
    VERTEX *v;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) ROMP_PF_continue;

    sse = mrisSampleAshburnerTriangleEnergy(mris, vno, parms, v->x, v->y, v->z);
    if (sse < 0) DiagBreak();

    sse_ashburner += sse;
    if (vno == Gdiag_no) printf("E_ash_triangle: vno = %d, E = %f\n", vno, sse);
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  sse_ashburner /= 2;
  return (sse_ashburner);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeThicknessMinimizationTerm(MRI_SURFACE *mris, double l_thick_min, INTEGRATION_PARMS *parms)
{
  int vno;

  if (FZERO(l_thick_min)) return (0.0);

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    float dE_de1, dE_de2, e1p, e1m, e2p, e2m;
    float e1x, e1y, e1z, e2x, e2y, e2z, norm, dx, dy, dz, E0, E1;
    VERTEX *v;
    double d_dist = D_DIST * mris->avg_vertex_dist;

    v = &mris->vertices[vno];
    if (v->ripflag) ROMP_PF_continue;

    if (vno == Gdiag_no) DiagBreak();

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;

    /*
      sample the coordinate functions along the tangent plane axes and
      compute the derivates using them.
    */
    e1p = mrisSampleMinimizationEnergy(mris, v, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleMinimizationEnergy(mris, v, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleMinimizationEnergy(mris, v, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleMinimizationEnergy(mris, v, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
    dE_de1 = (e1p - e1m) / (2 * d_dist);
    dE_de2 = (e2p - e2m) / (2 * d_dist);

    norm = sqrt(dE_de1 * dE_de1 + dE_de2 * dE_de2);
    if (norm > 1) {
      dE_de1 /= norm;
      dE_de2 /= norm;
    }

    dx = -l_thick_min * (dE_de1 * e1x + dE_de2 * e2x);
    dy = -l_thick_min * (dE_de1 * e1y + dE_de2 * e2y);
    dz = -l_thick_min * (dE_de1 * e1z + dE_de2 * e2z);
    E0 = mrisSampleMinimizationEnergy(mris, v, parms, v->x, v->y, v->z);
    E1 = mrisSampleMinimizationEnergy(
        mris, v, parms, v->x + parms->dt * dx, v->y + parms->dt * dy, v->z + parms->dt * dz);

    if (E1 > E0) {
      double E2;

      DiagBreak();
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      E2 = mrisSampleMinimizationEnergy(
          mris, v, parms, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz);
      if (E2 < E0) {
        dx *= -1;
        dy *= -1;
        dz *= -1;
      }
      else {
        dx = dy = dz = 0;
      }
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (Gdiag_no == vno) {
      printf(
          "l_thick_min: v %d: E0 %2.3f, E1 %2.3f, dE = (%2.2f, %2.2f), e1 = (%2.1f, %2.1f, %2.1f), e2 = (%2.1f, %2.1f, "
          "%2.1f), DX = (%2.2f, %2.2f, %2.2f)\n",

          vno,
          sqrt(E0),
          sqrt(E1),
          dE_de1,
          dE_de2,
          e1x,
          e1y,
          e1z,
          e2x,
          e2y,
          e2z,
          l_thick_min * (dE_de1 * e1x + dE_de2 * e2x),
          l_thick_min * (dE_de1 * e1y + dE_de2 * e2y),
          l_thick_min * (dE_de1 * e1z + dE_de2 * e2z));
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeThicknessParallelTerm(MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms)
{
  int vno;

  if (FZERO(l_thick_parallel)) return (0.0);

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time

  ROMP_PF_begin
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental)
  // endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    float dE_de1, dE_de2, e1p, e1m, e2p, e2m, dx, dy, dz;
    float e1x, e1y, e1z, e2x, e2y, e2z, norm, E0, E1;
    
    double d_dist = D_DIST * mris->avg_vertex_dist;

    VERTEX * const v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    if (v->ripflag) continue;

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;
    /*
      sample the coordinate functions along the tangent plane axes and
      compute the derivates using them.
    */
    e1p = mrisSampleParallelEnergy(mris, vno, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleParallelEnergy(mris, vno, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleParallelEnergy(mris, vno, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleParallelEnergy(mris, vno, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
    dE_de1 = (e1p - e1m) / (2 * d_dist);
    dE_de2 = (e2p - e2m) / (2 * d_dist);
    norm = sqrt(dE_de1 * dE_de1 + dE_de2 * dE_de2);
    if (norm > 1) {
      dE_de1 /= norm;
      dE_de2 /= norm;
    }

    dx = -l_thick_parallel * (dE_de1 * e1x + dE_de2 * e2x);
    dy = -l_thick_parallel * (dE_de1 * e1y + dE_de2 * e2y);
    dz = -l_thick_parallel * (dE_de1 * e1z + dE_de2 * e2z);
    E0 = mrisSampleParallelEnergy(mris, vno, parms, v->x, v->y, v->z);
    E1 = mrisSampleParallelEnergy(mris, vno, parms, v->x + parms->dt * dx, v->y + parms->dt * dy, v->z + parms->dt * dz);
#if 1
    if (E1 > E0) {
      double E2;
      DiagBreak();
      if (vno == Gdiag_no) DiagBreak();

      E2 =
          mrisSampleParallelEnergy(mris, vno, parms, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz);
      if (E2 < E0) {
        dx *= -1;
        dy *= -1;
        dz *= -1;
      }
      else {
        dx = dy = dz = 0;
      }
    }
#endif
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (Gdiag_no == vno) {
      printf(
          "l_thick_parallel: v %d: E0=%f, E1=%f, dE = (%2.3f, %2.3f), e1 = (%2.1f, %2.1f, %2.1f), e2 = (%2.1f, %2.1f, "
          "%2.1f), DX = (%2.3f, %2.3f, %2.3f)\n",

          vno,
          E0,
          E1,
          dE_de1,
          dE_de2,
          e1x,
          e1y,
          e1z,
          e2x,
          e2y,
          e2z,
          parms->dt * dx,
          parms->dt * dy,
          parms->dt * dz);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MAX_NORM 0.1
int mrisComputeThicknessNormalTerm(MRI_SURFACE *mris, double l_thick_normal, INTEGRATION_PARMS *parms)
{
  int vno;
  //  int     missed = 0 ;

  if (FZERO(l_thick_normal)) return (0.0);

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    float dE_de1, dE_de2, e1p, e1m, e2p, e2m, cx, cy, cz;
    float E0, E1, dx, dy, dz;
    float e1x, e1y, e1z, e2x, e2y, e2z, norm;
    VERTEX *v;
    double d_dist = D_DIST * mris->avg_vertex_dist;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;
    /*
      sample the coordinate functions along the tangent plane axes and
      compute the derivates using them.
    */
    e1p = mrisSampleNormalEnergy(mris, v, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleNormalEnergy(mris, v, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleNormalEnergy(mris, v, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleNormalEnergy(mris, v, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
    dE_de1 = (e1p - e1m) / (2 * d_dist);
    dE_de2 = (e2p - e2m) / (2 * d_dist);
    norm = sqrt(dE_de1 * dE_de1 + dE_de2 * dE_de2);
    if (norm > MAX_NORM) {
      dE_de1 = MAX_NORM * dE_de1 / norm;
      dE_de2 = MAX_NORM * dE_de2 / norm;
    }

    dx = -l_thick_normal * (dE_de1 * e1x + dE_de2 * e2x);
    dy = -l_thick_normal * (dE_de1 * e1y + dE_de2 * e2y);
    dz = -l_thick_normal * (dE_de1 * e1z + dE_de2 * e2z);
    cx = v->x + parms->dt * dx;
    cy = v->y + parms->dt * dy;
    cz = v->z + parms->dt * dz;
    E0 = mrisSampleNormalEnergy(mris, v, parms, v->x, v->y, v->z);
    E1 = mrisSampleNormalEnergy(mris, v, parms, cx, cy, cz);
    if (E1 > E0) {
      double E2;
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      DiagBreak();
      E2 = mrisSampleNormalEnergy(mris, v, parms, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz);
      if (E2 < E0) {
        dx *= -1;
        dy *= -1;
        dz *= -1;
      }
      else {
        dx = dy = dz = 0;
      }
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (Gdiag_no == vno) {
      float len, xw, yw, zw, xp, yp, zp;

      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);
      dx = xp - xw;
      dy = yp - yw;
      dz = zp - zw;
      len = sqrt(dx * dx + dy * dy + dz * dz);
      if (FZERO(len) == 0) {
        dx /= len;
        dy /= len;
        dz /= len;
      }

      printf(
          "l_thick_normal: v %d: E=%f-->%f (%f), dE = (%2.2f, %2.2f), e1 = (%2.1f, %2.1f, %2.1f), e2 = (%2.1f, %2.1f, "
          "%2.1f), DX = (%2.2f, %2.2f, %2.2f)\n",

          vno,
          E0,
          E1,
          E0 - E1,
          dE_de1,
          dE_de2,
          e1x,
          e1y,
          e1z,
          e2x,
          e2y,
          e2z,
          l_thick_normal * (dE_de1 * e1x + dE_de2 * e2x),
          l_thick_normal * (dE_de1 * e1y + dE_de2 * e2y),
          l_thick_normal * (dE_de1 * e1z + dE_de2 * e2z));
      printf("    N = (%2.2f, %2.2f, %2.2f), D = (%2.2f, %2.2f, %2.2f), dot= %f\n",
             v->wnx,
             v->wny,
             v->wnz,
             dx,
             dy,
             dz,
             v->wnx * dx + v->wny * dy + v->wnz * dz);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  return (NO_ERROR);
}

int mrisComputeThicknessSpringTerm(MRI_SURFACE *mris, double l_thick_spring, INTEGRATION_PARMS *parms)
{
  int vno, max_vno;
  float dE_de1, dE_de2, e1p, e1m, e2p, e2m, cx, cy, cz;
  float E0, E1, dx, dy, dz;
  float e1x, e1y, e1z, e2x, e2y, e2z, max_DE, norm;
  double d_dist = .05 * D_DIST * mris->avg_vertex_dist;
  
  //  int     missed = 0 ;

  if (FZERO(l_thick_spring)) return (0.0);

  max_DE = 0;
  max_vno = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();
    if (v->ripflag) continue;

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;
    /*
      sample the coordinate functions along the tangent plane axes and
      compute the derivates using them.
    */
    E0  = mrisSampleSpringEnergy(mris, vno, v->x, v->y, v->z, parms);
    e1p = mrisSampleSpringEnergy(mris, vno, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z, parms);
    e1m = mrisSampleSpringEnergy(mris, vno, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z, parms);
    e2p = mrisSampleSpringEnergy(mris, vno, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z, parms);
    e2m = mrisSampleSpringEnergy(mris, vno, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z, parms);
    dE_de1 = (e1p - e1m) / (2 * d_dist);
    dE_de2 = (e2p - e2m) / (2 * d_dist);
    if (e1p > E0 && e1m > E0)  // local max in this direction
      dE_de1 = 0;
    if (e2p > E0 && e2m > E0)  // local max in this direction
      dE_de2 = 0;
    norm = sqrt(dE_de1 * dE_de1 + dE_de2 * dE_de2);
    if (norm > max_DE) {
      max_vno = vno;
      max_DE = norm;
    }
    if (norm > MAX_NORM) {
      dE_de1 = MAX_NORM * dE_de1 / norm;
      dE_de2 = MAX_NORM * dE_de2 / norm;
    }

    dx = -l_thick_spring * (dE_de1 * e1x + dE_de2 * e2x);
    dy = -l_thick_spring * (dE_de1 * e1y + dE_de2 * e2y);
    dz = -l_thick_spring * (dE_de1 * e1z + dE_de2 * e2z);
    cx = v->x + parms->dt * dx;
    cy = v->y + parms->dt * dy;
    cz = v->z + parms->dt * dz;
    E1 = mrisSampleSpringEnergy(mris, vno, cx, cy, cz, parms);
    if (E1 > E0) {
      double E2;
      if (vno == Gdiag_no) DiagBreak();

      E2 = mrisSampleSpringEnergy(mris, vno, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz, parms);
      if (E2 < E0) {
        dx *= -1;
        dy *= -1;
        dz *= -1;
      }
      else {
        if (e1m < e2m && e1m < e1p && e1m < e2p && e1m < E0)  // e1m is best
        {
          dx = -d_dist * e1x;
          dy = -d_dist * e1y;
          dz = d_dist * e1z;
        }
        else if (e2m < e1p && e2m < e2p && e2m < E0)  // e2m is best
        {
          dx = -d_dist * e2x;
          dy = -d_dist * e2y;
          dz = -d_dist * e2z;
        }
        else if (e1p < e2p && e1p < E0)  // e1p is best
        {
          dx = d_dist * e1x;
          dy = d_dist * e1y;
          dz = d_dist * e1z;
        }
        else if (e2p < E0)  // e2p is best
        {
          dx = d_dist * e2x;
          dy = d_dist * e2y;
          dz = d_dist * e2z;
        }
        else {
          dx = dy = dz = 0;
        }
      }
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (Gdiag_no == vno) {
      float len, xw, yw, zw, xp, yp, zp;

      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);
      dx = xp - xw;
      dy = yp - yw;
      dz = zp - zw;
      len = sqrt(dx * dx + dy * dy + dz * dz);
      if (FZERO(len) == 0) {
        dx /= len;
        dy /= len;
        dz /= len;
      }

      printf(
          "l_thick_spring: v %d: E=%f-->%f (%f), dE = (%2.2f, %2.2f), e1 = (%2.1f, %2.1f, %2.1f), e2 = (%2.1f, %2.1f, "
          "%2.1f), DX = (%2.2f, %2.2f, %2.2f)\n",

          vno,
          E0,
          E1,
          E0 - E1,
          dE_de1,
          dE_de2,
          e1x,
          e1y,
          e1z,
          e2x,
          e2y,
          e2z,
          l_thick_spring * (dE_de1 * e1x + dE_de2 * e2x),
          l_thick_spring * (dE_de1 * e1y + dE_de2 * e2y),
          l_thick_spring * (dE_de1 * e1z + dE_de2 * e2z));
      printf("    N = (%2.2f, %2.2f, %2.2f), D = (%2.2f, %2.2f, %2.2f), dot= %f\n",
             v->wnx,
             v->wny,
             v->wnz,
             dx,
             dy,
             dz,
             v->wnx * dx + v->wny * dy + v->wnz * dz);
    }
  }
  return (NO_ERROR);
}

int mrisComputeAshburnerTriangleTerm(MRI_SURFACE *mris, double l_ashburner_triangle, INTEGRATION_PARMS *parms)
{
  int vno;
  //  int     missed = 0 ;

  if (FZERO(l_ashburner_triangle)) {
    return (0.0);
  }

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    float dE_de1, dE_de2, e1p, e1m, e2p, e2m, cx, cy, cz;
    float E0, E1, dx, dy, dz;
    float e1x, e1y, e1z, e2x, e2y, e2z, norm;
    VERTEX *v;
    double d_dist = D_DIST * mris->avg_vertex_dist;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) ROMP_PF_continue;

    e1x = v->e1x;
    e1y = v->e1y;
    e1z = v->e1z;
    e2x = v->e2x;
    e2y = v->e2y;
    e2z = v->e2z;
    /*
      sample the coordinate functions along the tangent plane axes and
      compute the derivates using them.
    */
    e1p = mrisSampleAshburnerTriangleEnergy(
        mris, vno, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleAshburnerTriangleEnergy(
        mris, vno, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleAshburnerTriangleEnergy(
        mris, vno, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleAshburnerTriangleEnergy(
        mris, vno, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
    dE_de1 = (e1p - e1m) / (2 * d_dist);
    dE_de2 = (e2p - e2m) / (2 * d_dist);
    norm = sqrt(dE_de1 * dE_de1 + dE_de2 * dE_de2);
    if (norm > MAX_NORM) {
      dE_de1 = MAX_NORM * dE_de1 / norm;
      dE_de2 = MAX_NORM * dE_de2 / norm;
    }

    dx = -l_ashburner_triangle * (dE_de1 * e1x + dE_de2 * e2x);
    dy = -l_ashburner_triangle * (dE_de1 * e1y + dE_de2 * e2y);
    dz = -l_ashburner_triangle * (dE_de1 * e1z + dE_de2 * e2z);
    cx = v->x + parms->dt * dx;
    cy = v->y + parms->dt * dy;
    cz = v->z + parms->dt * dz;
    E0 = mrisSampleAshburnerTriangleEnergy(mris, vno, parms, v->x, v->y, v->z);
    E1 = mrisSampleAshburnerTriangleEnergy(mris, vno, parms, cx, cy, cz);
    if (E1 > E0) {
      double E2;
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      DiagBreak();
      E2 = mrisSampleAshburnerTriangleEnergy(
          mris, vno, parms, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz);
      if (E2 < E0) {
        dx *= -1;
        dy *= -1;
        dz *= -1;
      }
      else {
        dx = dy = dz = 0;
      }
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (Gdiag_no == vno) {
      printf(
          "l_ashburner_triangle: v %d: E=%f-->%f (%f), dE = (%2.2f, %2.2f), e1 = (%2.1f, %2.1f, %2.1f), e2 = (%2.1f, "
          "%2.1f, %2.1f), DX = (%2.2f, %2.2f, %2.2f)\n",

          vno,
          E0,
          E1,
          E0 - E1,
          dE_de1,
          dE_de2,
          e1x,
          e1y,
          e1z,
          e2x,
          e2y,
          e2z,
          l_ashburner_triangle * (dE_de1 * e1x + dE_de2 * e2x),
          l_ashburner_triangle * (dE_de1 * e1y + dE_de2 * e2y),
          l_ashburner_triangle * (dE_de1 * e1z + dE_de2 * e2z));
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (NO_ERROR);
}

/*!
  \fn double mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht, MHT *mht_faces)
  \brief The repulsive term causes vertices to push away from each
  other based on the distance in 3D space (does not apply to nearest
  neighbors). This helps to prevent self-intersection. The force is
  inversely proportional to the distance to the 6th power (hidden
  parameter). Sets v->{dx,dy,dz}. 
  Hidden parameters:
    REPULSE_K - scaling term
    REPULSE_E - sets minimum distance
    4 - scaling term
*/
double mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht, MHT *mht_faces)
{
  int vno, num, min_vno, i, n;
  float dist, dx, dy, dz, x, y, z, min_d;
  double sse_repulse, v_sse;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  min_d = 1000.0;
  min_vno = 0;
  sse_repulse = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) 
      continue;

    x = v->x;
    y = v->y;
    z = v->z;

    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket)
      continue;
    
    MHB *bin;
    for (v_sse = 0.0, bin = bucket->bins, num = i = 0; i < bucket->nused; i++, bin++) {

      /* don't be repelled by myself */
      if (bin->fno == vno)
        continue; 

      /* don't be repelled by a neighbor */
      for (n = 0; n < vt->vtotal; n++){
        if (vt->v[n] == bin->fno) {
          break;
        }
      }
      if (n < vt->vtotal) 
        continue;

      VERTEX const * const vn = &mris->vertices[bin->fno];
      if (!vn->ripflag) {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        dist = sqrt(dx * dx + dy * dy + dz * dz) + REPULSE_E;

        if (vno == Gdiag_no) {
          if (dist - REPULSE_E < min_d) {
            min_vno = bin->fno;
            min_d = dist - REPULSE_E;
          }
        }

        dist = dist * dist * dist;
        dist *= dist; /* dist^6 */
	// dist = pow(dist,6.0); 
        v_sse += REPULSE_K / dist;
      }
    } // loop over bucket

    sse_repulse += v_sse; // does not divide by the number of bins

    if (vno == Gdiag_no && !FZERO(v_sse)) {
      printf("v %d: repulse sse:    min_dist=%2.4f, v_sse %2.4f\n", vno, min_d, v_sse);
    }
    
    MHTrelBucket(&bucket);
  }

  return (l_repulse * sse_repulse);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

double mrisComputeRepulsiveRatioEnergy(MRI_SURFACE *mris, double l_repulse)
{
  int vno, n;
  double sse_repulse, v_sse, dist, dx, dy, dz, x, y, z, canon_dist, cdx, cdy, cdz;

  if (FZERO(l_repulse))
    return (0.0);

  for (sse_repulse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (v->ripflag) 
      continue;

    x = v->x;
    y = v->y;
    z = v->z;
    for(v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (!vn->ripflag) {
        dx = x - vn->x;
        dy = y - vn->y;
        dz = z - vn->z;
        dist = sqrt(dx * dx + dy * dy + dz * dz);
        cdx = vn->cx - v->cx;
        cdy = vn->cy - v->cy;
        cdz = vn->cz - v->cz;
        canon_dist = sqrt(cdx * cdx + cdy * cdy + cdz * cdz) + REPULSE_E;
        dist /= canon_dist;
        dist += REPULSE_E;
#if 0
        v_sse += REPULSE_K / (dist*dist*dist*dist) ;
#else
        v_sse += REPULSE_K / (dist * dist);
#endif
      }
    }
    sse_repulse += v_sse;
  }
  return (l_repulse * sse_repulse);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeThicknessSmoothnessTerm(MRI_SURFACE *mris, double l_tsmooth, INTEGRATION_PARMS *parms)
{
  int vno, n, num;
  float dx, dy, dz, x, y, z, dn, d0, vx, vy, vz, delta;

  if (FZERO(l_tsmooth)) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    x = v->x;
    y = v->y;
    z = v->z;
    vx = v->x - v->origx;
    vy = v->y - v->origy;
    vz = v->z - v->origz;
    d0 = vx * vx + vy * vy + vz * vz;
    dx = dy = dz = 0.0;
    for (num = n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (!vn->ripflag) {
        dn = SQR(vn->x - vn->origx) + SQR(vn->origy - vn->y) + SQR(vn->origz - vn->z);
        delta = d0 - dn;
        dx -= delta * vx;
        dy -= delta * vy;
        dz -= delta * vz;
        num++;
      }
    }
    if (num) {
      dx /= num;
      dy /= num;
      dz /= num;
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
    if (vno == Gdiag_no) {
      fprintf(stdout, "v %d tsmooth term:        (%2.3f, %2.3f, %2.3f)\n", vno, dx, dy, dz);
    }
  }
  return (NO_ERROR);
}

/*!
  \fn int mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht, MHT *mht_faces)
  \brief The repulsive term causes vertices to push away from each
  other based on the distance in 3D space (does not apply to nearest
  neighbors). This helps to prevent self-intersection. The force is
  inversely proportional to the distance to the 7th power (hidden
  parameter). Sets v->{dx,dy,dz}. 
  Hidden parameters:
    REPULSE_K - scaling term
    REPULSE_E - sets minimum distance
    4 - scaling term
  Does not appear to use annotation
*/
int mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht, MHT *mht_faces)
{
  int vno, num, min_vno, i, n;
  float dist, dx, dy, dz, x, y, z, sx, sy, sz, min_d, min_scale, norm;
  double scale;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  min_d = 100000.0;
  min_scale = 1.0;
  min_vno = 0;
  // loop thru vertices
  for (vno = 0; vno < mris->nvertices; vno++) {
    // VERTEX_TOPOLOGY is a subset of VERTEX
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) 
      continue;

    if (vno == Gdiag_no)
      DiagBreak();

    x = v->x;
    y = v->y;
    z = v->z;

    // Get the list of vertices that are close in 3d space.
    // How close is close? Determined by bucket size?
    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if(!bucket)
      continue;

    // Go through list of vertices in bucket
    MHB *bin;
    sx = sy = sz = 0.0;
    for (bin = bucket->bins, num = i = 0; i < bucket->nused; i++, bin++) {

      /* don't be repelled by myself */
      if (bin->fno == vno) 
        continue; 

      /* don't be repelled by a neighbor */
      for (n = 0; n < vt->vtotal; n++){
        if (vt->v[n] == bin->fno) {
          break;
        }
      }
      if (n < vt->vtotal) 
        continue;

      VERTEX const * const vn = &mris->vertices[bin->fno];
      if (!vn->ripflag) {
	// Compute the distance between the two vertices
        dx = x - vn->x;
        dy = y - vn->y;
        dz = z - vn->z;
        dist = sqrt(dx * dx + dy * dy + dz * dz) + REPULSE_E;
	// REPULSE_E is a hidden parameter

	// Cost = K/pow(dist,6) (see mrisComputeRepulsiveEnergy())
        // dCost/dx = -dx*K/pow(dist,8) but it is incorrectly computed
        // here as dCost/dx = -dx*K/pow(dist,7). pow8 still not right
        // when E!=0. Maybe it does not matter much because you just
        // want to push them apart.
	// The multiplication by dx, dy, dz happens below. The
	// negative sign is not applied because it is the step that is
	// actually computed.
        scale = 4 * REPULSE_K / (dist * dist * dist * dist * dist * dist * dist); /* ^-7 */
	// REPULSE_K is a hidden parameter
	// 4 is a hidden parameter (?)

        if (vno == Gdiag_no) {
          if (dist - REPULSE_E < 0.75) {
            DiagBreak();
          }
          if (dist - REPULSE_E < min_d) {
            min_vno = bin->fno;
            min_d = dist - REPULSE_E;
            min_scale = scale;
          }
        }

	// Normalize dx, dy, dz
        norm = sqrt(dx * dx + dy * dy + dz * dz);
        if(FZERO(norm))  norm = 1.0;
        dx /= norm;
        dy /= norm;
        dz /= norm;

        if (!isfinite(dx) || !isfinite(dy) || !isfinite(dz)) 
          DiagBreak();

        sx += scale * dx;
        sy += scale * dy;
        sz += scale * dz;

        num++; // number of hits in the bucket

      } // not ripped

    } // loop over bucket

    if (num) {
      // "scale" here is a way to compute the mean (div by num) and
      // apply the weighting factor at the same time. Not to be
      // confused with "scale" above.  
      // NOTE: this dividing by num is not consistent with
      // mrisComputeRepulsiveEnergy().
      scale = l_repulse / (double)num;
      sx *= scale;
      sy *= scale;
      sz *= scale;
    }

    v->dx += sx;
    v->dy += sy;
    v->dz += sz;

    if ((vno == Gdiag_no) && min_d < 1) {
      fprintf(stdout, "v %d self repulse term:   (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
      fprintf(stdout, "min_dist @ %d = %2.2f, scale = %2.1f\n", min_vno, min_d, min_scale);
    }
    
    MHTrelBucket(&bucket);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeRepulsiveRatioTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int vno, num, min_vno, i;
  float dist, dx, dy, dz, x, y, z, sx, sy, sz, min_d, min_scale, canon_dist, cdx, cdy, cdz;
  double scale;
  VERTEX *v, *vn;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  min_d = 1000.0;
  min_scale = 1.0;
  min_vno = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    x = v->x;
    y = v->y;
    z = v->z;
    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) {
      continue;
    }
    sx = sy = sz = 0.0;
    MHB *bin;
    for (bin = bucket->bins, num = i = 0; i < bucket->nused; i++, bin++) {
      if (bin->fno == vno) {
        continue; /* don't be repelled by myself */
      }
      vn = &mris->vertices[bin->fno];
      if (!vn->ripflag) {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        dist = sqrt(dx * dx + dy * dy + dz * dz);
        cdx = vn->cx - v->cx;
        cdy = vn->cy - v->cy;
        cdz = vn->cz - v->cz;
        canon_dist = sqrt(cdx * cdx + cdy * cdy + cdz * cdz) + REPULSE_E;
        dist /= canon_dist;
        dist += REPULSE_E;
#if 0
        scale = -4*REPULSE_K / (dist*dist*dist*dist*dist) ;
#else
        scale = -4 * REPULSE_K / (dist * dist * dist);
#endif
        if (vno == Gdiag_no) {
          if (dist - REPULSE_E < min_d) {
            min_vno = bin->fno;
            min_d = dist - REPULSE_E;
            min_scale = scale;
          }
        }
        sx += scale * dx;
        sy += scale * dy;
        sz += scale * dz;
        num++;
      }
      MHTrelBucket(&bucket);
    }
    if (num) {
      scale = l_repulse / (double)num;
      sx *= scale;
      sy *= scale;
      sz *= scale;
    }
    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) {
      vn = &mris->vertices[min_vno];
      dx = x - vn->x;
      dy = y - vn->y;
      dz = z - vn->z;

      fprintf(stdout, "v %d repulse term:        (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
      fprintf(stdout, "min_dist @ %d = %2.2f, scale = %2.1f\n", min_vno, min_d, min_scale);
    }
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeSpringEnergy(MRI_SURFACE *mris)
{
  int vno, n;
  double area_scale, sse_spring, v_sse;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      v_sse += (v->dist[n] * v->dist[n]);
    }
    sse_spring += area_scale * v_sse;
  }
  return (sse_spring);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description


  Note: this function assumes that the mris surface has the original
  (i.e. after global rotational alignment)
  spherical coordinates in the TMP2_VERTICES
  ------------------------------------------------------*/
double mrisComputeLaplacianEnergy(MRI_SURFACE *mris)
{
  int vno, n;
  double area_scale, sse_lap, v_sse, dx, dy, dz, vx, vy, vz, vnx, vny, vnz, error;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  for (sse_lap = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    vx = v->x - v->tx2;
    vy = v->y - v->ty2;
    vz = v->z - v->tz2;
    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      vnx = vn->x - vn->tx2;
      vny = vn->y - vn->ty2;
      vnz = vn->z - vn->tz2;
      dx = vnx - vx;
      dy = vny - vy;
      dz = vnz - vz;
      error = dx * dx + dy * dy + dz * dz;
      v_sse += error;
    }
    sse_lap += area_scale * v_sse;
  }
  return (sse_lap);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris)
{
  int vno, n;
  double area_scale, sse_spring, v_sse;
  float dx, dy, dz, x, y, z, nc, dist_sq;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = FZERO(mris->total_area) ? 1.0 : mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    x = v->x;
    y = v->y;
    z = v->z;

    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - x;
      dy = vn->y - y;
      dz = vn->z - z;
      nc = dx * v->nx + dy * v->ny + dz * v->nz;
      dx -= nc * v->nx;
      dy -= nc * v->ny;
      dz -= nc * v->nz;
      dist_sq = dx * dx + dy * dy + dz * dz;
      v_sse += dist_sq;
    }
    sse_spring += area_scale * v_sse;
  }
  return (sse_spring);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
    compute nonlinear spring energy (stolen from BET, thanks Steve)
------------------------------------------------------*/
#define RMIN 1
#define RMAX 5

int mrisComputeNonlinearSpringTerm(MRI_SURFACE *mris, double l_nlspring, INTEGRATION_PARMS *parms)
{
  int vno, n;
  double area_scale, sse_spring, E, F, f, rmin, rmax;
  float dx, dy, dz, nc, r, lsq, mean_vdist;

  if (FZERO(parms->l_nlspring)) {
    return (NO_ERROR);
  }

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  mean_vdist = MRIScomputeVertexSpacingStats(mris, NULL, NULL, NULL, NULL, NULL, CURRENT_VERTICES);
  lsq = mean_vdist * mean_vdist;

  rmax = parms->rmax;
  rmin = parms->rmin;
  if (FZERO(rmin) || FZERO(rmax))
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisComputeNonlinearSpringTerm: rmin or rmax = 0!"));

  F = 6.0 / (1.0 / rmin - 1.0 / rmax);
  E = (1.0 / rmin + 1.0 / rmax) / 2;
  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - v->x;
      dy = vn->y - v->y;
      dz = vn->z - v->z;
      //      lsq = dx*dx + dy*dy + dz*dz ;
      nc = dx * v->nx + dy * v->ny + dz * v->nz;
      dx = nc * v->nx;
      dy = nc * v->ny;
      dz = nc * v->nz;  // sn
      r = lsq / fabs(2.0 * nc);
      if (r < rmin) {
        DiagBreak();
      }
      f = nc * (1 + tanh(F * (1.0 / r - E))) / 2.0;
      if (vno == Gdiag_no)
        printf("l_nlspring: f = %2.3f (r = %2.2f), dx = (%2.2f, %2.2f, %2.2f)\n",
               f,
               r,
               v->nx * f * l_nlspring,
               v->ny * f * l_nlspring,
               v->nz * f * l_nlspring);
      v->dx += v->nx * f * l_nlspring;
      v->dy += v->ny * f * l_nlspring;
      v->dz += v->nz * f * l_nlspring;
    }
  }
  return (NO_ERROR);
}

double mrisComputeNonlinearSpringEnergy(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno, n;
  double area_scale, sse_spring, E, F, f, rmin, rmax, ftotal;
  float dx, dy, dz, nc, r, lsq, mean_vdist;

  mean_vdist = MRIScomputeVertexSpacingStats(mris, NULL, NULL, NULL, NULL, NULL, CURRENT_VERTICES);
  lsq = mean_vdist * mean_vdist;
#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  rmin = parms->rmin;
  rmax = parms->rmax;
  if (FZERO(rmin) || FZERO(rmax))
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisComputeNonlinearSpringTerm: rmin or rmax = 0!"));

  F = 6.0 / (1.0 / rmin - 1.0 / rmax);
  E = (1.0 / rmin + 1.0 / rmax) / 2;
  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    for (ftotal = r = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - v->x;
      dy = vn->y - v->y;
      dz = vn->z - v->z;
      //      lsq = dx*dx + dy*dy + dz*dz ;
      nc = dx * v->nx + dy * v->ny + dz * v->nz;
      dx = nc * v->nx;
      dy = nc * v->ny;
      dz = nc * v->nz;  // sn
      r = lsq / fabs(2.0 * nc);
      f = (1 + tanh(F * (1.0 / r - E)));
      ftotal += f * f;
    }
    if (vno == Gdiag_no) {
      printf("E_nlspring: f = %2.3f\n", ftotal / vt->vnum);
    }
    sse_spring += area_scale * ftotal / vt->vnum;
  }
  return (sse_spring);
}

double mrisComputeSurfaceRepulsionEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int vno, max_vno, i;
  float dx, dy, dz, x, y, z, sx, sy, sz, norm[3], dot;
  float max_scale, max_dot;
  double scale, sse;
  VERTEX *v, *vn;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    x = v->x;
    y = v->y;
    z = v->z;
    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) {
      continue;
    }
    sx = sy = sz = 0.0;
    max_dot = max_scale = 0.0;
    max_vno = 0;
    MHB *bin = bucket->bins;
    for (i = 0; i < bucket->nused; i++, bin++) {
      vn = &mris->vertices[bin->fno];
      if (bin->fno == Gdiag_no) {
        DiagBreak();
      }
      if (vn->ripflag) {
        continue;
      }
      dx = x - vn->origx;
      dy = y - vn->origy;
      dz = z - vn->origz;
      mrisComputeOrigNormal(mris, bin->fno, norm);
      dot = dx * norm[0] + dy * norm[1] + dz * norm[2];
      if (dot > 1) {
        continue;
      }
      if (dot < 0 && vno == Gdiag_no) {
        DiagBreak();
      }
      if (dot > MAX_NEG_RATIO) {
        dot = MAX_NEG_RATIO;
      }
      else if (dot < -MAX_NEG_RATIO) {
        dot = -MAX_NEG_RATIO;
      }
#if 0
      scale = l_repulse / (1.0+exp(NEG_AREA_K*dot)) ;
#else
      scale = l_repulse * pow(1.0 - (double)dot, 4.0);
#endif
      if (scale > max_scale) {
        max_scale = scale;
        max_vno = bin->fno;
        max_dot = dot;
      }
      sx += (scale * v->nx);
      sy += (scale * v->ny);
      sz += (scale * v->nz);
    }

    sse += (sx * sx + sy * sy + sz * sz);
    if (vno == Gdiag_no) fprintf(stdout, "v %d inside repulse energy %2.3f\n", vno, (sx * sx + sy * sy + sz * sz));
    
    MHTrelBucket(&bucket);
  }
  return (sse);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  each face has 2 triangles defined by it:

  V0       d      V3
  o--------------o
  |              |
  | A0           |
  a |              | c
  |              |
  |           A1 |
  o--------------o
  V1      b        V2

  a = V1 - V0
  d = V3 - V0
  A0 = 0.5 (a x d) . n

  b = V1 - V2
  c = V3 - V2
  A1 = 0.5 (c x b) . n

  each face has 1 triangle defined by it:

  V0    b     V2
  o----------o
  |         /
  | A0    /
  a |     /
  |   /
  | /
  o
  V1

  a = V1 - V0
  b = V2 - V0
  A0 = 0.5 (a x b) . n

  ------------------------------------------------------*/
int mrisComputeAngleAreaTerms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int fno, ano;
  VERTEX *v0, *v1, *v2, *va, *vb, *vo;
  VECTOR *v_a, *v_b, *v_a_x_n, *v_b_x_n, *v_n, *v_tmp, *v_sum;
  FACE *face;
  float orig_area, area, l_parea, l_area, l_angle, delta, len, area_scale;

#if METRIC_SCALE
  if (mris->patch || (mris->status != MRIS_SPHERE && mris->status != MRIS_PARAMETERIZED_SPHERE)) {
    area_scale = 1.0f;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0f;
#endif

  l_angle = parms->l_angle;
  l_area = parms->l_area;
  l_parea = parms->l_parea;
  if (!FZERO(parms->l_nlarea)) {
    mrisComputeNonlinearAreaTerm(mris, parms);
  }

  if (FZERO(l_area) && FZERO(l_angle) && FZERO(l_parea) && FZERO(parms->l_pangle)) {
    return (NO_ERROR);
  }

  v_a = VectorAlloc(3, MATRIX_REAL);
  v_b = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL);

  v_tmp = VectorAlloc(3, MATRIX_REAL);
  v_sum = VectorAlloc(3, MATRIX_REAL);
  v_a_x_n = VectorAlloc(3, MATRIX_REAL);
  v_b_x_n = VectorAlloc(3, MATRIX_REAL);

  /* calculcate movement of each vertex caused by each triangle */
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    if (fno == Gdiag_no2) {
      DiagBreak();
    }
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
    VECTOR_LOAD(v_n, fNorm->nx, fNorm->ny, fNorm->nz);
    v0 = &mris->vertices[face->v[0]];
    v1 = &mris->vertices[face->v[1]];
    v2 = &mris->vertices[face->v[2]];
    VERTEX_EDGE(v_a, v0, v1);
    VERTEX_EDGE(v_b, v0, v2);
    orig_area = fNorm->orig_area;
    area = area_scale * face->area;
    delta = 0.0;
    if (!FZERO(l_parea)) {
      delta += l_parea * (area - orig_area);
    }

    if (!FZERO(l_area)) {
      if (area <= 0.0f) {
        delta += l_area * (area - orig_area);
      }
    }

    if (!FZERO(l_area) && (face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no)) {
      printf("face %d, orig area %2.2f, area, %2.2f, delta = %2.2f\n", fno, orig_area, area, delta);
      DiagBreak();
    }

    V3_CROSS_PRODUCT(v_a, v_n, v_a_x_n);
    V3_CROSS_PRODUCT(v_b, v_n, v_b_x_n);

    /* calculate movement of vertices in order, 0-3 */

    /* v0 */
    V3_SCALAR_MUL(v_a_x_n, -1.0f, v_sum);
    V3_ADD(v_sum, v_b_x_n, v_sum);
    V3_SCALAR_MUL(v_sum, delta, v_sum);
    v0->dx += V3_X(v_sum);
    v0->dy += V3_Y(v_sum);
    v0->dz += V3_Z(v_sum);
    if (face->v[0] == Gdiag_no && !FZERO(parms->l_area)) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[0], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }

    /* v1 */
    V3_SCALAR_MUL(v_b_x_n, -delta, v_sum);
    v1->dx += V3_X(v_sum);
    v1->dy += V3_Y(v_sum);
    v1->dz += V3_Z(v_sum);
    if (face->v[1] == Gdiag_no && !FZERO(parms->l_area)) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[1], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }

    /* v2 */
    V3_SCALAR_MUL(v_a_x_n, delta, v_sum);
    v2->dx += V3_X(v_sum);
    v2->dy += V3_Y(v_sum);
    v2->dz += V3_Z(v_sum);
    if (face->v[2] == Gdiag_no && !FZERO(parms->l_area)) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[2], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }

    /* now calculate the angle contributions */
    if (!FZERO(l_angle) || !FZERO(parms->l_pangle)) {
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        switch (ano) {
          default:
          case 0:
            vo = v0;
            va = v2;
            vb = v1;
            break;
          case 1:
            vo = v1;
            va = v0;
            vb = v2;
            break;
          case 2:
            vo = v2;
            va = v1;
            vb = v0;
            break;
        }
        delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f && ((mris->status == MRIS_PLANE) || (mris->status == MRIS_SPHERE) ||
                                         (mris->status == MRIS_PARAMETERIZED_SPHERE))) {
          delta = 0.0f;
        }
#endif

        // for pangle term don't penalize angles that are wider, just narrower ones to avoid pinching
        if (FZERO(parms->l_angle) && !FZERO(parms->l_pangle) && delta < 0) delta = 0.0;

        if (!FZERO(parms->l_angle + parms->l_pangle) &&
            ((face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no)))
          printf("face %d, ano %d, orig angle %2.1f, current angle %2.1f, delta %2.1f\n",
                 fno,
                 ano,
                 DEGREES(face->orig_angle[ano]),
                 DEGREES(face->angle[ano]),
                 DEGREES(delta));
        delta *= parms->l_angle;
        VERTEX_EDGE(v_a, vo, va);
        VERTEX_EDGE(v_b, vo, vb);

        /* this angle's contribution to va */
        V3_CROSS_PRODUCT(v_a, v_n, v_tmp);
        len = V3_DOT(v_a, v_a);
        if (!FZERO(len)) {
          V3_SCALAR_MUL(v_tmp, delta / len, v_tmp);
        }
        else {
          V3_SCALAR_MUL(v_tmp, 0.0f, v_tmp);
        }
        va->dx += V3_X(v_tmp);
        va->dy += V3_Y(v_tmp);
        va->dz += V3_Z(v_tmp);
        if (va - mris->vertices == Gdiag_no) {
          printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", Gdiag_no, V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
        }

        /* this angle's contribution to vb */
        V3_CROSS_PRODUCT(v_n, v_b, v_sum);
        len = V3_DOT(v_b, v_b);
        if (!FZERO(len)) {
          V3_SCALAR_MUL(v_sum, delta / len, v_sum);
        }
        else {
          V3_SCALAR_MUL(v_sum, 0.0f, v_sum);
        }
        vb->dx += V3_X(v_sum);
        vb->dy += V3_Y(v_sum);
        vb->dz += V3_Z(v_sum);
        if (vb - mris->vertices == Gdiag_no) {
          printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", Gdiag_no, V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
        }

        /* this angle's contribution to vo */
        V3_ADD(v_tmp, v_sum, v_sum);
        vo->dx -= V3_X(v_sum);
        vo->dy -= V3_Y(v_sum);
        vo->dz -= V3_Z(v_sum);
        if (vo - mris->vertices == Gdiag_no) {
          printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", Gdiag_no, V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
        }
      }
    }
  } /* done with all faces */

  VectorFree(&v_a);
  VectorFree(&v_b);
  VectorFree(&v_tmp);
  VectorFree(&v_sum);
  VectorFree(&v_n);

  VectorFree(&v_a_x_n);
  VectorFree(&v_b_x_n);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeNonlinearAreaTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int fno;
  VERTEX *v0, *v1, *v2;
  VECTOR *v_a, *v_b, *v_a_x_n, *v_b_x_n, *v_n, *v_tmp, *v_sum;
  FACE *face;
  double orig_area, area, delta, area_scale, scale, l_nlarea, ratio, nlscale;

#if METRIC_SCALE
  if (mris->patch || (mris->status != MRIS_SPHERE && mris->status != MRIS_PARAMETERIZED_SPHERE)) {
    area_scale = 1.0f;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0f;
#endif

  l_nlarea = parms->l_nlarea;

  if (FZERO(l_nlarea)) {
    return (NO_ERROR);
  }

  v_a = VectorAlloc(3, MATRIX_REAL);
  v_b = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL);

  v_tmp = VectorAlloc(3, MATRIX_REAL);
  v_sum = VectorAlloc(3, MATRIX_REAL);
  v_a_x_n = VectorAlloc(3, MATRIX_REAL);
  v_b_x_n = VectorAlloc(3, MATRIX_REAL);

  /* calculcate movement of each vertex caused by each triangle */
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    if (face->area < 0) {
      DiagBreak();
    }
    if (face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no) {
      DiagBreak();
    }
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
    VECTOR_LOAD(v_n, fNorm->nx, fNorm->ny, fNorm->nz);
    v0 = &mris->vertices[face->v[0]];
    v1 = &mris->vertices[face->v[1]];
    v2 = &mris->vertices[face->v[2]];
    VERTEX_EDGE(v_a, v0, v1);
    VERTEX_EDGE(v_b, v0, v2);
    orig_area = fNorm->orig_area;
    area = area_scale * face->area;
#if SCALE_NONLINEAR_AREA
    if (!FZERO(orig_area)) {
      ratio = area / orig_area;
    }
    else {
      ratio = 0.0f;
    }
#else
    ratio = area;
#endif

    if (ratio > MAX_NEG_RATIO) {
      ratio = MAX_NEG_RATIO;
    }
    else if (ratio < -MAX_NEG_RATIO) {
      ratio = -MAX_NEG_RATIO;
    }
#if 0
    scale = l_nlarea * (1 - (1/(1.0+exp(-NEG_AREA_K*ratio)))) ;
#else
    scale = l_nlarea / (1.0 + exp(NEG_AREA_K * ratio));
#endif
    delta = scale * (area - orig_area);
    nlscale = (mris->total_area / mris->nfaces) / ((area < 0 ? 0 : area) + 0.001);
    nlscale = 1;
    delta *= (nlscale * nlscale);

    if (face->v[0] == Gdiag_no || face->v[1] == Gdiag_no || face->v[2] == Gdiag_no) {
      printf("face %d, orig area %2.2f, area, %2.2f, delta = %2.2f\n", fno, orig_area, area, scale);
      DiagBreak();
    }

    V3_CROSS_PRODUCT(v_a, v_n, v_a_x_n);
    V3_CROSS_PRODUCT(v_b, v_n, v_b_x_n);

    /* calculate movement of vertices in order, 0-3 */

    /* v0 */
    V3_SCALAR_MUL(v_a_x_n, -1.0f, v_sum);
    V3_ADD(v_sum, v_b_x_n, v_sum);
    V3_SCALAR_MUL(v_sum, delta, v_sum);
    v0->dx += V3_X(v_sum);
    v0->dy += V3_Y(v_sum);
    v0->dz += V3_Z(v_sum);
    if (face->v[0] == Gdiag_no) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[0], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }

    /* v1 */
    V3_SCALAR_MUL(v_b_x_n, -delta, v_sum);
    v1->dx += V3_X(v_sum);
    v1->dy += V3_Y(v_sum);
    v1->dz += V3_Z(v_sum);
    if (face->v[1] == Gdiag_no) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[1], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }

    /* v2 */
    V3_SCALAR_MUL(v_a_x_n, delta, v_sum);
    v2->dx += V3_X(v_sum);
    v2->dy += V3_Y(v_sum);
    v2->dz += V3_Z(v_sum);
    if (face->v[2] == Gdiag_no) {
      printf("\tv %d, area/angle term = (%2.3f %2.3f, %2.3f)\n", face->v[2], V3_X(v_sum), V3_Y(v_sum), V3_Z(v_sum));
    }
  } /* done with all faces */

  VectorFree(&v_a);
  VectorFree(&v_b);
  VectorFree(&v_tmp);
  VectorFree(&v_sum);
  VectorFree(&v_n);

  VectorFree(&v_a_x_n);
  VectorFree(&v_b_x_n);
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeSurfaceRepulsionTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht)
{
  int vno, max_vno, i;
  float dx, dy, dz, x, y, z, sx, sy, sz, norm[3], dot;
  float max_scale, max_dot;
  double scale;
  VERTEX *v, *vn;

  if (FZERO(l_repulse)) {
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
    if (v->cropped)  // turn off this term for vertices that are intersecting so they don't drive their neighbors crazy
      continue;

    x = v->x;
    y = v->y;
    z = v->z;

    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) {
      continue;
    }
    sx = sy = sz = 0.0;
    max_dot = max_scale = 0.0;
    max_vno = 0;

    MHB *bin = bucket->bins;
    for (i = 0; i < bucket->nused; i++, bin++) {
      vn = &mris->vertices[bin->fno];
      if (bin->fno == Gdiag_no) {
        DiagBreak();
      }
      if (vn->ripflag) {
        continue;
      }
      dx = x - vn->origx;
      dy = y - vn->origy;
      dz = z - vn->origz;
      mrisComputeOrigNormal(mris, bin->fno, norm);
      dot = dx * norm[0] + dy * norm[1] + dz * norm[2];
      if (dot > 1) {
        continue;
      }
      if (dot < 0 && vno == Gdiag_no) {
        DiagBreak();
      }
      if (dot > MAX_NEG_RATIO) {
        dot = MAX_NEG_RATIO;
      }
      else if (dot < -MAX_NEG_RATIO) {
        dot = -MAX_NEG_RATIO;
      }
#if 0
      scale = l_repulse / (1.0+exp(NEG_AREA_K*dot)) ;
#else
      scale = l_repulse * pow(1.0 - (double)dot, 4.0);
#endif
      if (scale > max_scale) {
        max_scale = scale;
        max_vno = bin->fno;
        max_dot = dot;
      }
      sx += (scale * v->nx);
      sy += (scale * v->ny);
      sz += (scale * v->nz);
    }

    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) {
      vn = &mris->vertices[max_vno];
      dx = x - vn->x;
      dy = y - vn->y;
      dz = z - vn->z;

      fprintf(stdout, "v %d inside repulse term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
      fprintf(stdout, "max_scale @ %d = %2.2f, max dot = %2.2f\n", max_vno, max_scale, max_dot);
    }
    
    MHTrelBucket(&bucket);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int mrisComputeWhichSurfaceRepulsionTerm(
    MRI_SURFACE *mris, double l_repulse, MHT *mht, int which, float dot_thresh)
{
  int vno, max_vno, i;
  float dx, dy, dz, x, y, z, sx, sy, sz, norm[3], dot;
  float max_scale, max_dot;
  double scale, sgn;
  VERTEX *v, *vn;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  if (l_repulse < 0) {
    sgn = -1;
  }
  else {
    sgn = 1;
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    x = v->x;
    y = v->y;
    z = v->z;

    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) {
      continue;
    }
    
    sx = sy = sz = 0.0;
    max_dot = max_scale = 0.0;
    max_vno = 0;
    MHB *bin = bucket->bins;
    for (i = 0; i < bucket->nused; i++, bin++) {
      vn = &mris->vertices[bin->fno];
      if (bin->fno == Gdiag_no) {
        DiagBreak();
      }
      if (vn->ripflag) {
        continue;
      }
      switch (which) {
        default:
        case ORIGINAL_VERTICES:
          mrisComputeOrigNormal(mris, bin->fno, norm);
          dx = x - vn->origx;
          dy = y - vn->origy;
          dz = z - vn->origz;
          break;
        case WHITE_VERTICES:
          mrisComputeWhiteNormal(mris, bin->fno, norm);
          dx = x - vn->whitex;
          dy = y - vn->whitey;
          dz = z - vn->whitez;
          break;
        case PIAL_VERTICES:
          mrisComputePialNormal(mris, bin->fno, norm);
          dx = x - vn->pialx;
          dy = y - vn->pialy;
          dz = z - vn->pialz;
          break;
      }
      dot = dx * norm[0] + dy * norm[1] + dz * norm[2];
      if (sgn * dot > dot_thresh) {
        continue;
      }
      if (dot < 0 && vno == Gdiag_no) {
        DiagBreak();
      }
      if (dot > MAX_NEG_RATIO) {
        dot = MAX_NEG_RATIO;
      }
      else if (dot < -MAX_NEG_RATIO) {
        dot = -MAX_NEG_RATIO;
      }
#if 0
      scale = l_repulse / (1.0+exp(NEG_AREA_K*dot)) ;
#else
      scale = l_repulse * pow(1.0 - (double)dot, 15.0);
#endif
      if (scale > max_scale) {
        max_scale = scale;
        max_vno = bin->fno;
        max_dot = dot;
      }
      sx += (scale * v->nx);
      sy += (scale * v->ny);
      sz += (scale * v->nz);
    }

    v->dx += sgn * sx;
    v->dy += sgn * sy;
    v->dz += sgn * sz;
    if (vno == Gdiag_no) {
      vn = &mris->vertices[max_vno];
      dx = x - vn->x;
      dy = y - vn->y;
      dz = z - vn->z;

      fprintf(stdout, "v %d inside repulse term:  (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
      fprintf(stdout, "max_scale @ %d = %2.2f, max dot = %2.2f\n", max_vno, max_scale, max_dot);
    }
    MHTrelBucket(&bucket);
  }
  return (NO_ERROR);
}

/*!
  \fn double mrisComputeIntensityError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
  \brief Computes the sum of the squares of the value at a vertex minus the v->val.
   Ignores ripped vertices or any with v->val<0. Does not normalize by the number
   of vertices. Basically same computation as mrisRmsValError() but that func
   does normalize.
*/
double mrisComputeIntensityError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno,nhits;
  VERTEX *v;
  double val0, xw, yw, zw;
  double sse, del0;

  if (FZERO(parms->l_intensity))
    return (0.0f);

  nhits = 0;
  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0)
      continue;
    nhits++;
    // Sample mri_brain at vertex
    MRISvertexToVoxel(mris, v, parms->mri_brain, &xw, &yw, &zw);
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0);
    del0 = v->val - val0;
    sse += (del0 * del0);
  }
  //printf("mrisComputeIntensityError() %f %d\n",sse,nhits);
  return (sse);
}
/*! -----------------------------------------------------
  \fn static double mrisComputeTargetLocationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
  \brief Computes the distance squared between v->{xyz} and v->targ{xyz} and sums up over
  all unripped vertices. See also mrisRmsDistanceError(mris).
  ------------------------------------------------------*/
double mrisComputeTargetLocationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno, max_vno;
  VERTEX *v;
  double dx, dy, dz;
  double sse, mag, max_mag, last_mag;
  static double last_error[500000];

  if (FZERO(parms->l_location)) return (0.0f);

  last_mag = max_mag = 0;
  max_vno = -1;
  sse = 0.0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    dx = v->x - v->targx;
    dy = v->y - v->targy;
    dz = v->z - v->targz;

    mag = dx * dx + dy * dy + dz * dz;
    if (mag > 50) DiagBreak();

    if (!devFinite(mag)) DiagBreak();

    // Not sure what this bit of code is for since it does not affect
    // the output of the function. Maybe just a diagnosis.  This is
    // not thread safe because of the static, but maybe that does not
    // matter. 
    if (mag > last_error[vno]) {
      if (mag > max_mag) {
        last_mag = last_error[vno];
        max_mag = mag;
        max_vno = vno;
      }
    }
    last_error[vno] = mag;
    sse += mag;
  }

  if (last_mag > 0) DiagBreak();
  return (sse);
}

/*!
  \fn static double mrisRmsDistanceError(MRI_SURFACE *mris)
  Computes the RMS of the distance error (v->{xyz} - v->targ{xyz})
  over all unripped vertices. See also  mrisComputeTargetLocationError().
*/
double mrisRmsDistanceError(MRI_SURFACE *mris)
{
  INTEGRATION_PARMS parms;
  double rms;

  memset(&parms, 0, sizeof(parms));
  parms.l_location = 1;
  rms = mrisComputeTargetLocationError(mris, &parms);
  return (sqrt(rms / MRISvalidVertices(mris)));
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeIntensityGradientError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  float x, y, z;
  double mag0, xw, yw, zw, dx, dy, dz;
  double sse, del0;

  if (FZERO(parms->l_grad)) {
    return (0.0f);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

// MRIworldToVoxel(parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
#if 0
    MRISsurfaceRASToVoxelCached(mris, parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
#else
    MRISsurfaceRASToVoxelCached(mris, parms->mri_smooth, x, y, z, &xw, &yw, &zw);
#endif
    MRIsampleVolumeGradient(parms->mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag0 = sqrt(dx * dx + dy * dy + dz * dz);

    del0 = v->mean - mag0;
    sse += (del0 * del0);
  }

  return (sse);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeSphereError(MRI_SURFACE *mris, double l_sphere, double r0)
{
  int vno;
  double sse, x0, y0, z0;

  if (FZERO(l_sphere)) {
    return (0.0f);
  }

  x0 = (mris->xlo + mris->xhi) / 2.0f;
  y0 = (mris->ylo + mris->yhi) / 2.0f;
  z0 = (mris->zlo + mris->zhi) / 2.0f;

  // This code appears to have, but does not have, a numeric stability problem
  // Typical inputs are
  //  x:7.50467 y:-43.4641 z:-9.50775 r:45.1203 del:154.88
  // But typical results are
  //  sse:4.33023e+09
  // I tried summing in groups of 1024 numbers then summing the sums, and got the same answer in %g format
  //        /Bevin
  //
  sse = 0.0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX *v;
    double del, x, y, z, r;

    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = (double)v->x - x0;
    y = (double)v->y - y0;
    z = (double)v->z - z0;
    r = sqrt(x * x + y * y + z * z);

    del = r0 - r;
    sse += (del * del);
#if 0
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d sphere term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
#endif
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (sse);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size, int max_nbrs)
{
  int vno, n, nvertices;
  double dist_scale, pct, dist, odist, mean, mean_error, smean, total_mean_error, total_mean;

  MRIScomputeMetricProperties(mris);
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }

  total_mean = total_mean_error = mean = 0.0;
  for (pct = 0.0, nvertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      v->val = 0.0f;
      continue;
    }
    for (smean = mean_error = mean = 0.0, n = 0; n < vt->vtotal; n++) {
      nvertices++;
      dist = dist_scale * v->dist[n];
      odist = v->dist_orig[n];
#if 0
      mean += (dist - odist) * (dist - odist) ;
#else
      mean += odist;
#endif
      total_mean += odist;
      smean += dist - odist;
#define USE_FABS 1
#if USE_FABS
      mean_error += fabs(dist - odist);
      total_mean_error += fabs(dist - odist);
#else
      mean_error += (dist - odist) * (dist - odist);
      total_mean_error += (dist - odist) * (dist - odist);
#endif
      if (!FZERO(odist)) {
        pct += fabs(dist - odist) / odist;
      }
    }
    mean /= (double)vt->vtotal;
#if USE_FABS
    mean_error /= (double)vt->vtotal;
#else
    mean_error = sqrt(mean_error / (double)v->vtotal);
#endif
#if 0
    if (smean < 0.0f)
    {
      mean_error *= -1.0f ;
    }
#endif
    v->val = mean_error / mean;
  }

#if USE_FABS
  total_mean_error /= (double)nvertices;
#else
  total_mean_error = sqrt(total_mean_error / (double)nvertices);
#endif
  total_mean /= (double)nvertices;
  total_mean_error /= total_mean;
  fprintf(stdout, "mean dist = %2.3f, rms error = %2.2f%%\n", total_mean, 100.0 * total_mean_error);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static double mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris)
{
  double area_scale;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  double sse;

#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
  int trial; 
  double sse_trial0;
  for (trial = 0; trial < 2; trial++) {
#endif

  sse = 0;
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_REPRODUCIBLE
  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse  ROMP_PARTIALSUM(0)

#else
  int fno;
  
  ROMP_PF_begin     // mris_register
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
  #pragma omp parallel for if(trial==0) reduction(+ : sse)
#else
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse)
#endif
#endif
  for (fno = 0; fno < mris->nfaces; fno++) {
    ROMP_PFLB_begin

#endif
    
    double error, ratio;
    FACE *face;

    face = &mris->faces[fno];
    if (face->ripflag) {
      ROMP_PF_continue;
    }
#define SCALE_NONLINEAR_AREA 0
#if SCALE_NONLINEAR_AREA
    if (!FZERO(face->orig_area)) {
      ratio = area_scale * face->area / face->orig_area;
    }
    else {
      ratio = 0.0f;
    }
#else
    ratio = area_scale * face->area;
#endif
    if (ratio > MAX_NEG_RATIO) {
      ratio = MAX_NEG_RATIO;
    }
    else if (ratio < -MAX_NEG_RATIO) {
      ratio = -MAX_NEG_RATIO;
    }
#if 0
    error = (1.0 / NEG_AREA_K) * log(1.0+exp(-NEG_AREA_K*ratio)) ;
#else
    error = (log(1.0 + exp(NEG_AREA_K * ratio)) / NEG_AREA_K) - ratio;
#endif

    sse += error;
    if (!isfinite(sse) || !isfinite(error)) {
      ErrorExit(ERROR_BADPARM, "nlin area sse not finite at face %d!\n", fno);
    }
    
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_REPRODUCIBLE

    #undef sse
  #include "romp_for_end.h"
#else
  
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
    if (trial == 0) {
        sse_trial0 = sse;
    } else { 
        if (sse_trial0 != sse) {
            fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
               sse_trial0, sse, sse_trial0-sse);
        }
    }
  } // trial
#endif

  return (sse);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static double mrisComputeNonlinearDistanceSSE(MRI_SURFACE *mris)
{
  int vno, n, nvertices, max_v, max_n;
  double dist_scale, sse_dist, delta, v_sse, max_del, ratio;

#if METRIC_SCALE
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else if (mris->status == MRIS_PARAMETERIZED_SPHERE) {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
  else
    dist_scale = mris->neg_area < mris->total_area ? sqrt(mris->orig_area / (mris->total_area - mris->neg_area))
                                                   : sqrt(mris->orig_area / mris->total_area);
#else
  dist_scale = 1.0;
#endif
  max_del = -1.0;
  max_v = max_n = -1;
  for (sse_dist = 0.0, nvertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    nvertices++;
    for (v_sse = 0.0, n = 0; n < vt->vtotal; n++) {
      if (FZERO(v->dist_orig[n])) {
        continue;
      }
      ratio = dist_scale * v->dist[n] / v->dist_orig[n];
      delta = log(1 + exp(ratio));
      v_sse += delta;
      if (!isfinite(delta) || !isfinite(v_sse)) {
        DiagBreak();
      }
    }
    sse_dist += v_sse;
    if (!isfinite(sse_dist) || !isfinite(v_sse)) {
      DiagBreak();
    }
  }

  return (sse_dist);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Fit a 2-d quadratic to the surface locally and compute the SSE as
  the square of the constant term (the distance the quadratic fit surface
  is from going through the central vertex)
  ------------------------------------------------------*/
static double mrisComputeQuadraticCurvatureSSE(MRI_SURFACE *mris, double l_curv)            // BEVIN mris_make_surfaces 3
{
  if (FZERO(l_curv)) {
    return (NO_ERROR);
  }

  mrisComputeTangentPlanes(mris);
  
  typedef struct Reused {
    VECTOR * v_n  ;
    VECTOR * v_P  ;
    VECTOR * v_e1 ;
    VECTOR * v_e2 ;
    VECTOR * v_nbr;
  } Reused;
  
#ifdef HAVE_OPENMP
  int const maxThreads = omp_get_max_threads();
#else
  int const maxThreads = 1;
#endif
  Reused* reusedByThread = (Reused*)calloc(maxThreads, sizeof(Reused));
  
  double sse = 0.0;

#ifdef BEVIN_MRISCOMPUTEQUADRATICCURVATURESSE_REPRODUCIBLE

  #define ROMP_VARIABLE       vno 
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse ROMP_PARTIALSUM(0)

#else
  
  int vno;
  
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+:sse)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif

#ifdef HAVE_OPENMP
    int const tid = omp_get_thread_num();
#else
    int const tid = 0;
#endif

    Reused* reused = reusedByThread + tid;
    #define REUSE(NAME,DIM) \
        VECTOR* NAME = reused->NAME; if (!NAME) NAME = reused->NAME = VectorAlloc(DIM, MATRIX_REAL);
    REUSE(v_n   ,3)
    REUSE(v_P   ,5)
    REUSE(v_e1  ,3)
    REUSE(v_e2  ,3)
    REUSE(v_nbr ,3)
    #undef REUSE

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    VECTOR* v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL);    /* heights above TpS */
    VECTOR_LOAD(v_n,  v->nx,  v->ny,  v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);

    MATRIX* m_X = MatrixAlloc(vt->vtotal, 5, MATRIX_REAL); /* 2-d quadratic fit */
    int n;
    for (n = 0; n < vt->vtotal; n++) /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      float ui = V3_DOT(v_e1, v_nbr);
      float vi = V3_DOT(v_e2, v_nbr);

      *MATRIX_RELT(m_X, n + 1, 1) = ui * ui;
      *MATRIX_RELT(m_X, n + 1, 2) = vi * vi;
      *MATRIX_RELT(m_X, n + 1, 3) = ui;
      *MATRIX_RELT(m_X, n + 1, 4) = vi;
      *MATRIX_RELT(m_X, n + 1, 5) = 1;
    }

    MATRIX* m_X_inv = MatrixPseudoInverse(m_X, NULL);
    if (!m_X_inv) {
      MatrixFree(&m_X);
      VectorFree(&v_Y);
      continue;
    }
    
    v_P = MatrixMultiply(m_X_inv, v_Y, v_P);
    //oat a = VECTOR_ELT(v_P, 1);
    float e = VECTOR_ELT(v_P, 5);

    sse += e * e;
    if (vno == Gdiag_no) printf("v %d: e=%2.2f, curvature sse %2.2f\n", vno, e, e * e);

    MatrixFree(&m_X);
    MatrixFree(&m_X_inv);
    VectorFree(&v_Y);
    
#ifdef BEVIN_MRISCOMPUTEQUADRATICCURVATURESSE_REPRODUCIBLE
    #undef sse
  #include "romp_for_end.h"
#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

  { int tid;
    for (tid = 0; tid < maxThreads; tid++) {
      Reused* reused = reusedByThread + tid;
      if (reused->v_n)   VectorFree(&reused->v_n);
      if (reused->v_e1)  VectorFree(&reused->v_e1);
      if (reused->v_e2)  VectorFree(&reused->v_e2);
      if (reused->v_nbr) VectorFree(&reused->v_nbr);
      if (reused->v_P)   VectorFree(&reused->v_P);
  } }
  
  return (NO_ERROR);
}


static int mrisCreateLikelihoodHistograms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int x, y, z, wlabel, plabel;
  VECTOR *v_brain, *v_hires;
  MATRIX *m_hires_to_brain;
  MRI *mri_pial;
  double xv, yv, zv, val, dist;

  if (parms->mri_white == NULL)  // white isn't moving, so only have to do it once
  {
    parms->mri_white = MRIupsampleN(parms->mri_brain, NULL, 3);
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, WHITE_VERTICES);
    MRISfillInterior(mris, parms->mri_white->xsize, parms->mri_white);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
  }
  mri_pial = MRIclone(parms->mri_white, NULL);
  MRISfillInterior(mris, mri_pial->xsize, mri_pial);
  if (parms->mri_labels == NULL) {
    parms->mri_labels = MRIclone(parms->mri_white, NULL);
  }
  parms->mri_dist = MRIdistanceTransform(mri_pial, parms->mri_dist, 1, 5 / mri_pial->xsize, DTRANS_MODE_SIGNED, NULL);

  parms->h_wm = HISTOinit(parms->h_wm, 256, 0, 255);
  parms->h_gm = HISTOinit(parms->h_gm, 256, 0, 255);
  parms->h_nonbrain = HISTOinit(parms->h_nonbrain, 256, 0, 255);
  m_hires_to_brain = MRIgetVoxelToVoxelXform(parms->mri_labels, parms->mri_brain);

  v_brain = VectorAlloc(4, MATRIX_REAL);
  v_hires = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_brain, 4) = VECTOR_ELT(v_hires, 4) = 1.0;
  for (x = 0; x < mri_pial->width; x++) {
    V3_X(v_hires) = x;
    for (y = 0; y < mri_pial->height; y++) {
      V3_Y(v_hires) = y;
      for (z = 0; z < mri_pial->height; z++) {
        V3_Z(v_hires) = z;
        MatrixMultiply(m_hires_to_brain, v_hires, v_brain);
        xv = V3_X(v_brain);
        yv = V3_Y(v_brain);
        zv = V3_Z(v_brain);
        if (MRIindexNotInVolume(parms->mri_brain, xv, yv, zv)) {
          val = 0;
        }
        else {
          MRIsampleVolume(parms->mri_brain, xv, yv, zv, &val);
        }

        wlabel = MRIgetVoxVal(parms->mri_white, x, y, z, 0);
        plabel = MRIgetVoxVal(mri_pial, x, y, z, 0);
        dist = MRIgetVoxVal(parms->mri_dist, x, y, z, 0);
        if (dist > 3) {
          continue;  // don't consider the millions of voxels far from the surface
        }

        if (wlabel) {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_WHITE_INTERIOR);
          HISTOaddSample(parms->h_wm, val, 0, 255);
        }
        else if (plabel) {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_PIAL_INTERIOR);
          HISTOaddSample(parms->h_gm, val, 0, 255);
        }
        else {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_NONBRAIN);
          HISTOaddSample(parms->h_nonbrain, val, 0, 255);
        }
      }
    }
  }
  HISTOmakePDF(parms->h_nonbrain, parms->h_nonbrain);
  HISTOmakePDF(parms->h_wm, parms->h_wm);
  HISTOmakePDF(parms->h_gm, parms->h_gm);
  MatrixFree(&m_hires_to_brain);
  MatrixFree(&v_brain);
  MatrixFree(&v_hires);
  MRIfree(&mri_pial);
  return (NO_ERROR);
}


double mrisComputeHistoNegativeLikelihood(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double likelihood, entropy;
  int x, y, z, label;
  VECTOR *v_brain, *v_hires;
  MATRIX *m_brain_to_hires;
  double xv, yv, zv, val, pval;

  if (DZERO(parms->l_histo)) return (NO_ERROR);

  mrisCreateLikelihoodHistograms(mris, parms);

  v_brain = VectorAlloc(4, MATRIX_REAL);
  v_hires = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_brain, 4) = VECTOR_ELT(v_hires, 4) = 1.0;
  m_brain_to_hires = MRIgetVoxelToVoxelXform(parms->mri_brain, parms->mri_labels);
  for (likelihood = 0, x = 0; x < parms->mri_brain->width; x++) {
    V3_X(v_brain) = x;
    for (y = 0; y < parms->mri_brain->height; y++) {
      V3_Y(v_brain) = y;
      for (z = 0; z < parms->mri_brain->height; z++) {
        V3_Z(v_brain) = z;
        MatrixMultiply(m_brain_to_hires, v_brain, v_hires);
        xv = V3_X(v_hires);
        yv = V3_Y(v_hires);
        zv = V3_Z(v_hires);
        if (MRIindexNotInVolume(parms->mri_labels, xv, yv, zv)) {
          label = MRI_NONBRAIN;
        }
        else {
          label = MRIgetVoxVal(parms->mri_labels, nint(xv), nint(yv), nint(zv), 0);
        }

        val = MRIgetVoxVal(parms->mri_brain, x, y, z, 0);
        switch (label) {
          default:
          case MRI_NONBRAIN:
            pval = HISTOgetCount(parms->h_nonbrain, val);
            break;
          case MRI_WHITE_INTERIOR:
            pval = HISTOgetCount(parms->h_wm, val);
            break;
          case MRI_PIAL_INTERIOR:
            pval = HISTOgetCount(parms->h_gm, val);
            break;
        }
        likelihood += pval;
      }
    }
  }

  entropy = HISTOgetEntropy(parms->h_nonbrain) + HISTOgetEntropy(parms->h_gm);
  MatrixFree(&v_brain);
  MatrixFree(&v_hires);
  MatrixFree(&m_brain_to_hires);

  return (entropy);
}

int mrisComputeHistoTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  double xv, yv, zv, val, p, dx, dy, dz, d, d2, x, y, z, best_d, best_p;
  int num;

  if (DZERO(parms->l_histo)) return (NO_ERROR);

  mrisCreateLikelihoodHistograms(mris, parms);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    best_p = 0;
    best_d = 0;
    for (d = -2; d < .2; d += .1) {
      for (num = 0, p = 0, d2 = -2; d2 < .2; d2 += .1) {
        x = v->x + d2 * v->nx;
        y = v->y + d2 * v->ny;
        z = v->z + d2 * v->nz;
        MRISsurfaceRASToVoxelCached(mris, parms->mri_brain, x, y, z, &xv, &yv, &zv);
        MRIsampleVolume(parms->mri_brain, xv, yv, zv, &val);
        if (d2 < d) {
          p += HISTOgetCount(parms->h_gm, val);
        }
        else {
          p += HISTOgetCount(parms->h_nonbrain, val);
        }
        num++;
      }
      p /= num;
      if (p > best_p) {
        best_p = p;
        best_d = d;
      }
    }

    dx = v->nx * best_d * parms->l_histo;
    dy = v->ny * best_d * parms->l_histo;
    dz = v->nz * best_d * parms->l_histo;
    if (vno == Gdiag_no) {
      printf(
          "histoTerm for v %d: val = %2.1f, optimal distance=%2.3f, p = %f moving by (%2.3f, %2.3f, %2.3f), dot = "
          "%2.3f\n",
          vno,
          val,
          best_d,
          best_p,
          dx,
          dy,
          dz,
          dx * v->nx + dy * v->ny + dz * v->nz);

      DiagBreak();
      if (Gx < 0) {
        Gx = x;
        Gy = y;
        Gz = z;
      }
    }
    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
  }
  return (NO_ERROR);
}


double mrisComputeNegativeLogPosterior(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox)
{
  MRI *mri = parms->mri_brain;
  double sse = 0.0, ll, wm_frac, gm_frac, out_frac, Ig, Ic, pval;
  float vmin, vmax, val;
  HISTOGRAM *hin, *hout;
  int x, y, z, nvox, label;
  MRI *mri_ll = NULL;
  static double last_sse = 0.0;

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    mri_ll = MRIcloneDifferentType(mri, MRI_FLOAT);
    sprintf(fname, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t);
    MRIwrite(parms->mri_volume_fractions, fname);
  }

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }
  if (parms->hgm == NULL) {
    printf("creating intensity histograms\n");
    parms->hgm = hin = HISTOinit(parms->hgm, 256, (double)vmin, (double)vmax);
    parms->hout = hout = HISTOinit(parms->hout, 256, (double)vmin, (double)vmax);
    MRISclearMarks(mris);

    // build histogram estimates of PDFs of interior and exterior of ribbon
    for (x = 0; x < mri->width; x++)
      for (y = 0; y < mri->height; y++)
        for (z = 0; z < mri->depth; z++) {
          if (Gx == x && Gy == y && Gz == z) DiagBreak();
          if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
          val = MRIgetVoxVal(mri, x, y, z, 0);
          if (FZERO(val)) continue;

          wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
          gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1);
          out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2);
          if (parms->mri_aseg) {
            label = MRIgetVoxVal(parms->mri_aseg, x, y, z, 0);
            if (FZERO(gm_frac) && IS_CORTEX(label))  // aseg thinks it is but outside ribbon - ambiguous
              continue;
          }
          HISTOaddFractionalSample(hout, val, 0, 0, out_frac);
          HISTOaddFractionalSample(hin, val, 0, 0, gm_frac);
        }

    HISTOmakePDF(hin, hin);
    HISTOmakePDF(hout, hout);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      sprintf(fname, "hin.%s.%3.3d.plt", parms->base_name, parms->t);
      HISTOplot(hin, fname);
      sprintf(fname, "hout.%s.%3.3d.plt", parms->base_name, parms->t);
      HISTOplot(hout, fname);
    }
  }
  else  // use previously computed ones
  {
    hin = parms->hgm;
    hout = parms->hout;
  }

  for (sse = 0.0, nvox = x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        if (Gx == x && Gy == y && Gz == z) DiagBreak();
        if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (FZERO(val)) continue;
        wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
        if (wm_frac > 0) continue;

        gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1);
        out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2);

        if (FZERO(out_frac))  // all gm
          pval = HISTOgetCount(hin, val);
        else if (FZERO(gm_frac))
          pval = HISTOgetCount(hout, val);
        else  // partial volume voxel
        {
          for (pval = 0.0, Ig = 0; Ig <= 256; Ig++) {
            Ic = (val - gm_frac * Ig) / out_frac;
            pval += HISTOgetCount(hout, Ic) * HISTOgetCount(hin, Ig);
          }
        }

        if (pval > 1 || pval < 0) DiagBreak();
        if (DZERO(pval)) pval = 1e-20;
        ll = -log10(pval);
        sse += ll;
        nvox++;
        if (mri_ll) MRIsetVoxVal(mri_ll, x, y, z, 0, ll);
        if (Gx == x && Gy == y && Gz == z)
          printf("voxel(%d, %d, %d) = %d, vfracs = (%2.1f, %2.1f, %2.1f), ll = %2.1f\n",
                 x,
                 y,
                 z,
                 nint(val),
                 wm_frac,
                 gm_frac,
                 out_frac,
                 ll);
      }

  if (!FZERO(last_sse) && sse > last_sse) DiagBreak();
  if (mri_ll) {
    char fname[STRLEN];
    sprintf(fname, "%s.ll.%4.4d.mgz", parms->base_name, parms->t);
    printf("writing log likelihood volume to %s\n", fname);
    MRIwrite(mri_ll, fname);
    MRIfree(&mri_ll);
  }
  //  HISTOfree(&hin) ; HISTOfree(&hout) ;
  if (pnvox) *pnvox = nvox;
  if (Gdiag_no >= 0) printf("E_logPosterior, %3.3d: %.1f (nvox=%d)\n", parms->t, sse, nvox);

  last_sse = sse;  // diagnostics
  return (sse);
}


double mrisComputeNegativeLogPosterior2D(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox)
{
  MRI *mri = parms->mri_brain;
  MHT *mht;
  double sse = 0.0, ll, wm_frac, pval, dist, vdist, xs, ys, zs;
  float vmin, vmax, val;
  int x, y, z, nvox, label, xd, yd, zd, vno;
  VERTEX *v;
  MRI *mri_ll = NULL;
  static double last_sse = 0.0;
  HISTOGRAM2D *hs;
  MATRIX *m_vox2vox;
  VECTOR *v1, *v2;

  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);

  if (parms->mri_dtrans == NULL) {
    int nvox;
    MRI *mri_tmp;

    nvox = nint(
        MAX(MAX((mri->xsize / parms->resolution), (mri->ysize / parms->resolution)), (mri->zsize / parms->resolution)));
    mri_tmp = MRIcloneDifferentType(mri, MRI_FLOAT);
    parms->mri_dtrans = MRIupsampleN(mri_tmp, NULL, nvox);
    MRIfree(&mri_tmp);
    MRIScomputeDistanceToSurface(mris, parms->mri_dtrans, parms->resolution);
    if (parms->t > 0) DiagBreak();
    DiagBreak();
  }

  m_vox2vox = MRIgetVoxelToVoxelXform(mri, parms->mri_dtrans);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    mri_ll = MRIcloneDifferentType(mri, MRI_FLOAT);
    sprintf(fname, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t);
    MRIwrite(parms->mri_volume_fractions, fname);
  }

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }

  if (parms->h2d == NULL) {
    printf("creating 2D intensity histograms\n");
    parms->h2d = HISTO2Dinit(parms->h2d, 128, 101, (double)vmin, (double)vmax, -10, 10);
    MRISclearMarks(mris);

    // build histogram estimates of PDFs of interior and exterior of ribbon
    for (x = 0; x < mri->width; x++)
      for (y = 0; y < mri->height; y++)
        for (z = 0; z < mri->depth; z++) {
          if (Gx == x && Gy == y && Gz == z) DiagBreak();
          if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
          val = MRIgetVoxVal(mri, x, y, z, 0);
          wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
          if (wm_frac > 0) continue;

          V3_X(v1) = x;
          V3_Y(v1) = y;
          V3_Z(v1) = z;
          MatrixMultiply(m_vox2vox, v1, v2);
          xd = nint(V3_X(v2));
          yd = nint(V3_Y(v2));
          zd = nint(V3_Z(v2));
          if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd)) continue;
          dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
          if (FZERO(val) && (dist > 1 || dist < 0))  // don't allow 0s inside ribbon or too far out
            continue;
          if (FZERO(val)) continue;
          if (FZERO(val) && dist < 0) DiagBreak();

          // update distance with continuum measure
          MRIvoxelToSurfaceRAS(mri, x, y, z, &xs, &ys, &zs);
          MHTfindClosestVertexGeneric(mht, mris, xs, ys, zs, 10, 4, &v, &vno, &vdist);
          if (v != NULL)  // compute distance from surface in normal direction
          {
            double dx, dy, dz;

            dx = xs - v->x;
            dy = ys - v->y;
            dz = zs - v->z;
            dist = dx * v->nx + dy * v->ny + dz * v->nz;
          }

          if (parms->mri_aseg) {
            label = MRIgetVoxVal(parms->mri_aseg, x, y, z, 0);
            if (dist > 1 && IS_CORTEX(label))  // aseg thinks it is but outside ribbon - ambiguous
              continue;
          }
          HISTO2DaddSample(parms->h2d, val, dist, 0, 0, 0, 0);
        }
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];

      sprintf(fname, "h.%s.%3.3d.plt", parms->base_name, parms->t);
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(parms->h2d, fname);
    }

    //    hs = HISTO2DsoapBubbleZeros(parms->h2d, NULL, 100) ;
    hs = HISTO2DmakePDF(parms->h2d, NULL);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      sprintf(fname, "h.%s.%3.3d.pdf.plt", parms->base_name, parms->t);
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(hs, fname);
    }
    //    HISTO2DsmoothBins2(hs, parms->h2d, 5) ;
    if (parms->h2d_out) HISTO2Dfree(&parms->h2d_out);
    parms->h2d_out = HISTO2DsmoothAnisotropic(hs, NULL, .25, 1);
    HISTO2Dsmooth(hs, parms->h2d, 1);
    HISTO2DmakePDF(parms->h2d_out, hs);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      sprintf(fname, "h.%s.%3.3d.pdf.smooth.plt", parms->base_name, parms->t);
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(hs, fname);
    }
    HISTO2Dfree(&parms->h2d);
    parms->h2d = hs;

    if (Gx > 0) {
      int b1, b2;
      float val, c;

      x = Gx;
      y = Gy;
      z = Gz;
      V3_X(v1) = x;
      V3_Y(v1) = y;
      V3_Z(v1) = z;
      MatrixMultiply(m_vox2vox, v1, v2);
      xd = nint(V3_X(v2));
      yd = nint(V3_Y(v2));
      zd = nint(V3_Z(v2));
      if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd))
        dist = 10;
      else
        dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
      val = MRIgetVoxVal(mri, x, y, z, 0);
      b1 = HISTO2DfindBin1(parms->h2d, val);
      b2 = HISTO2DfindBin2(parms->h2d, dist);
      c = HISTO2DgetCount(parms->h2d, val, dist);
      printf("voxel (%d, %d, %d) = %2.0f, dist=%2.2f, bin=%d,%d, count=%f\n", Gx, Gy, Gz, val, dist, b1, b2, c);

      DiagBreak();
    }
  }
  else  // use previously computed ones
  {
  }

  for (sse = 0.0, nvox = x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        if (Gx == x && Gy == y && Gz == z) DiagBreak();
        if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
        val = MRIgetVoxVal(mri, x, y, z, 0);
        wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
        if (wm_frac > 0) continue;

        V3_X(v1) = x;
        V3_Y(v1) = y;
        V3_Z(v1) = z;
        MatrixMultiply(m_vox2vox, v1, v2);
        xd = nint(V3_X(v2));
        yd = nint(V3_Y(v2));
        zd = nint(V3_Z(v2));
        if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd))
          dist = 10;
        else
          dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
        if (FZERO(val) && dist > 1) continue;
        if (FZERO(val)) continue;
        if (FZERO(val)) DiagBreak();

        pval = HISTO2DgetCount(parms->h2d, val, dist);
        if (pval > 1 || pval < 0) DiagBreak();
        if (DZERO(pval)) pval = 1e-20;
        ll = -log10(pval);
        if (!isfinite(ll)) DiagBreak();
        sse += ll;
        nvox++;
        if (mri_ll) MRIsetVoxVal(mri_ll, x, y, z, 0, ll);
        if (Gx == x && Gy == y && Gz == z) {
          printf("voxel(%d, %d, %d) = %d, dist=%2.2f, ll = %2.1f, bins = %d, %d\n",
                 x,
                 y,
                 z,
                 nint(val),
                 dist,
                 ll,
                 HISTO2DfindBin1(parms->h2d, val),
                 HISTO2DfindBin2(parms->h2d, dist));
        }
      }

  if (!FZERO(last_sse) && sse > last_sse) DiagBreak();
  if (mri_ll) {
    char fname[STRLEN];
    sprintf(fname, "%s.ll.%4.4d.mgz", parms->base_name, parms->t);
    printf("writing log likelihood volume to %s\n", fname);
    MRIwrite(mri_ll, fname);
    MRIfree(&mri_ll);
  }
  //  HISTOfree(&hin) ; HISTOfree(&hout) ;
  if (pnvox) *pnvox = nvox;
  if (Gdiag_no >= 0) printf("E_logPosterior, %3.3d: %.1f (nvox=%d)\n", parms->t, sse, nvox);

  MHTfree(&mht);
  MatrixFree(&m_vox2vox);
  VectorFree(&v1);
  VectorFree(&v2);
  last_sse = sse;  // diagnostics
  return (sse);
}


static double vlst_loglikelihood(
    MRI_SURFACE *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM *hin, HISTOGRAM *hout)
{
  double ll = 0.0, dot, dx, dy, dz, pval, dist, Ig, Ic, gm_frac, out_frac;
  int i;
  float val;
  VERTEX *v;
  double xs, ys, zs;

  v = &mris->vertices[vno];
  xs = v->x + displacement * v->nx;
  ys = v->y + displacement * v->ny;
  zs = v->z + displacement * v->nz;
  for (i = 0; i < vl->nvox; i++) {
    dx = vl->xd[i] - xs;
    dy = vl->yd[i] - ys;
    dz = vl->zd[i] - zs;
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    dot = dx * v->nx + dy * v->ny + dz * v->nz;
    val = MRIgetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    if (dist < .5)  // distance to center<.5 --> distance to edge <1
    {
      if (dot > 0) {
        out_frac = dist + .5;
        gm_frac = 1 - out_frac;
      }
      else {
        gm_frac = dist + .5;
        out_frac = 1 - gm_frac;
      }
      for (pval = 0.0, Ig = 0; Ig <= 256; Ig++) {
        Ic = (val - gm_frac * Ig) / out_frac;
        pval += HISTOgetCount(hout, Ic) * HISTOgetCount(hin, Ig);
      }
    }
    else if (dot > 0)  // outside surface
      pval = HISTOgetCount(hout, val);
    else  // inside the surface
      pval = HISTOgetCount(hin, val);
    if (DZERO(pval)) pval = 1e-10;
    ll += -log(pval);
  }

  return (ll);
}

static double vlst_loglikelihood2D(
    MRI_SURFACE *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM2D *h, FILE *fp)
{
  double ll = 0.0, dot, dx, dy, dz, pval, dist;
  int i;
  float val;
  VERTEX *v;
  double xs, ys, zs;

  if (fp) fprintf(fp, "%f ", displacement);

  v = &mris->vertices[vno];
  xs = v->x + displacement * v->nx;
  ys = v->y + displacement * v->ny;
  zs = v->z + displacement * v->nz;
  for (i = 0; i < vl->nvox; i++) {
    dx = vl->xd[i] - xs;
    dy = vl->yd[i] - ys;
    dz = vl->zd[i] - zs;
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    dot = dx * v->nx + dy * v->ny + dz * v->nz;
    val = MRIgetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    pval = HISTO2DgetCount(h, val, dot);
    if (DZERO(pval)) pval = 1e-10;
    if (fp) fprintf(fp, "%d %2.2f %2.2f ", (int)val, dot, -log(pval));
    ll += -log(pval);
  }

  if (fp) fprintf(fp, "\n");
  return (ll);
}


#define MAX_VOXELS 1500
#define MAX_DISPLACEMENT 5
#define DISPLACEMENT_DELTA 0.1


int mrisComputePosterior2DTerm(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  static const double MAX_MOVEMENT = 1.0;

  MRI *mri = parms->mri_brain;
  FILE *fp, *fp2;
  double dist, xs, ys, zs, dn, best_dn, best_ll, ll;
  float vmin, vmax, val, wdist;
  int x, y, z, vno, n;
  MATRIX *m_vox2vox;
  MHT *mht;
  VOXEL_LIST **vl, **vl2;
  static MRI *mri_white_dist = NULL;

  if (FZERO(parms->l_map2d)) return (NO_ERROR);
  vl  = vlst_alloc(mris, MAX_VOXELS);
  vl2 = vlst_alloc(mris, MAX_VOXELS);

  if (mri_white_dist == NULL) {
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, WHITE_VERTICES);
    mri_white_dist = MRIScomputeDistanceToSurface(mris, NULL, 0.5);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
  }
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);

  m_vox2vox = MRIgetVoxelToVoxelXform(mri, mri_white_dist);
  VECTOR* v1 = VectorAlloc(4, MATRIX_REAL);
  VECTOR* v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_white_dist, "wd.mgz");

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }
  MRISclearMarks(mris);

  // find set of voxels that are closest to each vertex
  for (x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        val = MRIgetVoxVal(mri, x, y, z, 0);
        V3_X(v1) = x;
        V3_Y(v1) = y;
        V3_Z(v1) = z;
        MatrixMultiply(m_vox2vox, v1, v2);

        if (MRIindexNotInVolume(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2))) continue;
        wdist = MRIgetVoxVal(mri_white_dist, V3_X(v2), V3_Y(v2), V3_Z(v2), 0);
        if (wdist < 0) continue;

        // add this voxel to the list of voxels of the vertex it is closest to
        MRIvoxelToSurfaceRAS(mri, x, y, z, &xs, &ys, &zs);
        VERTEX * v;
        MHTfindClosestVertexGeneric(mht, mris, xs, ys, zs, 10, 4, &v, &vno, &dist);
        if (v == NULL) continue;
        if (FZERO(val) && dist > 1) continue;
        if (vno == Gdiag_no) DiagBreak();
        v->marked++;
        VLSTadd(vl[vno], x, y, z, xs, ys, zs);
        if (x == Gx && y == Gy && z == Gz) {
          printf("voxel (%d, %d, %d) --> vertex %d\n", x, y, z, vno);
          if (Gdiag_no < 0) {
            printf("setting Gdiag_no to %d\n", vno);
            Gdiag_no = vno;
          }
        }
      }

  if (Gdiag_no >= 0) printf("%d nearest voxels found to vertex %d\n", vl[Gdiag_no]->nvox, Gdiag_no);

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
    v->marked = vl2[vno]->nvox;
    if (vno == Gdiag_no) printf("%d total voxels found close to vertex %d after nbr adding\n", vl2[vno]->nvox, vno);
  }

  vlst_free(mris, &vl);
  vl = vl2;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri, parms->mri_dtrans);  // dist to pial surface
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
      sprintf(fname, "vno%d.%d.l.vals.dat", vno, parms->t);
      fp2 = fopen(fname, "w");
    }
    else
      fp = fp2 = NULL;
    best_ll = -1e10;
    best_dn = 0;
    if (fp2) {
      int i;
      float dot, dx, dy, dz, val;
      char fname[STRLEN];

      for (i = 0; i < vl[vno]->nvox; i++) {
        val = MRIgetVoxVal(mri, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], 0);
        dx = vl[vno]->xd[i] - v->x;
        dy = vl[vno]->yd[i] - v->y;
        dz = vl[vno]->zd[i] - v->z;
        dot = dx * v->nx + dy * v->ny + dz * v->nz;
        fprintf(fp2, "%d %d %d %f %f\n", vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], dot, val);
      }
      fclose(fp2);
      sprintf(fname, "vno%d.%d.l.good.dat", vno, parms->t);
      fp2 = fopen(fname, "w");
    }
    for (dn = -MAX_DISPLACEMENT; dn <= MAX_DISPLACEMENT; dn += DISPLACEMENT_DELTA) {
      ll = -vlst_loglikelihood2D(mris, mri, vno, dn, vl[vno], parms->h2d, fp2);
      if (devIsnan(ll)) DiagBreak();
      if (fp) fprintf(fp, "%f %f\n", dn, ll);
      if (ll > best_ll) {
        best_dn = dn;
        best_ll = ll;
      }
    }
    if (fp) fclose(fp);
    if (fp2) fclose(fp2);

    if (vno == Gdiag_no && parms->h2d_out != NULL)  // diags
    {
      char fname[STRLEN];
      double best_dn = 0, best_ll = -1e10;
      sprintf(fname, "vno%d.%d.l.bad.dat", vno, parms->t);
      fp = fopen(fname, "w");
      for (dn = -MAX_DISPLACEMENT; dn <= MAX_DISPLACEMENT; dn += DISPLACEMENT_DELTA) {
        ll = -vlst_loglikelihood2D(mris, mri, vno, dn, vl[vno], parms->h2d_out, fp);
        if (devIsnan(ll)) DiagBreak();
        if (ll > best_ll) {
          best_dn = dn;
          best_ll = ll;
        }
      }
      fclose(fp);
    }

    if (vno == Gdiag_no) {
      int i;
      char fname[STRLEN];
      double dx, dy, dz, dist, dot, p;

      sprintf(fname, "vno%d.%d.vox.l.dat", vno, parms->t);
      fp = fopen(fname, "w");

      for (i = 0; i < vl[vno]->nvox; i++) {
        val = MRIgetVoxVal(mri, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], 0);
        dx = vl[vno]->xd[i] - v->x;
        dy = vl[vno]->yd[i] - v->y;
        dz = vl[vno]->zd[i] - v->z;
        dist = dx * dx + dy * dy + dz * dz;
        dot = dx * v->nx + dy * v->ny + dz * v->nz;
        if (dot < 0) dist *= -1;
        p = HISTO2DgetCount(parms->h2d, val, dot);
        fprintf(fp, "%d %d %d %d %d %f %f %f\n", vno, i, vl[vno]->xi[i], vl[vno]->yi[i], vl[vno]->zi[i], val, dist, p);
      }
      fclose(fp);

      printf("l_map: vno %d, best displacement %2.3f, ll = %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",
             vno,
             best_dn,
             best_ll,
             best_dn * v->nx * parms->l_map2d,
             best_dn * v->ny * parms->l_map2d,
             best_dn * v->nz * parms->l_map2d);
      DiagBreak();
    }
#if 1
    if (fabs(best_dn) > MAX_MOVEMENT) {
      best_dn = MAX_MOVEMENT * best_dn / fabs(best_dn);
      if (vno == Gdiag_no) printf("cropping displacement to %2.3f\n", best_dn);
    }
#endif
    v->dx += best_dn * v->nx * parms->l_map2d;
    v->dy += best_dn * v->ny * parms->l_map2d;
    v->dz += best_dn * v->nz * parms->l_map2d;
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
      v->dx += dn * v->nx * parms->l_map2d;
      v->dy += dn * v->ny * parms->l_map2d;
      v->dz += dn * v->nz * parms->l_map2d;
      v->d = dn;
      if (vno == Gdiag_no)
        printf("l_map: vno %d, soap bubble displacement %2.3f, D = (%2.3f, %2.3f, %2.3f)\n",
               vno,
               dn,
               dn * v->nx * parms->l_map2d,
               dn * v->ny * parms->l_map2d,
               dn * v->nz * parms->l_map2d);
    }
  }

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN], path[STRLEN];

    FileNamePath(mris->fname, path);
    sprintf(fname,
            "%s/%s.%d.dist.mgz",
            path,
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : mris->hemisphere == BOTH_HEMISPHERES ? "both" : "rh",
            parms->t);
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


/*-----------------------------------------------------
  Description
  MRIScomputeSSE and MRIScomputeSSEExternal
  are used for the numerical integration.
  As such, they represent the exact error function being minimized, as opposed to computeError above.
  ------------------------------------------------------*/
#include "mrisurf_deform_computeSSE.h"

double MRIScomputeSSEExternal(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, double *ext_sse)
{
  double sse;

  if (gMRISexternalSSE) {
    sse = (*gMRISexternalSSE)(mris, parms);
  }
  else {
    sse = 0;
  }
  *ext_sse = sse;
  sse = MRIScomputeSSE(mris, parms); /* throw out ext_sse
                                        as it will be recomputed */

  return (sse);
}

double mrisComputeError(MRI_SURFACE *mris,
                               INTEGRATION_PARMS *parms,
                               float *parea_rms,
                               float *pangle_rms,
                               float *pcurv_rms,
                               float *pdist_rms,
                               float *pcorr_rms)
{
  double rms, sse_area, sse_angle, sse_curv, delta, area_scale, sse_dist, sse_corr;
  int ano, fno, ntriangles, total_neighbors;
  FACE *face;
  float nv;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  sse_angle = sse_area = 0.0;
  for (ntriangles = fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
    ntriangles++;
    delta = (double)(area_scale * face->area - fNorm->orig_area);
    sse_area += delta * delta;
    for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
      delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
      sse_angle += delta * delta;
    }
    if (!isfinite(sse_area) || !isfinite(sse_angle))
      ErrorExit(ERROR_BADPARM, "sse (%f, %f) is not finite at face %d!\n", sse_area, sse_angle, fno);
  }

  sse_corr = mrisComputeCorrelationError(mris, parms, 1);
  if (!DZERO(parms->l_dist)) {
    sse_dist = mrisComputeDistanceError(mris, parms);
  }
  else {
    sse_dist = 0;
  }

#if 0
  if (!FZERO(parms->l_spring))
  {
    sse_spring = mrisComputeSpringEnergy(mris) ;
  }
  if (!FZERO(parms->l_tspring))
  {
    sse_tspring = mrisComputeTangentialSpringEnergy(mris) ;
  }
#endif
  sse_curv = mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv);

  total_neighbors = mrisCountTotalNeighbors(mris);

  nv = (float)MRISvalidVertices(mris);
  if (mris->status != MRIS_PLANE) {
    *pcurv_rms = (float)sqrt(sse_curv / nv);
  }
  else {
    *pcurv_rms = 0.0f;
  }
  *pdist_rms = (float)sqrt(sse_dist / (double)total_neighbors);
  *parea_rms = (float)sqrt(sse_area / (double)ntriangles);
  *pangle_rms = (float)sqrt(sse_angle / (double)(ntriangles * ANGLES_PER_TRIANGLE));
  *pcorr_rms = (float)sqrt(sse_corr / (double)nv);

  rms = MRIScomputeSSE(mris, parms);

#if 0
  rms =
    *pdist_rms * parms->l_dist +
    *parea_rms * parms->l_area +
    *parea_rms * parms->l_parea +
    *pangle_rms * parms->l_angle +
    *pcorr_rms * (parms->l_corr+parms->l_pcorr) +
    sqrt(sse_spring/nv) * parms->l_spring ;
#endif
  return (rms);
}


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


double mrisComputeDistanceError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  if (!(mris->dist_alloced_flags & 1)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceError","should have computed distances already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeVertexDistances(mris);
    }
  }
  if (!(mris->dist_alloced_flags & 2)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceTerm","should have computed dist_origs already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeOriginalVertexDistances(mris);
    }
  }

  if (false) {
    fprintf(stdout, "%s:%d calling mrisCheckDistOrig\n", __FILE__, __LINE__);
    if (!mrisCheckDistOrig(mris)) 
      fprintf(stdout, "  failed mrisCheckDistOrig\n");
  }
  
  int max_v, max_n, err_cnt, max_errs;
  volatile int count_dist_orig_zeros = 0;
  double dist_scale, sse_dist, max_del;

#if METRIC_SCALE
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else if (mris->status == MRIS_PARAMETERIZED_SPHERE) {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
  else
    dist_scale = mris->neg_area < mris->total_area ? sqrt(mris->orig_area / (mris->total_area - mris->neg_area))
                                                   : sqrt(mris->orig_area / mris->total_area);
#else
  dist_scale = 1.0;
#endif
  max_del = -1.0;
  max_v = max_n = -1;

  err_cnt = 0;
  max_errs = 100;

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
  int trial;
  double sse_dist_trial0;
  for (trial = 0; trial < 2; trial++) {
#endif

  sse_dist = 0.0;
  
  int const acceptableNumberOfZeros = 
    ( mris->status == MRIS_PARAMETERIZED_SPHERE
    ||mris->status == MRIS_SPHERE)
    ? mris->nvertices * mris->avg_nbrs * 0.01       // only the direction from 000 has to match
    : mris->nvertices * mris->avg_nbrs * 0.001;     // xyz has to match

  const char* vertexRipflags = MRISexportVertexRipflags(mris);  
    // since we have to read them a lot, get them into 100KB in the L2 cache
    // rather than reading them in 6.4MB of cache lines
  
#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_REPRODUCIBLE

  #define ROMP_VARIABLE       vno 
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse_dist
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse_dist ROMP_PARTIALSUM(0)
    
#else
  int vno;
  
  ROMP_PF_begin         // mris_register

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
  #pragma omp parallel for if(trial==0) reduction(+ : sse_dist)
#else
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse_dist)
#endif
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif    

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) ROMP_PF_continue;

    if (vno == Gdiag_no) DiagBreak();

#if NO_NEG_DISTANCE_TERM
    if (v->neg) ROMP_PF_continue;
#endif

    double v_sse = 0.0;

    int n;
    for (n = 0; n < vt->vtotal; n++) {
      int const vn_vno = vt->v[n];
      if (vertexRipflags[vn_vno]) continue;
      
#if NO_NEG_DISTANCE_TERM
      if (mris->vertices[vn_vno].neg) continue;

#endif
      float const dist_orig_n = !v->dist_orig ? 0.0 : v->dist_orig[n];
      
      if (dist_orig_n >= UNFOUND_DIST) continue;

      if (DZERO(dist_orig_n) && (count_dist_orig_zeros++ > acceptableNumberOfZeros)) {
        fprintf(stderr, "v[%d]->dist_orig[%d] = %f!!!!, count_dist_orig_zeros:%d\n", vno, n, dist_orig_n, count_dist_orig_zeros);
        fflush(stderr);
        DiagBreak();
        if (++err_cnt > max_errs) {
          fprintf(stderr, ">>> dump of head zeroes \n");
          int dump_vno,dump_n,dump_count_dist_orig_zeros = 0;
          for (dump_vno = 0; dump_vno < mris->nvertices; dump_vno++) {
            VERTEX_TOPOLOGY const * const dump_vt = &mris->vertices_topology[dump_vno];
            VERTEX          const * const dump_v  = &mris->vertices         [dump_vno];
            if (dump_v->ripflag) continue;
            for (dump_n = 0; dump_n < dump_vt->vtotal; dump_n++) {
              int const dump_vn_vno = dump_vt->v[n];
              if (vertexRipflags[dump_vn_vno]) continue;
              float const dump_dist_orig_n = !dump_v->dist_orig ? 0.0 : dump_v->dist_orig[n];
              if (!DZERO(dump_dist_orig_n)) continue;
              dump_count_dist_orig_zeros++;
              fprintf(stderr, "v[%d]->dist_orig[%d] = %f!!!!, count_dist_orig_zeros:%d\n", dump_vno, dump_n, dump_dist_orig_n, dump_count_dist_orig_zeros);
              if (dump_count_dist_orig_zeros > 30) goto dump_done;
            }
          }
          dump_done:;
          ErrorExit(ERROR_BADLOOP, "mrisComputeDistanceError: Too many errors!\n");
        }
      }

      double delta = dist_scale * v->dist[n] - dist_orig_n;
      if (parms->vsmoothness)
        v_sse += (1.0 - parms->vsmoothness[vno]) * (delta * delta);
      else
        v_sse += delta * delta;

      if (!isfinite(delta) || !isfinite(v_sse)) DiagBreak();
    }
    if (v_sse > 10000) DiagBreak();

    if (parms->dist_error) parms->dist_error[vno] = v_sse;

    sse_dist += v_sse;
    if (!isfinite(sse_dist) || !isfinite(v_sse)) DiagBreak();
    
#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_REPRODUCIBLE
    #undef sse_dist 
  #include "romp_for_end.h"
#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
    if (trial == 0) {
        sse_dist_trial0 = sse_dist;
    } else { 
        if (sse_dist_trial0 != sse_dist) {
            fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
               sse_dist_trial0, sse_dist, sse_dist_trial0-sse_dist);
        }
    }
  } // trial
#endif

  freeAndNULL(vertexRipflags);  

  /*fprintf(stdout, "max_del = %f at v %d, n %d\n", max_del, max_v, max_n) ;*/
  return (sse_dist);
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
    if (!isfinite(scale)) {
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
    if (!isfinite(scale)) {
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

    if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
    }
    if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz)) {
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
      sphericalProjection(x, y, z, &v->x, &v->y, &v->z);
      orig_area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
      if (orig_area <= 0) {
        continue;
      }
      while (step < last_step) {
        eps = epsilon[step];
        sphericalProjection(x + eps * dx, y + eps * dy, z + eps * dz, &v->x, &v->y, &v->z);
        area = mrisComputeArea(mris, vt->f[n], (int)vt->n[n]);
        if (area > 0) {
          break; /* we can stop here */
        }
        step++;
      }
    }

    eps = epsilon[step];
    sphericalProjection(x + eps * dx, y + eps * dy, z + eps * dz, &v->x, &v->y, &v->z);
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

int MRISapplyGradient(MRI_SURFACE *mris, double dt)
{
  int vno, nvertices;

  nvertices = mris->nvertices;
  MRISstoreCurrentPositions(mris);
  if (mris->status == MRIS_RIGID_BODY) {
    MRISrotate(mris, mris, dt * mris->alpha, dt * mris->beta, dt * mris->gamma);
  }
  else {
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(assume_reproducible) schedule(static, 1)
#endif
    for (vno = 0; vno < nvertices; vno++) {
      ROMP_PFLB_begin
      
      VERTEX *v;
      v = &mris->vertices[vno];
      if (v->ripflag) ROMP_PF_continue;

      if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z))
        ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
      if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz))
        ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
      v->x += dt * v->dx;
      v->y += dt * v->dy;
      v->z += dt * v->dz;
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }
  return (NO_ERROR);
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

    if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
      ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
    }

    x = v->x;
    y = v->y;
    z = v->z;

    if (which_gradient) {
      if (!isfinite(v->odx) || !isfinite(v->ody) || !isfinite(v->odz)) {
        ErrorPrintf(ERROR_BADPARM, "vertex %d position is not finite!\n", vno);
      }

      dx = v->odx;
      dy = v->ody;
      dz = v->odz;
    }
    else {
      if (!isfinite(v->dx) || !isfinite(v->dy) || !isfinite(v->dz)) {
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
          sphericalProjection(x, y, z, &v->x, &v->y, &v->z);
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
            sphericalProjection(
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
        sphericalProjection(
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
          sphericalProjection(x, y, z, &v->x, &v->y, &v->z);
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
            sphericalProjection(
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
          sphericalProjection(x, y, z, &v->x, &v->y, &v->z);
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
          sphericalProjection(
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
        sphericalProjection(
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
          sphericalProjection(x, y, z, &v->x, &v->y, &v->z);
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
          sphericalProjection(
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
        sphericalProjection(
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
        VERTEX * v;
        MHTfindClosestVertexGeneric(mht, mris, xs, ys, zs, 10, 4, &v, &vno, &dist);
        if (v == NULL) continue;
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
    sprintf(fname,
            "%s/%s.%d.dist.mgz",
            path,
            mris->hemisphere == LEFT_HEMISPHERE ? "lh" : mris->hemisphere == BOTH_HEMISPHERES ? "both" : "rh",
            parms->t);
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

  MHT *mht_v_current = NULL, *mht_f_current;
  if (!FZERO(parms->l_repulse)) {
    mht_v_current = MHTcreateVertexTable(mris, CURRENT_VERTICES);
    mht_f_current = MHTcreateFaceTable(mris);
  }

  mrisComputeTargetLocationTerm(mris, parms->l_location, parms);
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  ------------------------------------------------------*/
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

    sprintf(fname, "%s.%s.out", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
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
  MRIS *mris_corrected;
  int *face_trans, *vertex_trans;
  FACE *f, *fdst;
  int vno, i, n, fno;

  face_trans = (int *)calloc(mris->nfaces, sizeof(int));
  vertex_trans = (int *)calloc(mris->nvertices, sizeof(int));
  memset(vertex_trans, -1, mris->nvertices * sizeof(int));
  memset(face_trans, -1, mris->nfaces * sizeof(int));
  // create a new surface
  mris_corrected = MRISalloc(mris->nvertices, mris->nfaces);
  // keep the extra info into the new one
  mris_corrected->useRealRAS = mris->useRealRAS;
  copyVolGeom(&mris->vg, &mris_corrected->vg);

  mris_corrected->type = MRIS_TRIANGULAR_SURFACE;

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
    for (n = vdstt->vnum = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].marked == 0) {
        vdstt->vnum++;
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
    vdst = MHTfindClosestVertex(mht, mris_total, v);
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
    v = MHTfindClosestVertex(mht_src, mris_src, vdst);
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
    VERTEX const * const v = MHTfindClosestVertex(mht_src, mris_src, vdst);
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
  MRIS *mris_corrected;

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
  mris_corrected = MRISalloc(kept_vertices, kept_faces);
  // keep the extra info into the new one
  mris_corrected->useRealRAS = mris->useRealRAS;
  copyVolGeom(&mris->vg, &mris_corrected->vg);

  mris_corrected->type = MRIS_TRIANGULAR_SURFACE;

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
    for (n = vdstt->vnum = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].ripflag == 0) {
        vdstt->vnum++;
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
int MRISremoveIntersections(MRI_SURFACE *mris)
{
  int n, num, writeit = 0, old_num, nbrs, min_int, no_progress = 0;

  n = 0;

  num = mrisMarkIntersections(mris);
  if (num == 0) return (NO_ERROR);
  printf("removing intersecting faces\n");
  min_int = old_num = mris->nvertices;
  MRISsaveVertexPositions(mris, TMP2_VERTICES);
  nbrs = 0;
  MRISclearMarks(mris);
  num = mrisMarkIntersections(mris);
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
    num = mrisMarkIntersections(mris);
  }

  if (num > min_int) {
    MRISrestoreVertexPositions(mris, TMP2_VERTICES);  // the least number of negative vertices
    num = mrisMarkIntersections(mris);
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

  memset(&parms, 0, sizeof(parms));
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

  printf(
      "MRISrepositionSurfaceToCoordinate(%d, %f, %f, %f, %d, %f, %x)\n", target_vno, tx, ty, tz, nsize, sigma, flags);
  memset(&parms, 0, sizeof(parms));
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

  memmove(&mris->vg, &mris1->vg, sizeof(mris1->vg));
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
  when the orig.nofix is input. Not sure why, but probably because
  the lengths of the edges are all either 1 or sqrt(2) thus
  creating some abiguity which is handled differently on different
  runs.
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
    vtxtnew->vnum = vtxtold->vnum; // number of neighboring vertices
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
  char *tmpName = strdup("/tmp/mris_decimateXXXXXX");
  int fd = mkstemp(tmpName);
  if(fd == -1) {
    printf("Error creating temporary file: %s\n",tmpName);
    return(NULL);
  }
  char tmp_fpath[STRLEN];
  FileNameAbsolute(tmpName, tmp_fpath);
  MRISwrite(mris, tmp_fpath);
  MRISfree(&mris);
  mris = MRISread(tmp_fpath);
  remove(tmp_fpath);

  
  return(mris);
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
    if (!isfinite(d)) {
      ErrorPrintf(ERROR_BADPARM, "point (%2.2f,%2.2f,%2.2f) cannot be projected on cylinder", x, y, z);
    }
    dx = d * x;
    dz = d * z;
    v->x = x + dx;
    v->z = z + dz;

    if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
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
  Perform a projection onto an sphere moving each
  point on the cortical surface to the closest spherical
  coordinate.
  ------------------------------------------------------*/
MRI_SURFACE *MRISprojectOntoSphere(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, double r)
{
  VERTEX *v;
  int vno;
  double x, y, z, d, dx, dy, dz, dist, total_dist, x2, y2, z2;

  if (FZERO(r)) {
    r = DEFAULT_RADIUS;
  }

  if (!mris_dst) {
    mris_dst = MRISclone(mris_src);
  }

  if ((mris_dst->status != MRIS_SPHERE) && (mris_dst->status != MRIS_PARAMETERIZED_SPHERE)) {
    MRIScenter(mris_dst, mris_dst);
  }

  mris_dst->radius = r;

  MRISfreeDistsButNotOrig(mris_dst);
  
  for (total_dist = vno = 0; vno < mris_dst->nvertices; vno++) {
    v = &mris_dst->vertices[vno];
    if (v->ripflag) /* shouldn't happen */
    {
      continue;
    }
    if (false && vno == 118009) {
      DiagBreak();
    }
    x = (double)v->x;
    y = (double)v->y;
    z = (double)v->z;

    x2 = x * x;
    y2 = y * y;
    z2 = z * z;
    dist = sqrt(x2 + y2 + z2);
    if (FZERO(dist)) {
      d = 0;
    }
    else {
      d = 1 - r / dist;
    }
    dx = d * x;
    dy = d * y;
    dz = d * z;
    v->x = x - dx;
    v->y = y - dy;
    v->z = z - dz;

    if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
      DiagBreak();
    }

    /*    if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON)*/
    {
      dist = sqrt((double)(dx * dx + dy * dy + dz * dz));
      total_dist += dist;
    }
#if 1
    x = (double)v->x;
    y = (double)v->y;
    z = (double)v->z;
    x2 = x * x;
    y2 = y * y;
    z2 = z * z;
    dist = sqrt(x2 + y2 + z2);
#endif
  }
  if ((Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
    fprintf(stdout, "sphere_project: total dist = %f\n", total_dist);
  }
  MRISupdateEllipsoidSurface(mris_dst);
  mris_dst->status = mris_src->status == MRIS_PARAMETERIZED_SPHERE ? MRIS_PARAMETERIZED_SPHERE : MRIS_SPHERE;
  return (mris_dst);
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
    if (!isfinite(d)) {
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

    if (!isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)) {
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


int MRISrigidBodyAlignLocal(MRI_SURFACE *mris, INTEGRATION_PARMS *old_parms)
{
  int old_status, steps;
  INTEGRATION_PARMS parms;

  /* dx,dy,dz interpreted as rotations in applyGradient when status is rigid */
  old_status = mris->status; /* okay, okay, this is a hack too... */
  mris->status = MRIS_RIGID_BODY;
  memset(&parms, 0, sizeof(parms));
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
  int n, old_status, steps;
  INTEGRATION_PARMS parms;

  /* dx,dy,dz interpreted as rotations in applyGradient when status is rigid */
  old_status = mris->status; /* okay, okay, this is a hack too... */
  mris->status = MRIS_RIGID_BODY;
  memset(&parms, 0, sizeof(parms));
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

  int const old_status = mris->status;
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

  double const ext_sse = gMRISexternalSSE ? (*gMRISexternalSSE)(mris, parms) : 0.0;

  if (use_new) {
  
    printf("Starting new MRISrigidBodyAlignGlobal_findMinSSE()\n");
    Timer new_timer;
    
    // This does not modify either mris or params until after the old code has executed
    //
    double ext_sse = (gMRISexternalSSE) ? (*gMRISexternalSSE)(mris, parms) : 0.0;
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
          center_sse_known ? (float)(center_sse + ext_sse) : -666.0);
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
            
            double sse = mrisComputeCorrelationErrorTraceable(mris, parms, 1, trace);
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
                (float)DEGREES(mina),  (float)DEGREES(minb), (float)DEGREES(ming),   (float)(min_sse + ext_sse),
                (float)DEGREES(alpha), (float)DEGREES(beta), (float)DEGREES(gamma),  (float)(    sse + ext_sse));
              fflush(stdout);
            }
          }  // gamma
          if (false && tracing) {
            fprintf(stdout, "\r  beta "
              "min @ (%2.2f, %2.2f, %2.2f) sse:%2.1f   try @ (%+2.2f, %+2.2f, *)",
              (float)DEGREES(mina),  (float)DEGREES(minb), (float)DEGREES(ming),   (float)(min_sse + ext_sse),
              (float)DEGREES(alpha), (float)DEGREES(beta));
          }
        }    // beta
        if (false && tracing) {
          fprintf(stdout, "\n");
        }
      }      // alpha
    
      int msec = old_timer.milliseconds();
      printf("  min @ (%2.2f, %2.2f, %2.2f) sse = %2.1f, elapsed since starting=%6.4f min\n",
        (float)DEGREES(mina), (float)DEGREES(minb), (float)DEGREES(ming),
        (float)(min_sse + ext_sse),
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
      (float)(new_sse + ext_sse));
  }
  if (tracingWithSnapshot) {
    fprintf(parms->fp,
      "rotating brain by (%2.2f, %2.2f, %2.2f), final sse: %2.2f\n",
      (float)DEGREES(new_mina),
      (float)DEGREES(new_minb),
      (float)DEGREES(new_ming),
      (float)(new_sse + ext_sse));
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
  int old_status = mris->status;

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

  // ath - temporarily removing (#626)
  // cheapAssert(mrisCanAttachFaceToVertices(mris, vno0, vno1, vno2));
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

