/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2019 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_compute_dxyz.h"

#include "mrisurf_project.h"
#include "mrisurf_sseTerms.h"

#include "mrisurf_base.h"

double LOCATION_MOVE_LEN = 0.25;
void mrisDxyzSetLocationMoveLen(double newval){
  LOCATION_MOVE_LEN = newval;
}

#define MAX_VOXELS          mrisurf_sse_MAX_VOXELS
#define MAX_DISPLACEMENT    mrisurf_sse_MAX_DISPLACEMENT 
#define DISPLACEMENT_DELTA  mrisurf_sse_DISPLACEMENT_DELTA
#define DEFAULT_STD         mrisurf_sse_DEFAULT_STD

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

        if (!std::isfinite(dx) || !std::isfinite(dy) || !std::isfinite(dz)) 
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
    vx = v->x - v->t2x;
    vy = v->y - v->t2y;
    vz = v->z - v->t2z;
    dx = dy = dz = 0.0f;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        vnx = vn->x - vn->t2x;
        vny = vn->y - vn->t2y;
        vnz = vn->z - vn->t2z;
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
int mrisComputeNormalizedSpringTerm(MRIS* mris, double const l_spring)
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


/**
  Computes the normal spring term for a given vertex.
*/
inline void vertexComputeNormalSpringTerm(MRIS* mris, int vno, float* dx, float* dy, float* dz, double l_spring)
{
  VERTEX * const vertex = &mris->vertices[vno];
  VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];

  if (vertex->ripflag) return;

  float x = vertex->x;
  float y = vertex->y;
  float z = vertex->z;

  float sx = 0, sy = 0, sz = 0;

  int n = 0;
  for (int m = 0; m < vertext->vnum; m++) {
    VERTEX const * const vn = &mris->vertices[vertext->v[m]];
    if (!vn->ripflag) {
      sx += vn->x - x;
      sy += vn->y - y;
      sz += vn->z - z;
      n++;
    }
  }

  if (n > 0) {
    sx /= n;
    sy /= n;
    sz /= n;
  }

  float nx = vertex->nx;
  float ny = vertex->ny;
  float nz = vertex->nz;

  // project onto normal
  float nc = sx * nx + sy * ny + sz * nz;

  // move in normal direction
  *dx = l_spring * nc * nx;
  *dy = l_spring * nc * ny;
  *dz = l_spring * nc * nz;
}


void mrisComputeNormalSpringTerm(MRIS *mris, double l_spring)
{
  if (FZERO(l_spring)) return;

  for (int vno = 0; vno < mris->nvertices; vno++) {
    float dx = 0, dy = 0, dz = 0;
    vertexComputeNormalSpringTerm(mris, vno, &dx, &dy, &dz, l_spring);

    VERTEX * const vertex = &mris->vertices[vno];
    vertex->dx += dx;
    vertex->dy += dy;
    vertex->dz += dz;
  }
}


/**
  Computes the tangential spring term for a given vertex.
*/
inline void vertexComputeTangentialSpringTerm(MRIS* mris, int vno, float* dx, float* dy, float* dz, double l_spring)
{
  VERTEX * const vertex = &mris->vertices[vno];
  VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];

  if (vertex->ripflag) return;
  if (vertex->border && !vertex->neg) return;

  float x = vertex->x;
  float y = vertex->y;
  float z = vertex->z;

  float sx = 0, sy = 0, sz = 0;

  int n = 0;
  for (int m = 0; m < vertext->vnum; m++) {
    VERTEX const * const vn = &mris->vertices[vertext->v[m]];
    if (!vn->ripflag) {
      sx += vn->x - x;
      sy += vn->y - y;
      sz += vn->z - z;
      n++;
    }
  }

  if (n > 0) {
    sx /= n;
    sy /= n;
    sz /= n;
  }

  float nx = vertex->nx;
  float ny = vertex->ny;
  float nz = vertex->nz;

  // project onto normal
  float nc = sx * nx + sy * ny + sz * nz;

  // remove normal component and scale
  *dx = l_spring * (sx - nc * nx);
  *dy = l_spring * (sy - nc * ny);
  *dz = l_spring * (sz - nc * nz);
}


void mrisComputeTangentialSpringTerm(MRIS *mris, double l_spring)
{
  if (FZERO(l_spring)) return;

  for (int vno = 0; vno < mris->nvertices; vno++) {
    float dx = 0, dy = 0, dz = 0;
    vertexComputeTangentialSpringTerm(mris, vno, &dx, &dy, &dz, l_spring);

    VERTEX * const vertex = &mris->vertices[vno];
    vertex->dx += dx;
    vertex->dy += dy;
    vertex->dz += dz;
  }
}


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

  step = mht->vres() / 2;
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
    target = MRISPfunctionVal(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no);
    std    = MRISPfunctionVal(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no + 1);
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
    up1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x + ux, y + uy, z + uz, fno);
    um1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x - ux, y - uy, z - uz, fno);
    vp1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x + vx, y + vy, z + vz, fno);
    vm1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x - vx, y - vy, z - vz, fno);
    du = (up1 - um1) / (2 * d_dist);
    dv = (vp1 - vm1) / (2 * d_dist);
    v->dx -= coef * (du * e1x + dv * e2x);
    v->dy -= coef * (du * e1y + dv * e2y);
    v->dz -= coef * (du * e1z + dv * e2z);

    mag = sqrt(v->dx * v->dx + v->dy * v->dy + v->dz * v->dz);
    if (mag > max_mag) {
      max_mag = mag;
    }
    if (!std::isfinite(v->dx) || !std::isfinite(v->dy) || !std::isfinite(v->dz)) {
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
    if (!std::isfinite(v->dx) || !std::isfinite(v->dy) || !std::isfinite(v->dz)) {
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
    target = MRISPfunctionVal(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no);
    std    = MRISPfunctionVal(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no + 1);
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
    am1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x - dx, y - dy, z - dz, 0);
    ap1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x + dx, y + dy, z + dz, 0);
    da = (ap1 - am1) / (2 * D_ANGLE);

    /* compute beta term - differential rotation around y axis */
    dx = -z * D_ANGLE;
    dy = 0;
    dz = x * D_ANGLE;
    bm1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x - dx, y - dy, z - dz, 0);
    bp1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x + dx, y + dy, z + dz, 0);
    db = (bp1 - bm1) / (2 * D_ANGLE);

    /* compute gamma term - differential rotation around x axis */
    dx = 0;
    dy = -z * D_ANGLE;
    dz = y * D_ANGLE;
    gm1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x - dx, y - dy, z - dz, 0);
    gp1 = MRISPfunctionVal(parms->mrisp_template, mris->radius, x + dx, y + dy, z + dz, 0);
    dg = (gp1 - gm1) / (2 * D_ANGLE);

    mris->gamma -= coef * dg; /* around x-axis */
    mris->beta -= coef * db;  /* around y-axis */
    mris->alpha -= coef * da; /* around z-axis */

    mag = sqrt(v->dx * v->dx + v->dy * v->dy + v->dz * v->dz);
    if (mag > max_mag) {
      max_mag = mag;
    }
    if (!std::isfinite(v->dx) || !std::isfinite(v->dy) || !std::isfinite(v->dz)) {
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
    double d_dist = D_DIST * mris->avg_vertex_dist;

    VERTEX * const v = &mris->vertices[vno];
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
    e1p = mrisSampleMinimizationEnergy(mris, vno, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleMinimizationEnergy(mris, vno, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleMinimizationEnergy(mris, vno, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleMinimizationEnergy(mris, vno, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
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
    E0 = mrisSampleMinimizationEnergy(mris, vno, parms, v->x, v->y, v->z);
    E1 = mrisSampleMinimizationEnergy(
        mris, vno, parms, v->x + parms->dt * dx, v->y + parms->dt * dy, v->z + parms->dt * dz);

    if (E1 > E0) {
      double E2;

      DiagBreak();
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      E2 = mrisSampleMinimizationEnergy(
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
    e1p = mrisSampleNormalEnergy(mris, vno, parms, v->x + d_dist * e1x, v->y + d_dist * e1y, v->z + d_dist * e1z);
    e1m = mrisSampleNormalEnergy(mris, vno, parms, v->x - d_dist * e1x, v->y - d_dist * e1y, v->z - d_dist * e1z);
    e2p = mrisSampleNormalEnergy(mris, vno, parms, v->x + d_dist * e2x, v->y + d_dist * e2y, v->z + d_dist * e2z);
    e2m = mrisSampleNormalEnergy(mris, vno, parms, v->x - d_dist * e2x, v->y - d_dist * e2y, v->z - d_dist * e2z);
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
    E0 = mrisSampleNormalEnergy(mris, vno, parms, v->x, v->y, v->z);
    E1 = mrisSampleNormalEnergy(mris, vno, parms, cx, cy, cz);
    if (E1 > E0) {
      double E2;
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      DiagBreak();
      E2 = mrisSampleNormalEnergy(mris, vno, parms, v->x - parms->dt * dx, v->y - parms->dt * dy, v->z - parms->dt * dz);
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
        MHTfindClosestVertexGeneric(mht, xs, ys, zs, 10, 4, &vno, &dist);
        if (vno < 0) continue;
        if (FZERO(val) && dist > 1) continue;
        if (vno == Gdiag_no) DiagBreak();
        VERTEX * v = &mris->vertices[vno];
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












