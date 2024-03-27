/**
 * @brief more surface processing utils
 *
 */
/*
 * Original Authors: Segonne and Greve
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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include "mrisutils.h"
#include "surfgrad.h"

#include "romp_support.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_sseTerms.h"
#include "mrisurf_vals.h"
#include "mrisurf_project.h"

#include "cma.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "fsenv.h"
#include "gca.h"
#include "gcamorph.h"
#include "icosahedron.h"
#include "macros.h"
#include "matrix.h"
#include "mri_identify.h"
#include "mrisegment.h"
#include "mrishash.h"
#include "proto.h"
#include "sig.h"
#include "stats.h"
#include "timer.h"
#include "tritri.h"
#include "resample.h"
#include "mri2.h"

#include "annotation.h"

#include "version.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
////                    USEFUL ROUTINES               ////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////


int MRIScanonicalToWorld(MRI_SURFACE *mris, double phi, double theta, double *pxw, double *pyw, double *pzw)
{
  double x, y, z, radius;

  radius = mris->radius;
  *pxw = x = radius * sin(phi) * cos(theta);
  *pyw = y = radius * sin(phi) * sin(theta);
  *pzw = z = radius * cos(phi);
  return (NO_ERROR);
}


// peel a volume from a surface (World coordinates)
// *** type == 0:  keep the inside
//     type == 1:  peel the outside surface and set the inside value to 'val'
//     type == 2: keep the outside
//     type == 3: peel the inside surface and set the outside value to 'val'
// *** NbVoxels corresponds to:
//               - the number of kept voxels if (type>=0)
//               - the number of removed voxels if (type<0)
MRI *MRISpeelVolume(MRIS *mris, MRI *mri_src, MRI *mri_dst, int type, unsigned char val, unsigned long *NbVoxels)
{
  int i, j, k, imnr;
  float x0, y0, z0, x1, y1, z1, x2, y2, z2, d0, d1, d2, dmax, u, v;
  float px, py, pz, px0, py0, pz0, px1, py1, pz1;
  int numu, numv, totalfilled, newfilled;
  double tx, ty, tz;
  unsigned char tmpval;
  unsigned long size;
  int width, height, depth;
  MRI *mri_buff;

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  if (mri_dst == NULL) mri_dst = MRIalloc(width, height, depth, mri_src->type);

  mri_buff = MRIalloc(width, height, depth, MRI_UCHAR);

  for (k = 0; k < mris->nfaces; k++) {
    x0 = mris->vertices[mris->faces[k].v[0]].x;
    y0 = mris->vertices[mris->faces[k].v[0]].y;
    z0 = mris->vertices[mris->faces[k].v[0]].z;
    x1 = mris->vertices[mris->faces[k].v[1]].x;
    y1 = mris->vertices[mris->faces[k].v[1]].y;
    z1 = mris->vertices[mris->faces[k].v[1]].z;
    x2 = mris->vertices[mris->faces[k].v[2]].x;
    y2 = mris->vertices[mris->faces[k].v[2]].y;
    z2 = mris->vertices[mris->faces[k].v[2]].z;
    d0 = sqrt(SQR(x1 - x0) + SQR(y1 - y0) + SQR(z1 - z0));
    d1 = sqrt(SQR(x2 - x1) + SQR(y2 - y1) + SQR(z2 - z1));
    d2 = sqrt(SQR(x0 - x2) + SQR(y0 - y2) + SQR(z0 - z2));
    dmax = (d0 >= d1 && d0 >= d2) ? d0 : (d1 >= d0 && d1 >= d2) ? d1 : d2;
    numu = ceil(2 * d0);
    numv = ceil(2 * dmax);

    for (v = 0; v <= numv; v++) {
      px0 = x0 + (x2 - x0) * v / numv;
      py0 = y0 + (y2 - y0) * v / numv;
      pz0 = z0 + (z2 - z0) * v / numv;
      px1 = x1 + (x2 - x1) * v / numv;
      py1 = y1 + (y2 - y1) * v / numv;
      pz1 = z1 + (z2 - z1) * v / numv;
      for (u = 0; u <= numu; u++) {
        px = px0 + (px1 - px0) * u / numu;
        py = py0 + (py1 - py0) * u / numu;
        pz = pz0 + (pz1 - pz0) * u / numu;

        // MRIworldToVoxel(mri_src,px,py,pz,&tx,&ty,&tz);
        MRIsurfaceRASToVoxel(mri_src, px, py, pz, &tx, &ty, &tz);
        imnr = (int)(tz + 0.5);
        j = (int)(ty + 0.5);
        i = (int)(tx + 0.5);

        if (i >= 0 && i < width && j >= 0 && j < height && imnr >= 0 && imnr < depth)
          MRIvox(mri_buff, i, j, imnr) = 255;
      }
    }
  }

  MRIvox(mri_buff, 1, 1, 1) = 64;
  totalfilled = newfilled = 1;
  while (newfilled > 0) {
    newfilled = 0;
    for (k = 1; k < depth - 1; k++)
      for (j = 1; j < height - 1; j++)
        for (i = 1; i < width - 1; i++)
          if (MRIvox(mri_buff, i, j, k) == 0)
            if (MRIvox(mri_buff, i, j, k - 1) == 64 || MRIvox(mri_buff, i, j - 1, k) == 64 ||
                MRIvox(mri_buff, i - 1, j, k) == 64) {
              MRIvox(mri_buff, i, j, k) = 64;
              newfilled++;
            }
    for (k = depth - 2; k >= 1; k--)
      for (j = height - 2; j >= 1; j--)
        for (i = width - 2; i >= 1; i--)
          if (MRIvox(mri_buff, i, j, k) == 0)
            if (MRIvox(mri_buff, i, j, k + 1) == 64 || MRIvox(mri_buff, i, j + 1, k) == 64 ||
                MRIvox(mri_buff, i + 1, j, k) == 64) {
              MRIvox(mri_buff, i, j, k) = 64;
              newfilled++;
            }
    totalfilled += newfilled;
  }

  size = 0;
  switch (type) {
    case 0:
      for (k = 1; k < depth - 1; k++)
        for (j = 1; j < height - 1; j++)
          for (i = 1; i < width - 1; i++) {
            if (MRIvox(mri_buff, i, j, k) == 64)
              MRIvox(mri_dst, i, j, k) = 0;
            else {
              tmpval = MRIvox(mri_src, i, j, k);
              MRIvox(mri_dst, i, j, k) = tmpval;
              size++;
            }
          }
      break;
    case 1:
      for (k = 1; k < depth - 1; k++)
        for (j = 1; j < height - 1; j++)
          for (i = 1; i < width - 1; i++) {
            if (MRIvox(mri_buff, i, j, k) == 64)
              MRIvox(mri_dst, i, j, k) = 0;
            else {
              MRIvox(mri_dst, i, j, k) = val;
              size++;
            }
          }
      break;
    case 2:
      for (k = 1; k < depth - 1; k++)
        for (j = 1; j < height - 1; j++)
          for (i = 1; i < width - 1; i++) {
            if (MRIvox(mri_buff, i, j, k) == 64) {
              tmpval = MRIvox(mri_src, i, j, k);
              MRIvox(mri_dst, i, j, k) = tmpval;
            }
            else {
              MRIvox(mri_dst, i, j, k) = 0;
              size++;
            }
          }
      break;
    case 3:
      for (k = 1; k < depth - 1; k++)
        for (j = 1; j < height - 1; j++)
          for (i = 1; i < width - 1; i++) {
            if (MRIvox(mri_buff, i, j, k) == 64)
              MRIvox(mri_dst, i, j, k) = val;
            else {
              MRIvox(mri_dst, i, j, k) = 0;
              size++;
            }
          }
      break;
  }
  if (NbVoxels) (*NbVoxels) = size;

  MRIfree(&mri_buff);
  return mri_dst;
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
////   ROUTINES FOR MATCHING A SURFACE TO A VOLUME LABEL //////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

static MRI *mriIsolateLabel(MRI *mri_seg, int label, MRI_REGION *bbox);
static double mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, float delta_t, MHT *mht, float max_mag);
static int mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno);
static int mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno);
static int mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno);
static int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv);
static int mrisComputeIntensityTerm(MRI_SURFACE *mris, double l_intensity, MRI *mri, double sigma);
static int mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring);
static int mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring);



MRI_REGION *MRIlocateRegion(MRI *mri, int label)
{
  int i, j, k;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  MRI_REGION *mri_region = (MRI_REGION *)malloc(sizeof(MRI_REGION));

  zmax = ymax = xmax = 0;
  zmin = ymin = xmin = 10000;

  for (k = 0; k < mri->depth; k++)
    for (j = 0; j < mri->height; j++)
      for (i = 0; i < mri->width; i++)
        if (MRIvox(mri, i, j, k) == label) {
          if (k < zmin) zmin = k;
          if (j < ymin) ymin = j;
          if (i < xmin) xmin = i;
          if (k > zmax) zmax = k;
          if (j > ymax) ymax = j;
          if (i > xmax) xmax = i;
        }

  mri_region->x = xmin;
  mri_region->y = ymin;
  mri_region->z = zmin;
  mri_region->dx = xmax - xmin;
  mri_region->dy = ymax - ymin;
  mri_region->dz = zmax - zmin;

  return mri_region;
}

static MRI *mriIsolateLabel(MRI *mri_seg, int label, MRI_REGION *bbox)
{
  int i, j, k;
  int xplusdx, yplusdy, zplusdz;
  MRI *mri = MRIalloc(mri_seg->width, mri_seg->height, mri_seg->depth, mri_seg->type);

  xplusdx = bbox->x + bbox->dx + 1;
  yplusdy = bbox->y + bbox->dy + 1;
  zplusdz = bbox->z + bbox->dz + 1;

  for (k = bbox->z; k < zplusdz; k++)
    for (j = bbox->y; j < yplusdy; j++)
      for (i = bbox->x; i < xplusdx; i++)
        if (MRIvox(mri_seg, i, j, k) == label) MRIvox(mri, i, j, k) = 1;

  return mri;
}


#define MIN_NBR_DIST (0.25)

static int mrisRemoveNormalGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX *v;
  float dot;

  v = &mris->vertices[vno];
  if (v->ripflag) return (NO_ERROR);

  dot = v->nx * v->odx + v->ny * v->ody + v->nz * v->odz;
  v->odx -= dot * v->nx;
  v->ody -= dot * v->ny;
  v->odz -= dot * v->nz;

  return (NO_ERROR);
}

static int mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];
  if (v->ripflag) return (NO_ERROR);

  float const x = v->x;
  float const y = v->y;
  float const z = v->z;

  int n;
  for (n = 0; n < vt->vnum; n++) {

    VERTEX const * const vn = &mris->vertices[vt->v[n]];
    float dx = vn->x - x;
    float dy = vn->y - y;
    float dz = vn->z - z;
    float const dist = sqrt(dx * dx + dy * dy + dz * dz);

    /* too close - take out gradient component in this dir. */
    if (dist <= MIN_NBR_DIST) {
      dx /= dist;
      dy /= dist;
      dz /= dist;
      float const dot = dx * v->odx + dy * v->ody + dz * v->odz;
      if (dot > 0.0) {
        v->odx -= dot * dx;
        v->ody -= dot * dy;
        v->odz -= dot * dz;
      }
    }
  }

  return (NO_ERROR);
}

static int mrisLimitGradientDistance(MRI_SURFACE *mris, MHT *mht, int vno)
{
  VERTEX *v = &mris->vertices[vno];

  mrisRemoveNeighborGradientComponent(mris, vno);
  if (MHTisVectorFilled(mht, vno, v->odx, v->ody, v->odz)) {
    mrisRemoveNormalGradientComponent(mris, vno);
    if (MHTisVectorFilled(mht, vno, v->odx, v->ody, v->odz)) {
      v->odx = v->ody = v->odz = 0.0;
      return (NO_ERROR);
    }
  }

  return (NO_ERROR);
}

static double mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, float delta_t, MHT *mht, float max_mag)
{
  static int direction = 1;
  
  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  int i;
  for (i = 0; i < mris->nvertices; i++) {
    int const vno =
      (direction < 0) ? (mris->nvertices - i - 1) : (i);
      
    //VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    
    if (v->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();
    
    v->odx = delta_t * v->dx + momentum * v->odx;
    v->ody = delta_t * v->dy + momentum * v->ody;
    v->odz = delta_t * v->dz + momentum * v->odz;
    
    double mag = sqrt(v->odx * v->odx + v->ody * v->ody + v->odz * v->odz);
    if (mag > max_mag) /* don't let step get too big */
    {
      mag = max_mag / mag;
      v->odx *= mag;
      v->ody *= mag;
      v->odz *= mag;
    }

    /* erase the faces this vertex is part of */

    if (mht) {
      MHTremoveAllFaces(mht, mris, vno);
      mrisLimitGradientDistance(mris, mht, vno);
    }
    
    MRISsetXYZ(mris, vno, 
      v->x + v->odx,
      v->y + v->ody,
      v->z + v->odz);
    
    if (mht) MHTaddAllFaces(mht, mris, vno);
  }

  direction *= -1;
  return (delta_t);
}

#define VERTEX_EDGE(vec, v0, v1) VECTOR_LOAD(vec, v1->x - v0->x, v1->y - v0->y, v1->z - v0->z)

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Fit a 1-d quadratic to the surface locally and move the
          vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static int mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv)
{
  MATRIX *m_R, *m_R_inv;
  VECTOR *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr;
  int vno, n;
  float ui, vi, rsq, a, b;

  if (FZERO(l_curv)) return (NO_ERROR);

  mrisComputeTangentPlanes(mris);
  v_n = VectorAlloc(3, MATRIX_REAL);
  v_A = VectorAlloc(2, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_nbr = VectorAlloc(3, MATRIX_REAL);
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;
    v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL);    /* heights above TpS */
    m_R = MatrixAlloc(vt->vtotal, 2, MATRIX_REAL); /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);
    for (n = 0; n < vt->vtotal; n++) /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      ui = V3_DOT(v_e1, v_nbr);
      vi = V3_DOT(v_e2, v_nbr);
      rsq = ui * ui + vi * vi;
      *MATRIX_RELT(m_R, n + 1, 1) = rsq;
      *MATRIX_RELT(m_R, n + 1, 2) = 1;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL);
    if (!m_R_inv) {
      MatrixFree(&m_R);
      VectorFree(&v_Y);
      continue;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A);
    a = VECTOR_ELT(v_A, 1);
    b = VECTOR_ELT(v_A, 2);
    b *= l_curv;
    v->dx += b * v->nx;
    v->dy += b * v->ny;
    v->dz += b * v->nz;

    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.1f\n",
              vno,
              b * v->nx,
              b * v->ny,
              b * v->nz,
              a,
              b);
    MatrixFree(&m_R);
    VectorFree(&v_Y);
    MatrixFree(&m_R_inv);
  }

  VectorFree(&v_n);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  VectorFree(&v_nbr);
  VectorFree(&v_A);
  return (NO_ERROR);
}

static int mrisComputeIntensityTerm(MRI_SURFACE *mris, double l_intensity, MRI *mri, double sigma)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz, dx, dy, dz;
  double val0, xw, yw, zw, del, val_outside, val_inside, delI, delV;
  // int k,ktotal ;

  if (FZERO(l_intensity)) return (NO_ERROR);

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) continue;
    if (vno == Gdiag_no) DiagBreak();

    x = v->x;
    y = v->y;
    z = v->z;

    // MRIworldToVoxel(mri, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri, xw, yw, zw, &val0);
    //    sigma = v->val2 ;

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

#if 1
    {
      double val;
      xw = x + nx;
      yw = y + ny;
      zw = z + nz;
      // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw);
      MRIsampleVolume(mri, xw, yw, zw, &val);
      val_outside = val;

      xw = x - nx;
      yw = y - ny;
      zw = z - nz;
      // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
      MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw);
      MRIsampleVolume(mri, xw, yw, zw, &val);
      val_inside = val;
    }
#else
    /* compute intensity gradient using smoothed volume */
    {
      double dist, val, step_size;
      int n;

      step_size = MIN(sigma / 2, 0.5);
      ktotal = 0.0;
      for (n = 0, val_outside = val_inside = 0.0, dist = step_size; dist <= 2 * sigma; dist += step_size, n++) {
        k = exp(-dist * dist / (2 * sigma * sigma));
        ktotal += k;
        xw = x + dist * nx;
        yw = y + dist * ny;
        zw = z + dist * nz;
        // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri, xw, yw, zw, &val);
        val_outside += k * val;

        xw = x - dist * nx;
        yw = y - dist * ny;
        zw = z - dist * nz;
        // MRIworldToVoxel(mri, xw, yw, zw, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xw, &yw, &zw);
        MRIsampleVolume(mri, xw, yw, zw, &val);
        val_inside += k * val;
      }
      val_inside /= (double)ktotal;
      val_outside /= (double)ktotal;
    }
#endif

    delV = v->val - val0;
    delI = (val_outside - val_inside) / 2.0;

    if (!FZERO(delI))
      delI /= fabs(delI);
    else if (delV < 0) /*we are inside*/
      delI = -1;
    else
      delI = 1; /* intensities tend to increase inwards */

    del = l_intensity * delV * delI;

    dx = nx * del;
    dy = ny * del;
    dz = nz * del;

    v->dx += dx;
    v->dy += dy;
    v->dz += dz;
  }

  return (NO_ERROR);
}

static int mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, nx, ny, nz, nc, x, y, z;

  if (FZERO(l_spring)) return (NO_ERROR);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vertex_topology = &mris->vertices_topology[vno];
    VERTEX                * const vertex          = &mris->vertices         [vno];
    if (vertex->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();

    nx = vertex->nx;
    ny = vertex->ny;
    nz = vertex->nz;
    x = vertex->x;
    y = vertex->y;
    z = vertex->z;

    sx = sy = sz = 0.0;
    n = 0;
    for (m = 0; m < vertex_topology->vnum; m++) {
      VERTEX const * const vn = &mris->vertices[vertex_topology->v[m]];
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
    // if (vno == Gdiag_no)
    // fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
    //        vno, sx, sy, sz) ;
  }

  return (NO_ERROR);
}

static int mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring)
{
  int vno, n, m;
  float sx, sy, sz, x, y, z, nc;

  if (FZERO(l_spring)) return (NO_ERROR);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();

    if (v->border && !v->neg) continue;

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
    sx -= l_spring * nc * v->nx;               /* remove  normal component */
    sy -= l_spring * nc * v->ny;
    sz -= l_spring * nc * v->nz;

    v->dx += sx;
    v->dy += sy;
    v->dz += sz;
    if (vno == Gdiag_no) fprintf(stdout, "v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n", vno, sx, sy, sz);
  }

  return (NO_ERROR);
}

#define MOMENTUM 0.8
#define MAX_ASYNCH_MM 0.3
#define MAX_REDUCTIONS 0
#define REDUCTION_PCT 0.5
#define AVERAGES_NBR 1
#define BASE_DT_SCALE 1.0

MRIS *MRISmatchSurfaceToLabel(
    MRIS *mris, MRI *mri_seg, int label, MRI_REGION *mri_region, INTEGRATION_PARMS *integration_parms, int connectivity)
{
  int bbox_indicator = 0, done, niterations, n, nreductions = 0;
  MRI_REGION *bbox;
  MRI *mri;
  INTEGRATION_PARMS *parms;
  int parms_indicator = 0;
  double sse, last_sse, rms, last_rms, base_dt, dt, delta_t = 0.0;
  double tol;
  MHT *mht = NULL;
  int avgs = AVERAGES_NBR;
  Timer then;
  int msec;

  then.reset();

  if (integration_parms == NULL) {
    parms = (INTEGRATION_PARMS *)calloc(1, sizeof(INTEGRATION_PARMS));
    parms_indicator = 1;
    parms->projection = NO_PROJECTION;
    parms->niterations = 5;
    parms->dt = 0.5f;
    parms->base_dt = BASE_DT_SCALE * parms->dt;
    parms->tol = 1e-3;
    parms->l_spring = 0.0f;
    parms->l_curv = 1.0;
    parms->l_intensity = 1.0f;
    parms->l_tspring = 1.0f;
    parms->l_nspring = 0.2f;
    parms->momentum = MOMENTUM;
    parms->dt_increase = 1.0 /* DT_INCREASE */;
    parms->dt_decrease = 0.50 /* DT_DECREASE*/;
    parms->error_ratio = 50.0 /*ERROR_RATIO */;
    /*  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;*/
    parms->l_surf_repulse = 1.0;
    parms->l_repulse = 0.0;  // ;
    parms->sigma = 2.0f;
    parms->mri_brain = mri_seg; /*necessary for using MRIcomputeSSE*/
  }
  else
    parms = integration_parms;

  niterations = parms->niterations;
  base_dt = parms->base_dt;
  tol = parms->tol;

  if (mri_region)
    bbox = mri_region;
  else {
    bbox = MRIlocateRegion(mri_seg, label);
    bbox_indicator = 1;
  }

  mri = mriIsolateLabel(mri_seg, label, bbox);

  mrisClearMomentum(mris);
  MRIScomputeMetricProperties(mris);
  MRISstoreMetricProperties(mris);
  MRIScomputeNormals(mris);

  switch (connectivity) {
    case 1:
      mrisSetVal(mris, 0.65);
      break;
    case 2:
      mrisSetVal(mris, 0.5);
      break;
    case 3:
      mrisSetVal(mris, 0.75);
      break;
    case 4:
      mrisSetVal(mris, 0.3);
      break;
    default:
      mrisSetVal(mris, 0.5);
      break;
  }

  last_rms = rms = mrisRmsValError(mris, mri);

  mris->noscale = 1;
  last_sse = sse = MRIScomputeSSE(mris, parms) / (double)mris->nvertices;

  if (1)  // Gdiag & DIAG_SHOW)
    fprintf(stdout, "%3.3d: dt: %2.4f, sse=%2.4f, rms=%2.4f\n", 0, 0.0f, (float)sse, (float)rms);

  for (n = 0; n < niterations; n++) {
    dt = base_dt;
    nreductions = 0;

    MHTfree(&mht);
    mht = MHTcreateFaceTable(mris);

    mrisClearGradient(mris);
    MRISstoreMetricProperties(mris);
    MRIScomputeNormals(mris);

    /*intensity based terms*/
    mrisComputeIntensityTerm(mris, parms->l_intensity, mri, parms->sigma);

    mrisAverageSignedGradients(mris, avgs);

    /*smoothness terms*/
    mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv);
    mrisComputeNormalSpringTerm(mris, parms->l_nspring);
    mrisComputeTangentialSpringTerm(mris, parms->l_tspring);

    //      sprintf(fname,"./zurf%d",n);
    // MRISwrite(mris,fname);

    do {
      MRISsaveVertexPositions(mris, TMP_VERTICES);
      delta_t = mrisAsynchronousTimeStepNew(mris, parms->momentum, dt, mht, MAX_ASYNCH_MM);

      MRIScomputeMetricProperties(mris);
      rms = mrisRmsValError(mris, mri);
      sse = MRIScomputeSSE(mris, parms) / (double)mris->nvertices;

      done = 1;
      if (rms > last_rms - tol) /*error increased - reduce step size*/
      {
        nreductions++;
        dt *= REDUCTION_PCT;
        if (0)  // Gdiag & DIAG_SHOW)
          fprintf(stdout,
                  "      sse=%2.1f, last_sse=%2.1f,\n"
                  "      ->  time setp reduction %d of %d to %2.3f...\n",
                  sse,
                  last_sse,
                  nreductions,
                  MAX_REDUCTIONS + 1,
                  dt);

        mrisClearMomentum(mris);
        if (rms > last_rms) /*error increased - reject step*/
        {
          MRISrestoreVertexPositions(mris, TMP_VERTICES);
          MRIScomputeMetricProperties(mris);
          done = (nreductions > MAX_REDUCTIONS);
        }
      }
    } while (!done);
    last_sse = sse;
    last_rms = rms;

    rms = mrisRmsValError(mris, mri);
    sse = MRIScomputeSSE(mris, parms) / (double)mris->nvertices;
    fprintf(stdout,
            "%3d, sse=%2.1f, last sse=%2.1f,\n"
            "     rms=%2.4f, last_rms=%2.4f\n",
            n,
            sse,
            last_sse,
            rms,
            last_rms);
  }

  msec = then.milliseconds();
  if (1)  // Gdiag & DIAG_SHOW)
    fprintf(stdout, "positioning took %2.2f minutes\n", (float)msec / (60 * 1000.0f));

  if (bbox_indicator) free(bbox);
  if (parms_indicator) free(parms);
  MRIfree(&mri);
  MHTfree(&mht);

  return mris;
}

// smooth a surface 'niter' times with a step (should be around 0.5)
void MRISsmoothSurface2(MRI_SURFACE *mris, int niter, float step, int avrg)
{

  if (step > 1) step = 1.0f;

  int iter;
  for (iter = 0; iter < niter; iter++) {

    MRIScomputeMetricProperties(mris);

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
      
      float x = 0, y = 0, z = 0;
      
      int m;
      for (m = 0; m < vt->vnum; m++) {
        x += mris->vertices[vt->v[m]].tx;
        y += mris->vertices[vt->v[m]].ty;
        z += mris->vertices[vt->v[m]].tz;
      }
      x /= vt->vnum;
      y /= vt->vnum;
      z /= vt->vnum;

      v->dx = step * (x - v->x);
      v->dy = step * (y - v->y);
      v->dz = step * (z - v->z);
    }
    
    mrisAverageSignedGradients(mris, avrg);

    MRISfreeDistsButNotOrig(mris);
      // MRISsetXYZ will invalidate all of these,
      // so make sure they are recomputed before being used again!
    
    for (k = 0; k < mris->nvertices; k++) {
      VERTEX * const v = &mris->vertices[k];
      MRISsetXYZ(mris, k,
        v->x + v->dx,
        v->y + v->dy,
        v->z + v->dz);
    }
  }
}

/*--------------------------------------------------------------------*/
MRIS *MRISloadSurfSubject(const char *subj, const char *hemi, const char *surfid, const char *SUBJECTS_DIR)
{
  MRIS *Surf;
  char fname[2000];

  if (SUBJECTS_DIR == NULL) {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR == NULL) {
      printf("ERROR: SUBJECTS_DIR not defined in environment.\n");
      return (NULL);
    }
  }

  sprintf(fname, "%s/%s/surf/%s.%s", SUBJECTS_DIR, subj, hemi, surfid);
  printf("  INFO: loading surface  %s\n", fname);
  fflush(stdout);
  Surf = MRISread(fname);
  if (Surf == NULL) {
    printf("ERROR: could not load registration surface\n");
    exit(1);
  }
  printf("nvertices = %d\n", Surf->nvertices);
  fflush(stdout);

  return (Surf);
}

/*-----------------------------------------------------------------
  MRISfdr2vwth() - computes the voxel-wise (or vertex-wise) threshold
    needed to realize the given False Discovery Rate (FDR) based on the
    values in the val field. The val field is copied to the val2 field,
    and then the val2 field is thresholded.

  fdr - false dicovery rate, between 0 and 1, eg: .05
  signid -
      0 = use all values regardless of sign
     +1 = use only positive values
     -1 = use only negative values
     If a vertex does not meet the sign criteria, its val2 is 0
  log10flag - interpret val field as -log10(p)
  maskflag - use the undefval field as a mask. If the undefval of
     a vertex is 1, then its val will be used to compute the threshold
     (if the val also meets the sign criteria). If the undefval is
     0, then val2 will be set to 0.
  vwth - voxel-wise threshold between 0 and 1. If log10flag is set,
     then vwth = -log10(vwth). Vertices with p values
     GREATER than vwth have val2=0. Note that this is the same
     as requiring -log10(p) > -log10(vwth).

  So, for the val2 to be set to something non-zero, the val must
  meet the sign, mask, and threshold criteria. If val meets all
  the criteria, then val2=val (ie, no log10 transforms). The val
  field itself is not altered.

  Return Values:
    0 - evertying is OK
    1 - no vertices met the mask and sign criteria

  Ref: http://www.sph.umich.edu/~nichols/FDR/FDR.m
  Thresholding of Statistical Maps in Functional Neuroimaging Using
  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,
  Thomas E. Nichols (2002).  NeuroImage 15:870-878.

  See also: fdr2vwth() in sig.c
  ---------------------------------------------------------------*/
int MRISfdr2vwth(MRIS *surf, double fdr, int signid, int log10flag, int maskflag, double *vwth)
{
  double *p = NULL, val, val2null;
  int vtxno, np;

  if (log10flag)
    val2null = 0;
  else
    val2null = 1;

  p = (double *)calloc(surf->nvertices, sizeof(double));
  np = 0;
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    if ((maskflag && !surf->vertices[vtxno].undefval) || surf->vertices[vtxno].ripflag) continue;
    val = surf->vertices[vtxno].val;
    if (signid == -1 && val > 0) continue;
    if (signid == +1 && val < 0) continue;
    val = fabs(val);
    if (log10flag) val = pow(10, -val);
    p[np] = val;
    np++;
  }

  // Check that something met the match criteria,
  // otherwise return an error
  if (np == 0) {
    printf("WARNING: MRISfdr2vwth(): no vertices met threshold\n");
    free(p);
    return (1);
  }

  *vwth = fdr2vwth(p, np, fdr);

  // Perform the thresholding
  for (vtxno = 0; vtxno < surf->nvertices; vtxno++) {
    if ((maskflag && !surf->vertices[vtxno].undefval) || surf->vertices[vtxno].ripflag) {
      // Set to null if masking and not in the mask
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }
    val = surf->vertices[vtxno].val;
    if (signid == -1 && val > 0) {
      // Set to null if wrong sign
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }
    if (signid == +1 && val < 0) {
      // Set to null if wrong sign
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }

    val = fabs(val);
    if (log10flag) val = pow(10, -val);

    if (val > *vwth) {
      // Set to null if greather than thresh
      surf->vertices[vtxno].val2 = val2null;
      continue;
    }

    // Otherwise, it meets all criteria, so
    // pass the original value through
    surf->vertices[vtxno].val2 = surf->vertices[vtxno].val;
  }

  // Change the vwth to log10 if needed
  if (log10flag) *vwth = -log10(*vwth);
  free(p);

  printf("MRISfdr2vwth(): np = %d, nv = %d, fdr = %g, vwth=%g\n", np, surf->nvertices, fdr, *vwth);

  return (0);
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*---------------------------------------------------------------*/
int MRISfwhm2niters(double fwhm, MRIS *surf)
{
  double avgvtxarea, gstd;
  int niters;

  MRIScomputeMetricProperties(surf);
  avgvtxarea = surf->total_area / surf->nvertices;

  if (surf->group_avg_surface_area > 0) {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    avgvtxarea *= (surf->group_avg_surface_area / surf->total_area);
    // if(getenv("FIX_VERTEX_AREA") != NULL)    {
    // printf("INFO: fwhm2niters: Fixing group surface area\n");
    // avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //}
    // else printf("INFO: fwhm2niters: NOT fixing group surface area\n");
  }

  gstd = fwhm / sqrt(log(256.0));
  // 1.14 is a fudge factor based on empirical fit of nearest neighbor
  niters = floor(1.14 * (4 * PI * (gstd * gstd)) / (7 * avgvtxarea) + 0.5);
  return (niters);
}

/*----------------------------------------------------------------------*/
double MRISniters2fwhm(int niters, MRIS *surf)
{
  double avgvtxarea, gstd, fwhm;

  MRIScomputeMetricProperties(surf);
  avgvtxarea = surf->total_area / surf->nvertices;

  if (surf->group_avg_surface_area > 0) {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    avgvtxarea *= (surf->group_avg_surface_area / surf->total_area);
    // if (getenv("FIX_VERTEX_AREA") != NULL)    {
    // printf("INFO: niters2fwhm: Fixing group surface area\n");
    // avgvtxarea *= (surf->group_avg_surface_area/surf->total_area);
    //}
    // else printf("INFO: fwhm2niters: NOT fixing group surface area\n");
  }

  gstd = sqrt(7 * avgvtxarea * niters / (1.14 * 4 * PI));
  fwhm = gstd * sqrt(log(256.0));
  return (fwhm);
}

/*---------------------------------------------------------------*/
int MRISfwhm2nitersSubj(double fwhm, const char *subject, const char *hemi, const char *surfname)
{
  char *SUBJECTS_DIR, surfpath[2000];
  MRIS *surf;
  int niters;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR == NULL) {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    return (-1);
  }
  sprintf(surfpath, "%s/%s/surf/%s.%s", SUBJECTS_DIR, subject, hemi, surfname);

  surf = MRISread(surfpath);
  if (surf == NULL) {
    printf("ERROR: could not read %s\n", surfpath);
    return (-1);
  }
  MRIScomputeMetricProperties(surf);

  niters = MRISfwhm2niters(fwhm, surf);

  return (niters);
}

/*!
  \fn MRI *MRISfwhmFromAR1Map(MRIS *surf, MRI *mask, MRI *ar1map)
  \brief Computes FWHM for each vertex
 */
MRI *MRISfwhmFromAR1Map(MRIS *surf, MRI *mask, MRI *ar1map)
{
  double ar1, fwhm, gstd, InterVertexDistAvg;
  MRIScomputeMetricProperties(surf);

  InterVertexDistAvg = surf->avg_vertex_dist;
  if (surf->group_avg_surface_area > 0)
    InterVertexDistAvg *= sqrt(surf->group_avg_surface_area / surf->total_area);

  MRI *fwhmmap = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  MRIcopyHeader(ar1map,fwhmmap);
  MRIcopyPulseParameters(ar1map,fwhmmap);

  int n;
  for(n=0; n < surf->nvertices; n++){
    if(mask && MRIgetVoxVal(mask,n,0,0,0)<0.5) continue;
    ar1 = MRIgetVoxVal(ar1map,n,0,0,0);
    // Cant just call MRISfwhmFromAR1() because it runs compute metric properties
    if(ar1 <= 0.0) continue;
    gstd = InterVertexDistAvg / sqrt(-4 * log(ar1));
    fwhm = gstd * sqrt(log(256.0));
    MRIsetVoxVal(fwhmmap,n,0,0,0,fwhm);
  }
  return(fwhmmap);
}


/*----------------------------------------------------------------------
  MRISfwhm() - estimates fwhm from global ar1 mean
  ----------------------------------------------------------------------*/
double MRISfwhmFromAR1(MRIS *surf, double ar1)
{
  double fwhm, gstd, InterVertexDistAvg;

  MRIScomputeMetricProperties(surf);
  InterVertexDistAvg = surf->avg_vertex_dist;
  if (surf->group_avg_surface_area > 0) {
    // This should be ok even if metric properties have been scaled ??
    // Always do this now (4/9/10)
    InterVertexDistAvg *= sqrt(surf->group_avg_surface_area / surf->total_area);
    // if (getenv("FIX_VERTEX_AREA") != NULL)    {
    // printf("INFO: fwhmFromAR1: Fixing group surface area\n");
    // InterVertexDistAvg *=
    //  sqrt(surf->group_avg_surface_area/surf->total_area);
    //}
    // else printf("INFO: fwhmFromAR1: NOT fixing group surface area\n");
  }

  if (ar1 > 0.0) {
    gstd = InterVertexDistAvg / sqrt(-4 * log(ar1));
    fwhm = gstd * sqrt(log(256.0));
  }
  else {
    printf("WARNING: ar1 = %g <= 0. Setting fwhm to 0.\n", ar1);
    fwhm = 0.0;
  }

  return (fwhm);
}

/*-------------------------------------------------------------------------
  MRISsmoothingArea() - computes the area coverted by the give number
  of iterations at the given vertex. Essentially it's the area of the
  neighborhood.
  -------------------------------------------------------------------------*/
double MRISsmoothingArea(MRIS *mris, int vtxno, int niters)
{
  MRI *mri;
  int n, nhits;
  double val, area;

  // alloc a surface mri, zeros for everybody
  mri = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT);

  // create a delta function at the target vertex
  MRIsetVoxVal(mri, vtxno, 0, 0, 0, 1);

  // smooth it by the number of iterations
  MRISsmoothMRI(mris, mri, niters, NULL, mri);

  // find the non-zero vertices. these are the vertices in the neighborhood
  // add up the number of vertices and area
  nhits = 0;
  area = 0.0;
  for (n = 0; n < mris->nvertices; n++) {
    val = MRIgetVoxVal(mri, n, 0, 0, 0);
    if (val > 0.0) {
      nhits++;
      if (mris->group_avg_vtxarea_loaded)
        area += mris->vertices[n].group_avg_area;
      else
        area += mris->vertices[n].area;
    }
  }
  MRIfree(&mri);
  printf("%6d  %3d   %4d %7.1f\n", vtxno, niters, nhits, area);
  fflush(stdout);

  return (area);
}

/*-------------------------------------------------------------------
  MRISseg2annot() - converts a segmentation in a surface-encoded
  MRI struct to an annotation in an MRIS struct. If ctab is null,
  it tries to use mris->ct. See also MRISannotIndex2Seg().
  -------------------------------------------------------------------*/
int MRISseg2annot(MRIS *mris, MRI *surfseg, COLOR_TABLE *ctab)
{
  int vtxno, segid, ano;

  if(ctab == NULL) {
    if(mris->ct == NULL) {
      printf("ERROR: MRISseg2annot: both ctab and mris->ct are NULL\n");
      return (1);
    }
  }
  else mris->ct = ctab;
  set_atable_from_ctable(mris->ct);

  int nhits = 0, nmisses=0;
  for (vtxno = 0; vtxno < mris->nvertices; vtxno++) {
    segid = MRIgetVoxVal(surfseg, vtxno, 0, 0, 0);
    ano = index_to_annotation(segid);
    if(ano != -1 && segid != 0) nhits ++;
    if(ano == -1 && segid != 0) nmisses ++;
    if(ano == -1) ano = 0;
    mris->vertices[vtxno].annotation = ano;
    if (vtxno == Gdiag_no)
      printf("%5d %2d %2d %s\n",vtxno,segid,ano,index_to_name(segid));
  }
  if(nhits == 0) printf("WARNING: MRISseg2annot: no voxels matching ctab\n");
  if(nmisses>0)  printf("WARNING: MRISseg2annot: found %d voxels with no entry in ctab\n",nmisses);

  return (0);
}

/*----------------------------------------------------------------
  MRISannotIndex2Seg() - creates an MRI struct where each voxel is the
  annotation index. The struct has nvertices columns, 1 row, 1 slice,
  and 1 frame. It should not be misconstrued as a volume. See also
  MRISseg2annot().
  ----------------------------------------------------------------*/
MRI *MRISannotIndex2Seg(MRIS *mris)
{
  MRI *seg;
  int vno, annot, annotid;

  annotid = 0;
  seg = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT);
  for (vno = 0; vno < mris->nvertices; vno++) {
    annot = mris->vertices[vno].annotation;
    if (mris->ct)
      CTABfindAnnotation(mris->ct, annot, &annotid);
    else
      annotid = annotation_to_index(annot);
    MRIsetVoxVal(seg, vno, 0, 0, 0, annotid);
  }
  return (seg);
}

/*!/
  \fn double MRISvolumeInSurf(MRIS *mris)
  \brief Computes the volume inside a surface. (Xiao)
*/
double MRISvolumeInSurf(MRIS *mris)
{
  int fno;
  double total_volume, face_area;
  VECTOR *v_a, *v_b, *v_n, *v_cen;
  VERTEX *v0, *v1, *v2;
  FACE *face;

  v_a = VectorAlloc(3, MATRIX_REAL);
  v_b = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL);   /* normal vector */
  v_cen = VectorAlloc(3, MATRIX_REAL); /* centroid vector */

  total_volume = 0;
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) continue;

    v0 = &mris->vertices[face->v[0]];
    v1 = &mris->vertices[face->v[1]];
    v2 = &mris->vertices[face->v[2]];

    VERTEX_EDGE(v_a, v0, v1);
    VERTEX_EDGE(v_b, v0, v2);

    /* face normal vector */
    V3_CROSS_PRODUCT(v_a, v_b, v_n);
    face_area = V3_LEN(v_n) * 0.5f;

    V3_NORMALIZE(v_n, v_n); /* make it a unit vector */

    /* compute face centroid */
    V3_X(v_cen) = (v0->x + v1->x + v2->x) / 3.0;
    V3_Y(v_cen) = (v0->y + v1->y + v2->y) / 3.0;
    V3_Z(v_cen) = (v0->z + v1->z + v2->z) / 3.0;

    total_volume += V3_DOT(v_cen, v_n) * face_area;
  }

  MatrixFree(&v_cen);
  MatrixFree(&v_a);
  MatrixFree(&v_b);
  MatrixFree(&v_n);

  total_volume /= 3.0;
  return (total_volume);
}
/*!
  \fn MRISvolumeTH3(MRIS *w, MRIS *p, MRI *vol)
  \brief Compute vertex-wise volume based on dividing each obliquely truncated
  trilateral pyramid into three tetrahedra. Based on Anderson M. Winkler's
  srf2vol matlab script.
  \params w is white surface, p is pial. Output is an MRI struct of size nvertices
 */
MRI *MRISvolumeTH3(MRIS *w, MRIS *p, MRI *vol, MRI *mask, double *totvol)
{
  MRI *facevol;
  int nthface, v1, v2, v3;
  int vno, nfacenbrs, nthfacenbr, nbrfaceno;
  MATRIX *Ap, *Bp, *Cp, *Aw, *Bw, *Cw;
  MATRIX *C1 = NULL, *C2 = NULL, *C3 = NULL;
  double T1, T2, T3, volsum, dval, vtxvolsum;
  VERTEX *vtx;

  if (!vol) vol = MRIalloc(w->nvertices, 1, 1, MRI_FLOAT);

  Ap = MatrixAlloc(3, 1, MATRIX_REAL);
  Bp = MatrixAlloc(3, 1, MATRIX_REAL);
  Cp = MatrixAlloc(3, 1, MATRIX_REAL);
  Aw = MatrixAlloc(3, 1, MATRIX_REAL);
  Bw = MatrixAlloc(3, 1, MATRIX_REAL);
  Cw = MatrixAlloc(3, 1, MATRIX_REAL);

  facevol = MRIalloc(w->nfaces, 1, 1, MRI_FLOAT);
  volsum = 0.0;
  for (nthface = 0; nthface < w->nfaces; nthface++) {
    MRIsetVoxVal(facevol, nthface, 0, 0, 0, 0.0);
    if (w->faces[nthface].ripflag) continue;

    v1 = w->faces[nthface].v[0];
    v2 = w->faces[nthface].v[1];
    v3 = w->faces[nthface].v[2];

    Ap->rptr[1][1] = p->vertices[v1].x;
    Ap->rptr[2][1] = p->vertices[v1].y;
    Ap->rptr[3][1] = p->vertices[v1].z;
    Bp->rptr[1][1] = p->vertices[v2].x;
    Bp->rptr[2][1] = p->vertices[v2].y;
    Bp->rptr[3][1] = p->vertices[v2].z;
    Cp->rptr[1][1] = p->vertices[v3].x;
    Cp->rptr[2][1] = p->vertices[v3].y;
    Cp->rptr[3][1] = p->vertices[v3].z;

    Aw->rptr[1][1] = w->vertices[v1].x;
    Aw->rptr[2][1] = w->vertices[v1].y;
    Aw->rptr[3][1] = w->vertices[v1].z;
    Bw->rptr[1][1] = w->vertices[v2].x;
    Bw->rptr[2][1] = w->vertices[v2].y;
    Bw->rptr[3][1] = w->vertices[v2].z;
    Cw->rptr[1][1] = w->vertices[v3].x;
    Cw->rptr[2][1] = w->vertices[v3].y;
    Cw->rptr[3][1] = w->vertices[v3].z;

    // As the Ap is the common vertex for all three, it can be used as the origin.
    Bp = MatrixSubtract(Bp, Ap, Bp);
    Cp = MatrixSubtract(Cp, Ap, Cp);
    Aw = MatrixSubtract(Aw, Ap, Aw);
    Bw = MatrixSubtract(Bw, Ap, Bw);
    Cw = MatrixSubtract(Cw, Ap, Cw);

    // The next lines compute the volume for each, using a scalar triple product:
    C1 = VectorCrossProduct(Bw, Cw, C1);
    T1 = fabs(V3_DOT(Aw, C1));

    C2 = VectorCrossProduct(Cp, Bw, C2);
    T2 = fabs(V3_DOT(Bp, C2));

    C3 = VectorCrossProduct(Cw, Bw, C3);
    T3 = fabs(V3_DOT(Cp, C3));

    dval = (T1 + T2 + T3) / 6.0;
    MRIsetVoxVal(facevol, nthface, 0, 0, 0, dval);
    volsum += dval;
  }
  printf("Total face volume %g\n", volsum);

  vtxvolsum = 0;
  for (vno = 0; vno < w->nvertices; vno++) {
    vtx = &w->vertices[vno];
    if (vtx->ripflag) continue;
    if (mask && MRIgetVoxVal(mask, vno, 0, 0, 0) < 0.5) continue;
    nfacenbrs = w->vertices_topology[vno].num;
    volsum = 0.0;
    for (nthfacenbr = 0; nthfacenbr < nfacenbrs; nthfacenbr++) {
      nbrfaceno = w->vertices_topology[vno].f[nthfacenbr];
      volsum += (MRIgetVoxVal(facevol, nbrfaceno, 0, 0, 0) / 3.0);
      // divide by 3 because each face participates in 3 vertices
    }
    MRIsetVoxVal(vol, vno, 0, 0, 0, volsum);
    vtxvolsum += volsum;
  }
  printf("Total vertex volume %g (mask=%d)\n", vtxvolsum, mask == NULL);
  *totvol = vtxvolsum;

  MRIfree(&facevol);
  MatrixFree(&Ap);
  MatrixFree(&Bp);
  MatrixFree(&Cp);
  MatrixFree(&Aw);
  MatrixFree(&Bw);
  MatrixFree(&Cw);
  MatrixFree(&C1);
  MatrixFree(&C2);
  MatrixFree(&C3);

  return (vol);
}

/*
  \fn LABEL *MRIScortexLabel(MRI_SURFACE *mris, MRI *mri_aseg, int min_vertices, int KeepHipAmyg)
  \brief This creates the label used for cortex.label. This label has
  two implications.  First, vertices outside of the label will be
  frozen in place (ie, white=pial) during surface placement.  Second,
  vertices outside of the label will be zeroed and/or masked out in
  many applications (eg, mris_preproc) so they will create holes in
  the group analysis. This label mostly defines the medial wall, but
  lesions ({Left,Right}_Lesion) will be excluded as well (but vertices
  near WMSAs are left intact because they are rarely next to cortex).
  Lesion vertices are not removed if they are in too small of a cluster.
  See also MRIScortexLabelDECC(). If KeepHipAmyg is set to 1, then vertices
  near hippocampus and amygdala are included in the label. In general,
  this is NOT desired as they are part of the medial wall. However, when
  placing the pial surface, it is beneficial to let these vertices float
  so that entorhinal and parahippo gyri are not affected. After the pial
  is placed, those vertices can be shifted back to the white surface.

  Note that if the v->ripflag or v->marked2 fields are set in vertices
  (e.g. from a .annot file), then these vertices will not be
  considered in the search for non-cortical vertices (that is, they
  will be labeled cortex).
*/
LABEL *MRIScortexLabel(MRI_SURFACE *mris, MRI *mri_aseg, int min_vertices, int KeepHipAmyg)  // BEVIN mris_make_surfaces 5
{
  LABEL *lcortex;
  int vno, label, nvox, total_vox, adjacent, x, y, z, target_label, l, base_label, left, right;
  VERTEX *v;
  double xv, yv, zv, val, xs, ys, zs, d;
  MRI_REGION box;

  // Control whether NucAcc is medial wall. Historically, it has been
  // if a NucAcc voxel was just outside the white.preaparc. NucAcc is
  // very close to the true medial wall, so sometimes these marked
  // vertices would get tangled up in the erode/dilate/cluster
  // operations on the medial wall (unlike putamen where there is a
  // clear distinction). In these cases, true cortex near the NucAcc
  // can get masked out. Changed so that by default NucAcc is now
  // NOT considered part of the medial wall. DNG 3/5/2020
  int NucAccIsMedialWall = 0; 
  if(getenv("FS_NUCACC_IS_MEDIAL_WALL") != NULL)
    sscanf(getenv("FS_NUCACC_IS_MEDIAL_WALL"),"%d",&NucAccIsMedialWall);

  printf(" Generating cortex label... RemoveHipAmgy=%d\n",KeepHipAmyg);
  printf(" NucAccIsMedialWall=%d\n",NucAccIsMedialWall);
  printf(" mris->useRealRAS=%d\n",mris->useRealRAS);

  mri_aseg = MRIcopy(mri_aseg, NULL);  // so we can mess with it

  // remove the posterior few mm of the ventricles to prevent
  // them poking into the calcarine
  int erase_mm = 3;
  for (l = 0; l < 2; l++) {
    if (l == 0)
      target_label = Left_Lateral_Ventricle;
    else
      target_label = Right_Lateral_Ventricle;
    MRIlabelBoundingBox(mri_aseg, target_label, &box);
    for (z = box.z + box.dz - (erase_mm + 1); z < box.z + box.dz; z++)
      for (y = 0; y < mri_aseg->height; y++)
        for (x = 0; x < mri_aseg->width; x++) {
          label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0);
          if (label == target_label) MRIsetVoxVal(mri_aseg, x, y, z, 0, 0);  // erase it
        }
  }

  // First, set all v->marked=1. Then removed vertices below
  MRISsetMarks(mris, 1);

  if(Gdiag_no > 0) printf("vno %d, Start: marked=%d\n", Gdiag_no, mris->vertices[Gdiag_no].marked);

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(serial)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

    v = &mris->vertices[vno];
    if (v->ripflag || v->marked2 > 0)  // already must be cortex
      continue;
    if (vno == Gdiag_no) DiagBreak();
    v->marked3 = 0;

    // don't sample inside here due to thin parahippocampal wm.
    // The other interior labels
    // will already have been ripped (and hence not marked, but all marks set to 1 above)
    base_label = 0;
    for (d = 0; d <= 2; d += 0.5) {
      xs = v->x + d * v->nx;
      ys = v->y + d * v->ny;
      zs = v->z + d * v->nz;
      if (mris->useRealRAS)
        MRIworldToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      else
        MRIsurfaceRASToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST);
      label = nint(val);
      if (FZERO(d)) base_label = label;

      if(vno == Gdiag_no) printf("vno %d d=%g label=%d marked=%d (%g %g %g) (%g %g %g) (%d %d %d)\n", 
    vno,d,label,v->marked,v->x,v->y,v->z,xs,ys,zs,(int)xv,(int)yv,(int)zv);

      if(label == Left_Lateral_Ventricle ||
          ( !KeepHipAmyg && (IS_HIPPO(label) || IS_AMYGDALA(label)) ) ||
          (IS_WM(label) && IS_WM(base_label) && (base_label != label)) ||  // crossed hemi staying in wm
          label == Right_Lateral_Ventricle || label == Third_Ventricle || 
          (NucAccIsMedialWall && (label == Left_Accumbens_area || label == Right_Accumbens_area)) ||
          label == Left_Caudate || label == Right_Caudate || IS_CC(label) ||
          label == Left_Pallidum || label == Right_Pallidum || 
          IS_LAT_VENT(label) || label == Third_Ventricle || label == Right_Thalamus ||
          label == Left_Thalamus || label == Brain_Stem || label == Left_VentralDC || label == Right_VentralDC) {
        if (label == Left_Putamen || label == Right_Putamen) DiagBreak();
        if (vno == Gdiag_no) DiagBreak();
        v->marked = 0;
        if (vno == Gdiag_no) printf("vno %d bad label %d  (%g %g %g)\n", vno,label,xv,yv,zv);
      }
      if(label == Left_Lesion || label == Right_Lesion) {
        v->marked = 0;
        v->marked3 = 1; // keep track of lesion
        if (vno == Gdiag_no) printf("vno %d is lesion\n", vno);
      }
    }

    // check to see if we are in region between callosum and thalamus on midline
    MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv);
    x = nint(xv);
    y = nint(yv);
    z = nint(zv);
    left = MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Left_Thalamus);
    right = MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Right_Thalamus);
    if (left && left >= right)  // near left thalamus
    {
      if (MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Left_Lateral_Ventricle) > 0) {
        if (vno == Gdiag_no) DiagBreak();
        v->marked = 0;
        if (vno == Gdiag_no) printf("vno %d near left thal\n", vno);
      }
    }
    else if (right && right >= left)  // near right thalamus
    {
      if (MRIlabelsInNbhd(mri_aseg, x, y, z, 2, Right_Lateral_Ventricle) > 0) {
        if (vno == Gdiag_no) DiagBreak();
        v->marked = 0;
        if (vno == Gdiag_no) printf("vno %d near right thal\n", vno);
      }
    }
    MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST);

    // putamen can be adjacent to insula in aseg, but shouldn't be inferior
    /* now check for putamen superior to this point. If there's a lot
       of it there, then we are in basal forebrain and not cortex. */
    for (adjacent = total_vox = nvox = 0, d = 0; d <= 10; d += 0.5, total_vox++) {
      xs = v->x;
      ys = v->y;
      zs = v->z + d;  // sample superiorly
      if (mris->useRealRAS)
        MRIworldToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      else
        MRIsurfaceRASToVoxel(mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST);
      label = nint(val);
      if (label == Left_Putamen || label == Right_Putamen) {
        nvox++;
        if (d < 1.5) adjacent = 1;
      }
    }
    if (adjacent && (double)nvox / (double)total_vox > 0.5){
      // more than 50% putamen
      v->marked = 0;
      if (vno == Gdiag_no) printf("vno %d more than 50%% putamen\n", vno);
    }

    if (v->marked)  // check to see if there is any cortical gm in the region in aseg
    {
      int whalf, lh, rh;
      whalf = 5;
      MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv);
      lh = MRIlabelsInNbhd(mri_aseg, xv, yv, zv, whalf, Left_Cerebral_Cortex);
      rh = MRIlabelsInNbhd(mri_aseg, xv, yv, zv, whalf, Right_Cerebral_Cortex);
      if (vno == Gdiag_no) DiagBreak();
      if (lh == 0 && rh == 0) {
        v->marked = 0;
        if (vno == Gdiag_no) printf("no cortical GM found in vicinity - removing %d  vertex from cortex\n", vno);
      }
    }

    ROMP_PFLB_end
  }
  ROMP_PF_end

  if(Gdiag_no > 0) printf("vno %d, Post Loop marked=%d\n", Gdiag_no, mris->vertices[Gdiag_no].marked);

  // remove small holes that shouldn't be non-cortex
  {
    LABEL **label_array;
    int nlabels, n, i;

    MRISinvertMarks(mris);  // marked->not cortex now
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "%s.marked.orig.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      MRISwriteMarked(mris, fname);
    }
    MRISdilateMarked(mris, 1);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "%s.marked.dilated.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      MRISwriteMarked(mris, fname);
    }
    MRISerodeMarked(mris, 1);  // do a first-order close on the non-cortex marks
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "%s.marked.closed.mgz", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      MRISwriteMarked(mris, fname);
    }
    MRISsegmentMarked(mris, &label_array, &nlabels, 0);
    printf("%d non-cortical segments detected\n", nlabels);
    if (min_vertices < 0)  // means only keep max segment
    {
      for (n = 0; n < nlabels; n++)
        if (label_array[n]->n_points > min_vertices) min_vertices = label_array[n]->n_points;
      printf("only using segment with %d vertices\n", min_vertices);
    }

    for (n = 0; n < nlabels; n++) {
      if (label_array[n]->n_points < min_vertices) {
        printf("erasing segment %d (vno[0] = %d)\n", n, label_array[n]->lv[0].vno);
        for (i = 0; i < label_array[n]->n_points; i++){
	  v = &(mris->vertices[label_array[n]->lv[i].vno]);
          if(v->marked3) continue; // don't unmark lesions
          v->marked = 0;  // mark it as cortex
	}
      }
      LabelFree(&label_array[n]);
    }
    free(label_array);
  }

  MRISinvertMarks(mris);  // marked --> is cortex again
  if (Gdiag_no >= 0) printf("v %d: is %sin cortex label\n", Gdiag_no, mris->vertices[Gdiag_no].marked ? "" : "NOT ");
  lcortex = LabelFromMarkedSurface(mris);

  MRIfree(&mri_aseg);  // a locally edited copy, not the original
  return (lcortex);
}

/*!
  \fn LABEL *MRIScortexLabelDECC(MRIS *mris, MRI *mri_aseg, int min_vertices, int KeepHipAmyg)
  \brief Creates a label of cortex that is better defined than that
  produced by MRIScortexLabel(). It runs MRIScortexLabel() first, then
  refines it by dilating (D) and eroding (E) then excluding all but
  the largest connected component (CC). This replaces what is in
  mris_make_surfaces when ndilate=nerode=4 and min_vertices=-1.
 */
LABEL *MRIScortexLabelDECC(MRIS *mris, MRI *mri_aseg, int ndilate, int nerode, int min_vertices, int KeepHipAmyg)
{
  LABEL *lcortex, **labels ;
  int   n, max_l, max_n, nlabels ;
  
  lcortex = MRIScortexLabel(mris, mri_aseg, min_vertices, KeepHipAmyg) ;
  LabelErode(lcortex, mris, nerode) ; // Erode the label by 4
  LabelDilate(lcortex, mris, ndilate, CURRENT_VERTICES) ;// Dilate the label by 4

  // Create clusters of the labels and take the biggest. This works on
  // the label to remove label islands. If there are wholes (eg, made
  // by lesions), then these are left intact. 
  MRISclearMarks(mris) ; // only clears v->marked
  LabelMark(lcortex, mris) ;
  MRISsegmentMarked(mris, &labels, &nlabels, 1) ; // 1=minlabelarea

  // Find the label with the most points in it
  max_n = 0 ;
  max_l = labels[0]->n_points ;
  for (n = 1 ; n < nlabels ; n++){
    if (labels[n]->n_points > max_l){
      max_l = labels[n]->n_points ;
      max_n = n ;
    }
  }

  // Unmark the vertices if they are not part of the biggest label 
  for (n = 0 ; n < nlabels ; n++){
    if (n != max_n){
      LabelUnmark(labels[n], mris) ;
    }
    LabelFree(&labels[n]) ;
  }
  LabelFree(&lcortex) ;
  
  // Get the final cortex label
  lcortex = LabelFromMarkedSurface(mris) ;
  return(lcortex);
 }

/*!
  \fn int MRISfindPath ( int *vert_vno, int num_vno, int max_path_length,
                         int *path, int *path_length, MRIS *mris ).
  \brief Finds a path that connects all the vertices listed in vert_vno
  \params vert_vno[] - array of vertex numbers to connect
  \params num_vno - number in array of vertex numbers to connect
  \params max_path_length - max number in path
  \params path - array of connected vertex numbers
  \params path_length - pointer to length of array of connected vertex numbers
  \params mris - surface
  This was copied from tksurfer.c find_path() and modified slightly.
*/
int MRISfindPath(int *vert_vno, int num_vno, int max_path_length, int *path, int *path_length, MRIS *mris)
{
  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char *check;
  float *dist;
  int *pred;
  char done;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;
  int flag2d = 0;  // for flattend surface?

  dist = (float *)calloc(mris->nvertices, sizeof(float));
  pred = (int *)calloc(mris->nvertices, sizeof(int));
  check = (char *)calloc(mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  for (cur_vert_vno = 0; cur_vert_vno < num_vno - 1; cur_vert_vno++) {
    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++) {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno + 1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices || dest_vno < 0 || dest_vno >= mris->nvertices) continue;

    if (src_vno == dest_vno) continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done) {
      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist) {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      VERTEX_TOPOLOGY const * const vt = &(mris->vertices_topology[closest_vno]);
      VERTEX          const * const v  = &(mris->vertices         [closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno) {
        done = TRUE;
      }
      else {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < vt->vnum; neighbor++) {
          neighbor_vno = vt->v[neighbor];
          VERTEX const * const u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
          if (flag2d)
            vu_z = 0;
          else
            vu_z = u->z - v->z;

          /* recalc the weight. */
          if (flag2d)
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) + ((v->y - u->y) * (v->y - u->y)));
          else
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) + ((v->y - u->y) * (v->y - u->y)) +
                           ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno]) {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0) {
        printf(".");
        fflush(stdout);
      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno && (*path_length) < max_path_length) {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  printf(" done\n");
  fflush(stdout);

  free(dist);
  free(pred);
  free(check);

  return (ERROR_NONE);
}

MRI *MRIScomputeFlattenedVolume(MRI_SURFACE *mris,
                                MRI *mri,
                                double res,
                                int nsamples,
                                int normalize,
                                MRI **pmri_vertices,
                                int smooth_iters,
                                double wm_dist,
                                double outside_dist)
{
  MRI *mri_flat, *mri_mask, *mri_counts, *mri_vno;
  int vno, width, height, u, v, w, fno, num;
  int uk, vk, ui, vi, whalf = 3, nv, wm_samples, outside_samples;
  double xmin, xmax, ymin, ymax, fdist, x, y, z, dx, dy, dz, norm, xf, yf, val, xv, yv, zv, oval, max_change;
  VERTEX *v0, *v1, *v2;
  FACE *face;
  MHT *mht;

  wm_samples = nint(wm_dist / res);
  outside_samples = nint(outside_dist / res);
  mht = MHTcreateFaceTable_Resolution(mris, FLATTENED_VERTICES, 2.0);
  ymax = xmax = -1e10;
  ymin = xmin = 1e10;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v0 = &mris->vertices[vno];
    if (v0->ripflag) continue;
    if (v0->fx < xmin) xmin = v0->fx;
    if (v0->fy < ymin) ymin = v0->fy;

    if (v0->fx > xmax) xmax = v0->fx;
    if (v0->fy > ymax) ymax = v0->fy;
  }

  width = ceil((xmax - xmin) / res);
  height = ceil((ymax - ymin) / res);
  mri_flat = MRIalloc(width, height, nsamples, MRI_FLOAT);
  mri_mask = MRIalloc(width, height, nsamples, MRI_UCHAR);
  mri_counts = MRIalloc(width, height, nsamples, MRI_INT);

  /*
    the first frame of mri_vno contains the # of vertices mapped to that (i,j) position, then
    the subsequent frames contain the vertex numbers
  */
  mri_vno = MRIalloc(width, height, nsamples, MRI_INT);
  MRIsetResolution(mri_flat, res, res, 3.0 / (float)nsamples);
  MRIsetResolution(mri_mask, res, res, 3.0 / (float)nsamples);
  MRIsetResolution(mri_vno, res, res, 3.0 / (float)nsamples);
  num = 0;
  whalf = ceil(2.0 / res);
  printf("using mask window size = %d\n", 2 * whalf + 1);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v0 = &mris->vertices[vno];
    u = ceil(v0->fx - xmin) / res;
    v = ceil(v0->fy - ymin) / res;
    for (uk = -whalf; uk <= whalf; uk++) {
      ui = mri_mask->xi[u + uk];
      for (vk = -whalf; vk <= whalf; vk++) {
        vi = mri_mask->yi[v + vk];
        if (MRIgetVoxVal(mri_mask, ui, vi, 0, 0) == 0) num++;
        MRIsetVoxVal(mri_mask, ui, vi, 0, 0, 1);
      }
    }
  }
  printf("%d voxels set in mask\n", num);

  for (num = u = 0; u < width; u++) {
    if (!(u % 100)) {
      printf("u = %d of %d\n", u, width);
      fflush(stdout);
    }
    for (v = 0; v < height; v++) {
      if (u == Gx && v == Gy) DiagBreak();
#if 1
      if (MRIgetVoxVal(mri_mask, u, v, 0, 0) == 0) continue;
#endif
      num++;
      xf = u * res + xmin;
      yf = v * res + ymin;
      MHTfindClosestFaceGeneric(mht, mris, xf, yf, 0, 1000, -1, -1, &face, &fno, &fdist);
      v0 = &mris->vertices[face->v[0]];
      v1 = &mris->vertices[face->v[1]];
      v2 = &mris->vertices[face->v[2]];
      if (v0->ripflag || v1->ripflag || v2->ripflag /* || fdist > 2*/) continue;
      if (v0 - mris->vertices == Gdiag_no || v1 - mris->vertices == Gdiag_no || v2 - mris->vertices == Gdiag_no)
        DiagBreak();

      dx = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->nx, v1->nx, v2->nx);
      dy = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->ny, v1->ny, v2->ny);
      dz = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->nz, v1->nz, v2->nz);
      norm = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
      if (FZERO(norm)) continue;
      dx /= norm;
      dy /= norm;
      dz /= norm;     // make it unit length
      if (normalize)  // divide the ribbon at this point into equally spaced samples
      {
        norm = (sqrt(SQR(v0->pialx - v0->whitex) + SQR(v0->pialy - v0->whitey) + SQR(v0->pialz - v0->whitez)) +
                sqrt(SQR(v1->pialx - v1->whitex) + SQR(v1->pialy - v1->whitey) + SQR(v1->pialz - v1->whitez)) +
                sqrt(SQR(v2->pialx - v2->whitex) + SQR(v2->pialy - v2->whitey) + SQR(v2->pialz - v2->whitez))) /
               3;
        norm += wm_dist + outside_dist;
        norm /= nsamples;  // divide average thickness into this many samples
      }
      else  // use uniform spacing
        norm = 1.0 / nsamples;

      dx *= norm;
      dy *= norm;
      dz *= norm;

      x = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->whitex, v1->whitex, v2->whitex);
      y = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->whitey, v1->whitey, v2->whitey);
      z = MRISsampleFace(mris, fno, FLATTENED_VERTICES, xf, yf, 0, v0->whitez, v1->whitez, v2->whitez);
      nv = MRIgetVoxVal(mri_vno, u, v, 0, 0);  // # of vertices mapped to this location
      vno = v0 - mris->vertices;
      for (w = 0; w < nv; w++) {
        int vno2;
        vno2 = (int)MRIgetVoxVal(mri_vno, u, v, w + 1, 0);
        if (vno2 == vno)  // already in the list
          break;
      }
      MRIsetVoxVal(mri_vno, u, v, nv + 1, 0, vno);
      if (w == nv) MRIsetVoxVal(mri_vno, u, v, 0, 0, nv + 1);  // increment # of vertices mapping here

      x -= dx * wm_dist;
      y -= dy * wm_dist;
      z -= dz * wm_dist;
      for (w = 0; w < nsamples; w++) {
        if (u == Gx && y == Gy && w == Gz) DiagBreak();
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
        MRIsampleVolume(mri, xv, yv, zv, &val);
        if (nint(xv) == Gx && nint(yv) == Gy && nint(zv) == Gz) DiagBreak();
        MRIsetVoxVal(mri_flat, u, v, w, 0, val);
        MRIsetVoxVal(mri_counts, u, v, w, 0, 1);
        if (v0 - mris->vertices == Gdiag_no || v1 - mris->vertices == Gdiag_no || v2 - mris->vertices == Gdiag_no)
          printf("(%2.1f %2.1f %2.1f) --> (%d, %d, %d) : %2.1f\n", xv, yv, zv, u, v, w, val);

        x += dx;
        y += dy;
        z += dz;
      }
      if (v0 - mris->vertices == Gdiag_no || v1 - mris->vertices == Gdiag_no || v2 - mris->vertices == Gdiag_no)
        DiagBreak();
    }
  }
  printf("%d voxel visited - %2.1f %% of total %d\n", num, 100.0 * num / (width * height), width * height);
  MRIremoveNaNs(mri_vno, mri_vno);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "vno.mgz");
    printf("writing vertex numbers to %s\n", fname);
    MRIwrite(mri_vno, fname);
  }
// fill in holes in the flatmap where no vertices mapped by averaging with neighboring filled locations
#define MAX_ITERS 50
  for (num = 0; num < MAX_ITERS; num++) {
    max_change = 0.0;
    for (u = 0; u < width; u++) {
      for (v = 0; v < height; v++) {
        for (w = 0; w < nsamples; w++) {
          if (MRIgetVoxVal(mri_flat, u, v, w, 0) > 30000) {
            fprintf(stderr, "voxel (%d, %d, %d) = %f\n", u, v, w, MRIgetVoxVal(mri_flat, u, v, w, 0));
            DiagBreak();
          }
          if (MRIgetVoxVal(mri_counts, u, v, w, 0) > 0)  // not a hole
            continue;
          for (val = 0.0, uk = -1; uk <= 1; uk++) {
            ui = mri_flat->xi[u + uk];
            for (vk = -1; vk <= 1; vk++) {
              vi = mri_flat->yi[v + vk];
              val += MRIgetVoxVal(mri_flat, ui, vi, w, 0);
            }
          }
          oval = MRIgetVoxVal(mri_flat, u, v, w, 0);
          val /= 9.0;  // # of voxels visited
          if (fabs(oval - val) > max_change) max_change = fabs(oval - val);
          MRIsetVoxVal(mri_flat, u, v, w, 0, val);
#if 1
          if (fabs(oval - val) < 1 && val > 50) MRIsetVoxVal(mri_counts, u, v, w, 0, 1);
#endif
        }
      }
    }
    if (max_change < 1) break;
    printf("%d of %d: max change %2.1f\n", num + 1, MAX_ITERS, max_change);
    fflush(stdout);
  }

  // now apply a bit of tangential smoothing
  for (num = 0; num < smooth_iters; num++) {
    for (u = 0; u < width; u++) {
      for (v = 0; v < height; v++) {
        for (w = 0; w < nsamples; w++) {
          for (val = 0.0, uk = -1; uk <= 1; uk++) {
            ui = mri_flat->xi[u + uk];
            for (vk = -1; vk <= 1; vk++) {
              vi = mri_flat->yi[v + vk];
              val += MRIgetVoxVal(mri_flat, ui, vi, w, 0);
            }
          }
          val /= 9.0;  // # of voxels visited
          MRIsetVoxVal(mri_flat, u, v, w, 0, val);
        }
      }
    }
  }

  MHTfree(&mht);
  MRIfree(&mri_mask);
  MRIfree(&mri_counts);
  if (pmri_vertices)
    *pmri_vertices = mri_vno;
  else
    MRIfree(&mri_vno);
  return (mri_flat);
}

/*!
  \fn MRI_SP *MRISmakeTemplate(int nsubjects, char **subjlist,
      int nhemis, char **hemilist, char *surfregname)
  \brief Creates a surface registration template. Overlaps with
    mris_make_template. Produces the same result as with -norot -aparc.
    Can be used on both lh and rh.
  \params number of subjects in subjlist
  \params subjlist - list of subjects
  \params nhemis - number of hemispheres in hemilist
  \params hemilist - list of hemispheres to use
  \params surfregname - name of surface registration
  Example 1:
    subjlist[0] = "s02.ghent";subjlist[1] = "s05.ghent";
    hemilist[0] = "lh";
    mrisp_template=MRISmakeTemplate(2, subjlist, 1, hemilist, "sphere.reg");
  Example 2:
    subjlist[0] = "s02.ghent";subjlist[1] = "s05.ghent";
    hemilist[0] = "rh";  hemilist[1] = "lh";
    mrisp_template=MRISmakeTemplate(2, subjlist, 2, hemilist, "sphere.left_right");
*/
MRI_SP *MRISmakeTemplate(int nsubjects, char **subjlist, int nhemis, char **hemilist, char *surfregname)
{
  static const char *surface_names[] = {"inflated", "smoothwm", "smoothwm"};
  static const char *curvature_names[] = {"inflated.H", "sulc", NULL};
  char tmpstr[2000];
  int images_per_surface = 3;
  int nsurfaces = sizeof(curvature_names) / sizeof(curvature_names[0]);
  int nparam_images = images_per_surface * nsurfaces;
  float scale = 1;
  const char *annot_name = "aparc";
  char *SUBJECTS_DIR, *hemi, *subject;
  INTEGRATION_PARMS parms;
  int which_norm = NORM_MEAN;
  int navgs = 0, nthhemi, sno, nthsubject;
  int nbrs = 3, err;
  MRI_SP *mrisp, /* *mrisp_aligned,*/ *mrisp_template;
  MRIS *mris;

  /* default template fields*/
  parms.nfields = 3;
  SetFieldLabel(&parms.fields[0], INFLATED_CURV_CORR_FRAME, 0, 0.0, 0.0, 0, which_norm);
  /* only use sulc for rigid registration */
  SetFieldLabel(&parms.fields[1], SULC_CORR_FRAME, 1, 1.0, 0.0, 0, which_norm);
  SetFieldLabel(&parms.fields[2], CURVATURE_CORR_FRAME, 2, 0.0, 0.0, 0, which_norm);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  mrisp_template = MRISPalloc(scale, nparam_images);
  for (nthsubject = 0; nthsubject < nsubjects; nthsubject++) {
    subject = subjlist[nthsubject];

    for (nthhemi = 0; nthhemi < nhemis; nthhemi++) {
      hemi = hemilist[nthhemi];
      printf("subject %s hemi %s\n", subject, hemi);
      sprintf(tmpstr, "%s/%s/surf/%s.%s", SUBJECTS_DIR, subject, hemi, surfregname);
      printf("   reading surface %s...\n", tmpstr);
      mris = MRISread(tmpstr);
      if (mris == NULL) {
        printf("ERROR: could not load %s\n", tmpstr);
        return (NULL);
      }

      err = MRISreadAnnotation(mris, annot_name);
      if (err) {
        printf("ERROR: could not load %s\n", annot_name);
        return (NULL);
      }

      MRISripMedialWall(mris);
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES);
      MRIScomputeMetricProperties(mris);
      MRISstoreMetricProperties(mris);

      for (sno = 0; sno < nsurfaces; sno++) {
        if (curvature_names[sno]) {
          /* read in precomputed curvature file */
          sprintf(tmpstr, "%s/%s/surf/%s.%s", SUBJECTS_DIR, subject, hemi, curvature_names[sno]);
          err = MRISreadCurvatureFile(mris, tmpstr);
          if (err) {
            printf("ERROR: could not load %s\n", tmpstr);
            return (NULL);
          }
          MRISaverageCurvatures(mris, navgs);
          MRISnormalizeCurvature(mris, which_norm);
        }
        else {
          sprintf(tmpstr, "%s/%s/surf/%s.%s", SUBJECTS_DIR, subject, hemi, surface_names[sno]);
          err = MRISreadVertexPositions(mris, tmpstr);
          if (err) {
            printf("ERROR: could not load %s\n", tmpstr);
            return (NULL);
          }
          if (nbrs > 1) MRISsetNeighborhoodSizeAndDist(mris, nbrs);
          MRIScomputeMetricProperties(mris);
          MRIScomputeSecondFundamentalForm(mris);
          MRISuseMeanCurvature(mris);
          MRISaverageCurvatures(mris, navgs);
          MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
          MRISnormalizeCurvature(mris, which_norm);
        }
        printf("  computing parameterization for surface %s...\n", tmpstr);
        mrisp = MRIStoParameterization(mris, NULL, scale, 0);
        MRISPcombine(mrisp, mrisp_template, sno * 3);
        MRISPfree(&mrisp);
      }
      MRISfree(&mris);
    }
  }
  return (mrisp_template);
}

/*
\fn int MRISsetPialUnknownToWhite(MRIS *white, MRIS *pial)
\brief Sets the vertex xyz of the pial equal to that of the white
  in vertices that have no annotation ("unknown"). Either the
  white or the pial must have an annotation. If both have an
  annotation, then white's is used.
*/
int MRISsetPialUnknownToWhite(const MRIS *white, MRIS *pial)
{
  int vtxno, UseWhite;
  int annot = 0, annotid = 0;

  if (white->ct == NULL && pial->ct == NULL) {
    printf("MRISsetPialUnknownToWhite(): neither white nor pial have an annotation\n");
    return (1);
  }

  if (white->ct != NULL)
    UseWhite = 1;
  else
    UseWhite = 0;
 
  MRISfreeDistsButNotOrig(pial);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate(annot, annotid)
#endif
  for (vtxno = 0; vtxno < white->nvertices; vtxno++) {
    ROMP_PFLB_begin
    
    // Convert annotation number to an entry number
    if (UseWhite) {
      annot = white->vertices[vtxno].annotation;
      CTABfindAnnotation(white->ct, annot, &annotid);
    }
    else {
      annot = pial->vertices[vtxno].annotation;
      CTABfindAnnotation(pial->ct, annot, &annotid);
    }
    if (annotid == -1 || white->vertices[vtxno].ripflag || pial->vertices[vtxno].ripflag) {
      MRISsetXYZ(pial,vtxno,
        white->vertices[vtxno].x,
        white->vertices[vtxno].y,
        white->vertices[vtxno].z);
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  MRIScomputeMetricProperties(pial);
  return (0);
}

/*!
  \fn int MRISshiftCRAS(MRIS *mris, int shift)
  \brief Shift surface xyz by scanner CRAS. This is not too useful but
  replicates some behavior in mris_convert that has been used by
  HCP. It results in the surface xyz being converted to (+1) or from
  (-1) scanner coords when the input is non-oblique. Better to use
  MRIStkr2Scanner() or MRISscanner2Tkr() which will have the same
  result for non-oblique but do the right thing for oblique.
 */
int MRISshiftCRAS(MRIS *mris, int shift)
{
  double dx, dy, dz;

  if (shift == +1) {
    mris->useRealRAS = 1;
    dx = +mris->vg.c_r;
    dy = +mris->vg.c_a;
    dz = +mris->vg.c_s;
    mris->xctr = mris->vg.c_r;
    mris->yctr = mris->vg.c_a;
    mris->zctr = mris->vg.c_s;
  }
  else {
    mris->useRealRAS = 0;
    dx = -mris->vg.c_r;
    dy = -mris->vg.c_a;
    dz = -mris->vg.c_s;
    mris->xctr = 0;
    mris->yctr = 0;
    mris->zctr = 0;
  }
  
  MRIStranslate(mris, dx, dy, dz);
    // The old code did not call mrisComputeSurfaceDimensions
    // leaving various values wrong
  
  return (0);
}

/*!
  \fn int MRIStkr2Scanner(MRIS *mris)
  \brief Convert surface xyz coords from tkregister space (ie, the
  native space that the surface xyz coords are generated in) to
  scanner space. The surface coords must be in tkr space, no checking
  is done to make sure, though useRealRAS must be 0. Use
  MRISscanner2Tkr() to reverse.  Sets mris->useRealRAS = 1;
 */
int MRIStkr2Scanner(MRIS *mris)
{
  if(mris->useRealRAS) return(0);
  MATRIX *M, *T, *Tinv, *S;
  S = vg_i_to_r(&mris->vg);
  T = TkrVox2RASfromVolGeom(&mris->vg);
  Tinv = MatrixInverse(T, NULL);
  M = MatrixMultiply(S, Tinv, NULL);
  MRISmatrixMultiply(mris, M);
  mris->useRealRAS = 1;
  MatrixFree(&S);
  MatrixFree(&T);
  MatrixFree(&Tinv);
  MatrixFree(&M);
  return (0);
}
/*!
  \fn int MRISscanner2Tkr(MRIS *mris)
  \brief Convert surface xyz coords from scanner space to tkregister
  space (ie, the native space that the surface xyz coords are
  generated in). The surface coords must be in scanner space, no
  checking is done to make sure. If useRealRAS=0, just returns without
  doing anything. Use MRIStkr2Scanner() to reverse.  Sets
  mris->useRealRAS = 0;
 */
int MRISscanner2Tkr(MRIS *mris)
{
  if(mris->useRealRAS == 0) return(0);
  MATRIX *Q, *T, *Sinv, *S;

  S = vg_i_to_r(&mris->vg);
  Sinv = MatrixInverse(S, NULL);
  T = TkrVox2RASfromVolGeom(&mris->vg);
  Q = MatrixMultiply(T, Sinv, NULL);
  MRISmatrixMultiply(mris, Q);
  mris->useRealRAS = 0;

  MatrixFree(&S);
  MatrixFree(&T);
  MatrixFree(&Sinv);
  MatrixFree(&Q);
  return (0);
}

/*!
  \fn int ComputeMRISvolumeTH3(char *subject, char *hemi, int DoMask, char *outfile)
  \brief Computes and saves the vertex-wise volume of cortex using TH3.
 */
int ComputeMRISvolumeTH3(char *subject, char *hemi, int DoMask, char *outfile)
{
  MRIS *w, *p;
  MRI *mrisvol;
  int err = 0;
  double totvol;
  FSENV *env;
  char fname[2000];
  LABEL *label;
  MRI *mask = NULL;

  env = FSENVgetenv();
  sprintf(fname, "%s/%s/surf/%s.white", env->SUBJECTS_DIR, subject, hemi);
  w = MRISread(fname);
  if (!w) return (1);
  sprintf(fname, "%s/%s/surf/%s.pial", env->SUBJECTS_DIR, subject, hemi);
  p = MRISread(fname);
  if (!p) return (1);

  if (DoMask) {
    sprintf(fname, "%s/%s/label/%s.cortex.label", env->SUBJECTS_DIR, subject, hemi);
    printf("masking with %s\n", fname);
    label = LabelRead(NULL, fname);
    if (label == NULL) return (1);
    mask = MRISlabel2Mask(w, label, NULL);
    if (mask == NULL) return (1);
    LabelFree(&label);
  }

  mrisvol = MRISvolumeTH3(w, p, NULL, mask, &totvol);
  printf("#@# %s %s %g\n", subject, hemi, totvol);

  if (IDextensionFromName(outfile)) {
    // output file has a known extention
    err = MRIwrite(mrisvol, outfile);
  }
  else {
    // otherwise assume it is a curv file
    MRIScopyMRI(w, mrisvol, 0, "curv");
    err = MRISwriteCurvature(w, outfile);
  }
  MRISfree(&w);
  MRISfree(&p);
  MRIfree(&mrisvol);
  if (mask) MRIfree(&mask);

  return (err);
}

/*!
  \fn int L2SaddPoint(LABEL2SURF *l2s, double col_or_vno, double row, double slice, int PointType, int Operation)

  \brief Adds (Operation==1) or removes (Operation!=1) a voxel or a
  vertex from a label based on its proximity to a surface. If
  PointType==0, then the point is treated like a voxel (col, row,
  slice). If PointType>0, then col_or_vno is treated as a vertex
  number (row and slice ignored), and the value of PointType is used
  to determine the surface number to use (surfno = PointType-1, so
  PointType==1 means the 1st surface).

  For label points derived from a voxel, if the CRS is within dmax to
  a surface, then it is assigned to the label of the closest
  surface. If it is further away, it is assigned to the label of the
  first surface (labels[0]) with vno = -1. XYZ of the CRS is
  determined from the template volume and a transform (can be NULL).
  Surface labels can be dilated along the surface based on the options
  set in the l2s structure.  The caller must have run L2Salloc() and
  L2Sinit() and set the appropriate variables in the L2S structure
  before running this function.  There is a label for each
  surface. Specifying multiple surfaces allows a label to be assigned,
  for example, to the lh but not the rh if both lh and rh surfaces are
  in the voxel. This function is not thread-safe. See also
  L2SaddVoxel() and L2Stest().

  Returns 1 if a point was added or removed and 0 if there was no change.

  Example usage:
  mri  = MRIread("template.nii.gz");
  surf = MRISread("lh.white");
  surf2 = MRISread("rh.white");
  lta = LTAread("register.dof6.lta"); // set to NULL for header reg
  l2s = L2Salloc(2, "bert"); // two surfaces, subject bert
  l2s->mri_template = mri;
  l2s->surfs[0] = surf;
  l2s->surfs[1] = surf2;
  l2s->dmax = 3; // voxel must be within 3mm of a vertexx
  l2s->hashres = 16;
  l2s->vol2surf = lta;
  l2s->nhopsmax = 10;
  L2Sinit(l2s);
  // start off with somelabel; any surface label points are assigned to
  // the 1st surface (surfno=0); volume label points are assigned to
  // the vollabel element
  L2SimportLabel(l2s, somelabel, 0);
  // add a surface vertex to label of the 1st surface
  L2SaddPoint(l2s, 110027, -1, -1, 1, 1);
  // remove some surface vertices from the label of the 2nd surface
  L2SaddPoint(l2s, 111010, -1, -1, 2, 0);
  L2SaddPoint(l2s,   5000, -1, -1, 2, 0);
  // add a voxel at CRS
  L2SaddPoint(l2s, col, row, slice, 0, 1); // add a voxel
  LabelWrite(l2s->labels[0],"./lh.label");
  LabelWrite(l2s->labels[1],"./rh.label");
  L2SaddPoint(l2s, col, row, slice, 0, 0); // remove a voxel
  L2Sfree(&l2s);
*/
int L2SaddPoint(LABEL2SURF *l2s, double col, double row, double slice, int PointType, int Operation)
{
  int n, nmin, vtxnominmin, vtxno, pointno;
  struct { float x,y,z; } v;
  static MATRIX *crs = NULL, *ras = NULL;
  float dminsurf;
  float dminmin = 0.0;
  LV *lv;
  LABEL *label;

  if (crs == NULL) {
    crs = MatrixAlloc(4, 1, MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }

  if (PointType != 0) {
    // point passed is a vertex from a surface
    vtxno = col;
    nmin = PointType - 1;
    label = l2s->labels[nmin];
    pointno = MRIgetVoxVal(l2s->masks[nmin], vtxno, 0, 0, 0);
    if (Operation == 1) {            // Add vertex to label
      if (pointno != 0) return (0);  // already there
      // If it gets here, then add the vertex
      if (l2s->debug) printf("Adding surf=%d vtxno=%d  np=%5d \n", nmin, vtxno, label->n_points);
      // Set value of the mask to the label point number + 1 (so 0=nolabel)
      MRIsetVoxVal(l2s->masks[nmin], vtxno, 0, 0, 0, label->n_points + 1);
      // Check whether need to alloc more points in label
      if (label->n_points >= label->max_points) LabelRealloc(label, nint(label->max_points * 1.5));
      // Finally, add this point
      lv = &(label->lv[label->n_points]);
      lv->vno = vtxno;
      lv->x = l2s->surfs[nmin]->vertices[vtxno].x;
      lv->y = l2s->surfs[nmin]->vertices[vtxno].y;
      lv->z = l2s->surfs[nmin]->vertices[vtxno].z;
      lv->stat = 0.0;
      // Incr the number of points in the label
      label->n_points++;
      return (1);
    }
    else {                           // Operation != 1, Remove vertex from label
      if (pointno == 0) return (0);  // not already there, cant remove it
      lv = &(label->lv[pointno - 1]);
      // If it gets here, then remove the vertex
      if (l2s->debug)
        printf(
            "Removing surf=%d vtxno=%d %g %g %g \n", nmin, vtxno, col, row, slice);
      lv->deleted = 1;
      MRIsetVoxVal(l2s->masks[nmin], lv->vno, 0, 0, 0, 0);
      return (1);
    }
    return (0);  // should never get here
  }

  // If it gets here, then point is from a volume
  crs->rptr[1][1] = col;
  crs->rptr[2][1] = row;
  crs->rptr[3][1] = slice;

  // Compute the TkReg RAS of this CRS in the surface volume space
  ras = MatrixMultiplyD(l2s->volcrs2surfxyz, crs, ras);

  // Load it into a vertex structure
  v.x = ras->rptr[1][1];
  v.y = ras->rptr[2][1];
  v.z = ras->rptr[3][1];
  if (l2s->debug > 1) printf("---------------\n%g %g %g   (%5.2f %5.2f %5.2f)\n", col, row, slice, v.x, v.y, v.z);

  // Go through each surface to find the surface and vertex in the surface that
  // the CRS is closest to.
  dminmin = 10e10;   // minimum distance over all surfaces
  vtxnominmin = -1;  // number of closest vertex
  nmin = -1;         // index of the surface with the closest vertex
  for (n = 0; n < l2s->nsurfs; n++) {
    vtxno = MHTfindClosestVertexNoXYZ(l2s->hashes[n], l2s->surfs[n], v.x,v.y,v.z, &dminsurf);
    if (vtxno >= 0) {
      if (l2s->debug > 1)
        printf("%3d %6d (%5.2f %5.2f %5.2f) %g\n",
               n,
               vtxno,
               l2s->surfs[n]->vertices[vtxno].x,
               l2s->surfs[n]->vertices[vtxno].y,
               l2s->surfs[n]->vertices[vtxno].z,
               dminsurf);
    }
    // should check whether there was a hash error?
    if (dminsurf < dminmin && dminsurf < l2s->dmax) {
      dminmin = dminsurf;
      vtxnominmin = vtxno;
      nmin = n;
    }
  }

  // If it does not meet the distance criteria, keep it as a volume point
  if (vtxnominmin == -1) {
    label = l2s->labels[0];  // add to the first surf label
    pointno = MRIgetVoxVal(l2s->volmask, round(col), round(row), round(slice), 0);
    if (Operation == 1) {            // Add voxel to label
      if (pointno != 0) return (0);  // already there
      // If it gets here, then add the voxel
      // Set value of the mask to the point number + 1
      MRIsetVoxVal(l2s->volmask, round(col), round(row), round(slice), 0, label->n_points + 1);
      // Check whether need to alloc more points in label
      if (label->n_points >= label->max_points) LabelRealloc(label, nint(label->max_points * 1.5));
      // Finally, add this point
      lv = &(label->lv[label->n_points]);
      lv->vno = -1;
      lv->x = ras->rptr[1][1];
      lv->y = ras->rptr[2][1];
      lv->z = ras->rptr[3][1];
      lv->stat = 0;
      // Incr the number of points in the label
      label->n_points++;
      return (1);  // return=1 because point has been added
    }
    else {                           // Operation != 1, Remove vertex from label
      if (pointno == 0) return (0);  // not already there, cant remove it
      lv = &(label->lv[pointno - 1]);
      lv->deleted = 1;
      MRIsetVoxVal(l2s->volmask, round(col), round(row), round(slice), 0, 0);
      return (1);  // return=1 because point has been removed
    }
  }

  // Select the label of the winning surface
  label = l2s->labels[nmin];

  // Below is a method to expand the label along the surface to
  // include other nearby vertices to the winning vertex. The hoplist
  // is a list of vertices for each ring of nearest neighbors. The 0th
  // hop is always the center vertex.
  SURFHOPLIST *shl;
  int nthnbr, nnbrs, nthhop, nthhoplast, nhits;
  int LabelChanged = 0;  // Return flag
  shl = SetSurfHopList(vtxnominmin, l2s->surfs[nmin], l2s->nhopsmax);
  nthhoplast = 0;
  for (nthhop = 0; nthhop < l2s->nhopsmax; nthhop++) {
    nnbrs = shl->nperhop[nthhop];  // number of neighbrs in the nthhop ring
    nhits = 0;

    // loop through the neighbors nthhop rings away
    for (nthnbr = 0; nthnbr < nnbrs; nthnbr++) {
      vtxno = shl->vtxlist[nthhop][nthnbr];

      if (l2s->DilateWithinVoxel) {
        // Require nearby vertices to be within the given voxel, where
        // "within" is defined as +/-0.5 of the voxel because the
        // voxel coords are at the center of the voxel
        v.x = l2s->surfs[nmin]->vertices[vtxno].tx;
        v.y = l2s->surfs[nmin]->vertices[vtxno].ty;
        v.z = l2s->surfs[nmin]->vertices[vtxno].tz;
        if (nthhop != 0) {
          // have to allow nthhop=0 because the closest vertex may be outside of the
          // voxel, depending upon the distance threshold
          // but could also use this to ignore dist thresh and force the source vtx
          // to be inside the voxel. but maybe better to do this before even starting
          // the hops?
          if (fabs(v.x - col) > 0.501 || fabs(v.y - row) > 0.501 || fabs(v.z - slice) > 0.501) continue;
        }
      }

      nthhoplast = nthhop;  // keep track of the most recent hop with a hit
      nhits++;              // keep track of the number of vertices hit in this hop

      // Get the 1-based label point number. Using 1-based so that 0
      // can code for no label at that point (yet)
      pointno = MRIgetVoxVal(l2s->masks[nmin], vtxno, 0, 0, 0);

      if (Operation == 1) {          // Add vertex to label
        if (pointno != 0) continue;  // already there
        // If it gets here, then add the vertex
        if (l2s->debug)
          printf("Adding surf=%d hop=%d vtxno=%d  np=%5d   %g %g %g   (%5.2f %5.2f %5.2f)\n",
                 nmin,
                 nthhop,
                 vtxno,
                 label->n_points,
                 col,
                 row,
                 slice,
                 v.x,
                 v.y,
                 v.z);
        // Set value of the mask to the label point number + 1 (so 0=nolabel)
        MRIsetVoxVal(l2s->masks[nmin], vtxno, 0, 0, 0, label->n_points + 1);

        // Check whether need to alloc more points in label
        if (label->n_points >= label->max_points) LabelRealloc(label, nint(label->max_points * 1.5));

        // Finally, add this point
        lv = &(label->lv[label->n_points]);
        lv->vno = vtxno;
        lv->x = l2s->surfs[nmin]->vertices[vtxno].x;
        lv->y = l2s->surfs[nmin]->vertices[vtxno].y;
        lv->z = l2s->surfs[nmin]->vertices[vtxno].z;
        lv->stat = dminmin;
        // Incr the number of points in the label
        label->n_points++;
        LabelChanged = 1;
      }
      else {                         // Operation != 1, Remove vertex from label
        if (pointno == 0) continue;  // not already there, cant remove it
        lv = &(label->lv[pointno - 1]);
        // If it gets here, then remove the vertex
        if (l2s->debug)
          printf("Removing surf=%d vtxno=%d %g %g %g   (%5.2f %5.2f %5.2f)\n",
                 nmin,
                 vtxno,
                 col,
                 row,
                 slice,
                 v.x,
                 v.y,
                 v.z);
        lv->deleted = 1;
        MRIsetVoxVal(l2s->masks[nmin], lv->vno, 0, 0, 0, 0);
        LabelChanged = 1;
      }
    }

    // If it has gone through a ring and not gotten a hit (nhits==0), then future hits
    // will not be contiguous with the current label. This can happen if the surface
    // folds back into the voxel. This can easily happen in fMRI-sized voxels. It is not
    // clear how this should be handled, so I put l2s->DilateContiguous into the structure
    // to allow or forbid disconinuous labels.
    if (l2s->DilateContiguous && nhits == 0) break;
  }
  SurfHopListFree(&shl);

  // check whether it might be a good idea to increase the number of hops
  if (nthhoplast == l2s->nhopsmax - 1)
    if (l2s->debug) printf("WARNING: hop saturation nthsurf=%d, cvtxno = %d\n", nmin, vtxnominmin);

  return (LabelChanged);
}

/*!
  \fn int L2Sinit(LABEL2SURF *l2s)
  \brief Initializes the L2S structure. The caller must have run
  L2Salloc() and set the surfs and vol2surf structures. If a header
  registration is good enough, then set vol2surf=NULL. Note that
  vol2surf can go in either direction; this function will determine
  which way to go from the template volume and surface vg geometries
  and handle it appropriately. This also computes the vertex coords in
  CRS space and stores them in the tx, ty, tz elements of the vertex,
  so the caller needs make sure that this is not overwriting something
  needed and to make sure not to overwrite it. Also computes a value
  for nhopsmax; the value is a reasonable guess based on the voxel
  size.
*/
int L2Sinit(LABEL2SURF *l2s)
{
  int n;

  // initialize hashes
  for (n = 0; n < l2s->nsurfs; n++) {
    l2s->hashes[n] = MHTcreateVertexTable_Resolution(l2s->surfs[n], CURRENT_VERTICES, l2s->hashres);
    if (l2s->hashes[n] == NULL) {
      printf("ERROR: L2Sinit(): MHTcreateVertexTable_Resolution() failed\n");
      return (-1);
    }
    l2s->masks[n] = MRIalloc(l2s->surfs[n]->nvertices, 1, 1, MRI_INT);
  }
  l2s->volmask =
      MRIallocSequence(l2s->mri_template->width, l2s->mri_template->height, l2s->mri_template->depth, MRI_INT, 1);
  MRIcopyHeader(l2s->mri_template, l2s->volmask);
  MRIcopyPulseParameters(l2s->mri_template, l2s->volmask);

  // compute the matrix that maps the template volume CRS to surface RAS
  // volcrs2surfxyz = K*inv(Vs)*R*Vv
  MATRIX *K, *Vs, *invVs, *Vv, *R;
  LTA *lta;
  K = TkrVox2RASfromVolGeom(&l2s->surfs[0]->vg);  // vox2tkras of surface
  Vs = vg_i_to_r(&l2s->surfs[0]->vg);             // vox2scanneras of surface
  Vv = MRIxfmCRS2XYZ(l2s->mri_template, 0);       // vox2scanneras of template volume
  invVs = MatrixInverse(Vs, NULL);
  if (invVs == NULL) {
    printf("ERROR: L2Sinit(): surf[0]->vg is null\n");
    return (-1);
  }

  if (l2s->vol2surf == NULL)
    R = MatrixIdentity(4, NULL);
  else {
    // A registration LTA has been passed, make sure it is consisent
    // with the input geometries and that it points in the right
    // direction.  Note that these function use the global
    // vg_isEqual_Threshold variable which should be set to something
    // small, like 10^-4
    int DoInvert = 0;
    VOL_GEOM vgvol;
    getVolGeom(l2s->mri_template, &vgvol);
    if (!vg_isEqual(&vgvol, &(l2s->vol2surf->xforms[0].src))) {
      // The src does not match the template, so try the dst
      if (!vg_isEqual(&l2s->surfs[0]->vg, &(l2s->vol2surf->xforms[0].src))) {
        printf("ERROR: L2Sinit(): neither registration vgs match template %g\n", vg_isEqual_Threshold);
        return (-1);
      }
      // Verify that the template matches the dst
      if (!vg_isEqual(&vgvol, &(l2s->vol2surf->xforms[0].dst))) {
        printf("ERROR: L2Sinit(): registration does not match volume vg %g\n", vg_isEqual_Threshold);
        return (-1);
      }
      DoInvert = 1;
    }
    else {
      // The source matches, but does the target?
      if (!vg_isEqual(&l2s->surfs[0]->vg, &(l2s->vol2surf->xforms[0].dst))) {
        printf("ERROR: L2Sinit(): registration does not match surface vg %g\n", vg_isEqual_Threshold);
        return (-1);
      }
    }
    // Copy the LTA
    lta = LTAcopy(l2s->vol2surf, NULL);
    // Make sure the type is RAS2RAS
    LTAchangeType(lta, LINEAR_RAS_TO_RAS);
    if (DoInvert) {
      if (l2s->debug) printf("L2Sinit(): inverting reg\n");
      R = MatrixInverse(lta->xforms[0].m_L, NULL);
    }
    else
      R = MatrixCopy(lta->xforms[0].m_L, NULL);
    LTAfree(&lta);
  }
  // Now finally compute it
  l2s->volcrs2surfxyz = MatrixMultiplyD(K, invVs, NULL);
  MatrixMultiplyD(l2s->volcrs2surfxyz, R, l2s->volcrs2surfxyz);
  MatrixMultiplyD(l2s->volcrs2surfxyz, Vv, l2s->volcrs2surfxyz);
  l2s->surfxyz2volcrs = MatrixInverse(l2s->volcrs2surfxyz, NULL);
  if (l2s->debug) {
    printf("L2Sinit(): \n");
    printf("K = [\n");
    MatrixPrint(stdout, K);
    printf("];\n");
    printf("invVs = [\n");
    MatrixPrint(stdout, invVs);
    printf("];\n");
    printf("R = [\n");
    MatrixPrint(stdout, R);
    printf("];\n");
    printf("Vv = [\n");
    MatrixPrint(stdout, Vv);
    printf("];\n");
    printf("volcrs2surfxyz = [\n");
    MatrixPrint(stdout, l2s->volcrs2surfxyz);
    printf("];\n");
  }
  MatrixFree(&K);
  MatrixFree(&Vs);
  MatrixFree(&Vv);
  MatrixFree(&invVs);
  MatrixFree(&R);

  // Now precompute the CRS coords of each vertex and save in the
  // t{xyz} part of the vertex structure. This makes constraining the
  // label dilation to within the voxel easier
  int vtxno;
  MATRIX *ras, *crs = NULL;
  ras = MatrixAlloc(4, 1, MATRIX_REAL);
  ras->rptr[4][1] = 1;
  for (n = 0; n < l2s->nsurfs; n++) {
    for (vtxno = 0; vtxno < l2s->surfs[n]->nvertices; vtxno++) {
      ras->rptr[1][1] = l2s->surfs[n]->vertices[vtxno].x;
      ras->rptr[2][1] = l2s->surfs[n]->vertices[vtxno].y;
      ras->rptr[3][1] = l2s->surfs[n]->vertices[vtxno].z;
      crs = MatrixMultiplyD(l2s->surfxyz2volcrs, ras, crs);
      l2s->surfs[n]->vertices[vtxno].tx = crs->rptr[1][1];
      l2s->surfs[n]->vertices[vtxno].ty = crs->rptr[2][1];
      l2s->surfs[n]->vertices[vtxno].tz = crs->rptr[3][1];
    }
  }
  MatrixFree(&ras);
  MatrixFree(&crs);

  // Take a stab at automatically computing the maximum number of hops. The basic
  // idea is to scale the number based upon the voxel size. The 3 is just something
  // I guessed.
  if (l2s->nhopsmax == -1)
    l2s->nhopsmax = 3 * round(sqrt(pow(l2s->mri_template->xsize, 2.0) + pow(l2s->mri_template->ysize, 2) +
                                   pow(l2s->mri_template->zsize, 2)));

  return (0);
}

/*!
  \fn LABEL2SURF *L2Salloc(int nsurfs, char *subject)
  \brief Allocs the relevant elements of the LABEL2SURF structure. Does not
  alloc the surfaces. subject can be null (this info is just passed to the
  labels).
*/
LABEL2SURF *L2Salloc(int nsurfs, const char *subject)
{
  int n;
  LABEL2SURF *l2s;
  l2s = (LABEL2SURF *)calloc(sizeof(LABEL2SURF), 1);
  l2s->surfs = (MRIS **)calloc(sizeof(MRIS *), nsurfs);
  l2s->hashes = (MHT **)calloc(sizeof(MHT *), nsurfs);
  l2s->nsurfs = nsurfs;
  l2s->labels = (LABEL **)calloc(sizeof(LABEL *), nsurfs);
  for (n = 0; n < nsurfs; n++) l2s->labels[n] = LabelAlloc(100, subject, NULL);
  l2s->masks = (MRI **)calloc(sizeof(MRI *), nsurfs);
  l2s->nhopsmax = -1;          // use -1 to keep track of whether it has been set or not
  l2s->DilateWithinVoxel = 1;  // turn on by default
  l2s->DilateContiguous = 0;   // turn off by default
  l2s->vollabel = LabelAlloc(100, subject, NULL);
  return (l2s);
}
/*!
  \fn int L2Sfree(LABEL2SURF **pl2s)
  \brief Frees the relevant elements of the LABEL2SURF structure. Does not
  free the surfaces.
*/
int L2Sfree(LABEL2SURF **pl2s)
{
  int n;
  LABEL2SURF *l2s = *pl2s;

  free(l2s->surfs);  // only free the pointer to the pointers
  for (n = 0; n < l2s->nsurfs; n++) {
    LabelFree(&l2s->labels[n]);
    MHTfree(&l2s->hashes[n]);
    MRIfree(&l2s->masks[n]);
  }
  MRIfree(&l2s->volmask);
  LabelFree(&l2s->vollabel);
  free(l2s->hashes);
  free(l2s->masks);
  MatrixFree(&l2s->volcrs2surfxyz);
  free(l2s);
  return (0);
}

/*!
  \fn int L2SaddVoxel(LABEL2SURF *l2s, double col, double row, double slice, int nsegs, int Operation)
  \brief Divides the voxel into nsegs parts and then runs L2SaddPoint(). This is an
  alternative way to include vertices inside the voxel that are not necessarily the closest
  vertex. The problem is that it adds any vertex regardless of how far away it is along the
  surface from the center vertex. L2SaddPoint() has a built-in mechanism for solving this
  problem, and it may be better. See L2SaddPoint() for argument descriptions.
*/
int L2SaddVoxel(LABEL2SURF *l2s, double col, double row, double slice, int nsegs, int Operation)
{
  double c, r, s, dseg;
  int ret, kc, kr, ks;

  if (nsegs == 1) {
    ret = L2SaddPoint(l2s, col, row, slice, 0, Operation);
    return (ret);
  }

  dseg = 1.0 / (nsegs - 1);
  // using nsegs as the upper limit (instead of nsegs+1) means that
  // the "end" side of the voxel is not included (but the "start" side
  // is).
  for (kc = 0; kc < nsegs; kc++) {
    c = col + kc * dseg - 0.5;
    for (kr = 0; kr < nsegs; kr++) {
      r = row + kr * dseg - 0.5;
      for (ks = 0; ks < nsegs; ks++) {
        s = slice + ks * dseg - 0.5;
        ret = L2SaddPoint(l2s, c, r, s, 0, Operation);
        if (ret < 0) return (ret);
      }
    }
  }
  // Make sure to add the center voxel
  ret = L2SaddPoint(l2s, col, row, slice, 0, Operation);

  return (ret);
}

/*!
  \fn int L2SimportLabel(LABEL2SURF *l2s, LABEL *label, int surfno)
  \brief Imports the given label into the L2S structure. The given
  label may be volume, surface, or mixed. If a label point is a
  volume-based label (ie, vno<0), then it is added to the vollabel
  element of L2S, otherwise it is added to the label of the surfno
  surface. The label points are simply copied; no attempt is made to
  convert the coordinates or check whether any of the label points
  might be there already.
 */
int L2SimportLabel(LABEL2SURF *l2s, LABEL *label, int surfno)
{
  int n;
  LV *lv;

  for (n = 0; n < label->n_points; n++) {
    lv = &(label->lv[n]);
    if (lv->vno < 0)
      LabelAddPoint(l2s->vollabel, lv);
    else
      LabelAddPoint(l2s->labels[surfno], lv);
  }

  return (0);
}

/*!
  \fn int L2Stest(char *subject)
  \brief Runs some tests on L2S routines. These are not exhaustive.
  subject defaults to bert if NULL.  Note: fails with average subjects
  because surf vg not set. Returns 0 if all tests passed or 1 if any
  test failed.
*/
int L2Stest(const char *subject)
{
  char *SUBJECTS_DIR, tmpstr[2000];
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  MRI *mri, *apas;
  MRIS *surf;
  int vno, c, r, s, ok, err, n;
  LV *lv;
  LABEL2SURF *l2s;
  LABEL *label;

  // Note: might fail with average subject because vg not set
  if (subject == NULL) subject = "bert";
  printf("L2Stest: subject %s\n", subject);

  sprintf(tmpstr, "%s/%s/mri/orig.mgz", SUBJECTS_DIR, subject);
  mri = MRIread(tmpstr);
  if (mri == NULL) return (1);

  sprintf(tmpstr, "%s/%s/mri/aparc+aseg.mgz", SUBJECTS_DIR, subject);
  apas = MRIread(tmpstr);
  if (apas == NULL) return (1);

  sprintf(tmpstr, "%s/%s/surf/lh.white", SUBJECTS_DIR, subject);
  surf = MRISread(tmpstr);

  sprintf(tmpstr, "%s/%s/label/lh.aparc.annot", SUBJECTS_DIR, subject);
  MRISreadAnnotation(surf, tmpstr);

  l2s = L2Salloc(1, "");
  l2s->mri_template = mri;
  l2s->surfs[0] = surf;
  l2s->dmax = 3;
  l2s->hashres = 16;
  l2s->vol2surf = NULL;
  l2s->debug = 0;
  L2Sinit(l2s);

  // Add a surface vertex
  vno = round(surf->nvertices / 2);
  printf("L2Stest: vno = %d\n", vno);
  L2SaddPoint(l2s, vno, -1, -1, 1, 1);
  if (l2s->labels[0]->n_points != 1) {
    printf("L2Stest: Failed to add surf vertex\n");
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  lv = &(l2s->labels[0]->lv[0]);
  if (lv->vno != vno) {
    printf("L2Stest: Failed to add correct surf vertex (%d)\n", lv->vno);
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  // Remove a surface vertex
  L2SaddPoint(l2s, vno, -1, -1, 1, 0);
  err = 0;
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (lv->vno == vno && !lv->deleted) {
      err = 1;
      break;
    }
  }
  if (err) {
    printf("L2Stest: Failed to add correct surf vertex (%d)\n", lv->vno);
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  l2s->labels[0]->n_points = 0;  // reset label

  // Add a voxel near a surface vertex.
  // A little circular because tx,ty,tz are computed by L2S
  c = round(surf->vertices[vno].tx);
  r = round(surf->vertices[vno].ty);
  s = round(surf->vertices[vno].tz);
  L2SaddPoint(l2s, c, r, s, 0, 1);
  ok = 0;
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (lv->vno == vno) {
      ok = 1;
      break;
    }
  }
  if (!ok) {
    printf("L2Stest: Failed to add correct vox-based surf vertex (%d)\n", lv->vno);
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  // Remove voxel near a surface vertex
  L2SaddPoint(l2s, c, r, s, 0, 0);
  err = 0;
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (lv->vno == vno && !lv->deleted) {
      err = 1;
      break;
    }
  }
  if (err) {
    printf("L2Stest: Failed to add correct surf vertex (%d)\n", lv->vno);
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  l2s->labels[0]->n_points = 0;  // reset label

  // Add a voxel away from the surface
  c = 0;  // 0 0 0 should be guaranteed to not be near surf
  r = 0;
  s = 0;
  printf("L2Stest: crs = (%d, %d, %d)\n", c, r, s);

  L2SaddPoint(l2s, c, r, s, 0, 1);
  if (l2s->labels[0]->n_points != 1) {
    printf("L2Stest: Failed to add non-surf voxel\n");
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (lv->vno != -1) {
      printf("L2Stest: Added non-surf voxel, but interpreted as a surf-based\n");
      printf("break %s:%d\n", __FILE__, __LINE__);
      return (1);
    }
  }
  // Remove voxel away from the surface
  L2SaddPoint(l2s, c, r, s, 0, 0);
  err = 0;
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (!lv->deleted) {
      err = 1;
      break;
    }
  }
  if (err) {
    printf("L2Stest: Failed to remove non-surf-based voxel\n");
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  l2s->labels[0]->n_points = 0;  // reset label

  // Import a label
  label = annotation2label(1, surf);  // banks of the stss
  printf("L2Stest: importing label with np = %d\n", label->n_points);
  L2SimportLabel(l2s, label, 0);
  if (label->n_points != l2s->labels[0]->n_points) {
    printf("L2Stest: Failed to import label np = %d %d\n", label->n_points, l2s->labels[0]->n_points);
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  err = 0;
  for (n = 0; n < l2s->labels[0]->n_points; n++) {
    lv = &(l2s->labels[0]->lv[n]);
    if (lv->vno != label->lv[n].vno) {
      err = 1;
      break;
    }
  }
  if (err) {
    printf("L2Stest: Failed to correctly import label\n");
    printf("break %s:%d\n", __FILE__, __LINE__);
    return (1);
  }
  l2s->labels[0]->n_points = 0;  // reset label

  L2Sfree(&l2s);

  printf("L2Stest: passed tests, subject = %s\n", subject);
  return (0);
}

/*!
  \fn MRIS *MakeAverageSurf(AVERAGE_SURFACE_PARAMS *asp)
  \brief Computes an average surface given parameters in asp.  This
  routine differs from the current (7/23/18) version
  mris_make_average_surface in that it maps directly to the nearest
  vertex of the icosahedron rather than going through the
  parameterized surface (MRISP). This prevents the really big vertex
  at the poles. But otherwise, the surface appears to be near
  idential, though the actual placement of the vertices might not be
  the same.
  \example
    asp = MRISaverageSurfaceParamAlloc(40);
    asp->icoorder = 7;
    asp->hemi = "lh";
    asp->surfname = "white";
    asp->surfregname = "sphere.reg";
    asp->xform_name = "talairach.xfm";
    asp->subectlist[0] = "subject0";
    asp->subectlist[1] = "subject1";
    ...
 */
MRIS *MakeAverageSurf(AVERAGE_SURFACE_PARAMS *asp)
{
  MRIS *surf, *targsurfreg=NULL, *surfreg, *stpair[2];
  int nth;
  char *subject;
  FSENV *fsenv;
  char tmpstr[2000];
  MATRIX *XFM;
  GCA_MORPH *gcam=NULL;
  MRI *SrcXYZ, *TargXYZ, *TargXYZSum=NULL;
  double average_surface_area;

  fsenv = FSENVgetenv();

  printf("MakeAverageSurf(): %d %s %s %s %s %d %d\n",asp->nsubjects,asp->hemi,
	 asp->surfname,asp->surfregname,asp->xform_name,asp->ReverseMapFlag,
	 asp->UseHash);

  if(asp->targsurfreg != NULL){
    // If the actual target surface is specified, then use that (not tested yet)
    targsurfreg = MRISclone(asp->targsurfreg);
    printf("Using passed target surface\n");
  }
  else if(asp->icoorder > 0){
    // Otherwise, use the icosahedron order
    sprintf(tmpstr,"%s/lib/bem/ic%d.tri",fsenv->FREESURFER_HOME,asp->icoorder);
    printf("Loading %s as target surface\n",tmpstr);
    //targsurfreg = MRISread(tmpstr);
    targsurfreg = ReadIcoByOrder(asp->icoorder, 100.0);
    if(targsurfreg==NULL) return(NULL);
  }
  else if(asp->targsubject){
    // Otherwise, use the target surbject
    sprintf(tmpstr,"%s/%s/surf/%s.%s",fsenv->SUBJECTS_DIR,asp->targsubject,asp->hemi,asp->surfregname);
    printf("Loading %s as target surface\n",tmpstr);
    targsurfreg = MRISread(tmpstr);
    if(targsurfreg==NULL) return(NULL);
  }
  if(targsurfreg == NULL){
    printf("ERROR: no target specified\n");
    return(NULL);
  }

  average_surface_area = 0;
  for(nth=0; nth < asp->nsubjects; nth++){
    subject = asp->subjectlist[nth];
    printf("  %d/%d %s ---------------\n",nth+1,asp->nsubjects,subject);

    sprintf(tmpstr,"%s/%s/surf/%s.%s",fsenv->SUBJECTS_DIR,subject,asp->hemi,asp->surfregname);
    surfreg = MRISread(tmpstr);
    if(surfreg==NULL) return(NULL);

    sprintf(tmpstr,"%s/%s/surf/%s.%s",fsenv->SUBJECTS_DIR,subject,asp->hemi,asp->surfname);
    surf = MRISread(tmpstr);
    if(surf==NULL) return(NULL);

    if(asp->xform_name){
      char *ext = fio_extension(asp->xform_name);
      sprintf(tmpstr, "%s/%s/mri/transforms/%s", fsenv->SUBJECTS_DIR, subject,asp->xform_name);

      if(!strcmp(asp->xform_name,"talairach.xfm")) {
	printf("  Applying linear transform %s\n",asp->xform_name);
	XFM = DevolveXFMWithSubjectsDir(subject, NULL, asp->xform_name, fsenv->SUBJECTS_DIR);
	if(XFM == NULL) return(NULL);
	MRISmatrixMultiply(surf, XFM);
	MatrixFree(&XFM);
      } 
      else if(!strcmp(ext,"m3z") || !strcmp(ext,"mgz") || !strcmp(ext,"mgh")){
	printf("  Applying GCA Morph %s\n",asp->xform_name);
	gcam = GCAMread(tmpstr);
	if(gcam == NULL) return(NULL);
	GCAMmorphSurf(surf, gcam);
	GCAMfree(&gcam);
      } 
      else if(!strcmp(ext,"lta")){
	printf("  Applying linear transform %s\n",tmpstr);
	LTA *lta = LTAread(tmpstr);
	if(lta==NULL) exit(1);
	int err = MRISltaMultiply(surf, lta);
	if(err) return(NULL);
	LTAfree(&lta);
      }
      else {
	printf("ERROR: don't know what to do with %s\n",tmpstr);
	return(NULL);
      }
    }
    if(nth==0) copyVolGeom(&surf->vg,&targsurfreg->vg);

    // Could remove intersections here after applying xform

    // Copy the surface XYZ into an MRI structure
    SrcXYZ = MRIcopyMRIS(NULL, surf, 2, "z"); // start at z to autoalloc
    MRIcopyMRIS(SrcXYZ, surf, 0, "x");
    MRIcopyMRIS(SrcXYZ, surf, 1, "y");

    // Apply the surface registration to the XYZ
    stpair[0] = surfreg;
    stpair[1] = targsurfreg;
    TargXYZ = MRISapplyReg(SrcXYZ, stpair, 2, asp->ReverseMapFlag, 0, asp->UseHash);
    if(TargXYZ == NULL) return(NULL);

    // Accumulate
    if(nth == 0) TargXYZSum = MRIcopy(TargXYZ,TargXYZSum);
    else         TargXYZSum = MRIadd(TargXYZSum,TargXYZ,TargXYZSum);

    // keep track of the total area
    average_surface_area += surf->total_area ;

    MRIfree(&SrcXYZ);
    MRIfree(&TargXYZ);
    MRISfree(&surf);
    MRISfree(&surfreg);
  }

  // Compute average and copy back into the target surface XYZ
  TargXYZSum = MRImultiplyConst(TargXYZSum,1.0/asp->nsubjects,TargXYZSum);
  MRIScopyMRI(targsurfreg,TargXYZSum,0,"x");
  MRIScopyMRI(targsurfreg,TargXYZSum,1,"y");
  MRIScopyMRI(targsurfreg,TargXYZSum,2,"z");

  LTA *DestLTA = NULL;
  if(asp->DestLTA) DestLTA = asp->DestLTA;
  if(asp->Conform && !mriConformed(&targsurfreg->vg)){
    VOL_GEOM vgdest;
    copyVolGeom(&targsurfreg->vg,&vgdest);
    conformGeom(&vgdest,0,1,0);
    DestLTA = TransformRegDat2LTA(&vgdest,&targsurfreg->vg,NULL);
    LTAprint(stdout, DestLTA);
  }
  if(DestLTA){
    // This can be good when the xform does not take it all the way to the
    // target (eg, conformed) space.
    printf("Applying LTA to output (Conform=%d)\n",asp->Conform);
    int err = MRISltaMultiply(targsurfreg, DestLTA);
    if(err) return(NULL);
  }
  if(asp->Conform && !mriConformed(&targsurfreg->vg)) LTAfree(&DestLTA);

  average_surface_area /= (double)asp->nsubjects ;
  MRIScomputeMetricProperties(targsurfreg);
  printf("setting group surface area to be %2.1f cm^2 (scale=%2.2f)\n",
	 average_surface_area/100.0,sqrt(average_surface_area/targsurfreg->total_area)) ;
  targsurfreg->group_avg_surface_area = average_surface_area ;
  MRIScomputeMetricProperties(targsurfreg);

  FSENVfree(&fsenv);
  MRIfree(&TargXYZSum);

  printf("MakeAverageSurf(): done\n");
  return(targsurfreg);
}

/*!
  \fn AVERAGE_SURFACE_PARAMS *MRISaverageSurfaceParamAlloc(int nsubjects)
  \brief Allocates the subjectlist char string in the PARAMS
  struct. Also sets the icoorder to -1, ReverseMapFlag = 1, and
  UseHash = 1.
 */
AVERAGE_SURFACE_PARAMS *MRISaverageSurfaceParamAlloc(int nsubjects)
{
  AVERAGE_SURFACE_PARAMS *asp;
  asp = (AVERAGE_SURFACE_PARAMS *)calloc(1, sizeof(AVERAGE_SURFACE_PARAMS));
  asp->nsubjects = nsubjects;
  asp->subjectlist = (char **) calloc(sizeof(char*),nsubjects);
  asp->icoorder = -1;
  asp->ReverseMapFlag = 1;
  asp->UseHash = 1;
  return(asp);
}

/*!
  \fn int MRISaverageSurfaceParamFree(AVERAGE_SURFACE_PARAMS **pasp)
  \brief Frees the subjectlist and the ASP prointer.
 */
int MRISaverageSurfaceParamFree(AVERAGE_SURFACE_PARAMS **pasp)
{
  AVERAGE_SURFACE_PARAMS *asp = *pasp;
  free(asp->subjectlist);
  free(*pasp);
  *pasp = NULL;
  return(0);
}

/*!
  \fn int MRISeulerNoSeg(MRI_SURFACE *mris, MRI *surfseg, int segno, int *pnvertices, int *pnfaces, int *pnedges, int *pv0)
  \brief Computes the euler number for the set of vertices in which
  surfseg==segno. If *pnvertices, *pnfaces, *pnedges, *pv0 are
  non-null, then returns those values. Alters vertex->marked. v0 is a
  vertex in the seg.  Only considers vertices if they belong to a
  triangle where all the corners are in the segmentation.
*/
int MRISeulerNoSeg(MRI_SURFACE *mris, MRI *surfseg, int segno, int *pnvertices, int *pnfaces, int *pnedges, int *pv0)
{
  int eno, nfaces, nedges, nvertices, vno, fno, vnb, i, dno;
  int n, nhits,v0;

  // Only consider vertices if they belong to a triangle where all the corners
  // are in the segmentation. 
  nfaces = 0;
  for(fno = 0; fno < mris->nfaces; fno++){
    if(mris->faces[fno].ripflag) continue;
    nhits = 0;
    for(n=0; n < 3; n++){
      vno = mris->faces[fno].v[n];
      if(mris->vertices[vno].ripflag) continue;
      if(MRIgetVoxVal(surfseg,vno,0,0,0) != segno) break;
      nhits++;
    }
    if(nhits != 3) continue;
    for(n=0; n < 3; n++){
      vno = mris->faces[fno].v[n];
      mris->vertices[vno].marked = 1;
    }
    nfaces++;
  }

  // Count the number of vertices and get a vertex in the label
  v0 = -1;
  nvertices = 0;
  for(vno = 0; vno < mris->nvertices; vno++){
    if(!mris->vertices[vno].marked) continue;
    if(v0<0) v0 = vno;
    nvertices++;
  }

  // Count up the edges. This cannot be parallized.
  nedges = 0;
  for (vno = 0; vno < mris->nvertices; vno++){
    if(!mris->vertices[vno].marked) continue;
    VERTEX_TOPOLOGY const *v1 = &mris->vertices_topology[vno];
    for (i = 0; i < v1->vnum; i++) {
      vnb = v1->v[i];
      if(! mris->vertices[vnb].marked) continue;
      if(vnb <= vno) continue; /* already counted */
      nedges++;
    }
  }
  
  eno = nvertices - nedges + nfaces; // euler number
  dno = abs(2 - eno) + abs(2 * nedges - 3 * nfaces);
  //if(nvertices>0 && eno == 1)
  //printf("label=%d v0=%d nv=%d, ne=%d nf=%d eno=%d dno=%d\n",segno,v0,nvertices,nedges,nfaces,eno,dno);
  //printf("%3d %4d %5d %5d %5d %5d\n",segno,nvertices,nedges,nfaces,eno,dno);

  // Clear the marks
  for(vno = 0; vno < mris->nvertices; vno++)
    mris->vertices[vno].marked = 0;

  if(pnvertices != NULL) *pnvertices = nvertices;
  if(pnfaces != NULL) *pnfaces = nfaces;
  if(pnedges != NULL) *pnedges = nedges;
  if(pv0 != NULL)     *pv0 = v0;

  return(eno);
}


/*!
  \fn double *MRIStriangleAreaStats(MRIS *surf, double *stats)
  \brief Computes stats (nfaces, mean, stddev, min, max) over
  all faces that do not have a ripped vertex. If mask is non-null,
  then the mask value of a vertex must be greater than 0.5 to
  be included in the list. Runs MRIScomputeMetricProperties().
 */
double *MRIStriangleAreaStats(MRIS *surf, MRI *mask, double *stats)
{
  int fno, nfaces;
  double *area;
  MRIScomputeMetricProperties(surf);

  area = (double*)calloc(sizeof(double),surf->nfaces);
  nfaces = 0;
  for(fno=0; fno < surf->nfaces; fno++){
    FACE *f = &(surf->faces[fno]);
    int nthv, skip;
    skip = 0;
    for(nthv = 0; nthv < 3; nthv++){
      int vno = f->v[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      if(v->ripflag) skip = 1;
      if(mask && MRIgetVoxVal(mask,vno,0,0,0) < 0.5) skip = 1;
    }
    if(skip) continue;
    area[nfaces] = f->area;
    nfaces ++;
    //printf("%g\n",f->area);
  }
  stats = DListStats(area, nfaces, stats);
  free(area);

  return(stats);
}
/*!
  \fn double *MRISedgeStats(MRIS *surf, int metricid, MRI *mask, double *stats)
  \brief Computes stats (nedges, mean, stddev, min, max) over
  all edges that do not have a ripped vertex. If mask is non-null,
  then the mask value of a vertex must be greater than 0.5 to
  be included in the list. metricid: 0=length, 1=dot, 2=angle.
  Will create the edge structure if not already there. Runs
  MRIScomputeMetricProperties() and MRISedgeMetric(surf).
 */
double *MRISedgeStats(MRIS *surf, int metricid, MRI *mask, double *stats)
{
  int edgeno, nedges, nthv;
  MRI_EDGE *e;
  double *metric;

  if(surf->edges == NULL){
    MRISedges(surf);
  }
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,0);

  metric = (double*)calloc(sizeof(double),surf->nedges);
  nedges = 0;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    e = &(surf->edges[edgeno]);
    int skip = 0;
    for(nthv=0; nthv < 4; nthv++){
      int vno = e->vtxno[nthv];
      VERTEX  * const v = &(surf->vertices[vno]);
      if(v->ripflag) skip = 1;
      if(mask && MRIgetVoxVal(mask,vno,0,0,0) < 0.5) skip = 1;
    }
    if(skip) continue;
    switch(metricid){
    case 0: metric[nedges] = e->len; break;
    case 1: metric[nedges] = e->dot; break;
    case 2: metric[nedges] = e->angle; break;
    default:
      printf("ERROR: MRISedgeStats() metricid %d unrecognized\n",metricid);
      return(NULL);
    }
    //printf("%lf\n",metric[nedges]);
    nedges++;
  }

  stats = DListStats(metric, nedges, stats);
  free(metric);
  return(stats);
}
/*!
  \fn int MRISedgePrint(FILE *fp, MRIS *surf)
  \brief Print edge metrics to a stream. Runs
  MRIScomputeMetricProperties() and MRISedgeMetric()
 */
int MRISedgePrint(FILE *fp, MRIS *surf)
{
  int edgeno;
  MRI_EDGE *e;
  MRIScomputeMetricProperties(surf);
  MRISedgeMetric(surf,0);
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    e = &(surf->edges[edgeno]);
    double cost = (1.0-e->dot)*(1.0-e->dot);
    fprintf(fp,"%6d %6d %6d %10.8f %10.8f %8.4f %12.8f\n",edgeno,
	    e->vtxno[0],e->vtxno[1],e->len,e->dot,e->angle,cost);
  }
  fflush(fp);
  return(0);
}

/*!
  \fn int MRISedgeWrite(char *filename, MRIS *surf)
  \brief Print edge metrics to a file. Runs
  MRIScomputeMetricProperties() and MRISedgeMetric()
  via MRISedgePrint(). 
 */
int MRISedgeWrite(char *filename, MRIS *surf)
{
  FILE *fp;
  if(surf->edges == NULL){
    MRISedges(surf);
  }
  MRIScomputeMetricProperties(surf) ;
  MRISedges(surf);
  MRISfaceNormalGrad(surf, 0); //0=DoGrad
  MRISedgeMetric(surf,0);

  fp = fopen(filename,"w");
  if(fp == NULL) return(1);
  MRISedgePrint(fp,surf);
  return(0);
}


/*!
\fn MRIS *MRISupsampleCentroid(MRIS *insurf, int nupsamples)
\brief Upsample the surface by putting new vertices at the centroid of
each triangle. For one upsample, the new surface will have nvertices =
insurf->nvertices + insurf->nfaces and nfaces = 3*insurf->nfaces. Only
the minimal set of parameters are set in the new surface so functions
like MRIScomputeMetricProperties() might not do the right thing if run
directly on the output (though it is run at the end of the
function). It would be nice to have a function that updates all the
parameters given the minimum, but this does not appear to exist. But
you can save it to file and read it back. This function tends to make
long triangles (bad quality) with vertices non-uniformly distributed.
Probably better to use MRISupsampleSplit().
*/
MRIS *MRISupsampleCentroid(MRIS *srcsurf, int nupsamples)
{
  if(srcsurf == NULL) return(NULL);
  MRIS *upsurf;

  if(nupsamples > 1){
    upsurf=srcsurf;
    for(int n=0; n < nupsamples; n++){
      printf("MRISupsampleCentroid() upsample level = %d ======================\n",n);
      MRIS *tmpsurf = MRISupsampleCentroid(upsurf, 1);
      if(tmpsurf == NULL) {
	printf("ERROR: MRISupsampleCentroid() return is NULL\n");
	return(NULL);
      }
      if(n>0) MRISfree(&upsurf);
      upsurf = tmpsurf;
    }
    return(upsurf);
  }

  int nvertices = srcsurf->nvertices + srcsurf->nfaces;
  int nfaces = 3*srcsurf->nfaces;
  upsurf = MRISclone(srcsurf);
  if(upsurf == NULL) return(NULL);
  MRISreallocVerticesAndFaces(upsurf,nvertices, nfaces);
  printf("MRISupsampleCentroid(): source nvertices = %d nfaces = %d, new nv = %d nf = %d\n",
	 srcsurf->nvertices,srcsurf->nfaces,upsurf->nvertices,upsurf->nfaces);fflush(stdout);

  int nthface, newfaceno=srcsurf->nfaces;
  for(nthface=0; nthface < srcsurf->nfaces; nthface++){
    FACE *face = &(srcsurf->faces[nthface]); // can be source or up

    // Put a new vertex at the centroid of this face
    int upvtxno = srcsurf->nvertices+nthface; // index of new vertex
    VERTEX *vup = &(upsurf->vertices[upvtxno]); // must be upsurf

    // Compute the centroid
    int k,nbrvtxno;
    double xsum=0, ysum=0, zsum=0;
    for(k=0; k < 3; k++){
      nbrvtxno = face->v[k];
      VERTEX *vk = &(srcsurf->vertices[nbrvtxno]); // up or source here
      xsum += vk->x;
      ysum += vk->y;
      zsum += vk->z;
    }
    vup->x = xsum/3.0;
    vup->y = ysum/3.0;
    vup->z = zsum/3.0;

    // Create new faces for the three new triangles. Note that the
    // order of the vertices in the face is important as it determines
    // whether the norm points in or out
    FACE *newface;
    // For the first face, replace the values in the already existing face
    newface = &(upsurf->faces[nthface]); // nthface here
    newface->v[0] = face->v[0];
    newface->v[1] = face->v[1];
    newface->v[2] = upvtxno;

    newface = &(upsurf->faces[newfaceno]); // newfaceno here
    newface->v[0] = face->v[1];
    newface->v[1] = face->v[2];
    newface->v[2] = upvtxno;
    newfaceno++;

    newface = &(upsurf->faces[newfaceno]); // newfaceno here
    newface->v[0] = face->v[2];
    newface->v[1] = face->v[0];
    newface->v[2] = upvtxno;
    newfaceno++;

  }
  MRIScomputeMetricProperties(upsurf); // not sure if this actually works
  return(upsurf);
}

/*!
\fn MRIS *MRISupsampleSplit(MRIS *srcsurf, int nupsamples, int SortType)
\brief Upsample the surface by adding a vertex at the middle of an edge.
This splits the adjacent triangles. Not every edge will be split, because,
once a triangle is split, none of the adjacent edges can be used for spliting. 
For this reason, the order of splitting matters. The order can be controlled
with SortType where the edges are ordered by length and/or the size
of an adjoining triangle:
  case 0: Not sorting
  case 1: Length then Max Area, Descending (high to low)
  case 2: Length then Max Area, Ascending  (low to high)
  case 3: Max Area then Length, Descending (high to low)
  case 4: Max Area then Length, Ascending  (low to high)
Eg, use SortType=1 so split the largest edges first.

For one upsample, the new surface will have a maximum nvertices =
srcsurf->nvertices + srcsurf->nedges and nfaces = 2*srcsurf->nfaces.
The function does not set all the parameters in MRI_SURFACE structure,
so it writes the surface to a temporary file and then reads it back
in. Yes, silly. The temp file will have a name like
MRISupsample.tmpsurf.XXXXXX and will be in FREESURFER_TMP_DIR is set
or /tmp otherwise.
*/
MRIS *MRISupsampleSplit(MRIS *srcsurf, int nupsamples, int SortType)
{
  if(srcsurf == NULL) return(NULL);
  MRIS *upsurf;

  MRISedges(srcsurf);
  MRISfaceMetric(srcsurf,0);
  MRISedgeMetric(srcsurf,0);

  MRI_EDGE *edges;
  if(SortType == 0) edges = srcsurf->edges;
  else              edges = MRISedgeSort(srcsurf, SortType, NULL);

  if(nupsamples > 1){
    upsurf=srcsurf;
    for(int n=0; n < nupsamples; n++){
      printf("MRISupsampleSplit() upsample level = %d SortType=%d ======================\n",n,SortType);
      MRIS *tmpsurf = MRISupsampleSplit(upsurf, 1, SortType);
      if(tmpsurf == NULL) {
	printf("ERROR: MRISupsampleSplit() return is NULL\n");
	return(NULL);
      }
      if(n>0) MRISfree(&upsurf);
      upsurf = tmpsurf;
    }
    return(upsurf);
  }

  int nvertices = srcsurf->nvertices + srcsurf->nedges; // max number
  int nfaces = 2*srcsurf->nfaces; // max
  upsurf = MRISclone(srcsurf);
  if(upsurf == NULL) return(NULL);
  MRISreallocVerticesAndFaces(upsurf,nvertices, nfaces);
  printf("MRISupsampleSplit(): source nvertices = %d nfaces = %d, nedges=%d new nv = %d nf = %d\n",
	 srcsurf->nvertices,srcsurf->nfaces,srcsurf->nedges,upsurf->nvertices,upsurf->nfaces);fflush(stdout);

  int *FaceNotSplitable = (int*) calloc(srcsurf->nfaces,sizeof(int));
  int nthedge, upvtxno = srcsurf->nvertices, upfaceno=srcsurf->nfaces, db=0;
  for(nthedge = 0; nthedge < srcsurf->nedges; nthedge++){
    MRI_EDGE *e = &(edges[nthedge]);
    if(db) printf("nthedge = %d eno=%d len=%6.2f maxarea=%6.1f v=(%d,%d) f=(%d,%d) upv=%d upf=%d=========================\n",
		  nthedge,e->edgeno,e->len,e->maxarea,e->vtxno[0],e->vtxno[1],e->faceno[0],e->faceno[1],upvtxno,upfaceno);
    if(FaceNotSplitable[e->faceno[0]] || FaceNotSplitable[e->faceno[1]]) continue;

    VERTEX *vsrc0 = &(srcsurf->vertices[e->vtxno[0]]);
    VERTEX *vsrc1 = &(srcsurf->vertices[e->vtxno[1]]);
    VERTEX *vup = &(upsurf->vertices[upvtxno]);
    vup->x = (vsrc0->x+vsrc1->x)/2.0;
    vup->y = (vsrc0->y+vsrc1->y)/2.0;
    vup->z = (vsrc0->z+vsrc1->z)/2.0;

    // Create new faces for the two new triangles. Note that the order
    // of the vertices in the face is important as it determines
    // whether the norm points in or out. For the first face, replace
    // the values in the already existing face.
    for(int k=0; k<2; k++){
      FACE *face = &(srcsurf->faces[e->faceno[k]]);
      int splitedgeno;
      for(splitedgeno=0; splitedgeno < 3; splitedgeno++){
	int m = splitedgeno + 1;
	if(m==3) m=0;
	if( (face->v[splitedgeno] == e->vtxno[0] && face->v[m] == e->vtxno[1]) ||
	    (face->v[splitedgeno] == e->vtxno[1] && face->v[m] == e->vtxno[0]) ){
	  break;
	}
      }
      if(db)printf("  Spliting %dth faceno %d along edge %d\n",k,e->faceno[k],splitedgeno);
      FACE *newface;
      if(splitedgeno==0){// spliting 0-1
	newface = &(upsurf->faces[e->faceno[k]]); // e->faceno[k] here
	newface->v[0] = face->v[0];
	newface->v[1] = upvtxno;
	newface->v[2] = face->v[2];
	if(db){printf("    edge %d  (%d,%d,%d)\n",nthedge,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
	newface = &(upsurf->faces[upfaceno]); // upfaceno here
	newface->v[0] = upvtxno;
	newface->v[1] = face->v[1];
	newface->v[2] = face->v[2];
	if(db){printf("    newface %d  (%d,%d,%d)\n",upfaceno,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
      }
      if(splitedgeno==1){// splitting 1-2
	newface = &(upsurf->faces[e->faceno[k]]); // e->faceno[k] here
	newface->v[0] = face->v[0];
	newface->v[1] = face->v[1];
	newface->v[2] = upvtxno;
	if(db){printf("    edge %d  (%d,%d,%d)\n",nthedge,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
	newface = &(upsurf->faces[upfaceno]); // upfaceno here
	newface->v[0] = upvtxno;
	newface->v[1] = face->v[2];
	newface->v[2] = face->v[0];
	if(db){printf("    newface %d  (%d,%d,%d)\n",upfaceno,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
      }
      if(splitedgeno==2){// splitting 2-0
	newface = &(upsurf->faces[e->faceno[k]]); // e->faceno[k] here
	newface->v[0] = face->v[0];
	newface->v[1] = face->v[1];
	newface->v[2] = upvtxno;
	if(db){printf("    edge %d  (%d,%d,%d)\n",nthedge,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
	newface = &(upsurf->faces[upfaceno]); // upfaceno here
	newface->v[0] = face->v[1];
	newface->v[1] = face->v[2];
	newface->v[2] = upvtxno;
	if(db){printf("    newface %d  (%d,%d,%d)\n",upfaceno,newface->v[0],newface->v[1],newface->v[2]); fflush(stdout);}
      }
      upfaceno++;
      FaceNotSplitable[e->faceno[k]]=1;
    }// kth face in edge
    upvtxno++; // only one vertex added per edge
  }// nthedge
  upsurf->nvertices = upvtxno;
  upsurf->nfaces = upfaceno;
  printf("new surf nv = %d nf = %d\n",upsurf->nvertices,upsurf->nfaces);
  free(FaceNotSplitable);
  if(SortType != 0) free(edges);

  // OK, this is a little crazy, but save it and read it back in so that
  // everything that should be set gets set
  char tmpfile[2000];
  char *tmpdir = getenv("FREESURFER_TMP_DIR");
  if(tmpdir==NULL) sprintf(tmpfile,"/tmp/MRISupsample.tmpsurf.XXXXXX");
  else             sprintf(tmpfile,  "%s/MRISupsample.tmpsurf.XXXXXX",tmpdir);
  mkstemp(tmpfile);
  printf("MRISupsampleSplit(): writing/reading temporary surf to %s\n",tmpfile);
  int err = 0;
  err = MRISwrite(upsurf,tmpfile);
  if(err) {
    printf("MRISupsampleSplit(): could not save temp surf to %s\n",tmpfile);
    return(NULL);
  }
  MRISfree(&upsurf);
  upsurf = MRISread(tmpfile);
  if(upsurf==NULL){
    printf("MRISupsampleSplit(): could not read temp surf from %s\n",tmpfile);
    return(NULL);
  }
  unlink(tmpfile);

  return(upsurf);
}

/* This function only supports mris_sphere options -Q, -P, -A, -seed.
 * mris_sphere calls in recon-all:
 * mris_sphere -q -p 6 -a 128 -seed 1234 ../surf/lh.inflated.nofix ../surf/lh.qsphere.nofix 
 * mris_sphere -q -p 6 -a 128 -seed 1234 ../surf/rh.inflated.nofix ../surf/rh.qsphere.nofix
 *
 * These static variables are defined in mris_sphere.cpp:
 *
 * #define BASE_DT_SCALE     1.0
 * static float base_dt_scale = BASE_DT_SCALE ;
 * static int nbrs = 2 ;
 * static int do_inflate = 0 ;
 * static double disturb = 0 ;
 * static int   max_passes = 1 ;
 * static int   randomly_project = 0 ;
 * static int   talairach = 0 ;
 * static float scale = 1.0 ;
 * static int quick = 0 ;
 * static int load = 0 ;
 * static float inflate_area  = 0.0f ;
 * static float inflate_dt = 0.9 ;
 * static float inflate_tol  = 1.0f ;
 * static float inflate_nlarea  = 0.0f ;
 * static float inflate_spring  = 0.0f ;
 * static float inflate_tspring = 0 ;
 * 
 * static float ralpha = 0.0 ;
 * static float rbeta = 0.0 ;
 * static float rgamma = 0.0 ;
 * 
 * static double target_radius = DEFAULT_RADIUS ;
 * #if 0
 * static int   inflate_avgs = 64 ;
 * static int   inflate_iterations = 50 ;
 * static float l_convex = 1.0 ;
 * static float l_spring_norm = .1 ;
 * static float l_sphere = 0.25 ;
 * #else
 * static int   inflate_avgs = 0 ;
 * static int   inflate_iterations = 1000 ;
 * static float l_convex = 1.0 ;
 * static float l_expand = 0.0 ;
 * static float l_spring_norm = 1.0 ;
 * static float l_sphere = 0.025 ;
 * #endif
 * 
 * static const char *orig_name = "smoothwm" ;
 * static int smooth_avgs = 0 ;
 * 
 * static char *xform_fname = NULL ;
 * static char *vol_fname = NULL ;
 * 
 * static int remove_negative = 1 ;
 * char *rusage_file=NULL;
 */
MRIS *MRISQuickSphericalInflate(int max_passes, int n_averages, long seed, MRIS *inSurf, const char *outSurf)
{
  double target_radius = DEFAULT_RADIUS ;  // include/mrisurf.h:#define DEFAULT_RADIUS  100.0f
  float base_dt_scale = 1.0;
  float scale = 1.0;

  float inflate_area  = 0.0f ;
  float inflate_dt = 0.9 ;
  float inflate_tol  = 1.0f ;
  float inflate_nlarea  = 0.0f ;
  float inflate_spring  = 0.0f ;
  float inflate_tspring = 0 ;

  int   inflate_avgs = 0 ;
  //int   inflate_iterations = 1000 ;
  float l_convex = 1.0 ;
  float l_expand = 0.0 ;
  float l_spring_norm = 1.0 ;
  float l_sphere = 0.025 ;

  // these are the settings for -Q
  //int remove_negative = 0;
  //int quick = 1;
  //int do_inflate = 1;
  int inflate_iterations = 300;
  //int max_passes = 3;
  int nbrs = 1;

  Timer then ;
  INTEGRATION_PARMS  parms;

  then.reset();

  // default settings
  parms.dt = .05 ;
  parms.projection = PROJECT_ELLIPSOID ;
  parms.tol = .5 /*1e-1*/ ;
  parms.n_averages = 1024 ;
  parms.min_averages = 0 ;
  parms.l_angle = 0.0 /* L_ANGLE */ ;
  parms.l_area = 0.0 /* L_AREA */ ;
  parms.l_neg = 0.0 ;
  parms.l_dist = 1.0 ;
  parms.l_spring = 0.0 ;
  parms.l_area = 1.0 ;
  parms.l_boundary = 0.0 ;
  parms.l_curv = 0.0 ;
  parms.niterations = 25 ;
  parms.write_iterations = 1000 ;
  parms.a = parms.b = parms.c = 0.0f ;  /* ellipsoid parameters */
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.99 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.momentum = 0.9 ;
  parms.desired_rms_height = -1.0 ;
  parms.base_name[0] = 0 ;
  parms.Hdesired = 0.0 ;   /* a flat surface */
  parms.nbhd_size = 7 ;
  parms.max_nbrs = 8 ;

  // setting for quick spherical unfolding: mris_sphere -Q option
  fprintf(stderr, "doing quick spherical unfolding.\n");
  parms.l_spring = parms.l_dist = parms.l_parea = parms.l_area = 0.0 ;
  parms.l_nlarea = 1.0 ;
  parms.tol = 1e-1 ;
  parms.n_averages = 32 ;

  // mris_sphere -P 
  fprintf(stderr, "limitting unfolding to %d passes\n", max_passes);

  // mris_sphere -A
  fprintf(stderr, "using n_averages = %d\n", n_averages) ;
  parms.n_averages = n_averages;

  // mris_sphere -seed
  setRandomSeed(seed) ;
  fprintf(stderr,"setting seed for random number genererator to %ld\n", seed) ;

  parms.scale = scale ;
  parms.base_dt = base_dt_scale * parms.dt ;

  std::cout << "version: " << getVersion() << std::endl;

#ifdef HAVE_OPENMP
  std::cout << "available threads: " << omp_get_max_threads() << std::endl;
#endif

  if (parms.base_name[0] == 0 && outSurf != NULL)
  {
    char fname[STRLEN];
    FileNameOnly(outSurf, fname) ;
    char *cp = strchr(fname, '.') ;
    if (cp)
    {
      strcpy(parms.base_name, cp+1) ;
    }
    else
    {
      strcpy(parms.base_name, "sphere") ;
    }
  }

  float max_dim = MAX(abs(inSurf->xhi-inSurf->xlo), abs(inSurf->yhi-inSurf->ylo)) ;
  max_dim = MAX(max_dim,abs(inSurf->zhi-inSurf->zlo)) ;

  if (max_dim > .75*DEFAULT_RADIUS)
  {
    float ratio = .75*DEFAULT_RADIUS / (max_dim) ;
    printf("scaling brain by %2.3f...\n", ratio) ;
    MRISscaleBrain(inSurf, inSurf, ratio) ;
  }

    fflush(stdout); fflush(stderr);
    fprintf(stderr, "inflating...\n");
    fflush(stdout); fflush(stderr);

  // if (!load && do_inflate) 
  {
    INTEGRATION_PARMS inflation_parms ;

    MRIScenter(inSurf, inSurf) ;
    strcpy(inflation_parms.base_name, parms.base_name) ;
    inflation_parms.write_iterations = parms.write_iterations ;
    inflation_parms.niterations = inflate_iterations ;
    inflation_parms.l_spring_norm = l_spring_norm ;
    inflation_parms.l_spring = inflate_spring ;
    inflation_parms.l_nlarea = inflate_nlarea ;
    inflation_parms.l_area = inflate_area ;
    inflation_parms.n_averages = inflate_avgs ;
    inflation_parms.l_expand = l_expand ;
    inflation_parms.l_tspring = inflate_tspring ;
    inflation_parms.l_sphere = l_sphere ;
    inflation_parms.l_convex = l_convex ;
#define SCALE_UP 2
    inflation_parms.a = SCALE_UP*DEFAULT_RADIUS ;
    inflation_parms.tol = inflate_tol ;
    inflation_parms.integration_type = INTEGRATE_MOMENTUM ;
    inflation_parms.momentum = 0.9 ;
    inflation_parms.dt = inflate_dt ;

    /* store the inflated positions in the v->c? field so that they can
      be used in the repulsive term.
    */
    /*    inflation_parms.l_repulse_ratio = .1 ;*/
    MRISsaveVertexPositions(inSurf, CANONICAL_VERTICES) ;
    if (l_expand > 0)
    {
      MRISexpandSurface(inSurf, target_radius/2, &inflation_parms, 0, 1) ;
      l_expand = parms.l_expand = 0 ;
    }
    MRIScenter(inSurf, inSurf) ;
    inSurf->x0 = inSurf->xctr ;
    inSurf->y0 = inSurf->yctr ;
    inSurf->z0 = inSurf->zctr ;
    MRISinflateToSphere(inSurf, &inflation_parms) ;
    if (inflation_parms.l_expand > 0)
    {
      inflation_parms.l_expand = 0 ;
      inflation_parms.niterations += (inflate_iterations*.1) ;
      MRISinflateToSphere(inSurf, &inflation_parms) ;
    }
    MRISscaleBrain(inSurf, inSurf, target_radius/(DEFAULT_RADIUS*SCALE_UP)) ;
    parms.start_t = inflation_parms.start_t ;
    MRISresetNeighborhoodSize(inSurf, nbrs) ;
  }

  fs::debug() << "should but doesn't set orig xyz here";

  fprintf(stderr, "projecting onto sphere...\n");
  MRISprojectOntoSphere(inSurf, inSurf, target_radius) ;

  fprintf(stderr,"surface projected - minimizing metric distortion...\n");
  MRISsetNeighborhoodSize(inSurf, nbrs) ;

  int const countNegativeFaces   = MRIScountNegativeFaces(inSurf);
  int const allowedNegativeFaces = nint(.8*inSurf->nfaces);
  if (countNegativeFaces > allowedNegativeFaces)
  {
    printf("!!!!!!!!!  everted surface detected (countNegativeFaces:%d > allowedNegativeFaces:%d) - correcting !!!!!!!!!!!!!!\n", countNegativeFaces, allowedNegativeFaces) ;
    MRISevertSurface(inSurf) ;
  }

  MRISprintTessellationStats(inSurf, stderr) ;
  MRISquickSphere(inSurf, &parms, max_passes) ;

  // if (!load)
  if (outSurf != NULL)
  {
    fflush(stdout); fflush(stderr);
    fprintf(stderr, "writing spherical brain to %s\n", outSurf) ;
    fflush(stdout); fflush(stderr);
    MRISwrite(inSurf, outSurf) ;
  }

  fflush(stdout); fflush(stderr);

  int msec = then.milliseconds() ;
  fflush(stdout); fflush(stderr);
  fprintf(stderr, "spherical transformation took %2.4f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));
  fflush(stdout); fflush(stderr);

  // Output formatted so it can be easily grepped
#ifdef HAVE_OPENMP
  int n_omp_threads = omp_get_max_threads();
  printf("FSRUNTIME@ mris_sphere %7.4f hours %d threads\n",msec/(1000.0*60.0*60.0),n_omp_threads);
#else
  printf("FSRUNTIME@ mris_sphere %7.4f hours %d threads\n",msec/(1000.0*60.0*60.0),1);
#endif
  printf("#VMPC# mris_sphere VmPeak  %d\n",GetVmPeak());
  printf("MRISQuickSphericalInflate done\n");

  return inSurf;
}

/*!
  \fn MRIS *MRISaverageSurfaces(std::vector<MRIS*> surfs, MRIS *avgsurf)
  \brief Compute the average coordinates of the given surfaces
 */
MRIS *MRISaverageSurfaces(std::vector<MRIS*> surfs, MRIS *avgsurf)
{
  if(avgsurf == NULL) avgsurf = MRISclone(surfs[0]);
  else {
    MRIS *surf = surfs[0];
    if(avgsurf->nvertices != surf->nvertices){
      printf("ERROR: MRISaverageSurfaces(): surface 0 nvertices mismatch %d %d\n",avgsurf->nvertices,surf->nvertices);
      return(NULL);
    }
    for(int vno=0; vno < avgsurf->nvertices; vno++){
      VERTEX *va = &(avgsurf->vertices[vno]);
      VERTEX *vn = &(surf->vertices[vno]);
      va->x = vn->x;
      va->y = vn->y;
      va->z = vn->z;
    }
  }
  for(int n=1; n < surfs.size(); n++){
    MRIS *surf = surfs[n];
    if(avgsurf->nvertices != surf->nvertices){
      printf("ERROR: MRISaverageSurfaces(): surface %d nvertices mismatch %d %d\n",n,avgsurf->nvertices,surf->nvertices);
      return(NULL);
    }
    for(int vno=0; vno < avgsurf->nvertices; vno++){
      VERTEX *va = &(avgsurf->vertices[vno]);
      VERTEX *vn = &(surf->vertices[vno]);
      va->x += vn->x;
      va->y += vn->y;
      va->z += vn->z;
    }
  }
  for(int vno=0; vno < avgsurf->nvertices; vno++){
    VERTEX *va = &(avgsurf->vertices[vno]);
    va->x /= surfs.size();
    va->y /= surfs.size();
    va->z /= surfs.size();
  }

  return(avgsurf);
}
