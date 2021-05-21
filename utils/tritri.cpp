/*
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

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * int tri_tri_intersect(double V0[3],double V1[3],double V2[3],
 *                         double U0[3],double U1[3],double U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
#include <math.h>
 */

#include <math.h>
#include "macros.h"
#include "tritri.h"

/* if USE_EPSILON_TEST is true then we do a check:
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define USE_EPSILON_TEST TRUE
#define EPSILON 0.000001
#define BEPSILON (0.000001 * 10)

/* sort so that a<=b */
#define SORT(a, b) \
  if (a > b) {     \
    double c;      \
    c = a;         \
    a = b;         \
    b = c;         \
  }

#define ISECT(VV0, VV1, VV2, D0, D1, D2, isect0, isect1) \
  isect0 = VV0 + (VV1 - VV0) * D0 / (D0 - D1);           \
  isect1 = VV0 + (VV2 - VV0) * D0 / (D0 - D2);

#define COMPUTE_INTERVALS(VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2, isect0, isect1) \
  if (D0D1 > 0.0f) {                                                             \
    /* here we know that D0D2<=0.0 */                                            \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */   \
    ISECT(VV2, VV0, VV1, D2, D0, D1, isect0, isect1);                            \
  }                                                                              \
  else if (D0D2 > 0.0f) {                                                        \
    /* here we know that d0d1<=0.0 */                                            \
    ISECT(VV1, VV0, VV2, D1, D0, D2, isect0, isect1);                            \
  }                                                                              \
  else if (D1 * D2 > 0.0f || D0 != 0.0f) {                                       \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */                            \
    ISECT(VV0, VV1, VV2, D0, D1, D2, isect0, isect1);                            \
  }                                                                              \
  else if (D1 != 0.0f) {                                                         \
    ISECT(VV1, VV0, VV2, D1, D0, D2, isect0, isect1);                            \
  }                                                                              \
  else if (D2 != 0.0f) {                                                         \
    ISECT(VV2, VV0, VV1, D2, D0, D1, isect0, isect1);                            \
  }                                                                              \
  else {                                                                         \
    /* triangles are coplanar */                                                 \
    return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);                         \
  }

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */
#define EDGE_EDGE_TEST(V0, U0, U1)                                  \
  Bx = U0[i0] - U1[i0];                                             \
  By = U0[i1] - U1[i1];                                             \
  Cx = V0[i0] - U0[i0];                                             \
  Cy = V0[i1] - U0[i1];                                             \
  f = Ay * Bx - Ax * By;                                            \
  d = By * Cx - Bx * Cy;                                            \
  if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) { \
    e = Ax * Cy - Ay * Cx;                                          \
    if (f > 0) {                                                    \
      if (e >= 0 && e <= f) return 1;                               \
    }                                                               \
    else {                                                          \
      if (e <= 0 && e >= f) return 1;                               \
    }                                                               \
  }

#define EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2) \
  {                                                \
    double Ax, Ay, Bx, By, Cx, Cy, e, d, f;        \
    Ax = V1[i0] - V0[i0];                          \
    Ay = V1[i1] - V0[i1];                          \
    /* test edge U0,U1 against V0,V1 */            \
    EDGE_EDGE_TEST(V0, U0, U1);                    \
    /* test edge U1,U2 against V0,V1 */            \
    EDGE_EDGE_TEST(V0, U1, U2);                    \
    /* test edge U2,U1 against V0,V1 */            \
    EDGE_EDGE_TEST(V0, U2, U0);                    \
  }

#define POINT_IN_TRI(V0, U0, U1, U2)          \
  {                                           \
    double a, b, c, d0, d1, d2;               \
    /* is T1 completly inside T2? */          \
    /* check if V0 is inside tri(U0,U1,U2) */ \
    a = U1[i1] - U0[i1];                      \
    b = -(U1[i0] - U0[i0]);                   \
    c = -a * U0[i0] - b * U0[i1];             \
    d0 = a * V0[i0] + b * V0[i1] + c;         \
                                              \
    a = U2[i1] - U1[i1];                      \
    b = -(U2[i0] - U1[i0]);                   \
    c = -a * U1[i0] - b * U1[i1];             \
    d1 = a * V0[i0] + b * V0[i1] + c;         \
                                              \
    a = U0[i1] - U2[i1];                      \
    b = -(U0[i0] - U2[i0]);                   \
    c = -a * U2[i0] - b * U2[i1];             \
    d2 = a * V0[i0] + b * V0[i1] + c;         \
    if (d0 * d1 > 0.0) {                      \
      if (d0 * d2 > 0.0) return 1;            \
    }                                         \
  }

int coplanar_tri_tri(double N[3], double V0[3], double V1[3], double V2[3], double U0[3], double U1[3], double U2[3]);

int coplanar_tri_tri(double N[3], double V0[3], double V1[3], double V2[3], double U0[3], double U1[3], double U2[3])
{
  double A[3];
  short i0, i1;
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  A[0] = fabs(N[0]);
  A[1] = fabs(N[1]);
  A[2] = fabs(N[2]);
  if (A[0] > A[1]) {
    if (A[0] > A[2]) {
      i0 = 1; /* A[0] is greatest */
      i1 = 2;
    }
    else {
      i0 = 0; /* A[2] is greatest */
      i1 = 1;
    }
  }
  else /* A[0]<=A[1] */
  {
    if (A[2] > A[1]) {
      i0 = 0; /* A[2] is greatest */
      i1 = 1;
    }
    else {
      i0 = 0; /* A[1] is greatest */
      i1 = 2;
    }
  }

  /* test all edges of triangle 1 against the edges of triangle 2 */
  EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
  EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
  EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);

  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  POINT_IN_TRI(V0, U0, U1, U2);
  POINT_IN_TRI(U0, V0, V1, V2);

  return 0;
}

int tri_tri_intersect(double V0[3], double V1[3], double V2[3], double U0[3], double U1[3], double U2[3])
{
  double E1[3], E2[3];
  double N1[3], N2[3], d1, d2;
  double du0, du1, du2, dv0, dv1, dv2, fdu0, fdu1, fdu2, fdv0, fdv1, fdv2;
  double D[3];
  double isect1[2], isect2[2];
  double du0du1, du0du2, dv0dv1, dv0dv2;
  short index;
  double vp0, vp1, vp2;
  double up0, up1, up2;
  double b, c, max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1, V1, V0);
  SUB(E2, V2, V0);
  CROSS(N1, E1, E2);
  d1 = -DOT(N1, V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = DOT(N1, U0) + d1;
  du1 = DOT(N1, U1) + d1;
  du2 = DOT(N1, U2) + d1;
  fdu0 = fabs(du0);
  fdu1 = fabs(du1);
  fdu2 = fabs(du2);
  /* coplanarity robustness check */

  /* same sign on all of them + not equal 0 ? */
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;
  if ((fdu0 > EPSILON || fdu1 > EPSILON || fdu2 > EPSILON) && (du0du1 > 0.0f && du0du2 > 0.0f))
    return 0; /* no intersection occurs */
#if USE_EPSILON_TEST == TRUE
  if (fdu0 < EPSILON) du0 = 0.0;
  if (fdu1 < EPSILON) du1 = 0.0;
  if (fdu2 < EPSILON) du2 = 0.0;
#endif
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;
  if ((fdu0 < BEPSILON) && (fdu1 < BEPSILON) && (fdu2 < EPSILON)) du0du1 = du0du2 = du0 = du1 = du2 = 0.0;

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1, U1, U0);
  SUB(E2, U2, U0);
  CROSS(N2, E1, E2);
  d2 = -DOT(N2, U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0 = DOT(N2, V0) + d2;
  dv1 = DOT(N2, V1) + d2;
  dv2 = DOT(N2, V2) + d2;
  fdv0 = fabs(dv0);
  fdv1 = fabs(dv1);
  fdv2 = fabs(dv2);

  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;

  if ((fdv0 > EPSILON || fdv1 > EPSILON || fdv2 > EPSILON) && (dv0dv1 > 0.0f && dv0dv2 > 0.0f))
    return 0; /* no intersection occurs */
#if USE_EPSILON_TEST == TRUE
  if (fabs(dv0) < EPSILON) dv0 = 0.0;
  if (fabs(dv1) < EPSILON) dv1 = 0.0;
  if (fabs(dv2) < EPSILON) dv2 = 0.0;
#endif

  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;
  if ((fdv0 < BEPSILON) && (fdv1 < BEPSILON) && (fdv2 < EPSILON)) dv0dv1 = dv0dv2 = dv0 = dv1 = dv2 = 0.0;

  /* compute direction of intersection line */
  CROSS(D, N1, N2);

  /* compute and index to the largest component of D */
  max = fabs(D[0]);
  index = 0;
  b = fabs(D[1]);
  c = fabs(D[2]);
  if (b > max) max = b, index = 1;
  if (c > max) max = c, index = 2;

  /* this is the simplified projection onto L*/
  vp0 = V0[index];
  vp1 = V1[index];
  vp2 = V2[index];

  up0 = U0[index];
  up1 = U1[index];
  up2 = U2[index];

  /* compute interval for triangle 1 */
  COMPUTE_INTERVALS(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, isect1[0], isect1[1]);

  /* compute interval for triangle 2 */
  COMPUTE_INTERVALS(up0, up1, up2, du0, du1, du2, du0du1, du0du2, isect2[0], isect2[1]);

  SORT(isect1[0], isect1[1]);
  SORT(isect2[0], isect2[1]);

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;
  return 1;
}

#define MAT_COL_SET(m, c, v)           \
  *MATRIX_RELT(m, 1, c) = (float)v[0]; \
  *MATRIX_RELT(m, 2, c) = (float)v[1]; \
  *MATRIX_RELT(m, 3, c) = (float)v[2];

#define MAT_COL_GET(v, m, c)    \
  v[0] = *MATRIX_RELT(m, 1, c); \
  v[1] = *MATRIX_RELT(m, 2, c); \
  v[2] = *MATRIX_RELT(m, 3, c);

#define VEC_GET(v, vec)      \
  v[0] = VECTOR_ELT(vec, 1); \
  v[1] = VECTOR_ELT(vec, 2); \
  v[2] = VECTOR_ELT(vec, 3);

#define VEC_SET(vec, v)      \
  VECTOR_ELT(vec, 1) = v[0]; \
  VECTOR_ELT(vec, 2) = v[1]; \
  VECTOR_ELT(vec, 3) = v[2];

/*
  find the intersection of a ray and a triangle. If none, return 0.
  If found, return the location in T.
  P is the origin of the ray. U0, U1, U2 contains the
  coordinates of the vertices, and D is the direction of the ray.
*/
#include "matrix.h"
int triangle_ray_intersect(double orig_pt[3], double dir[3], double U0[3], double U1[3], double U2[3], double int_pt[3])
{
  double basis1[3], basis2[3], tmp[3], L0[3], L1[3], len, t, dot, norm_proj[3], *V0, *V1, *V2, desc[3], Point[3];
  // douoble a, b;
  MATRIX *m_U, *m_inv;
  VECTOR *v_p, *v_r;
  int i;

  m_U = MatrixAlloc(3, 3, MATRIX_REAL);
  v_p = VectorAlloc(3, MATRIX_REAL);
  v_r = VectorAlloc(3, MATRIX_REAL);

  /* first compute basis vectors for plane defined by triangle */

  /* first basis vector is U1 - U0 */
  SUB(basis1, U1, U0);

  /* compute second basis vector by crossing legs */
  SUB(basis2, U2, U0);
  CROSS(tmp, basis1, basis2);
  CROSS(basis2, basis1, tmp); /* now basis2 is orthogonal to basis1 */

  /* normalize their lengths */
  len = VLEN(basis1);
  if (FZERO(len)) return (0);
  basis1[0] /= len;
  basis1[1] /= len;
  basis1[2] /= len;
  len = VLEN(basis2);
  if (FZERO(len)) return (0);
  basis2[0] /= len;
  basis2[1] /= len;
  basis2[2] /= len;

  /*
     build matrix: 1st two cols are basis vectors and third is
     negative of direction vector. Inverting this will yield solution.
  */
  MAT_COL_SET(m_U, 1, basis1);
  MAT_COL_SET(m_U, 2, basis2);
  MAT_COL_SET(m_U, 3, -dir);
  m_inv = MatrixSVDInverse(m_U, NULL);
  SUB(Point, orig_pt, U0);
  VEC_SET(v_p, Point); /* pt in U0 coord sys */
  MatrixMultiply(m_inv, v_p, v_r);

  /* a and b are coordinate of point in plane, t is parameterization of ray */
  // a = VECTOR_ELT(v_r, 1);
  // b = VECTOR_ELT(v_r, 2);
  t = VECTOR_ELT(v_r, 3);
  MatrixFree(&m_U);
  VectorFree(&v_p);
  VectorFree(&v_r);
  MatrixFree(&m_inv);

  /* coordinates of interesection point */
  int_pt[0] = orig_pt[0] + t * dir[0];
  int_pt[1] = orig_pt[1] + t * dir[1];
  int_pt[2] = orig_pt[2] + t * dir[2];

  /* now determine whether the point is in the triagle by seeing if it is
     on the 'right' halfplane defined by each triangle leg
  */

  for (i = 0; i < 3; i++) {
    /*
       build a coordinate system with V0 as the origin, then construct
       the vector connecting V2 with it's normal projection onto V0->V1.
       This will be a descriminant vector for dividing the plane by the
       V0->V1 line. A positive dot product with the desc. vector indicates
       that the point is on the positive side of the plane and therefore
       may be contained within the triangle. Doing this for each of the
       legs in sequence gives a test for being inside the triangle.
       */

    switch (i) {
      default:
      case 0:
        V0 = U0;
        V1 = U1;
        V2 = U2;
        break;
      case 1:
        V0 = U1;
        V1 = U2;
        V2 = U0;
        break;
      case 2:
        V0 = U2;
        V1 = U0;
        V2 = U1;
        break;
    }
    SUB(L0, V1, V0);
    SUB(L1, V2, V0);

    /* compute normal projection onto base of triangle */
    len = VLEN(L0);
    L0[0] /= len;
    L0[1] /= len;
    L0[2] /= len;
    dot = DOT(L0, L1);
    SCALAR_MUL(norm_proj, dot, L0);

    /* build descriminant vector */
    SUB(desc, L1, norm_proj);

    /*
       transform point in question into local coordinate system and build
       the vector from the point in question to the normal projection point.
       The dot product of this vector with the descrimant vector will then
       indicate which side of the V0->V1 line the point is on.
       */
    SUB(Point, int_pt, V0);
    SUB(Point, Point, norm_proj);
    dot = DOT(desc, Point);
    if (dot < 0 && !DZERO(dot)) return (0);
  }
  return (1);
}

// from http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
// intersect_RayTriangle(): intersect a ray with a 3D triangle
//    Input:  a ray n, and a triangle T
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 = disjoint (no intersect)
//             1 = intersect in unique point I1
//             2 = are in the same plane
int intersect_RayTriangle(double ray[2][3], double V0[3], double V1[3], double V2[3], double I[3])
{
  double u[3], v[3], n[3];     // triangle vectors
  double dir[3], w0[3], w[3];  // ray vectors
  float r, a, b;               // params to calc ray-plane intersect
  float uu, uv, vv, wu, wv, D;
  float s, t;

  // get triangle edge vectors and plane normal
  SUB(u, V1, V0);
  SUB(v, V2, V0);
  CROSS(n, u, v);    // cross product
  if (VLEN(n) == 0)  // triangle is degenerate
    return -1;       // do not deal with this case

  SUB(dir, ray[1], ray[0]);  // ray direction vector
  SUB(w0, ray[0], V0);
  a = -DOT(n, w0);
  b = DOT(n, dir);
#define SMALL_NUM 1e-8
  if (fabs(b) < SMALL_NUM) {  // ray is parallel to triangle plane
    if (a == 0)               // ray lies in triangle plane
      return 2;
    else
      return 0;  // ray disjoint from plane
  }

  // get intersect point of ray with triangle plane
  r = a / b;
  if (r < 0.0)  // ray goes away from triangle
    return 0;   // => no intersect
  // for a segment, also test if (r > 1.0) => no intersect

  SCALAR_MUL(dir, r, dir);
  ADD(I, ray[0], dir);

  // is I inside T?
  uu = DOT(u, u);
  uv = DOT(u, v);
  vv = DOT(v, v);
  SUB(w, I, V0);
  wu = DOT(w, u);
  wv = DOT(w, v);
  D = uv * uv - uu * vv;

  // get and test parametric coords
  s = (uv * wv - vv * wu) / D;
  if (s < 0.0 || s > 1.0)  // I is outside T
    return 0;
  t = (uv * wu - uu * wv) / D;
  if (t < 0.0 || (s + t) > 1.0)  // I is outside T
    return 0;

  return 1;  // I is in T
}
int project_point_to_plane(
    double point[3], double V0[3], double V1[3], double V2[3], double proj[3], double pe1[3], double pe2[3])
{
  double e1[3], e2[3], norm, d1, d2, tmp[3];

  // first basis vector
  SUB(e1, V1, V0);
  norm = VLEN(e1);
  SCALAR_MUL(e1, 1.0 / norm, e1);
  if (DZERO(norm)) {
    ADD(point, V0, V1);
    ADD(proj, point, V2);
    SCALAR_MUL(proj, 1.0 / 3.0, proj);
    return (-1);
  }

  // project 1st basis vector out of 2nd to orthonormalize it
  SUB(e2, V2, V0);
  norm = VLEN(e2);
  SCALAR_MUL(e2, 1.0 / norm, e2);
  if (DZERO(norm)) {
    ADD(point, V0, V1);
    ADD(proj, point, V2);
    SCALAR_MUL(proj, 1.0 / 3.0, proj);
    return (-2);
  }

  d1 = DOT(e1, e2);
  SCALAR_MUL(tmp, d1, e1);
  SUB(e2, e2, tmp);
  norm = VLEN(e2);
  if (DZERO(norm)) {
    ADD(point, V0, V1);
    ADD(proj, point, V2);
    SCALAR_MUL(proj, 1.0 / 3.0, proj);
    return (-1);
  }
  SCALAR_MUL(e2, 1.0 / norm, e2);

  if (pe1)  // return coordinate basis
  {
    memmove(pe1, e1, sizeof(e1));
    memmove(pe2, e2, sizeof(e2));
  }

  SUB(point, point, V0);
  d1 = DOT(e1, point);
  d2 = DOT(e2, point);
  SCALAR_MUL(e1, d1, e1);
  SCALAR_MUL(e2, d2, e2);
  ADD(proj, e1, e2);
  ADD(proj, V0, proj);
  return (0);
}
