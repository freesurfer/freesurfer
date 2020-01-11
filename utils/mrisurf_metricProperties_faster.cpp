#include <algorithm>

#include "mrisurf.h"
#include "mrisurf_metricProperties.h"
#include "timer.h"
#include "fnv_hash.h"
#include "romp_support.h"


// define some very straightfoward 3D vector operations - nothing fancy
#define vec_get(array, index) (&array[index * 3])
#define vec_length(v) sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
#define vec_dot(a, b) ( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] )
#define vec_add(v, d) { \
  v[0] += d[0]; \
  v[1] += d[1]; \
  v[2] += d[2]; }
#define vec_mult(v, d) { \
  v[0] *= d; \
  v[1] *= d; \
  v[2] *= d; }
#define vec_divide(v, d) { \
  v[0] /= d; \
  v[1] /= d; \
  v[2] /= d; }
#define vec_edge(a, b, v) { \
  v[0] = b[0] - a[0]; \
  v[1] = b[1] - a[1]; \
  v[2] = b[2] - a[2]; }
#define vec_cross(a, b, c) { \
  c[0] = a[1] * b[2] - a[2] * b[1]; \
  c[1] = a[2] * b[0] - a[0] * b[2]; \
  c[2] = a[0] * b[1] - a[1] * b[0]; }


class MetricPropertiesCalculator
{
public:
  MetricPropertiesCalculator(MRIS *surf) : mris(surf) { init(); }
  ~MetricPropertiesCalculator();

  void compute();
  void transfer();

  // if enabled, face normals will by weighted by their angle when
  // computing vertex normals (this increases runtime slightly)
  bool weight_norms_by_angle = false;

  // max possible distance to move when reorganizing overlapping vertices
  const double rand_dist = 0.001;

private:
  void init();

  // these steps are separated purely for readability, they should never
  // be run out of order or outside of the compute() function (hence why
  // they're private)
  void computeNormalsAreasAndAngles();
  void computeVertexDistances();
  void computeDimensions();
  void computeEuclideanVertexDistances();
  void computeSphericalVertexDistances();
  void computeAverageVertexDistance();
  void orient();
  void orientPlane();
  void orientEllipsoid();

  // utility to move overlapping vertices
  void moveVertexRandomly(int vno, int trial, int *random_counter);

  // keeps track of whether the vertex positions have been messed
  // with (this happens when two vertices in a face are overlapping)
  bool vertices_have_moved = false;

  // options
  float vertex_area_fix_value = 3.0;
  bool surf_is_3d = false;

  // vertex info arrays
  int nvertices;
  float *vertices;
  float *face_angles;
  float *vertex_areas;
  float *vertex_normals;
  bool *vertex_ripped;

  // vertex topology arrays
  int *vertex_nfaces;
  int **vertex_faces;
  int *vertex_vnums;
  int **vertex_neighbors;
  int *vertex_vtotals;
  float **vertex_vdistances;

  // face info arrays
  int nfaces;
  int *faces;
  float *face_areas;
  float *face_normals;
  bool *face_ripped;

  // stored mris structure
  MRIS* mris;
};


MetricPropertiesCalculator::~MetricPropertiesCalculator()
{
  // TODO: using vectors shouldn't slow things down significantlly, let's
  // look into using those so we don't have to worry about allocating / deallocating
  delete [] vertices;
  delete [] faces;
  delete [] face_areas;
  delete [] face_normals;
  delete [] face_angles;
  delete [] face_ripped;
  delete [] vertex_areas;
  delete [] vertex_normals;
  delete [] vertex_ripped;
  delete [] vertex_vdistances;
  delete [] vertex_neighbors;
  delete [] vertex_vtotals;
  delete [] vertex_vnums;
  delete [] vertex_faces;
  delete [] vertex_nfaces;
}


void MetricPropertiesCalculator::init()
{
  nfaces = mris->nfaces;
  nvertices = mris->nvertices;

  // configure some initial settings
  surf_is_3d = !(mris->status == MRIS_PLANE || mris->status == MRIS_CUT);
  vertex_area_fix_value = MRISgetFixVertexAreaValue() ? 3.0 : 2.0;

  // transfer vertex positions to xyz array
  vertices = new float[nvertices * 3];
  vertex_ripped = new bool[nvertices];
  std::fill(vertex_ripped, vertex_ripped + nvertices, false);

  vertex_vdistances = new float*[nvertices];
  vertex_neighbors = new int*[nvertices];
  vertex_vtotals = new int[nvertices];
  vertex_vnums = new int[nvertices];
  vertex_faces = new int*[nvertices];
  vertex_nfaces = new int[nvertices];

  #pragma omp parallel for
  for (int vno = 0 ; vno < nvertices ; vno++) {
    VERTEX *vertex = &mris->vertices[vno];

    float *v = vec_get(vertices, vno);
    v[0] = vertex->x;
    v[1] = vertex->y;
    v[2] = vertex->z;

    vertex_ripped[vno] = vertex->ripflag;
    
    MRISmakeDist(mris, vno);

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    vertex_vdistances[vno] = vertex->dist;
    vertex_neighbors[vno] = vt->v;
    vertex_vtotals[vno] = vt->vtotal;
    vertex_vnums[vno] = vt->vnum;
    vertex_faces[vno] = vt->f;
    vertex_nfaces[vno] = vt->num;
  }

  // transfer face vertex indices to N x 3 array
  faces = new int[nfaces * 3];
  face_ripped = new bool[nfaces];
  std::fill(face_ripped, face_ripped + nfaces, false);

  #pragma omp parallel for
  for (int fno = 0 ; fno < nfaces ; fno++) {
    FACE *face = &mris->faces[fno];

    int *f = vec_get(faces, fno);
    f[0] = face->v[0];
    f[1] = face->v[1];
    f[2] = face->v[2];
    face_ripped[fno] = face->ripflag;

    if (face->ripflag) {
      for (int n = 0; n < 3; n++) {
        #pragma omp critical
        { mris->vertices[face->v[n]].border = true; }
      }
    }
  }

  // allocate face metrics
  face_areas = new float[nfaces];
  face_normals = new float[nfaces * 3];
  face_angles = new float[nfaces * 3];

  // allocate vertex metrics
  vertex_areas = new float[nvertices];
  vertex_normals = new float[nvertices * 3];
}


/*
  Transfers any remaining computed values back into the MRIS structure.
*/
void MetricPropertiesCalculator::transfer()
{
  // start by transferring the vertex metrics
  #pragma omp parallel for
  for (int vno = 0 ; vno < nvertices ; vno++) {
    if (vertex_ripped[vno]) continue;
    VERTEX *vertex = &mris->vertices[vno];

    // vertex positions (maybe)
    if (vertices_have_moved) {
      float *pos = vec_get(vertices, vno);
      vertex->x = pos[0];
      vertex->y = pos[1];
      vertex->z = pos[2];
    }

    // vertex normals
    float *norm = vec_get(vertex_normals, vno);
    vertex->nx = norm[0];
    vertex->ny = norm[1];
    vertex->nz = norm[2];

    // vertex areas
    vertex->area = vertex_areas[vno];
    if (vertex->origarea < 0) vertex->origarea = abs(vertex_areas[vno]);
  }

  // now transfer the face metrics
  #pragma omp parallel for
  for (int fno = 0 ; fno < nfaces ; fno++) {
    if (face_ripped[fno]) continue;
    FACE *face = &mris->faces[fno];

    // face normals
    float *norm = vec_get(face_normals, fno);
    setFaceNorm(mris, fno, norm[0], norm[1], norm[2]);

    // face angles
    float *angles = vec_get(face_angles, fno);
    face->angle[0] = angles[0];
    face->angle[1] = angles[1];
    face->angle[2] = angles[2];

    // face areas
    face->area = face_areas[fno];
  }

  // set dist_nsize (this is replicating functionality from mrisComputeVertexDistances)
  mris->dist_nsize = mris->nsize;
}


/*
  Computes all metric properties. Order matters here!
*/
void MetricPropertiesCalculator::compute()
{
  // compute normals, areas, and face angles 
  computeNormalsAreasAndAngles();

  // compute surface bounding box
  computeDimensions();

  // compute distances between vertices and calculate average
  computeVertexDistances();

  // compute the average vertex distance
  computeAverageVertexDistance();

  // compute an average vertex area before orienting, since it might
  // change the total area (if there are negative faces)
  mris->avg_vertex_area = mris->total_area / nvertices;

  // fix any negatively-oriented face normals for spheres and planes 
  orient();

  // compute theoretical total area if sphere
  if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_RIGID_BODY || mris->status == MRIS_SPHERE) {
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0;
  }
}


/*
  Computes all metric properties. Order matters here!
*/
void MetricPropertiesCalculator::computeNormalsAreasAndAngles()
{
  // might take a few tries if we have a really ugly mesh
  for (int trial = 0 ; trial < 5 ; trial++) {

    // reset any cumulative values to zero
    double total_area = 0;
    std::fill(vertex_areas, vertex_areas + nvertices, 0);
    std::fill(vertex_normals, vertex_normals + (nvertices * 3), 0);

    // loop through faces and their calculate areas, normals, and angles
    #pragma omp parallel for reduction(+:total_area)
    for (int fno = 0 ; fno < nfaces ; fno++) {

      // skip ripped faces
      if (face_ripped[fno]) continue;

      // get triangle vertex numbers
      int *vnos = vec_get(faces, fno);
      float *pos0 = vec_get(vertices, vnos[0]);
      float *pos1 = vec_get(vertices, vnos[1]);
      float *pos2 = vec_get(vertices, vnos[2]);

      float a[3], b[3], c[3];

      // a is edge vector v0 -> v1
      // b is edge vector v0 -> v2
      vec_edge(pos0, pos1, a);
      vec_edge(pos0, pos2, b);

      // face norm = a x b
      float *fnorm = vec_get(face_normals, fno);
      vec_cross(a, b, fnorm);

      // this is actually area * 2 for now
      float area = vec_length(fnorm);

      // unitize face normal
      if (area > 0.0) vec_divide(fnorm, area);

      // halve to get true triangular area
      area /= 2.0;
      face_areas[fno] = area;
      total_area += area;

      // compute the angle in each face corner
      float *angles = vec_get(face_angles, fno);
      for (int corner = 0 ; corner < 3 ; corner++) {

        // find edges defining each corner
        switch (corner) {
          case 0:
            vec_edge(pos0, pos2, a);
            vec_edge(pos0, pos1, b);
            break;
          case 1:
            vec_edge(pos1, pos0, a);
            vec_edge(pos1, pos2, b);
            break;
          case 2:
            vec_edge(pos2, pos1, a);
            vec_edge(pos2, pos0, b);
            break;
        }

        // compute and set angle from edges
        vec_cross(b, a, c);
        float d1 = vec_dot(c, fnorm);
        float d2 = vec_dot(a, b);
        float angle = fastApproxAtan2f(d1, d2);
        angles[corner] = angle;
      }
    }

    // transfer to MRIS
    mris->total_area = total_area;

    // as we loop through the vertices below, let's check if any have invalid normals
    bool vertices_must_move = false;

    // for each vertex, sum surrounding the face normals and areas
    #pragma omp parallel for
    for (int vno = 0 ; vno < nvertices ; vno++) {
      if (vertex_ripped[vno]) continue;

      // prepare vertex norm array
      float *norm = vec_get(vertex_normals, vno);

      int numfaces = vertex_nfaces[vno];
      int *neighbor_faces = vertex_faces[vno];
      for (int nf = 0 ; nf < numfaces ; nf++) {
        int fno = neighbor_faces[nf];

        // sum face areas
        vertex_areas[vno] += face_areas[fno];

        // sum face normals
        float *face_norm = vec_get(face_normals, fno);
        if (!weight_norms_by_angle) {
          vec_add(norm, face_norm);
        } else {
          // We have to relocate the index of the vertex in each face to
          // find the corresponding corner angle. This introduces a roughly
          // 5% slowdown for all of computeMetricProperties.
          int *vnos = vec_get(faces, fno);
          float *angles = vec_get(face_angles, fno);
          int index = std::distance(vnos, std::find(vnos, vnos + 3, vno));
          float angle = angles[index];
          norm[0] += face_norm[0] * angle;
          norm[1] += face_norm[1] * angle;
          norm[2] += face_norm[2] * angle;
        }
      }

      // divide total face area, since each face is added to three vertices
      vertex_areas[vno] /= vertex_area_fix_value;

      // unitize vertex normal and ensure it has a valid length; if not,
      // lets move the vertices a bit
      float length = vec_length(norm);
      if (length > 0) {
        vec_divide(norm, length);
      } else {
        // mark vertex_areas as a nifty way to keep track of bad vertices
        vertices_must_move = true;
        vertex_areas[vno] = -1;
      }
    }

    // all done if vertices aren't overlapping
    if (!vertices_must_move) return;

    // shift overlapping vertices (and their neighbors)
    // this can't be safely parallelized, but it probably doesn't matter as it only gets reached in rare cases
    for (int vno = 0 ; vno < nvertices ; vno++) {
      // bad vertices were marked by setting the area to -1
      if (vertex_areas[vno] < 0) {
        int random_counter = 0;
        moveVertexRandomly(vno, trial, &random_counter);
        const int * neighbors = vertex_neighbors[vno];
        for (int n = 0; n < vertex_vnums[vno]; n++) moveVertexRandomly(neighbors[n], trial, &random_counter);
      }
    }

    vertices_have_moved = true;
  }
}


/*
  Shifts a vertex position in a random but deterministic direction.
*/
void MetricPropertiesCalculator::moveVertexRandomly(int vno, int trial, int *random_counter)
{
  float random[3];
  random[0] = fnv_hash(trial, vno, random_counter, -rand_dist, rand_dist);
  random[1] = fnv_hash(trial, vno, random_counter, -rand_dist, rand_dist);
  random[2] = surf_is_3d ? fnv_hash(trial, vno, random_counter, -rand_dist, rand_dist) : 0.0;
  float *pos = vec_get(vertices, vno);
  vec_add(pos, random);
}


/*
  Computes the bounding box min, max, and center coordinates defined by the vertex positions.
*/
void MetricPropertiesCalculator::computeDimensions()
{
  float xhi = vertices[0];
  float yhi = vertices[1];
  float zhi = vertices[2];
  float xlo = xhi;
  float ylo = yhi;
  float zlo = zhi;

  for (int vno = 1; vno < nvertices; vno++) {
    float const *pos = vec_get(vertices, vno);
    float const x = pos[0];
    float const y = pos[1];
    float const z = pos[2];
    if (x > xhi) xhi = x;
    if (x < xlo) xlo = x;
    if (y > yhi) yhi = y;
    if (y < ylo) ylo = y;
    if (z > zhi) zhi = z;
    if (z < zlo) zlo = z;
  }

  mris->xctr = 0.5f * (float)((double)xlo + (double)xhi);
  mris->yctr = 0.5f * (float)((double)ylo + (double)yhi);
  mris->zctr = 0.5f * (float)((double)zlo + (double)zhi);
  mris->xhi = xhi;
  mris->yhi = yhi;
  mris->zhi = zhi;
  mris->xlo = xlo;
  mris->ylo = ylo;
  mris->zlo = zlo;
}


/*
  Computes distances between a vertex and all vertices in its
  neighborhood. Distances are euclidean unless surface is a sphere.
*/
void MetricPropertiesCalculator::computeVertexDistances()
{
  if (mris->status == MRIS_SPHERE || mris->status == MRIS_PARAMETERIZED_SPHERE) {
    computeSphericalVertexDistances();
  } else {
    computeEuclideanVertexDistances();
  }
}


/*
  Computes the euclidean distances between a vertex and its neighbors.
*/
void MetricPropertiesCalculator::computeEuclideanVertexDistances()
{
  #pragma omp parallel for
  for (int vno = 0; vno < nvertices; vno++) {
    if (vertex_ripped[vno]) continue;

    float * const distances = vertex_vdistances[vno];
    const int * neighbors = vertex_neighbors[vno];

    int vtotal = vertex_vtotals[vno];

    float *pos = vec_get(vertices, vno);

    float diff[3];
    for (int n = 0; n < vtotal; n++) {
      float *npos = vec_get(vertices, neighbors[n]);
      diff[0] = pos[0] - npos[0];
      diff[1] = pos[1] - npos[1];
      diff[2] = pos[2] - npos[2];
      distances[n] = vec_length(diff);
    }
  }
}


/*
  Computes the arced distances between a vertex and its neighbors along the surface.
*/
void MetricPropertiesCalculator::computeSphericalVertexDistances()
{
  // sped up by computing the normalized vectors for all the vertices once rather than repeatedly
  // and storing them along with their ripflag in a structure!
  struct PerVertexInfo {
    char  ripped;
    char  rippedOrZeroLength;
    XYZ   xyz_unnormalized;
    float xyz_length;
    XYZ   xyz_normalized;
  };

  PerVertexInfo *vertex_infos =  new PerVertexInfo[nvertices];

  // in parallel, load the normalized xyz so its only computed once per vertex
  #pragma omp parallel for
  for (int vno = 0; vno < nvertices; vno++) {

    char const ripped = (char)vertex_ripped[vno];
    vertex_infos[vno].ripped = ripped;
    vertex_infos[vno].rippedOrZeroLength = ripped;
    if (ripped) continue;

    // length xyz_length
    float *pos = vec_get(vertices, vno);
    vertex_infos[vno].xyz_unnormalized.x = pos[0];
    vertex_infos[vno].xyz_unnormalized.y = pos[1];
    vertex_infos[vno].xyz_unnormalized.z = pos[2];

    // length 1 along radius vector
    XYZ_NORMALIZED_LOAD(&vertex_infos[vno].xyz_normalized, &vertex_infos[vno].xyz_length, pos[0], pos[1], pos[2]);
    if (FZERO(vertex_infos[vno].xyz_length)) vertex_infos[vno].rippedOrZeroLength = (char)true;
  }

  #pragma omp parallel for
  for (int vno = 0; vno < nvertices; vno++) {

    PerVertexInfo* vinfo = &vertex_infos[vno];

    // non-ripflag with a radius of zero still need to get their distances set to zero
    if (vinfo->ripped) continue;

    float const radius = vinfo->xyz_length;

    float * const distances = vertex_vdistances[vno];
    const int * neighbors = vertex_neighbors[vno];

    int vtotal = vertex_vtotals[vno];

    for (int n = 0; n < vtotal; n++) {
      PerVertexInfo* ninfo = &vertex_infos[neighbors[n]];
      float d = 0.0;

      if (!ninfo->rippedOrZeroLength) {
        // radians - so 2pi around the circumference
        float angle = fabs(XYZApproxAngle_knownLength(
              &vinfo->xyz_normalized,
              ninfo->xyz_unnormalized.x,
              ninfo->xyz_unnormalized.y,
              ninfo->xyz_unnormalized.z,
              ninfo->xyz_length));

        d = angle * radius;
      }

      distances[n] = d;
    }
  }

  delete [] vertex_infos;
}


/*
  Computes the average (and std) distance between vertices and their immediate neighbors.
*/
void MetricPropertiesCalculator::computeAverageVertexDistance()
{
  double sum = 0;
  double sum2 = 0;
  double N = 0;

  // ROMP: begin partial sum loop
  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             nvertices
  #define ROMP_SUMREDUCTION0  sum
  #define ROMP_SUMREDUCTION1  sum2
  #define ROMP_SUMREDUCTION2  N
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    #define sum  ROMP_PARTIALSUM(0)
    #define sum2 ROMP_PARTIALSUM(1)
    #define N    ROMP_PARTIALSUM(2)

    if (vertex_ripped[vno]) continue;

    float * const distances = vertex_vdistances[vno];
    const int * neighbors = vertex_neighbors[vno];

    int vnum = vertex_vnums[vno];

    for (int n = 0; n < vnum; n++) {
      if (vertex_ripped[neighbors[n]]) continue;
      double d = distances[n];
      sum += d;
      sum2 += (d * d);
      N += 1;
    }

    #undef sum
    #undef sum2
    #undef N
  #include "romp_for_end.h"
  // ROMP: end partial sum loop

  // NOTE: this is a poor algorithm for computing the std dev because of how the floating point
  // errors accumulate but it seems to work for us because the double has enough accuracy to sum
  // the few hundred thousand small (but not too small) floats that we have
  double avg = sum / N;
  double std = sqrt(N * (sum2 / N - avg * avg) / (N - 1));

  // transfer the computed distances
  mris->avg_vertex_dist = avg;
  mris->std_vertex_dist = std;
}


/*
  Fixes any negatively-oriented face normals for spheres and planes and sums
  the total negative area.
*/
void MetricPropertiesCalculator::orient()
{
  switch (mris->status) {
    case MRIS_RIGID_BODY:
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
    case MRIS_ELLIPSOID:
    case MRIS_SPHERICAL_PATCH:
      orientEllipsoid();
      break;
    case MRIS_PLANE:
      orientPlane();
      break;
    default:
      // do nothing for a regular surfaces
      break;
  }
}


/*
  Looks for vertex and face normals oriented negatively in the z plane.
*/
void MetricPropertiesCalculator::orientPlane()
{
  double total_area = 0.0;
  double neg_area = 0.0;
  double neg_orig_area = 0.0;

  // ROMP: begin partial sum loop
  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             nfaces
  #define ROMP_SUMREDUCTION0  total_area
  #define ROMP_SUMREDUCTION1  neg_area
  #define ROMP_SUMREDUCTION2  neg_orig_area
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    #define total_area      ROMP_PARTIALSUM(0)
    #define neg_area        ROMP_PARTIALSUM(1)
    #define neg_orig_area   ROMP_PARTIALSUM(2)
  
    if (face_ripped[fno]) ROMP_PFLB_continue;

    // if the unit normal is pointing downwards in the plane,
    // then the area should be negative, and the normal and angles
    // should be reversed

    float *norm = vec_get(face_normals, fno);
    if (norm[2] < 0.0f) {
      // not in same direction, area < 0 and reverse n
      face_areas[fno] *= -1.0f;
      vec_mult(norm, -1.0f);
      float *angles = vec_get(face_angles, fno);
      vec_mult(angles, -1.0f);
    }

    float area = face_areas[fno];
    if (area >= 0) {
      total_area += area;
    } else {
      neg_area -= area;
      // NOTE: I don't think neg_orig_area is ever even used
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
      neg_orig_area += fNorm->orig_area;
    }

    #undef total_area
    #undef neg_area
    #undef neg_orig_area
  #include "romp_for_end.h"
  // ROMP: end partial sum loop

  // transfer values back to MRIS
  mris->total_area = total_area;
  mris->neg_area = neg_area;
  mris->neg_orig_area = neg_orig_area;

  for (int vno = 0; vno < nvertices; vno++) {

    if (vertex_ripped[vno]) continue;

    VERTEX *vertex = &mris->vertices[vno];
    float *vnorm = vec_get(vertex_normals, vno);
    if (vnorm[2] < 0.0f) {
      vnorm[2] *= -1.0f;
      vertex->neg = 1;
    } else {
      vertex->neg = 0;
    }

    int *neighbor_faces = vertex_faces[vno];
    int nfaces = vertex_nfaces[vno];
    float area = 0.0f;
    for (int n = 0; n < nfaces; n++) area += face_areas[neighbor_faces[n]];
    
    // divide total face area, since each face is added to three vertices
    vertex_areas[vno] = area / vertex_area_fix_value;
  }
}


/*
  Looks for face normals that point inward in the sphere.
*/
void MetricPropertiesCalculator::orientEllipsoid()
{
  double total_area = 0.0;
  double neg_area = 0.0;
  double neg_orig_area = 0.0;

  // ROMP: begin partial sum loop
  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             nfaces
  #define ROMP_SUMREDUCTION0  total_area
  #define ROMP_SUMREDUCTION1  neg_area
  #define ROMP_SUMREDUCTION2  neg_orig_area
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    #define total_area      ROMP_PARTIALSUM(0)
    #define neg_area        ROMP_PARTIALSUM(1)
    #define neg_orig_area   ROMP_PARTIALSUM(2)
    
    if (face_ripped[fno]) ROMP_PFLB_continue

    // give the area an orientation: if the unit normal is pointing
    // inwards on the ellipsoid then the area should be negative
    
    int *vnos = vec_get(faces, fno);
    float *pos0 = vec_get(vertices, vnos[0]);
    float *pos1 = vec_get(vertices, vnos[1]);
    float *pos2 = vec_get(vertices, vnos[2]);
    float *norm = vec_get(face_normals, fno);

    // no need to divide by 3 since we only use the magnitude of the dot product
    float const xc = (pos0[0] + pos1[0] + pos2[0]) * norm[0];
    float const yc = (pos0[1] + pos1[1] + pos2[1]) * norm[1];
    float const zc = (pos0[2] + pos1[2] + pos2[2]) * norm[2];
    float const dot = xc + yc + zc;
    
    if (dot < 0.0f) {
      // normal faces inward, so reverse face area, normal, and angles
      face_areas[fno] *= -1.0f;
      vec_mult(norm, -1.0f);
      float *angles = vec_get(face_angles, fno);
      vec_mult(angles, -1.0f);
    }
    
    float area = face_areas[fno];
    if (area >= 0.0f) {
      total_area += area;
    } else {
      neg_area -= area;
      // NOTE: I don't think neg_orig_area is ever even used
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
      neg_orig_area += fNorm->orig_area;
    }

    #undef total_area
    #undef neg_area
    #undef neg_orig_area
  #include "romp_for_end.h"
  // ROMP: end partial sum loop

  // transfer values back to MRIS
  mris->total_area = total_area;
  mris->neg_area = neg_area;
  mris->neg_orig_area = neg_orig_area;
}


/*
  An even faster method for computing metric properties that also supports
  vertex normals calculated by weighted face-angles.
*/
void MRIScomputeMetricPropertiesFaster(MRIS* mris)
{
  MetricPropertiesCalculator calculator = MetricPropertiesCalculator(mris);
  calculator.compute();
  calculator.transfer();
}
