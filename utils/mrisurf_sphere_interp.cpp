#include "mrisurf_sphere_interp.h"
#include "mrishash_internals.h"


/*
  Constructs a interpolator cache from a spherical surface and its
  internal vertex curv values.
*/
SphericalInterpolator::SphericalInterpolator(MRIS *surf) : mris(surf)
{
  // build the hash table for fast lookup
  mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, mris->avg_vertex_dist / 2);

  // cache the radius for spherical to cartesian translations
  radius = MRISaverageRadius(mris);

  // store vertex positions and curv values as seperate lists for faster lookup
  vertices = std::vector<Vec3>(mris->nvertices);
  overlay = std::vector<float>(mris->nvertices);
  for (int vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX *v = &mris->vertices[vno];
    vertices[vno] = Vec3(v->x, v->y, v->z);
    overlay[vno] = v->curv;
  }
}


/*
Set the interpolator overlay values.
*/
void SphericalInterpolator::setOverlay(const float *array)
{
  for (int vno = 0 ; vno < mris->nvertices ; vno++) overlay[vno] = array[vno];
}


/*
  Searches for the face containing the spherical point defined by (phi, theta) and
  returns the (barycentric) interpolated value at that point.
*/
float SphericalInterpolator::interp(double phi, double theta)
{
  // translate from spherical to cartesian xyz
  double x = radius * sin(phi) * cos(theta);
  double y = radius * sin(phi) * sin(theta);
  double z = radius * cos(phi);

  // translate xyz point to hash coordinate
  float res = mht->vres();
  double hashx = (x / res) + TABLE_CENTER;
  double hashy = (y / res) + TABLE_CENTER;
  double hashz = (z / res) + TABLE_CENTER;

  // get hash index
  int ix = int(hashx);
  int iy = int(hashy);
  int iz = int(hashz);

  // get mod for rounding
  double modx = hashx - double(ix);
  double mody = hashy - double(iy);
  double modz = hashz - double(iz);

  // compute offset depending on rounding
  int offsetx = (modx <= 0.5) ? -1 : 0;
  int offsety = (mody <= 0.5) ? -1 : 0;
  int offsetz = (modz <= 0.5) ? -1 : 0;

  // initialize checklist for central 27 voxels
  unsigned char central27[3][3][3];
  memset(central27, 0, sizeof(central27));

  // return value
  float value = 0.0;

  // search home bucket and the surrounding 7 buckets
  for (int xvi = 0; xvi <= 1; xvi++) {
    int xv = xvi + offsetx;
    for (int yvi = 0; yvi <= 1; yvi++) {
      int yv = yvi + offsety;
      for (int zvi = 0; zvi <= 1; zvi++) {
        int zv = zvi + offsetz;
        if (searchBucket(ix + xv, iy + yv, iz + zv, x, y, z, &value)) return value;
        central27[xv + 1][yv + 1][zv + 1] = 1;
      }
    }
  }

  // search the rest of the central 27 buckets
  for (int xv = -1; xv <= 1; xv++) {
    for (int yv = -1; yv <= 1; yv++) {
      for (int zv = -1; zv <= 1; zv++) {
        // skip ones already done
        if (!central27[xv + 1][yv + 1][zv + 1]) {
          if (searchBucket(ix + xv, iy + yv, iz + zv, x, y, z, &value)) return value;
          central27[xv + 1][yv + 1][zv + 1] = 1;
        }
      }
    }
  }

  // search further-away buckets
  for (int RVox = 2; RVox <= 6; RVox++) {
    // jump from one side to the other across the empty middle
    int WallJump = RVox + RVox;
    for (int xv = -RVox; xv <= RVox; xv++) {
      for (int yv = -RVox; yv <= RVox; yv++) {
        bool isWall = ((xv == -RVox) || (xv == RVox)) || ((yv == -RVox) || (yv == RVox));
        for (int zv = -RVox; zv <= RVox; zv = isWall ? zv + 1 : zv + WallJump) {
          if (searchBucket(ix + xv, iy + yv, iz + zv, x, y, z, &value)) return value;
        }
      }
    }
  }

  fs::error() << "Could not locate a face corresponding to phi=" << phi << ", theta=" << theta
              << ". Are you sure this surface is a sphere?";
  return 0.0;
}


/*
  Searches the bucket defined by (bx, by, bz) for a face that intersects with the
  infinite ray defined by (x, y, z). If found, interpolates the face value at the point of
  intersection and returns true.
*/
bool SphericalInterpolator::searchBucket(int bx, int by, int bz, float x, float y, float z, float *value)
{
  // checkout the bucket
  MHBT* bucket = MHTacqBucketAtVoxIx(mht, bx, by, bz);
  if (!bucket) return false;

  bool found = false;
  MHB *bin = bucket->bins;
  for (int faceix = 0; faceix < bucket->nused; faceix++, bin++) {
    // check for ray-face intersection and interpolate if found
    if (testRayIntersection(bin->fno, x, y, z, value)) {
      found = true;
      break;
    }
  }

  // be sure to release the bucket
  MHTrelBucket(&bucket);
  return found;
}


/*
  Tests for intersection (using Möller–Trumbore algorithm) between face fno and the infinite
  ray defined by (x, y, z). If found, interpolates the face value at the point of
  intersection and returns true.
*/
bool SphericalInterpolator::testRayIntersection(int fno, float x, float y, float z, float *value)
{
  FACE *face = &mris->faces[fno];
  Vec3 vertex0 = vertices[face->v[0]];
  Vec3 vertex1 = vertices[face->v[1]];
  Vec3 vertex2 = vertices[face->v[2]];
  Vec3 ray = Vec3(x, y, z);

  Vec3 edge1 = vertex1 - vertex0;
  Vec3 edge2 = vertex2 - vertex0;

  Vec3 h = cross(ray, edge2);

  float det = dot(edge1, h);
  if (fabs(det) < 1e-6) return false;

  float f = 1.0 / det;
  Vec3 s = Vec3(0, 0, 0) - vertex0;
  float u = f * dot(s, h);
  if (u < 0.0 || u > 1.0) return false;
  
  Vec3 q = cross(s, edge1);
  float v = f * dot(ray, q);
  if (v < 0.0 || u + v > 1.0) return false;

  float l = 1.0f - v - u;

  // if the ray intersects, do interpolation here as we've already computed barycentric coordinates
  *value = overlay[face->v[0]] * l +
           overlay[face->v[1]] * u +
           overlay[face->v[2]] * v;

  return true;
}