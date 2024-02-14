#pragma once

#include <vector>

#include "mrisurf.h"
#include "mrishash.h"


/*
  A simple 3D vector class for easy operations while testing for ray-face intersections.
  This should be substituted for something for legit once we adopt a c++ matrix library.
*/
class Vec3
{
public:

  Vec3() {
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
  }

  Vec3(const double x, const double y, const double z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  inline double& operator[](int index) { return v[index]; }
  inline const double& operator[](int index) const { return v[index]; }

  inline Vec3& operator=(const Vec3 &vect) {
    v[0] = vect[0];
    v[1] = vect[1];
    v[2] = vect[2];
    return *this;
  }

  inline Vec3 operator+(const Vec3 &vect) const { return Vec3(v[0] + vect[0], v[1] + vect[1], v[2] + vect[2]); }
  inline Vec3 operator-(const Vec3 &vect) const { return Vec3(v[0] - vect[0], v[1] - vect[1], v[2] - vect[2]); }
  inline double operator*(const Vec3 &vect) const { return v[0] * vect[0] + v[1] * vect[1] + v[2] * vect[2]; }
  inline Vec3 operator*(const double scalar) const { return Vec3(v[0] * scalar, v[1] * scalar, v[2] * scalar); }

private:
  double v[3];
};

inline Vec3 operator*(const double scalar, const Vec3 vect) { return vect * scalar; }

inline float dot(const Vec3 &u, const Vec3 &v) { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; }

inline Vec3 cross(const Vec3 &u, const Vec3 &v) {
  return Vec3(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]);
}


/*
  A utility for interpolating overlay values on sphere at given (spherical) coordinates.
*/
class SphericalInterpolator
{
public:
  SphericalInterpolator(MRIS *surf, int which=CURRENT_VERTICES);
  ~SphericalInterpolator() { MHTfree(&mht); }

  float interp(double phi, double theta);
  void setOverlay(const float *array);
  bool nearestneighbor = false;
  bool testRayIntersection(int fno, float x, float y, float z, float *value, double *w=NULL, bool interp=true);

private:
  MRIS *mris;
  MRIS_HASH_TABLE *mht;

  double radius;
  std::vector<Vec3> vertices;
  std::vector<float> overlay;

  bool searchBucket(int bx, int by, int bz, float x, float y, float z, float *value);
};

int MHTfindClosestFaceSph(MRIS *surf, MRIS_HASH_TABLE *mht, SphericalInterpolator *si, 
			  double *cxyz, double *w=NULL, int debug=0);
