#pragma once

#include "mrisurf.h"
#include "mrisurf_sphere_interp.h"


/*
  Very simple 2D array utility. We really need to settle on a
  universal matrix library...
*/
template <typename T>
class ImageArray
{
public:

  ImageArray() {};

  ImageArray(int w, int h, T fill = 0) {
    xl = w;
    yl = h;
    arr = std::vector<T>(w * h, fill);
  };

  T& item(int x, int y) { return arr[x * yl + y]; };

  int xlen() { return xl; };
  int ylen() { return yl; };

private:
  int xl;
  int yl;
  std::vector<T> arr;
};


/*
  Projector class to (more) easily map spherical overlays and
  parameterization images. This is an initial attempt to cleanup the
  parameterization code a bit.
*/
class SphericalProjector
{
public:

  enum InterpMethod { Barycentric, Nearest }; 

  SphericalProjector(MRIS *surf, MRI_SP *param);
  ~SphericalProjector();

  void parameterizeOverlay(const float* overlay, int frameno, InterpMethod interp);
  void sampleParameterization(float* overlay, int frameno, InterpMethod interp);

private:
  MRIS *original;
  MRIS *mris;
  SphericalInterpolator *interpolator;

  MRI_SP *mrisp;
  int udim;
  int vdim;

  std::vector<int> vertex_u;
  std::vector<int> vertex_v;
  std::vector<float> vertex_uf;
  std::vector<float> vertex_vf;

  ImageArray<int> hits;
  ImageArray<int> nearest;
  ImageArray<float> distance;
};
