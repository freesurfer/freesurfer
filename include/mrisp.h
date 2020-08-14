#pragma once
/**
 * @brief utilities for spherical parameterization of an MRI_SURFACE
 *
 * mrisp = MRI_SURFACE parameterization contains utilities for writing various
 * fields over the surface (e.g. curvature, coordinate functions) into a
 * spherical (longitude/colatitude) parameterization in the form of a 2D
 * image. Also for Gaussian blurring and other utilities.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf.h"
#include "mrisurf_sphere_interp.h"


/*
  Helper class to (more) easily project overlays into an parameterization image.
  This is an initial attempt to cleanup the parameterization code a bit.
*/
class BarycentricSphericalProjector
{
public:
  BarycentricSphericalProjector(MRIS *mris, MRI_SP *param);
  ~BarycentricSphericalProjector();
  void projectOverlay(const float* overlay, int frameno = 0);

private:
  std::vector<int> vertex_u;
  std::vector<int> vertex_v;
  int **filled;
  float **distances;
  int u_max_index, v_max_index;

  MRIS *original = nullptr;
  MRIS *mris = nullptr;
  MRI_SP *mrisp = nullptr;

  SphericalInterpolator *interpolator = nullptr;
};
