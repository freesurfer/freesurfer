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
