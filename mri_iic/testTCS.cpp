/**
 * @brief test program of the TensorCubicSmoothing class
 *
 */
/*
 * Original Author: Benjamin Lewin
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

#include <vnl/vnl_matrix.h>

#include "tensorCubicSmoothing.h"
#include "matrix3d.h"

const int ROWS   = 400;
const int COLS   = 400;
const int SLICES = 400;
const int NUM_BX = 15;
const int NUM_BY = 15;
const int NUM_BZ = 15;

int main()
{

  vnl_matrix<float> Bx(ROWS, NUM_BX);
  vnl_matrix<float> By(COLS, NUM_BY);
  vnl_matrix<float> Bz(SLICES, NUM_BZ);

  Matrix3d data(ROWS, COLS, SLICES);

  TensorCubicSmoothing test;

  test.doBasisMultiplications(data, Bx, By, Bz);

  return 0;
}
