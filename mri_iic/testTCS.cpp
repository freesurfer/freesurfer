//testTCS.cpp
//
//Created 8/7/14
//By: Benjamin Lewin
//

#include "tensorCubicSmoothing.h"
#include "matrix3d.h"
#include <vnl/vnl_matrix.h>

const int ROWS   = 3;
const int COLS   = 4;
const int SLICES = 5;
const int NUM_BX = 4;
const int NUM_BY = 5;
const int NUM_BZ = 3;

int main()
{

  vnl_matrix<float> Bx(ROWS, NUM_BX);
  vnl_matrix<float> By(COLS, NUM_BY);
  vnl_matrix<float> Bz(SLICES, NUM_BZ);

  Matrix3d data(ROWS, COLS, SLICES);

  TensorCubicSmoothing test;

  test.doSneakySeperabilityMultiplications(data, Bx, By, Bz);

  return 0;
}


