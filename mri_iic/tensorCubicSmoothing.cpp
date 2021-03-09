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

#include <iostream>
#include <math.h>

#include "tensorCubicSmoothing.h"

//using namespace ctl;

TensorCubicSmoothing::TensorCubicSmoothing()
{
  //Nothing Needed
}
TensorCubicSmoothing::~TensorCubicSmoothing()
{
  //Nothing Needed
}

//Copy Constructor
TensorCubicSmoothing::TensorCubicSmoothing(const TensorCubicSmoothing& other)
{
  AtWA = other.AtWA;
  AtWr = other.AtWr;
  coefficients = other.coefficients;
}
// Assignment Overload
TensorCubicSmoothing& TensorCubicSmoothing::operator=(const TensorCubicSmoothing& other)
{
  AtWA = other.AtWA;
  AtWr = other.AtWr;
  coefficients = other.coefficients;
  return *this;
}

// takes in separable basis function in three dimensions and weights
// and constructs the matrix AtWA by exploiting the seperability
int TensorCubicSmoothing::constructAtWA(const vnl_matrix<float> &B2x,
                                        const vnl_matrix<float> &B2y,
                                        const vnl_matrix<float> &B2z,
                                        const Matrix3d &W,
                                        const vnl_vector<int> &indexMap)
{
  int X = B2x.cols();
  int Y = B2y.cols();
  int Z = B2z.cols();
  int size = (int)sqrt(X * Y * Z);
  Matrix3d temp = doBasisMultiplications(W, B2x, B2y, B2z);
  AtWA.set_size(size, size);

  //move to temporary 1 dimensional form
  vnl_vector<float> intermediate(size * size);
  int count = 0;
  for(int i = 0; i < Z; i++)
    for(int j = 0; j < Y; j++)
      for(int k = 0; k < X; k++)
      {
        intermediate(count++) = temp.getVal(j, k, i);
      }

  //put in 2d array using index map
  count = 0;
  for(int i = 0; i < size; i++)
    for(int j = 0; j < size; j++)
    {
      AtWA(j, i) = intermediate(indexMap(count++));
    }

  return 0;
}

// takes in separable basis function in three dimensions and weights
// and constructs the matrix AtWr by exploiting the seperability
int TensorCubicSmoothing::constructAtWr(const vnl_matrix<float> &Bx,
                                        const vnl_matrix<float> &By,
                                        const vnl_matrix<float> &Bz,
                                        const Matrix3d          &W,
                                        const Matrix3d          &r)
{

  int X = Bx.cols();
  int Y = By.cols();
  int Z = Bz.cols();
  Matrix3d Wr = W;

  // multiply each element of W by corresponding in r
  for(int i = 0; i < Z; i++)
    for(int j = 0; j < Y; j++)
      for(int k = 0; k < X; k++)
      {
        Wr.setVal(j, k, i, Wr.getVal(j, k, i) * r.getVal(j, k, i));
      }

  Matrix3d temp = doBasisMultiplications(Wr, Bx, By, Bz);

  AtWr.set_size(X * Y * Z);
  int count = 0;
  for(int i = 0; i < Z; i++)
    for(int j = 0; j < Y; j++)
      for(int k = 0; k < X; k++)
      {
        AtWr(count++) = temp.getVal(j, k, i);
      }


  //map into 1d vector
  return 0;
}

// solves the non-linear system by means of least squares given regularization parameters
// TODO is P 3 or 2 dimensions?
int TensorCubicSmoothing::solve(const Matrix3d &P, float lambda)
{
  //STUB
  return 0;
}
// solves the non-linear system by means of least squares without regularization
int TensorCubicSmoothing::solve()
{
  //STUB
  return 0;
}

// returns the solved coefficients
void TensorCubicSmoothing::getCoefficients(vnl_vector<float> &c) const
{
  c = coefficients;
}
// returns the smoothed data
void TensorCubicSmoothing::expandCoefficients(Matrix3d &d,
    const vnl_matrix<float> &Bx,
    const vnl_matrix<float> &By,
    const vnl_matrix<float> &Bz) const
{
  //STUB
}

// multiplies the 3d data matrix by basis functions in each direction
Matrix3d TensorCubicSmoothing::doBasisMultiplications(const Matrix3d &data,
    const vnl_matrix<float> &Bx,
    const vnl_matrix<float> &By,
    const vnl_matrix<float> &Bz)
{
  int rows = data.getHeight();
  //int cols = data.getWidth();
  int slices = data.getDepth();
  int numBx = Bx.columns();
  int numBy = By.columns();
  int numBz = Bz.columns();

  vnl_matrix<float> temp1(rows, numBy * slices);

  // multiplies each slice by By and places it in a temporary array shaped
  // such that it can be multiplied by Bxt
  std::cerr << "data columns: " << (data.getSlice(0)).columns() << "\nBy rows: " << By.rows() << '\n';
  for(int i = 0; i < slices; i++)
  {
    std::cerr << i <<"th slice being attempted\n";
    temp1.set_columns(i * numBy, data.getSlice(i) * By);
    std::cerr << "slice complete\n";
  }
  std::cerr << "Done By\n";

  vnl_matrix<float> Bxt = Bx.transpose();
  // multiply by Bxt
  vnl_matrix<float> temp2 = Bxt * temp1;

  temp1.set_size(numBx * numBy, slices);
  // reshape matrix for multiplication in z direction, this takes more effort
  // on our part because the data is not adjacent in our data structure
  for(int i = 0; i < slices; i++)
  {
    int count = 0;
    for(int j = 0; j < numBy; j++)
      for(int k = 0; k < numBx; k++)
      {
        temp1(count++, i) = temp2(k, j + i * numBy);
      }
  }
// multiply by Bz
  temp2 = temp1 * Bz;

  Matrix3d result(numBy, numBx, numBz);

  // reshape back into a 3 dimensional array
  for(int i = 0; i < numBz; i++)
  {
    int count = 0;
    for(int j = 0; j < numBx; j++)
      for(int k = 0; k < numBy; k++)
      {
        result.setVal(k, j, i, temp2(count++, i));
      }
  }
  return result;
}
