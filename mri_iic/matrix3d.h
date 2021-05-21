/**
 * @brief wrapper for 3d vnl_matrix operations
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

#ifndef MATRIX3D_INCLUDED
#define MATRIX3D_INCLUDED

#include <vnl/vnl_matrix.h>

class Matrix3d
{
public:
  // default constructor, creates matrix with all dimentions 0 and NULL data
  Matrix3d();
  // constructor that creates the matrix with specified dimensions
  Matrix3d(int rows, int cols, int slices);
  // destructor, deletes the array of matricies and all matricies within
  ~Matrix3d();
  // copy constructor, creates exact replica of other
  Matrix3d(const Matrix3d& other);
  // copies information from other matrix, overwriting old data
  Matrix3d& operator=(const Matrix3d& other);

  // returns the value at the given coordinate
  float getVal(int row, int col, int slice) const;
  // sets the given coordinate to the value provided
  void setVal(int row, int col, int slice, float const &val);

  // returns the width of the matrix (number of columns)
  int getWidth() const;
  // returns the height of the matrix (number of rows)
  int getHeight() const;
  // returns the depth of the matrix (number of slices)
  int getDepth() const;

  vnl_matrix<float> const& getSlice(int sliceNum) const;
  // DEBUG FUNCTIONS

  //prints matrix slices to terminal
  void print();

private:

  //Array of pointers to 2d matricies
  vnl_matrix<float> **data;

  int width;
  int height;
  int depth;
};

#endif
