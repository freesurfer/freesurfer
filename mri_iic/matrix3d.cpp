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

#include <vnl/vnl_matrix.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "matrix3d.h"

// default constructor, creates matrix with all dimentions 0 and NULL data
Matrix3d::Matrix3d()
{
  data = NULL;
  width = 0;
  height = 0;
  depth = 0;
}

// constructor that creates the matrix with specified dimensions
Matrix3d::Matrix3d(int rows, int cols, int slices)
{
  data = new vnl_matrix<float>* [slices];
  for(int i = 0; i < slices; i++)
  {
    data[i] = new vnl_matrix<float> (rows, cols);
  }
  width = cols;
  height = rows;
  depth = slices;
}

// destructor, deletes the array of matricies and all matricies within
Matrix3d::~Matrix3d()
{
  for(int i = 0; i < depth; i++)
  {
    std::cerr << i << '\n';
    delete data[i];
  }
  delete [] data;
  data = NULL;
}

// copy constructor, creates exact replica of other
Matrix3d::Matrix3d(const Matrix3d& other)
{
  width = other.width;
  height = other.height;
  depth = other.depth;

  data = new vnl_matrix<float>* [depth];
  for(int i = 0; i < depth; i++)
  {
    data[i] = new vnl_matrix<float> (height, width);
    data[i] = other.data[i];
  }
}

// copies information from other matrix, overwriting old data
Matrix3d& Matrix3d::operator=(const Matrix3d& other)
{
  if (this != &other)
  {
    // delete old matrix
    if (data != NULL)
    {
      for(int i = 0; i < depth; i++)
      {
        delete data[i];
      }
      delete [] data;
    }
    //find new dimensions
    width = other.width;
    height = other.height;
    depth = other.depth;

    //create new matrix and copy data
    data = new vnl_matrix<float>* [depth];
    for(int i = 0; i < depth; i++)
    {
      data[i] = new vnl_matrix<float> (height, width);
      data[i] = other.data[i];
    }
  }
  return *this;
}

// returns the value at the given coordinate
float Matrix3d::getVal(int row, int col, int slice) const
{
  if(row < 0 || row >= height || col < 0 || col >= width || slice < 0 || slice >= depth)
  {
    std::cerr << "ERROR: index out of bounds\n"; //Maybe Throw exeption instead
    exit(1);
  }
  return data[slice]->get(row, col);
}

// sets the given coordinate to the value provided
void Matrix3d::setVal(int row, int col, int slice, float const &val)
{
  if(row < 0 || row >= height || col < 0 || col >= width || slice < 0 || slice >= depth)
  {
    std::cerr << "ERROR: index out of bounds\n"; //Maybe Throw exeption instead
    exit(1);
  }

  data[slice]->put(row, col, val);
}

// returns the width of the matrix (number of columns)
int Matrix3d::getWidth() const
{
  return width;
}

// returns the height of the matrix (number of rows)
int Matrix3d::getHeight() const
{
  return height;
}

// returns the depth of the matrix (number of slices)
int Matrix3d::getDepth() const
{
  return depth;
}

vnl_matrix<float> const& Matrix3d::getSlice(int sliceNum) const
{
  return *data[sliceNum];
}

/**********************DEBUGGING***************************/
//prints matrix slices to terminal
void Matrix3d::print()
{
  std::cout << "PRINTING MATRIX:\n\n";
  for (int i = 0; i < depth; i++)
  {
    std::cout << "Slice: " << i << std::endl;
    for (int j = 0; j < height; j++)
    {
      for (int k = 0; k < width; k++)
      {
        std::cout << std::setw(5)    << std::setprecision(3)
                  << getVal(j, k, i) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}
