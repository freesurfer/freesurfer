//testMatrix.cpp
//
//Created 7/23/2014
//By: Benjamin Lewin
//
//A test file for the class Matrix3d
//


#include <iostream>
#include "matrix3d.h"


int main()
{
  int rows = 0;
  int cols = 0;
  int slices = 0;
  std::cout << "Enter matrix dimensions:\n" << "Rows: ";
  std::cin >> rows;
  std::cout << "Columns: ";
  std::cin >> cols;
  std::cout << "Slices: ";
  std::cin >> slices;

  Matrix3d matrix0(rows, cols, slices);
  Matrix3d matrix;
//  Matrix3d mempty;
//  matrix = mempty;
  matrix = matrix0;

  std::cout << "\nDIMENSIONS:\n";
  std::cout << "Height: " << matrix.getHeight() << std::endl;
  std::cout << "Width: " << matrix.getWidth() << std::endl;
  std::cout << "Depth: " << matrix.getDepth() << std::endl;

  //matrix.setVal(1,1,2,5.595);
  //std::cout << matrix.getVal(0,1,2) <<std::endl;

  matrix.print();
}

