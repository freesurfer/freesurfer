#ifndef FSMATRIX_H
#define FSMATRIX_H

//#include <iostream>
//#include <vector>
//using namespace std;

//#include <armadillo>
//#include "matrix.h"

template <class T, class MatClass>
class FSMatrix
{
private:
  // Data members of FSMatrix
  MatClass *mat;
  
public:
  int n_rows;
  int n_cols;
  int n_elem;

  FSMatrix();
  FSMatrix(int rows, int cols);
  FSMatrix(const FSMatrix& oth);

  // copy assignment
  FSMatrix& operator=(const FSMatrix& oth);

  // obtain the raw memory pointer to element data
  T* memptr();

  // getter/setter methods
  T  elem(int row, int col) const;       // Returns the element at row r and column c
  T& elem(int row, int col);             // Sets the element at row r and column c
  
  // transpose
  FSMatrix t();

  // inverse
  FSMatrix inv();

  // determinant
  float det();
  
  void zeros(int rows, int cols);
  void ones(int rows, int cols);
  void fill(T k);
  //void rand48(int rows, int cols);  // .randu(), .randn()

  void print();
  int save(const char *fname);
  int load(const char *fname);

  // matrix operations
  // FSMatrix * FSMatrix  
  FSMatrix operator*(const FSMatrix& oth);
  // FSMatrix * scalar  
  FSMatrix operator*(T scalar);
  // FSMatrix + FSMatrix
  FSMatrix operator+(const FSMatrix& oth);
  // FSMatrix + scalar
  FSMatrix operator+(T scalar);
  // FSMatrix - FSMatrix
  FSMatrix operator-(const FSMatrix& oth);
};

#endif // FSMARITX_H

