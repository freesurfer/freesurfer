#ifndef FSMATRIX_ARMA_H
#define FSMATRIX_ARMA_H

#include <iostream>
using namespace std;

#include <armadillo>

#include "FSMatrix.h"


/* To be implemented:
 * is_zero(), is_empty(), is_identity()
 * resize(), reshape()
 */

template <class T>
class FSMatrix<T, arma::Mat<T>>;

// partial specialization for armadillo matrix
template <class T>
class FSMatrix<T, arma::Mat<T>>
{
private:
  arma::Mat<T> *mat;

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

// default constructor
template <class T>
FSMatrix<T, arma::Mat<T>>::FSMatrix()
{
  // Initialization of data members
  //printf("[DEBUG] FSMatrix<T, arma::Mat<T>>::FSMatrix()\n");
  mat = new arma::Mat<T>();
    
  n_rows = 0;
  n_cols = 0;
  n_elem = 0;
}

// constructor
template <class T>
FSMatrix<T, arma::Mat<T>>::FSMatrix(int rows, int cols)
{
  // Initialization of data members
  //printf("[DEBUG] FSMatrix<T, arma::Mat<T>>::FSMatrix(%d, %d)\n", rows, cols);
  mat = new arma::Mat<T>(rows, cols);
    
  n_rows = rows;
  n_cols = cols;
  n_elem = n_rows * n_cols;
}

// copy constructor
template <class T>
FSMatrix<T, arma::Mat<T>>::FSMatrix(const FSMatrix& oth)
{
  //printf("[DEBUG] copy constructor FSMatrix<T, arma::Mat<T>>::FSMatrix(const FSMatrix& oth)\n");
  
  mat = new arma::Mat<T>();
  *mat = (*oth.mat);
  n_rows = oth.n_rows;
  n_cols = oth.n_cols;
  n_elem = oth.n_elem;
}

// copy assignment
template <class T>
FSMatrix<T, arma::Mat<T>>& FSMatrix<T, arma::Mat<T>>::operator=(const FSMatrix& oth)
{
  //printf("[DEBUG] copy assignment FSMatrix<T, arma::Mat<T>>:::operator=(const FSMatrix& oth)\n");
  // Guard self assignment
  if (this == &oth)
    return *this;
    
  *mat = (*oth.mat);
  n_rows = oth.n_rows;
  n_cols = oth.n_cols;
  n_elem = oth.n_elem;

  return *this;
}

// obtain the raw memory pointer to element data
template <class T>
T* FSMatrix<T, arma::Mat<T>>::memptr()
{
  return mat->memptr();
}

// getter/setter methods
// Returns the element at row r and column c
template <class T>
T FSMatrix<T, arma::Mat<T>>::elem(int row, int col) const
{
  return mat->at(row, col);
}


// Sets the element at row r and column c
template <class T>
T& FSMatrix<T, arma::Mat<T>>::elem(int row, int col)
{
  return mat->at(row, col);
}

// transpose
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix <T, arma::Mat<T>>::t()
{
  FSMatrix<T, arma::Mat<T>> trans(n_cols, n_rows);
  trans.mat = mat->t();
  return trans;
}

// inverse
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix <T, arma::Mat<T>>::inv()
{
  return mat->inv();
}

// determinant
template <class T>
float FSMatrix<T, arma::Mat<T>>::det()
{
  return mat->det();
}

// set all elements to zero
template <class T>
void FSMatrix <T, arma::Mat<T>>::zeros(int rows, int cols)
  {
    mat->zeros(rows, cols);
    n_rows = rows;
    n_cols = cols;
  }

// set all elements to one
template <class T>
void FSMatrix <T, arma::Mat<T>>::ones(int rows, int cols)
{
  mat->ones(rows, cols);    
  n_rows = rows;
  n_cols = cols;
}

  // set all elements to be equal to k
template <class T>
void FSMatrix<T, arma::Mat<T>>::fill(T k)
{
  mat->fill(k);
}

#if 0
// armadillo has .randu(rows, cols) and .randn(rows, cols)
// ???
template <class T>
void FSMatrix <T, arma::Mat<T>>::rand48(int rows, int cols)
{
  // ???
  n_rows = rows;
  n_cols = cols;
}
#endif

// methods for matrix operations
// FSMatrix * FSMatrix
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix<T, arma::Mat<T>>::operator*(const FSMatrix<T, arma::Mat<T>>& oth)
{
  FSMatrix<T, arma::Mat<T>> res;
  
  *res.mat = (*mat) * (*oth.mat);
  res.n_rows = res.mat->n_rows;
  res.n_cols = res.mat->n_cols;
  res.n_elem = res.mat->n_elem;

  return res;
}

// FSMatrix * scalar
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix<T, arma::Mat<T>>::operator*(T scalar)
{
  FSMatrix<T, arma::Mat<T>> res;
  
  *res.mat = (*mat) * scalar;  
  res.n_rows = res.mat->n_rows;
  res.n_cols = res.mat->n_cols;
  res.n_elem = res.mat->n_elem;

  return res;
}

// FSMatrix + FSMatrix
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix<T, arma::Mat<T>>::operator+(const FSMatrix<T, arma::Mat<T>>& oth)
{
  FSMatrix<T, arma::Mat<T>> res;

  *res.mat = (*mat) + (*oth.mat);
  res.n_rows = res.mat->n_rows;
  res.n_cols = res.mat->n_cols;
  res.n_elem = res.mat->n_elem;

  return res;
}

// FSMatrix + scalar
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix<T, arma::Mat<T>>::operator+(T scalar)
{
  FSMatrix<T, arma::Mat<T>> res;
  
  *res.mat = (*mat) + scalar;  
  res.n_rows = res.mat->n_rows;
  res.n_cols = res.mat->n_cols;
  res.n_elem = res.mat->n_elem;

  return res;
}

// FSMatrix - FSMatrix
template <class T>
FSMatrix<T, arma::Mat<T>> FSMatrix<T, arma::Mat<T>>::operator-(const FSMatrix<T, arma::Mat<T>>& oth)
{
  FSMatrix<T, arma::Mat<T>> res;

  *res.mat = (*mat) - (*oth.mat);
  res.n_rows = res.mat->n_rows;
  res.n_cols = res.mat->n_cols;
  res.n_elem = res.mat->n_elem;

  return res;
}

// print elements to the stdout
template <class T>
void FSMatrix<T, arma::Mat<T>>::print()
{
  mat->print();
}

// store matrix in the specified file
template <class T>
int FSMatrix<T, arma::Mat<T>>::save(const char *fname)
{
  mat->save(fname, arma::arma_ascii);
  return 0;
}

// retrieve matrix from the specified file
template <class T>
int FSMatrix<T, arma::Mat<T>>::load(const char *fname)
{
  mat->load(fname, arma::arma_ascii);
  return 0;
}  

#endif // FSMARITX_ARMA_H

