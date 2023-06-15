#ifndef FSMATRIX_MATRIX_H
#define FSMATRIX_MATRIX_H

#include <iostream>
using namespace std;

#include "FSMatrix.h"
#include "matrix.h"

/* To be implemented:
 * is_zero(), is_empty(), is_identity()
 * resize(), reshape()
 */

template <>
class FSMatrix <float, MATRIX>;

// specialization for <float, MATRIX>
// ??? another class for <complex_float, MATRIX> ???
template <>
class FSMatrix <float, MATRIX>
{
private:
  MATRIX *mat;  // MATRIX indexing starts at 1,
                // index will be adjusted internally so same interface for caller

public:
  int n_rows;
  int n_cols;
  int n_elem;

  FSMatrix();
  FSMatrix(int rows, int cols);
  FSMatrix(const FSMatrix& oth);

  // copy assignment
  FSMatrix& operator=(const FSMatrix& oth);
  
  ~FSMatrix();

  // obtain the raw memory pointer to element data
  float* memptr();
  
  // getter/setter methods
  float  elem(int row, int col) const;       // Returns the element at row r and column c
  float& elem(int row, int col);             // Sets the element at row r and column c
  
  // transpose
  FSMatrix t();

  // inverse
  FSMatrix inv();

  // determinant
  float det();

  void zeros(int rows, int cols);
  void ones(int rows, int cols);
  void fill(float k);
  //void rand48(int rows, int cols);

  void print(FILE *fp=stdout);
  int save(const char *fname);
  int load(const char *fname);

  // matrix operations
  // FSMatrix * FSMatrix
  FSMatrix operator*(const FSMatrix& oth);
  // FSMatrix * scalar
  FSMatrix operator*(const float scalar);
  // FSMatrix + FSMatrix
  FSMatrix operator+(const FSMatrix& oth);
  // FSMatrix + scalar
  FSMatrix operator+(const float scalar);
  // FSMatrix - FSMatrix
  FSMatrix operator-(const FSMatrix& oth);
};


// default constructor
FSMatrix<float, MATRIX>::FSMatrix()
{
  // Initialization of data members
  //printf("[DEBUG] FSMatrix<float, MATRIX>::FSMatrix()\n");
    
  n_rows = 0;
  n_cols = 0;
  n_elem = 0;

  mat = NULL;
}


// constructor
FSMatrix<float, MATRIX>::FSMatrix(int rows, int cols)
{
  // Initialization of data members
  //printf("[DEBUG] FSMatrix<float, MATRIX>::FSMatrix(%d, %d)\n", rows, cols);
  mat = MatrixAlloc(rows, cols, MATRIX_REAL);
    
  n_rows = rows;
  n_cols = cols;
  n_elem = n_rows * n_cols;
}


// copy constructor
FSMatrix<float, MATRIX>::FSMatrix(const FSMatrix& oth)
{
  // Initialization of data members
  //printf("[DEBUG] copy constructor FSMatrix<float, MATRIX>::FSMatrix(const FSMatrix& oth)\n");

  mat = MatrixCopy(oth.mat, NULL);
    
  n_rows = oth.n_rows;
  n_cols = oth.n_cols;
  n_elem = oth.n_elem;
}


// copy assignment
FSMatrix<float, MATRIX>& FSMatrix<float, MATRIX>::operator=(const FSMatrix& oth)
{
  // Initialization of data members
  //printf("[DEBUG] copy assignment FSMatrix<float, MATRIX>::operator=(const FSMatrix& oth)\n");

  // Guard self assignment
  if (this == &oth)
    return *this;

  if (mat != NULL)
    MatrixFree(&mat);

  mat = MatrixCopy(oth.mat, NULL);
    
  n_rows = oth.n_rows;
  n_cols = oth.n_cols;
  n_elem = oth.n_elem;

  return *this;  
}


// destructor
FSMatrix<float, MATRIX>::~FSMatrix()
{
  //printf("[DEBUG] FSMatrix<float, MATRIX>::~FSMatrix()\n");
  if (mat != NULL)
  {
    MatrixFree(&mat);
    mat = NULL;
  }
}


// obtain the raw memory pointer to element data
float* FSMatrix <float, MATRIX>::memptr()
{
  return mat->data;
}


// getter/setter methods
// Returns the element at row r and column c
float FSMatrix<float, MATRIX>::elem(int row, int col) const
{
  // MATRIX starts from index 1
  return *MATRIX_RELT(mat, row+1, col+1);
}


// Sets the element at row r and column c
float& FSMatrix<float, MATRIX>::elem(int row, int col)
{
  // MATRIX starts from index 1
  return *MATRIX_RELT(mat, row+1, col+1);
}


// transpose
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::t()
{
  FSMatrix trans(n_cols, n_rows);
  trans.mat = MatrixTranspose(mat, trans.mat);
  return trans;
}

// inverse
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::inv()
{
  FSMatrix trans(n_cols, n_rows);
  trans.mat = MatrixInverse(mat, trans.mat);
  return trans;
}

// determinant
float FSMatrix <float, MATRIX>::det()
{
  return MatrixDeterminant(mat);
}

// set all elements to zero
void FSMatrix <float, MATRIX>::zeros(int rows, int cols)
{
  n_rows = rows;
  n_cols = cols;
  mat =  MatrixConstVal(0, n_rows, n_cols, NULL); 
}

// set all elements to one
void FSMatrix <float, MATRIX>::ones(int rows, int cols)
{
  n_rows = rows;
  n_cols = cols;
  mat =  MatrixConstVal(1, n_rows, n_cols, NULL);
}

// set all elements to be equal to k
void FSMatrix <float, MATRIX>::fill(float k)
{
  mat =  MatrixConstVal(k, n_rows, n_cols, mat);
}

#if 0
// armadillo has .randu(rows, cols) and .randn(rows, cols) 
void FSMatrix <float, MATRIX>::rand48(int rows, int cols)
{
  n_rows = rows;
  n_cols = cols;
  mat = MatrixDRand48(n_rows, n_cols, mat);
}
#endif
  
// methods for matrix operations
// FSMatrix * FSMatrix
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::operator*(const FSMatrix& oth)
{
  FSMatrix res;
  
  res.mat = MatrixMultiplyD(mat, oth.mat, NULL);
  res.n_rows = (res.mat)->rows;
  res.n_cols = (res.mat)->cols;
  res.n_elem = n_rows * n_cols;
  
  return res;
}

// FSMatrix * scalar
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::operator*(const float scalar)
{
  FSMatrix res;
    
  res.mat = MatrixScalarMul(mat, scalar, NULL);
  res.n_rows = (res.mat)->rows;
  res.n_cols = (res.mat)->cols;
  res.n_elem = n_rows * n_cols;
    
  return res;
}

// FSMatrix + FSMatrix
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::operator+(const FSMatrix& oth)
{
  FSMatrix res;
  
  res.mat = MatrixAdd(mat, oth.mat, NULL);
  res.n_rows = (res.mat)->rows;
  res.n_cols = (res.mat)->cols;
  res.n_elem = n_rows * n_cols;
  
  return res;
}

// FSMatrix + scalar
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::operator+(const float scalar)
{
  FSMatrix res;
    
  res.mat = MatrixScalarAdd(mat, scalar, NULL);
  res.n_rows = (res.mat)->rows;
  res.n_cols = (res.mat)->cols;
  res.n_elem = n_rows * n_cols;
    
  return res;
}

// FSMatrix - FSMatrix
FSMatrix<float, MATRIX> FSMatrix<float, MATRIX>::operator-(const FSMatrix& oth)
{
  FSMatrix res;
  
  res.mat = MatrixSubtract(mat, oth.mat, NULL);
  res.n_rows = (res.mat)->rows;
  res.n_cols = (res.mat)->cols;
  res.n_elem = n_rows * n_cols;

  return res;
}

// print
void FSMatrix<float, MATRIX>::print(FILE *fp)
{
  MatrixPrint(fp, mat);
}

// store matrix in the specified file
int FSMatrix<float, MATRIX>::save(const char *fname)
{
  // MatrixWrite()/MatlabWrite(), MatrixWriteTxt(), MatrixWriteInto()
  // ??? what are the difference ???
  return MatrixAsciiWrite(fname, mat);  // call MatrixAsciiWriteInto() internally
}

// retrieve matrix from the specified file
int FSMatrix<float, MATRIX>::load(const char *fname)
{
  // MatrixRead()/MatlabRead(), MatrixReadTxt(), MatrixReadFrom(), MatrixAsciiReadRaw()
  // ??? what are the difference ???
  mat = MatrixAsciiRead(fname, mat);  // call MatrixAsciiReadFrom() internally
  return (mat == NULL) ? 1 : 0;
}

#endif // FSMARITX_MATRIX_H

