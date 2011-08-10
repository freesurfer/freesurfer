#ifndef _SBL_MATRIX_UTIL_H_
#define _SBL_MATRIX_UTIL_H_
#include <sbl/core/File.h>
#include <sbl/core/Array.h>
#include <sbl/core/Pointer.h>
#include <sbl/math/Matrix.h>
#include <sbl/math/Vector.h>
namespace sbl {


/*! \file MatrixUtil.h
    \brief The MatrixUtil module provides various functions for creating, manipulating,
    and analyzing matrices (represented by Matrix class instances).
*/


//-------------------------------------------
// MATRIX FILE I/O
//-------------------------------------------


/// write matrix (and dimensions) to file
void writeMatrixF( File &file, const MatrixF &mat );
void writeMatrixF( const String &fileName, const MatrixF &mat );
void writeSymMatrixF( const String &fileName, const MatrixF &mat );
void writeSparseMatrixF( const String &fileName, const MatrixF &mat );


/// read matrix (and dimensions) from file
aptr<MatrixF> readMatrixF( File &file );
aptr<MatrixF> readMatrixF( const String &fileName );
aptr<MatrixF> readSymMatrixF( const String &fileName );
aptr<MatrixF> readSparseMatrixF( const String &fileName );


/// load/save matrix from/to CSV text file
void exportMatrixF( const String &fileName, const MatrixF &mat );
aptr<MatrixF> importMatrixF( const String &fileName );


//-------------------------------------------
// ROW/COL MATRIX OPERATIONS
//-------------------------------------------


/// returns the mean of the rows (or cols) (the mean along each column (or row))
VectorF rowMean( const MatrixF &m );
VectorF colMean( const MatrixF &m );


/// returns the sum of the columns (the sum along each row)
template <typename T> Vector<T> rowSum( const Matrix<T> &m );


/// returns the sum of the rows (the sum along each col)
template <typename T> Vector<T> colSum( const Matrix<T> &m );


/// subtract vector from each row of matrix
aptr<MatrixF> subtractRow( const MatrixF &m, const VectorF &v );


//-------------------------------------------
// MATRIX MULTIPLICATION
//-------------------------------------------


/// returns matrix-vector product
VectorF multiply( const MatrixF &m, const VectorF &v );


/// returns matrix-vector product
VectorF multiplyXTY( const MatrixF &m, const VectorF &v );


/// multiplies each element in matrix with scalar
template <typename T> void multiply( const Matrix<T> &mIn, T v, Matrix<T> &mOut );


/// returns matrix-matrix product
aptr<MatrixF> multiply( const MatrixF &m1, const MatrixF &m2 );


/// returns matrix product
aptr<MatrixF> multiplyXTX( const MatrixF &m );


/// returns matrix-matrix product
aptr<MatrixF> multiplyXYT( const MatrixF &m1, const MatrixF &m2 );


//-------------------------------------------
// OTHER MATRIX MATH
//-------------------------------------------


/// computes col by col covariance, E[ x^T x ]
/// assumes mean already subtracted from x
aptr<MatrixF> innerCovarance( const MatrixF &x );


/// add a scalar to the diagonal elements of the matrix (assumes square)
void addToDiag( MatrixF &m, float val );


/// the mean element-wise difference 
float meanDiff( const MatrixF &m1, const MatrixF &m2 );


/// display stats about the difference between the matrices
void dispMatrixCompare( const MatrixF &m1, const MatrixF &m2 );


/// returns a distance matrix between rows of x
aptr<MatrixF> distSqdMatrix( const MatrixF &x );


/// solve system A x = b for x
VectorF solveEquation( const MatrixF &a, const VectorF &b );


/// compute eigenvectors and eigenvalues of symmetric matrix; returns eigenvectors as matrix columns
aptr<MatrixF> eigenSymmetric( const MatrixF &m, VectorF &eigenVals );


//-------------------------------------------
// MATRIX CONVERSION AND CREATION
//-------------------------------------------


/// returns identity matrix of size by size
aptr<MatrixF> identityF( int size );


/// each input vector becomes a row in the returned matrix
aptr<MatrixF> toMatrixF( const Array<VectorF> &vectorArray );


/// create a matrix of random values, sampled uniformly between min and max
aptr<MatrixF> randomMatrixF( int rows, int cols, float min, float max );


/// unroll the matrix into a vector (scanning across rows);
/// returns aptr because vector may be large
aptr<VectorF> toVector( MatrixF &m );


/// create a new matrix containing each row for which rowMask is non-zero and each column for which colMask is non-zero
aptr<MatrixF> subset( const MatrixF &m, const VectorI &rowMask, const VectorI &colMask );


//-------------------------------------------
// TEST COMMANDS
//-------------------------------------------


/// test efficiency of return value vs return pointer
void testReturn();


/// test memory allocation failures
void testMemoryAlloc();


} // end namespace sbl
#endif // _SBL_MATRIX_UTIL_H_

