// Licensed under MIT license; see license.txt.

#include <sbl/math/MatrixUtil.h>
#include <sbl/core/Display.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/system/Timer.h>
#ifdef USE_OPENCV 
    #include <opencv/cxcore.h>
#endif 
namespace sbl {


//-------------------------------------------
// MATRIX FILE I/O
//-------------------------------------------


/// write matrix (and dimensions) to file
void writeMatrixF( File &file, const MatrixF &mat ) {
    int rows = mat.rows(), cols = mat.cols();
    file.writeInt( rows );
    file.writeInt( cols );
    for (int i = 0; i < rows; i++) 
        for (int j = 0; j < cols; j++) 
            file.writeFloat( mat( i, j ) );
}


/// write symmetric matrix (and dimensions) to file
void writeSymMatrixF( File &file, const MatrixF &mat ) {
    int rows = mat.rows(), cols = mat.cols();
    assertDebug( rows == cols );
    file.writeInt( rows );
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j <= i; j++) {
            if (fAbs( mat( i, j ) - mat( j, i )) > 0.0000001f) {
                disp( 1, "i: %d, j: %d, rows: %d, cols: %d", i, j, rows, cols );
                disp( 1, "mat( i, j ): %g, mat( j, i ): %g", mat( i, j ), mat( j, i ));
                fatalError( "matrix not symmetric" );
            }
            file.writeFloat( mat( i, j ) );
        }
    }
}


/// write matrix (and dimensions) to file
void writeSparseMatrixF( File &file, const MatrixF &mat ) {
    int rows = mat.rows(), cols = mat.cols();
    file.writeInt( rows );
    file.writeInt( cols );
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (mat( i, j )) {
                file.writeInt( j );
                file.writeFloat( mat( i, j ) );
            }
        }
        file.writeInt( -1 );
    }
}


/// write matrix (and dimensions) to file
void writeMatrixF( const String &fileName, const MatrixF &mat ) {
    File file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess())
        writeMatrixF( file, mat );
}


/// write matrix (and dimensions) to file
void writeSymMatrixF( const String &fileName, const MatrixF &mat ) {
    File file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess())
        writeSymMatrixF( file, mat );
}


/// write matrix (and dimensions) to file
void writeSparseMatrixF( const String &fileName, const MatrixF &mat ) {
    File file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess())
        writeSparseMatrixF( file, mat );
}


/// read matrix (and dimensions) from file
aptr<MatrixF> readMatrixF( File &file ) {
    aptr<MatrixF> mat;
    int rows = file.readInt();
    int cols = file.readInt();
    if (rows < 0 || cols < 0) {
        warning( "invalid matrix" );
    } else {
        mat.reset( new MatrixF( rows, cols ) );
        for (int i = 0; i < rows; i++) 
            for (int j = 0; j < cols; j++) 
                mat->data( i, j ) = file.readFloat();
    }
    return mat;
}


/// read symmetric matrix (and dimensions) from file
aptr<MatrixF> readSymMatrixF( File &file ) {
    aptr<MatrixF> mat;
    int rows = file.readInt();
    int cols = rows;
    if (rows < 0 || cols < 0) {
        warning( "invalid matrix" );
    } else {
        mat.reset( new MatrixF( rows, cols, false ) );
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j <= i; j++) {
                float m = file.readFloat();
                mat->data( i, j ) = m;
                mat->data( j, i ) = m;
            }
        }
    }
    return mat;
}


/// read matrix (and dimensions) from file
aptr<MatrixF> readSparseMatrixF( File &file ) {
    aptr<MatrixF> mat;
    int rows = file.readInt();
    int cols = file.readInt();
    if (rows < 0 || cols < 0) {
        warning( "invalid matrix" );
    } else {
        mat.reset( new MatrixF( rows, cols ) );
        mat->clear( 0.0f );
        for (int i = 0; i < rows; i++) {
            while (1) {
                int j = file.readInt();
                if (j == -1)
                    break;
                mat->data( i, j ) = file.readFloat();
            }
        }
    }
    return mat;
}


/// read matrix (and dimensions) from file
aptr<MatrixF> readMatrixF( const String &fileName ) {
    aptr<MatrixF> mat;
    File file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess())
        mat = readMatrixF( file );
    return mat ;
}


/// read matrix (and dimensions) from file
aptr<MatrixF> readSymMatrixF( const String &fileName ) {
    aptr<MatrixF> mat;
    File file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess())
        mat = readSymMatrixF( file );
    return mat;
}


/// read matrix (and dimensions) from file
aptr<MatrixF> readSparseMatrixF( const String &fileName ) {
    aptr<MatrixF> mat;
    File file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess())
        mat = readSparseMatrixF( file );
    return mat;
}


/// save matrix to CSV file
void exportMatrixF( const String &fileName, const MatrixF &mat ) {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {
        int rows = mat.rows(), cols = mat.cols();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                file.writeF( "%f",  mat( i, j ) );
                if (j + 1 < cols)
                    file.writeF( ", " );
            }
            file.writeF( "\n" );
        }
    }
}


//-------------------------------------------
// ROW/COL MATRIX OPERATIONS
//-------------------------------------------


/// returns the mean of the rows (the mean along each column)
VectorF rowMean( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    VectorF result( cols );
    result.clear( 0 ); 
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) 
            result[ j ] += m( i, j );
    float factor = 1.0f / (float) rows;
    for (int j = 0; j < cols; j++)
        result[ j ] *= factor;
    return result;
}


/// returns the mean of the cols (the mean along each row)
VectorF colMean( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    VectorF result( rows );
    result.clear( 0 ); 
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) 
            result[ i ] += m( i, j );
    float factor = 1.0f / (float) cols;
    for (int j = 0; j < rows; j++)
        result[ j ] *= factor;
    return result;
}


/// returns the sum of the columns (the sum along each row)
template <typename T> Vector<T> rowSum( const Matrix<T> &m ) {
    int rows = m.rows(), cols = m.cols();
    Vector<T> result( rows );
    for (int i = 0; i < rows; i++) {
        T sum = 0;
        for (int j = 0; j < cols; j++) {
            sum += m( i, j );
        }
        result[ i ] = sum;
    }
    return result;
}
template Vector<int> rowSum( const Matrix<int> &m );
template Vector<float> rowSum( const Matrix<float> &m );
template Vector<double> rowSum( const Matrix<double> &m );


/// returns the sum of the columns (the sum along each row)
template <typename T> Vector<T> colSum( const Matrix<T> &m ) {
    int rows = m.rows(), cols = m.cols();
    Vector<T> result( rows );
    for (int i = 0; i < rows; i++) {
        T sum = 0;
        for (int j = 0; j < cols; j++) {
            sum += m( i, j );
        }
        result[ i ] = sum;
    }
    return result;
}
template Vector<int> colSum( const Matrix<int> &m );
template Vector<float> colSum( const Matrix<float> &m );
template Vector<double> colSum( const Matrix<double> &m );


/// subtract vector from each row of matrix
aptr<MatrixF> subtractRow( const MatrixF &m, const VectorF &v ) {
    assertDebug( m.cols() == v.length() );
    int rows = m.rows(), cols = m.cols();
    aptr<MatrixF> result( new MatrixF( rows, cols ) );
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result->data( i, j ) = m( i, j ) - v[ j ];
    return result;
}


//-------------------------------------------
// MATRIX MULTIPLICATION
//-------------------------------------------


/// returns matrix-vector product; treats vector as column vector
VectorF multiply( const MatrixF &m, const VectorF &v ) {
    int rows = m.rows(), cols = m.cols();
    assertDebug( v.length() == cols );
    VectorF prod( rows );
    for (int i = 0; i < rows; i++) {
        float sum = 0;
        for (int j = 0; j < cols; j++) 
            sum += m( i, j ) * v[ j ];
        prod[ i ] = sum;
    }
    return prod;
}


/// returns matrix-vector product
VectorF multiplyXTY( const MatrixF &m, const VectorF &v ) {
    int rows = m.rows(), cols = m.cols(), len = v.length();
    assertAlways( len == rows );
    VectorF prod( cols );
    for (int i = 0; i < cols; i++) {
        float sum = 0;
        for (int j = 0; j < rows; j++) 
            sum += m( j, i ) * v[ j ];
        prod[ i ] = sum;
    }
    return prod;
}


/// multiplies each element in matrix with scalar
template <typename T> void multiply( const Matrix<T> &mIn, T v, Matrix<T> &mOut ) {
    int rows = mIn.rows(), cols = mOut.cols();
    assertDebug( mOut.rows() == rows && mOut.cols() == cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mOut( i, j ) = v * mIn( i, j );
        }
    }
}
template void multiply( const Matrix<int> &mIn, int v, Matrix<int> &mOut );
template void multiply( const Matrix<float> &mIn, float v, Matrix<float> &mOut );
template void multiply( const Matrix<double> &mIn, double v, Matrix<double> &mOut );


/// returns matrix-matrix product
aptr<MatrixF> multiply( const MatrixF &m1, const MatrixF &m2 ) {
    int rows1 = m1.rows(), cols1 = m1.cols();
    int cols2 = m2.cols();
    assertDebug( cols1 == m2.rows() );
    aptr<MatrixF> m( new MatrixF(rows1, cols2) );
    for (int i = 0; i < rows1; i++) 
        for (int j = 0; j < cols2; j++) {
            float sum = 0;
            for (int k = 0; k < cols1; k++) {
                sum += m1( i, k ) * m2( k, j );
            }
            m->data( i, j ) = sum;
        }
    return m;
}


/// returns matrix product
aptr<MatrixF> multiplyXTX( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    aptr<MatrixF> prod( new MatrixF( cols, cols ) );
    for (int i = 0; i < cols; i++) 
        for (int j = 0; j < cols; j++) {
            float sum = 0;
            for (int k = 0; k < rows; k++) 
                sum += m( k, i ) * m( k, j );
            prod->data( i, j ) = sum;
        }
    return prod;
}


/// returns matrix-matrix product
aptr<MatrixF> multiplyXTY( const MatrixF &m1, const MatrixF &m2 ) {
    int rows1 = m1.cols(), cols1 = m1.rows(); // rows/cols of transposed m1
    int rows2 = m2.rows(), cols2 = m2.cols();
    assertAlways( cols1 != rows2 );
    aptr<MatrixF> m( new MatrixF( rows1, cols2 ) );
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            float sum = 0;
            for (int k = 0; k < cols1; k++) {
                sum += m1( k, i ) * m2( k, j ); // transpose indices of m1
            }
            m->data( i, j ) = sum;
        }
    }
    return m;
}


/// returns matrix-matrix product
aptr<MatrixF> multiplyXYT( const MatrixF &m1, const MatrixF &m2 ) {
    int rows1 = m1.rows(), cols1 = m1.cols();
    int rows2 = m2.rows();
    assertDebug( cols1 == m2.cols() );
    aptr<MatrixF> m( new MatrixF( rows1, rows2 ) );
    for (int i = 0; i < rows1; i++) 
        for (int j = 0; j < rows2; j++) {
            float sum = 0;
            for (int k = 0; k < cols1; k++) {
                sum += m1( i, k ) * m2( j, k );
            }
            m->data( i, j ) = sum;
        }
    return m;
}


//-------------------------------------------
// OTHER MATRIX MATH
//-------------------------------------------


/// returns min elem in matrix
float min( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    assertDebug( rows && cols );
    float min = m( 0, 0 );
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            float val = m( i, j );
            if (val < min) 
                min = val;
        }
    return min;
}


/// returns max elem in matrix
float max( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    assertDebug( rows && cols );
    float max = m( 0, 0 );
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float val = m( i, j );
            if (val > max) 
                max = val;
        }
    }
    return max;
}


/// returns max elem in matrix
float sum( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    float sum = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) 
            sum += m( i, j );
    return sum;
}


/// transpose the matrix
aptr<MatrixF> transpose( const MatrixF &m ) {
    int rows = m.rows(), cols = m.cols();
    aptr<MatrixF> result( new MatrixF( cols, rows ) );
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) 
            result->data( j, i ) = m( i, j );
    return result;
}


/// computes col by col covariance, E[ x^T x ];
/// assumes mean already subtracted from x
aptr<MatrixF> innerCovarance( const MatrixF &x ) {
    int rows = x.rows(), cols = x.cols();
    aptr<MatrixF> xt = transpose( x );
    aptr<MatrixF> cov( new MatrixF( cols, cols ) );
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j <= i; j++) {
            float *xti = xt->dataRow( i );
            float *xtj = xt->dataRow( j );
            float sum = 0;
            for (int k = 0; k < rows; k++) 
                sum += xti[ k ] * xtj[ k ];
            cov->data( i, j ) = cov->data( j, i ) = sum / (float) rows;
        }
    }
    return cov;
}


/// add a scalar to the diagonal elements of the matrix (assumes square)
void addToDiag( MatrixF &m, float val ) {
    assertDebug( m.rows() == m.cols() );
    int size = m.rows();
    for (int i = 0; i < size; i++)
        m.data( i, i ) += val;
}


/// the mean element-wise difference 
float meanDiff( const MatrixF &m1, const MatrixF &m2 ) {
    int rows = m1.rows(), cols = m1.cols();
    assertDebug( rows == m2.rows() && cols == m2.cols() );
    float sumDiff = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float diff = m1( i, j ) - m2( i, j );
            if (diff < 0) 
                diff = -diff;
            sumDiff += diff;
        }
    }
    return sumDiff / (float) (rows * cols);
}


/// display stats about the difference between the matrices
void dispMatrixCompare( const MatrixF &m1, const MatrixF &m2 ) {
    int rows = m1.rows(), cols = m1.cols();
    assertDebug( rows == m2.rows() && cols == m2.cols() );
    float sumDiff = 0, maxDiff = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            float diff = m1( i, j ) - m2( i, j );
            if (diff < 0) 
                diff = -diff;
            sumDiff += diff;
            if (diff > maxDiff)
                maxDiff = diff;
        }
    }
    float meanDiff = sumDiff / (float) (rows * cols);
    disp( 1, "max diff: %1.4f, mean diff: %1.4f", maxDiff, meanDiff );
}


/// returns a distance matrix between rows of x
aptr<MatrixF> distSqdMatrix( const MatrixF &x ) {
    int pointCount = x.rows();
    int dimCount = x.cols();
    aptr<MatrixF> dSqd( new MatrixF( pointCount, pointCount, false ) );
    for (int i = 0; i < pointCount; i++) {
        for (int j = 0; j < i; j++) {
            float d = distSqd( x.dataRow( i ), x.dataRow( j ), dimCount );
            dSqd->data( i, j ) = d;
            dSqd->data( j, i ) = d;
        }
        dSqd->data( i, i ) = 0;
    }
    return dSqd;
}


/// solve system A x = b for x
VectorF solveEquation( const MatrixF &a, const VectorF &b ) {
    int size = b.length();
    assertAlways( size == a.rows() && size == a.cols() );
    VectorF x( size );
#ifdef USE_OPENCV
    cv::Mat aMat( a.rows(), a.cols(), CV_32F, (void *) a.dataPtr() );
    cv::Mat bMat( b.length(), 1, CV_32F, (void *) b.dataPtr() );
    cv::Mat xMat( size, 1, CV_32F );
    cv::solve( aMat, bMat, xMat, cv::DECOMP_CHOLESKY );
    for (int i = 0; i < size; i++) 
        x[ i ] = xMat.at<float>( i, 0 );
#endif
    return x;
}


/// compute eigenvectors and eigenvalues of symmetric matrix; returns eigenvectors as matrix columns
aptr<MatrixF> eigenSymmetric( const MatrixF &m, VectorF &eigenVals ) {
    int size = m.rows();
    assertAlways( m.cols() == size );
    assertAlways( eigenVals.length() == size );
    aptr<MatrixF> eigenVects;
#ifdef USE_OPENCV

    // prepare data
    // note: these don't need to be deallocated by us
    CvMat *dataMat = cvCreateMat( size, size, CV_32F ); 
    for (int j = 0; j < size; j++)
        for (int i = 0; i < size; i++)
            cvmSet( dataMat, i, j, m.data( i, j ) );
    CvMat *vectMat = cvCreateMat( size, size, CV_32F );
    CvMat *valueMat = cvCreateMat( size, 1, CV_32F );

    // get eigenvalues and eigenvectors
    cvEigenVV( dataMat, vectMat, valueMat );

    // copy results 
    eigenVects.reset( new MatrixF( size, size ) );
    for (int i = 0; i < size; i++) {
        eigenVals.data( i ) = (float) cvmGet( valueMat, i, 1 );
        for (int j = 0; j < size; j++)
            eigenVects->data( i, j ) = (float) cvmGet( vectMat, j, i ); // note transpose
    }
    cvReleaseMat( &dataMat ); 
    cvReleaseMat( &vectMat ); 
    cvReleaseMat( &valueMat );
#endif
    return eigenVects;
}


//-------------------------------------------
// MATRIX CONVERSION AND CREATION
//-------------------------------------------


/// returns identity matrix of size by size
aptr<MatrixF> identityF( int size ) {
    assertDebug( size > 0 );
    aptr<MatrixF> m( new MatrixF( size, size ) );
    m->clear(0);
    for (int i = 0; i < size; i++)
        m->data( i, i ) = 1;
    return m;
}


/// each input vector becomes a row in the returned matrix
aptr<MatrixF> toMatrixF( const Array<VectorF> &vectorArray ) {
    int rows = vectorArray.count();
    int cols = vectorArray[ 0 ].length();
    assertDebug( rows && cols );
    aptr<MatrixF> m( new MatrixF( rows, cols ) );
    for (int i = 0; i < rows; i++) {
        const VectorF &v = vectorArray[ i ];
        for (int j = 0; j < cols; j++)
            m->data( i, j ) = v[ j ];
    }
    return m;
}



/// create a matrix of random values, sampled uniformly between min and max
aptr<MatrixF> randomMatrixF( int rows, int cols, float min, float max ) {
    aptr<MatrixF> m( new MatrixF( rows, cols ) );
    for (int i = 0; i < rows; i++) 
        for (int j = 0; j < cols; j++)
            m->data( i, j ) = randomFloat( min, max );
    return m;
}


/// unroll the matrix into a vector (scanning across rows);
/// returns aptr because vector may be large
aptr<VectorF> toVector( MatrixF &m ) {
    int rows = m.rows();
    int cols = m.cols();
    aptr<VectorF> v( new VectorF( rows * cols ) );
    int index = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) 
            v->data( index++ ) = m( i, j );
    return v;
}


/// create a new matrix containing each row for which rowMask is non-zero and each column for which colMask is non-zero
aptr<MatrixF> subset( const MatrixF &m, const VectorI &rowMask, const VectorI &colMask ) {
    int rows = m.rows(), cols = m.cols();
    int subRows = nonZeroCount( rowMask );
    int subCols = nonZeroCount( colMask );
    aptr<MatrixF> result( new MatrixF( subRows, subCols ) );
    int rowIndex = 0;
    for (int i = 0; i < rows; i++) {
        if (rowMask[ i ]) {
            int colIndex = 0;
            for (int j = 0; j < cols; j++) {
                if (colMask[ j ]) {
                    result->data( rowIndex, colIndex ) = m.data( i, j );
                    colIndex++;
                }
            }
            rowIndex++;
        }
    }
    return result;
}


//-------------------------------------------
// TEST COMMANDS
//-------------------------------------------


// test return by value
MatrixF returnByValue() {
    MatrixF m( 10000, 10000 );
    m.clear( 0 );
    return m;
}


// test return pointer
MatrixF *returnPtr() {
    MatrixF *m = new MatrixF( 10000, 10000 );
    m->clear( 0 );
    return m;
}


/// test efficiency of return value vs return pointer
void testReturn() {
    for (int i = 0; i < 5; i++) {

        // return pointer
        Timer timer2;
        MatrixF *m2 = returnPtr();
        disp( 1, "returnPtr: %f", timer2.timeSum() );

        // return by value
        Timer timer1;
        MatrixF m1 = returnByValue();
        disp( 1, "returnByValue: %f", timer1.timeSum() );

        // perform ops
        m1.clear( 1 );
        m2->clear( 1 );
        delete m2;
    }
}


/// test memory allocation failures
void testMemoryAlloc() {
    int size = 1024 * 1024;
    int bytes = 0;
    Array<VectorF> arr;
    for (int i = 0; i < 2000; i++) {
        VectorF *v = new VectorF( size / 4 );
        v->clear( 0 );
        arr.append( v );
        disp( 1, "bytes: %d", bytes );
        bytes += size;
    }
}


} // end namespace sbl

