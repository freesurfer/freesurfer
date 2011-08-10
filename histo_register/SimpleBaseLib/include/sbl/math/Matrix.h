#ifndef _SBL_MATRIX_H_
#define _SBL_MATRIX_H_
#include <sbl/core/Display.h>
#include <sbl/math/Vector.h>
#include <string.h> // for memcpy
namespace sbl {


/// min/max macros conflict with out vector min/max methods
#undef min
#undef max


//-------------------------------------------
// TEMPLATE MATRIX CLASS
//-------------------------------------------


/// The Matrix class represents a standard numeric matrix.
template <typename T> class Matrix {
public:

    /// constructor / destructor
    Matrix( int rows, int cols, bool contiguous = true );
    Matrix( const Matrix<T> &matrix );
    ~Matrix();

    /// matrix dimensions
    inline int rows() const { return m_rows; }
    inline int cols() const { return m_cols; }

    /// access to matrix data
    inline T operator()( int i, int j ) const { assertDebug( i >= 0 && i < m_rows && j >= 0 && j < m_cols ); return m_data[ i ][ j ]; }
    inline T &operator()( int i, int j ) { assertDebug( i >= 0 && i < m_rows && j >= 0 && j < m_cols ); return m_data[ i ][ j ]; }
    inline T data( int i, int j ) const { assertDebug( i >= 0 && i < m_rows && j >= 0 && j < m_cols ); return m_data[ i ][ j ]; }
    inline T &data( int i, int j ) { assertDebug( i >= 0 && i < m_rows && j >= 0 && j < m_cols ); return m_data[ i ][ j ]; }
    inline T **const dataPtr() const { return m_data; } 
    inline T **dataPtr() { return m_data; }
    inline T *dataVectPtr() { return m_dataVect; }
    inline const T *dataVectPtr() const { return m_dataVect; }
    inline const T *dataRow( int i ) const { assertDebug( i >= 0 && i < m_rows ); return m_data[ i ]; }
    inline T *dataRow( int i ) { assertDebug( i >= 0 && i < m_rows ); return m_data[ i ]; }

    /// set all items to the given value
    void clear( T val );

    /// extract a row or column 
    Vector<T> row( int row ) const;
    Vector<T> col( int col ) const;

    /// matrix element statistics
    T min() const;
    T max() const;
    T mean() const;
    T sum() const;

    /// number of bytes used by the matrix
    // fix(later): does not handle contiguity
    inline int memUsed() const { return sizeof(Matrix<T>) + sizeof(float) * m_rows * m_cols; }

private:

    /// size of the matrix
    int m_rows;
    int m_cols;

    /// if true, matrix is allocated as a single block of memory; otherwise allocated per row
    bool m_contiguous;

    /// the matrix data
    T **m_data;
    T *m_dataVect;

    // disable assignment operator
    Matrix &operator=( const Matrix &m );
};


// basic constructor
template <typename T> Matrix<T>::Matrix( int rows, int cols, bool contiguous ) {
    assertDebug( rows >= 0 && cols >= 1 );
    m_rows = rows;
    m_cols = cols;
    m_contiguous = contiguous;

    // allocate as single block with pointer to each row
    if (m_contiguous) {
        m_dataVect = new T[ m_rows * m_cols ];
        assertDebug( m_dataVect );
        m_data = new T*[ m_rows ];
        assertDebug( m_data );
        for (int i = 0; i < m_rows; i++)
            m_data[ i ] = m_dataVect + i * cols;

    // allocate each row seperately
    } else {
        m_dataVect = NULL;
        m_data = new T*[ m_rows ];
        assertDebug( m_data );
        for (int i = 0; i < m_rows; i++) {
            m_data[ i ] = new T[ m_cols ];
            assertDebug( m_data[ i ] );
        }
    }
}


// copy constructor
template <typename T> Matrix<T>::Matrix( const Matrix<T> &matrix ) {
    m_rows = matrix.rows();
    m_cols = matrix.cols();
    m_contiguous = true; // fix(later): create non-contiguous copies if source is non-contiguous
    m_dataVect = new T[ m_rows * m_cols ];
    assertDebug( m_dataVect );
    m_data = new T*[ m_rows ];
    assertDebug( m_data );
    for (int i = 0; i < m_rows; i++) {
        m_data[i] = m_dataVect + i * m_cols;
        memcpy( m_data[i], matrix.m_data[i], m_cols * sizeof( T ) );
    }
}


// deallocate matrix data
template <typename T> Matrix<T>::~Matrix() {
    if (m_contiguous) {
        delete [] m_dataVect;
    } else {
        for (int i = 0; i < m_rows; i++) 
            delete [] m_data[ i ];
    }
    delete [] m_data;
}


/// set all items to the given value
template <typename T> void Matrix<T>::clear( T val ) {
    for (int i = 0; i < m_rows; i++)
        for (int j = 0; j < m_cols; j++)
            m_data[ i ][ j ] = val;
}


/// extract a row 
template <typename T> Vector<T> Matrix<T>::row( int row ) const {
    assertDebug( row >= 0 && row < m_rows );
    Vector<T> vector( m_cols );
    memcpy( vector.dataPtr(), m_data[ row ], m_cols * sizeof( T ));
    return vector;
}


/// extract a column
template <typename T> Vector<T> Matrix<T>::col( int col ) const {
    assertDebug( col >= 0 && col < m_cols );
    Vector<T> vector( m_rows );
    for (int i = 0; i < m_rows; i++) 
        vector.data( i ) = m_data[ i ][ col ];
    return vector;
}


/// minumum element value
template <typename T> T Matrix<T>::min() const {
    if (m_rows == 0 || m_cols == 0)
        return 0;
    T min = m_data[ 0 ][ 0 ];
    for (int i = 0; i < m_rows; i++) { // fix(faster): could make faster if contiguous
        for (int j = 0; j < m_cols; j++) {
            T val = m_data[ i ][ j ];
            if (val < min) 
                min = val;
        }
    }
    return min;
}


/// maximum element value
template <typename T> T Matrix<T>::max() const {
    if (m_rows == 0 || m_cols == 0)
        return 0;
    T max = m_data[ 0 ][ 0 ];
    for (int i = 0; i < m_rows; i++) { // fix(faster): could make faster if contiguous
        for (int j = 0; j < m_cols; j++) {
            T val = m_data[ i ][ j ];
            if (max < val) 
                max = val;
        }
    }
    return max;
}


/// mean element value
template <typename T> T Matrix<T>::mean() const {
    if (m_rows && m_cols) {
        double sum = 0;
        for (int i = 0; i < m_rows; i++) { // fix(faster): could make faster if contiguous
            for (int j = 0; j < m_cols; j++) {
                sum += m_data[ i ][ j ];
            }
        }
        return (T) (sum / (double) (m_rows * m_cols));
    } else {
        return 0;
    }
}


/// sum of element values
template <typename T> T Matrix<T>::sum() const {
    double sum = 0;
    for (int i = 0; i < m_rows; i++) { // fix(faster): could make faster if contiguous
        for (int j = 0; j < m_cols; j++) {
            sum += m_data[ i ][ j ];
        }
    }
    return (T) sum;
}


// common matrix types
typedef Matrix<unsigned char> MatrixU;
typedef Matrix<int> MatrixI;
typedef Matrix<float> MatrixF;
typedef Matrix<double> MatrixD;


} // end namespace sbl
#endif // _SBL_MATRIX_H_

