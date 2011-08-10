#ifndef _SBL_VECTOR_H_
#define _SBL_VECTOR_H_
#include <sbl/core/Display.h>
#include <string.h> // for memcpy
#include <stdlib.h> // for qsort
namespace sbl {


/// min/max macros conflict with our vector min/max methods
#undef min
#undef max


//-------------------------------------------
// TEMPLATE VECTOR CLASS
//-------------------------------------------


/// The Vector class represents a resizeable numeric vector.
template<typename T> class Vector {
public:

    // basic constructor / destructor
    explicit inline Vector( int length = 0 ) { alloc( length ); }
    Vector( const Vector &vector );
    Vector( int length, const T *data ) { alloc( length ); for (int i = 0; i < length; i++) m_data[ i ] = data[ i ]; }
    Vector( T v1, T v2, T v3, T v4, T v5 ); // e.g. for unit tests
    inline ~Vector() { delete [] m_data; }

    //-------------------------------------------
    // ACCESS VECTOR DATA
    //-------------------------------------------

    /// the current length of the vector
    int inline size() const { return m_length; }
    int inline length() const { return m_length; }

    /// access vector data by index
    inline T operator[]( int i ) const { assertDebug( i >= 0 && i < m_length ); return m_data[ i ]; }
    inline T &operator[]( int i ) { assertDebug( i >= 0 && i < m_length ); return m_data[ i ]; }
    inline const T &data( int i ) const { assertDebug( i >= 0 && i < m_length ); return m_data[ i ]; }
    inline T &data( int i ) { assertDebug( i >= 0 && i < m_length ); return m_data[ i ]; }

    /// directly access data pointer
    inline const T *dataPtr() const { return m_data; }
    inline T *dataPtr() { return m_data; }

    /// retrieve last element in vector
    inline T endValue() const { assertDebug( m_length ); return m_data[ m_length - 1 ]; }

    //-------------------------------------------
    // MODIFY / MANIPULATE VECTOR DATA
    //-------------------------------------------

    /// append a single value, increasing the vector length by one
    void append( T val );

    /// append a vector of values
    inline void append( const Vector &vector ) { for (int i = 0; i < vector.length(); i++) append( vector[ i ] ); }

    /// set all elements to the given value
    inline void clear( T val ) { for (int i = 0; i < m_length; i++) m_data[ i ] = val; }

    /// resets vector length; clobbers any existing data; leaves data uninitialized
    inline void setLength( int length ) { delete [] m_data; alloc( length ); }

    /// assignment operator
    Vector &operator=( const Vector &vector );

    /// sort the vector elements (ascending) in place
    void sort();

    //-------------------------------------------
    // GET VECTOR INFO
    //-------------------------------------------

    /// compare with another vector
    bool operator==( const Vector &vector ) const;
    inline bool operator!=( const Vector &vector ) const { return !operator==( vector ); }

    /// true if contains the given value
    bool contains( T val ) const;

    /// number of bytes used by this object
    inline int memUsed() const { return sizeof(Vector) + sizeof(T) * m_allocLength; }

    /// vector element statistics
    T min() const;
    T max() const;
    T mean() const;
    T sum() const;

private:

    /// set vector length
    void alloc( int length );

    /// the data
    T *m_data;

    /// the number of defined elements
    int m_length;

    /// the allocated length (>= m_length)
    int m_allocLength;
};


/// append a single value, increasing the vector length by one
template <typename T> void Vector<T>::append( T val ) {
    if (m_length == m_allocLength) {
        m_allocLength *= 2;
        T *newData = new T[ m_allocLength ];
        assertDebug( newData );
        memcpy( newData, m_data, m_length * sizeof(T) );
        delete [] m_data;
        m_data = newData;
    }
    m_data[ m_length ] = val;
    m_length++;
}


/// compare with another vector
template <typename T> bool Vector<T>::operator==( const Vector<T> &vector ) const {
    if (m_length != vector.length())
        return false;
    for (int i = 0; i < m_length; i++) 
        if (m_data[ i ] != vector[ i ])
            return false;
    return true;
}


/// true if contains the given value
template <typename T> bool Vector<T>::contains( T val ) const {
    for (int i = 0; i < m_length; i++) 
        if (m_data[ i ] == val)
            return true;
    return false;
}


/// the minimum element value
template <typename T> T Vector<T>::min() const {
    if (m_length == 0)
        return 0;
    T min = m_data[ 0 ];
    for (int i = 0; i < m_length; i++) {
        T val = m_data[ i ];
        if (val < min) 
            min = val;
    }
    return min;
}


/// the maximum element value
template <typename T> T Vector<T>::max() const {
    if (m_length == 0)
        return 0;
    T max = m_data[ 0 ];
    for (int i = 0; i < m_length; i++) {
        T val = m_data[ i ];
        if (max < val) 
            max = val;
    }
    return max;
}


/// the mean element value
template <typename T> T Vector<T>::mean() const {
    if (m_length) {
        double sum = 0;
        for (int i = 0; i < m_length; i++) 
            sum += m_data[ i ];
        return (T) (sum / (double) m_length);
    } else {
        return 0;
    }
}


/// the sum of element values
template <typename T> T Vector<T>::sum() const {
    double sum = 0;
    for (int i = 0; i < m_length; i++) 
        sum += m_data[ i ];
    return (T) sum;
}


/// comparison used by Vector::sort
template <typename T> int compare( const void *v1, const void *v2 ) {
    T i = *((T *) v1);
    T j = *((T *) v2);
    if (i > j)
        return 1;
    if (i < j)
        return -1;
    return 0;
}


/// sort the vector elements (ascending) in place
template <typename T> void Vector<T>::sort() {
    qsort( m_data, m_length, sizeof(T), compare<T> );
}


/// basic copy constructor
template <typename T> Vector<T>::Vector( const Vector<T> &vector ) {
    alloc( vector.length() );
    memcpy( m_data, vector.m_data, m_length * sizeof(T) );
}


/// a wrapper constructor for creating concise unit tests
template <typename T> Vector<T>::Vector( T v1, T v2, T v3, T v4, T v5 ) {
    alloc( 5 );
    m_data[ 0 ] = v1;
    m_data[ 1 ] = v2;
    m_data[ 2 ] = v3;
    m_data[ 3 ] = v4;
    m_data[ 4 ] = v5;
}


/// set vector length
template <typename T> void Vector<T>::alloc( int length ) {
    assertDebug( length >= 0 );
    m_length = length;
    m_allocLength = length;
    if (m_allocLength == 0) m_allocLength = 10;
    m_data = new T[ m_allocLength ];
    assertDebug( m_data );
}


/// basic assignment operator
template <typename T> Vector<T> &Vector<T>::operator=( const Vector<T> &vector ) {
    if (this != &vector) {
        setLength( vector.length() );
        assertDebug( vector.length() == m_length );
        memcpy( m_data, vector.m_data, m_length * sizeof(T) );
    }
    return *this;
}


// common vector types
typedef Vector<unsigned char> VectorU;
typedef Vector<int> VectorI;
typedef Vector<float> VectorF;
typedef Vector<double> VectorD;


} // end namespace sbl
#endif // _SBL_VECTOR_H_

