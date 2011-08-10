#ifndef _SBL_TENSOR_H_
#define _SBL_TENSOR_H_
#include <sbl/core/Display.h>
namespace sbl {


/// The Tensor class represents a multi-dimensional matrix.
/// The representation is recursive, facilitating recursive algorithms to operate on the data.
/// The Tensor1F class is the base of the recursion (a 1D vector of floats).
class Tensor1F {
public:

    // basic constructor / destructor
    Tensor1F( int size = 1 ) { m_data = new float[ size ]; m_size = size; }
    ~Tensor1F() { delete [] m_data; }

    /// the number of dimensions 
    inline int dimCount() const { return 1; }

    /// the size of tensor at the current level (the current dimension)
    inline int size() const { return m_size; }

    /// the size at a given level 
    inline int size( int level ) const { assertAlways( level == 0 ); return m_size; }

    /// set the size at a given level
    inline void setSize( int level, int size ) {
        assertAlways( level == 0 );
        delete [] m_data;
        m_data = new float[ size ];
        m_size = size;
    }

    /// access a sub-tensor from the current level
    inline const float &operator[]( int i ) const { assertDebug( i >= 0 && i < m_size ); return m_data[ i ]; };
    inline float &operator[]( int i ) { assertDebug( i >= 0 && i < m_size ); return m_data[ i ]; };

    /// access an element given the position within each dimension
    inline const float &elem( const int *ind ) const { assertDebug( ind[ 0 ] >= 0 && ind[ 0 ] < m_size ); return m_data[ ind[ 0 ] ]; }
    inline float &elem( const int *ind ) { assertDebug( ind[ 0 ] >= 0 && ind[ 0 ] < m_size ); return m_data[ ind[ 0 ] ]; }

    /// clean the entire tensor to a specific value
    inline void operator=( float value ) {    for (int i = 0; i < m_size; i++) m_data[ i ] = value; }

    /// multiply every element by a scalar
    inline void operator*=( float value ) { for (int i = 0; i < m_size; i++) m_data[ i ] *= value; }

    /// copy the tensor from another (of the same size)
    inline void operator=( const Tensor1F &t ) { for (int i = 0; i < m_size; i++) m_data[ i ] = t.m_data[ i ]; }

private:

    // tensor data at this level
    float *m_data;
    int m_size;

    // disable copy constructor 
    Tensor1F( const Tensor1F &x );
};


/// The Tensor class represents a multi-dimensional matrix.
/// The representation is recursive, facilitating recursive algorithms to operate on the data.
// fix(clean): make this a sub-class of Tensor1F to reduce code replication? (might have large costs of space and/or speed)
template <typename T> class Tensor {
public:

    // basic constructor / destructor
    Tensor( int size = 1 ) { m_data = new T[ size ]; m_size = size; }
    ~Tensor() { delete [] m_data; }

    /// the number of dimensions 
    inline int dimCount() const { return 1 + m_data[ 0 ].dimCount(); }

    /// the size of tensor at the current level (the current dimension)
    inline int size() const { return m_size; }

    /// the size at a given level 
    inline int size( int level ) const {
        if (level)
            return m_data[ 0 ].size( level - 1 );
        return m_size;
//        assertDebug( level > 0 ); return m_data[ 0 ].size( level - 1 ); 
    }

    /// set the size at a given level
    inline void setSize( int level, int size ) {
/*        assertAlways( level > 0 );
        for (int i = 0; i < m_size; i++) 
            m_data[ i ].setSize( level - 1, size );
*/
        if (level == 0) {
            delete [] m_data;
            m_data = new T[ size ];
            m_size = size;
        } else {
            for (int i = 0; i < m_size; i++) 
                m_data[ i ].setSize( level - 1, size );
        }
    }

    /// access a sub-tensor from the current level
    inline const T &operator[]( int i ) const { assertDebug( i >= 0 && i < m_size ); return m_data[ i ]; };
    inline T &operator[]( int i ) { assertDebug( i >= 0 && i < m_size ); return m_data[ i ]; };

    /// access an element given the position within each dimension
    inline const float &elem( const int *ind ) const { assertDebug( ind[ 0 ] >= 0 && ind[ 0 ] < m_size ); return m_data[ ind[ 0 ] ].elem( ind + 1 ); }
    inline float &elem( const int *ind ) { assertDebug( ind[ 0 ] >= 0 && ind[ 0 ] < m_size ); return m_data[ ind[ 0 ] ].elem( ind + 1 ); }

    /// clean the entire tensor to a specific value
    inline void operator=( float value ) { for (int i = 0; i < m_size; i++) m_data[ i ] = value; }

    /// multiply every element by a scalar
    inline void operator*=( float value ) { for (int i = 0; i < m_size; i++) m_data[ i ] *= value; }

    /// copy the tensor from another (of the same size)
    inline void operator=( const Tensor<T> &t ) { for (int i = 0; i < m_size; i++) m_data[ i ] = t.m_data[ i ]; }

private:

    // tensor data at this level
    T *m_data;
    int m_size;

    // disable copy constructor 
    Tensor( const Tensor &x );
};


// common tensor types
typedef Tensor<Tensor1F> Tensor2F;
typedef Tensor<Tensor2F> Tensor3F;
typedef Tensor<Tensor3F> Tensor4F;
typedef Tensor<Tensor4F> Tensor5F;
typedef Tensor<Tensor5F> Tensor6F;
typedef Tensor<Tensor6F> Tensor7F;
typedef Tensor<Tensor7F> Tensor8F;


} // end namespace sbl
#endif // _SBL_TENSOR_H_

