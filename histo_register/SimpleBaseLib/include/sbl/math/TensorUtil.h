#ifndef _SBL_TENSOR_UTIL_H_
#define _TENSOR_UTIL_H_
#include <sbl/math/Tensor.h>
namespace sbl {


/*! \file TensorUtil.h
    \brief The TensorUtil module provides functions that operate on tensors 
    (represented using the Tensor class).
*/


//-------------------------------------------
// TENSOR FILTER UTILS
//-------------------------------------------


/// compute sum of all elements within radius of tensor element indexed by ind
float elementSum( const Tensor1F &tensor, const int *ind, int radius );


/// compute sum of all elements within radius of tensor element indexed by ind
template<typename T> float elementSum( const Tensor<T> &tensor, const int *ind, int radius ) {
    int minIndex = ind[ 0 ] - radius;
    int maxIndex = ind[ 0 ] + radius;
    if (minIndex < 0)
        minIndex = 0;
    if (maxIndex > tensor.size() - 1)
        maxIndex = tensor.size() - 1;
    float sum = 0;
    for (int i = minIndex; i <= maxIndex; i++) 
        sum += elementSum( tensor[ i ], ind + 1, radius );
    return sum / (float) (maxIndex - minIndex + 1);
}


/// apply a box filter to the tensor
template <typename T> void applyBlur( int *ind, int level, const T &source, T &dest, int blurRadius ) {
    if (level < source.dimCount()) {
        int size = source.size( level );
        for (int i = 0; i < size; i++) {
            ind[ level ] = i;
            applyBlur( ind, level + 1, source, dest, blurRadius );
        }
    } else {
        dest.elem( ind ) = elementSum( source, ind, blurRadius );
    }
}


/// apply a box filter to the tensor
template <typename T> void blurBox( const T &source, T &dest, int blurSize ) {
    int blurRadius = (blurSize - 1) / 2;
    int *ind = new int[ source.dimCount() ];
    for (int i = 0; i < source.dimCount(); i++)
        ind[ i ] = 0;
    applyBlur( ind, 0, source, dest, blurRadius );
}


//-------------------------------------------
// TENSOR HISTOGRAM UTILS
//-------------------------------------------


/// compute negative entropy of probability distribution represented with a tensor histogram
double negativeEntropy( const Tensor1F &tensor );


/// compute negative entropy of probability distribution represented with a tensor histogram
template<typename T> double negativeEntropy( const Tensor<T> &tensor ) {
    double sum = 0;
    for (int i = 0; i < tensor.size(); i++) 
        sum += negativeEntropy( tensor[ i ] );
    return sum;
}


/// returns sum of all values for with the given index for the given level
double marginalSum( const Tensor1F &tensor, int level, int index );


/// returns sum of all values for with the given index for the given level
template<typename T> double marginalSum( const Tensor<T> &tensor, int level, int index ) {
    double sum = 0;
    if (level == 0) {
        sum = marginalSum( tensor[ index ], level - 1, index );
    } else {
        for (int i = 0; i < tensor.size(); i++)
            sum += marginalSum( tensor[ i ], level - 1, index );
    }
    return sum;
}


} // end namespace sbl
#endif // _SBL_TENSOR_UTIL_H_

