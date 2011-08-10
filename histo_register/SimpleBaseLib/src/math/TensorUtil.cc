// Licensed under MIT license; see license.txt.

#include <sbl/math/TensorUtil.h>
#include <math.h>
namespace sbl {


//-------------------------------------------
// TENSOR FILTER UTILS
//-------------------------------------------


/// compute sum of all elements within radius of tensor element indexed by ind
float elementSum( const Tensor1F &tensor, const int *ind, int radius ) {
    int minIndex = ind[ 0 ] - radius;
    int maxIndex = ind[ 0 ] + radius;
    if (minIndex < 0)
        minIndex = 0;
    if (maxIndex > tensor.size() - 1)
        maxIndex = tensor.size() - 1;
    float sum = 0;
    for (int i = minIndex; i <= maxIndex; i++) 
        sum += tensor[ i ];
    return sum / (float) (maxIndex - minIndex + 1);
}


//-------------------------------------------
// TENSOR HISTOGRAM UTILS
//-------------------------------------------


/// compute negative entropy of probability distribution represented with a tensor histogram
double negativeEntropy( const Tensor1F &tensor ) {
    double sum = 0;
    for (int i = 0; i < tensor.size(); i++) {
        double v = tensor[ i ];
        if (v > 1e-10)
            sum += v * log( v );
    }
    return sum;
}


/// returns sum of all values for with the given index for the given level
double marginalSum( const Tensor1F &tensor, int level, int index ) {
    double sum = 0;
    if (level == 0) {
        sum = tensor[ index ];
    } else {
        for (int i = 0; i < tensor.size(); i++)
            sum += tensor[ i ];
    }
    return sum;
}


} // end namespace sbl

