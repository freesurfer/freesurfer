// Licensed under MIT license; see license.txt.

#include <sbl/math/VectorUtil.h>
#include <sbl/core/Display.h>
#include <sbl/core/UnitTest.h>
#include <sbl/core/Command.h>
#include <sbl/math/MathUtil.h>
#include <stdlib.h> // for qsort
#include <math.h>
namespace sbl {


//-------------------------------------------
// VECTOR MATH OPS
//-------------------------------------------


/// vector dot product
float dot( const VectorF &v1, const VectorF &v2 ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    float sumProd = 0;
    for (int i = 0; i < len; i++) 
        sumProd += v1[ i ] * v2[ i ];
    return sumProd;
}


/// make vector have unit length
void normalize( VectorF &v, float scale ) {
    int len = v.length();
    double norm = 0;
    for (int i = 0; i < len; i++) { 
        double val = v[ i ];
        norm += val * val;
    }
    norm = sqrt( norm );
    if (norm > 1e-8) {
        double factor = scale / norm;
        for (int i = 0; i < len; i++)
            v[ i ] = (float) ((double) v[ i ] * factor);
    } 
}


/// multiply a vector by a scalar
void multiply( const VectorF &src, float val, VectorF &dest ) {
    assertDebug( src.length() == dest.length() );
    int len = src.length();
    for (int i = 0; i < len; i++)
        dest[ i ] = src[ i ] * val;
}


/// add scalar to vector
void add( const VectorF &src, float val, VectorF &dest ) {
    assertDebug( src.length() == dest.length() );
    int len = src.length();
    for (int i = 0; i < len; i++)
        dest[ i ] = src[ i ] + val;
}


/// add scalar to vector
void add( const VectorI &src, int val, VectorI &dest ) {
    assertDebug( src.length() == dest.length() );
    int len = src.length();
    for (int i = 0; i < len; i++)
        dest[ i ] = src[ i ] + val;
}


/// add vector to vector
void add( const VectorF &v1, const VectorF &v2, VectorF &dest ) {
    assertDebug( v1.length() == v2.length() && v1.length() == dest.length() );
    int len = dest.length();
    for (int i = 0; i < len; i++)
        dest[ i ] = v1[ i ] + v2[ i ];
}


/// subtract vector from vector
VectorF subtract( const VectorF &v1, const VectorF &v2 ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    VectorF result( len  );
    for (int i = 0; i < len; i++)
        result[ i ] = v1[ i ] - v2[ i ];
    return result;
}


/// scale and shift vector so that min = 0, max = 1 
// (assuming not all elements are equal)
void normalizeBounds( VectorF &v ) {
    float min = v.min();
    float max = v.max();
    if (max - min > 1e-8) {
        float factor = 1.0f / (max - min);
        for (int i = 0; i < v.length(); i++) {
            v[ i ] = (v[ i ] - min) * factor;
        }
    } else {
//        warning( "normalizeBounds failed; min: %f, max: %f", min, max );
    }
}


/// clamp all values to be within given range
void clamp( VectorF &v, float min, float max ) {
    for (int i = 0; i < v.length(); i++) {
        if (v[ i ] < min)
            v[ i ] = min;
        if (v[ i ] > max)
            v[ i ] = max;
    }
}


//-------------------------------------------
// VECTOR COMPARISON
//-------------------------------------------


/// squared euclidean distance between vectors (sum of squared differences)
int distSqd( const VectorU &v1, const VectorU &v2 ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    int sum = 0;
    for (int i = 0; i < len; i++) {
        int diff = v1[ i ] - v2[ i ];
        sum += diff * diff;
    }
    return sum;
}


/// squared euclidean distance between vectors (sum of squared differences)
float distSqd( const VectorF &v1, const VectorF &v2 ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    float sum = 0;
    for (int i = 0; i < len; i++) {
        float diff = v1[ i ] - v2[ i ];
        sum += diff * diff;
    }
    return sum;
}


/// squared euclidean distance between vectors (sum of squared differences)
float distSqd( const float *p1, const float *p2, int len ) {
    assertDebug( p1 && p2 );
    float distSqd = 0;
    for (int i = 0; i < len; i++) {
        float diff = p1[ i ] - p2[ i ];
        distSqd += diff * diff;
    }
    return distSqd;
}


/// squared euclidean distance between vectors
float distSqd( const VectorF &v1, const VectorF &v2, float maxDistSqd ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    float sum = 0;
    for (int i = 0; i < len; i++) {
        float diff = v1[ i ] - v2[ i ];
        sum += diff * diff;
        if (sum > maxDistSqd)
            return maxDistSqd;
    }
    return sum;
}


/// sum of absolute value of differences (element-wise)
float sumAbsDiff( const VectorF &v1, const VectorF &v2 ) {
    assertDebug( v1.length() == v2.length() );
    float sum = 0;
    int len = v1.length();
    for (int i = 0; i < len; i++) {
        float diff = v1[ i ] - v2[ i ];
        if (diff < 0)
            diff = -diff;
        sum += diff;
    }
    return sum;
}


/// returns 1 if vectors are identical (or parallel), 0 if orthogonal 
float cosineComparison( const VectorF &v1, const VectorF &v2 ) {
    assertDebug( v1.length() == v2.length() );
    int len = v1.length();
    float sumProd = 0, v1sumSq = 0, v2sumSq = 0;
    for (int i = 0; i < len; i++) {
        float v1val = v1[ i ];
        float v2val = v2[ i ];
        sumProd += v1val * v2val;
        v1sumSq += v1val * v1val;
        v2sumSq += v2val * v2val;
    }
    return sumProd / (sqrtf( v1sumSq ) * sqrtf( v2sumSq ));
}


//-------------------------------------------
// VECTOR SORTING
//-------------------------------------------


/// storage used by SortIndex
/// fix(later): this is not thread safe
const VectorD *g_sortVectD = NULL;
const VectorF *g_sortVectF = NULL;
const VectorI *g_sortVectI = NULL;


/// comparison used by SortIndex
int refCompareDouble( const void *v1, const void *v2 ) {
    int index1 = *((int *) v1);
    int index2 = *((int *) v2);
    double val1 = g_sortVectD->data( index1 );
    double val2 = g_sortVectD->data( index2 );
    if (val1 > val2)
        return 1;
    if (val1 < val2)
        return -1;
    return 0;
}


/// comparison used by SortIndex
int refCompareFloat( const void *v1, const void *v2 ) {
    int index1 = *((int *) v1);
    int index2 = *((int *) v2);
    float val1 = g_sortVectF->data( index1 );
    float val2 = g_sortVectF->data( index2 );
    if (val1 > val2)
        return 1;
    if (val1 < val2)
        return -1;
    return 0;
}


/// comparison used by SortIndex
int refCompareInt( const void *v1, const void *v2 ) {
    int index1 = *((int *) v1);
    int index2 = *((int *) v2);
    int val1 = g_sortVectI->data( index1 );
    int val2 = g_sortVectI->data( index2 );
    if (val1 > val2)
        return 1;
    if (val1 < val2)
        return -1;
    return 0;
}


/// returns indices of sorted elements (small to large)
VectorI sortIndex( const VectorD &v ) {
    g_sortVectD = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareDouble );
    return index;
}


/// returns indices of sorted elements (small to large)
VectorI sortIndex( const VectorF &v ) {
    g_sortVectF = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareFloat );
    return index;
}


/// returns indices of sorted elements (small to large)
VectorI sortIndex( const VectorI &v ) {
    g_sortVectI = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareInt );
    return index;
}


// fix(later): use reverse compare?
/// returns indices of elements, reverse sorted
VectorI reverseSortIndex( const VectorD &v ) {
    g_sortVectD = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareDouble );
    
    // reverse the ordering
    VectorI revIndex( len );
    for (int i = 0; i < len; i++)
        revIndex[ i ] = index[ len - i - 1 ];
    return revIndex;
}


// fix(later): use reverse compare?
/// returns indices of elements, reverse sorted
VectorI reverseSortIndex( const VectorF &v ) {
    g_sortVectF = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareFloat );
    
    // reverse the ordering
    VectorI revIndex( len );
    for (int i = 0; i < len; i++)
        revIndex[ i ] = index[ len - i - 1 ];
    return revIndex;
}


// fix(later): use reverse compare?
/// returns indices of elements, reverse sorted
VectorI reverseSortIndex( const VectorI &v ) {
    g_sortVectI = &v;
    int len = v.length();
    VectorI index( len );
    for (int i = 0; i < len; i++) 
        index[ i ] = i;
    qsort( index.dataPtr(), len, sizeof(int), refCompareInt );
    
    // reverse the ordering
    VectorI revIndex( len );
    for (int i = 0; i < len; i++)
        revIndex[ i ] = index[ len - i - 1 ];
    return revIndex;
}


/// reverse sorts a vector (from high to low)
void reverseSort( VectorF &v ) {
    VectorF sort( v );
    sort.sort();
    int len = v.length();
    for (int i = 0; i < len; i++) 
        v[ i ] = sort[ len - i - 1 ];
}


//-------------------------------------------
// VECTOR CREATION / CONVERSION
//-------------------------------------------


/// convert vector to float vector
template <typename T> VectorF toFloat( const Vector<T> &v ) {
    VectorF result( v.length() );
    for (int i = 0; i < v.length(); i++) 
        result[ i ] = (float) v[ i ];
    return result;
}
template VectorF toFloat( const Vector<int> &v );
template VectorF toFloat( const Vector<double> &v );


/// convert vector to double vector
template <typename T> VectorD toDouble( const Vector<T> &v ) {
    VectorD result( v.length() );
    for (int i = 0; i < v.length(); i++) 
        result[ i ] = (double) v[ i ];
    return result;
}
template VectorD toDouble( const Vector<float> &v );
template VectorD toDouble( const Vector<double> &v );


/// create sequence of integers in [minItem, maxItem]
VectorI sequenceI( int minItem, int maxItem ) {
    int length = maxItem - minItem + 1;
    assertDebug( minItem <= maxItem );
    VectorI result( length );
    for (int i = 0; i < length; i++) 
        result[ i ] = i + minItem;
    return result;
}


/// reverse the order of the elements in the vector
VectorF reverse( const VectorF &v ) {
    int len = v.length(); 
    VectorF result( len );
    for (int i = 0; i < len; i++) 
        result[ i ] = v[ len - i - 1 ];
    return result;
}


/// re-order the elements using the given permutation vector
template <typename T> Vector<T> reorder( const Vector<T> &v, const VectorI &order ) {
    int len = v.length();
    Vector<T> result( len );
    for (int i = 0; i < len; i++) 
        result[ i ] = v[ order[ i ] ];
    return result;
}
template Vector<int> reorder( const Vector<int> &v, const VectorI &order );
template Vector<float> reorder( const Vector<float> &v, const VectorI &order );


/// extract a contiguous subset of a vector
template <typename T> Vector<T> subset( const Vector<T> &v, int minIndex, int maxIndex ) {
    Vector<T> result( maxIndex - minIndex + 1 );
    for (int i = minIndex; i <= maxIndex; i++) 
        result[ i - minIndex ] = v[ i ];
    return result;
}
template Vector<int> subset( const Vector<int> &v, int minIndex, int maxIndex );
template Vector<float> subset( const Vector<float> &v, int minIndex, int maxIndex );
template Vector<double> subset( const Vector<double> &v, int minIndex, int maxIndex );


/// extract a subset consisting of each item where mask is non-zero
VectorF subset( const VectorF &v, const VectorI &mask ) {
    assertDebug( v.length() == mask.length() );
    VectorF result;
    for (int i = 0; i < v.length(); i++) {
        if (mask[ i ])
            result.append( v[ i ] );
    }
    return result;
}


/// extract a subset by stepping over the vector
VectorF stepSubset( const VectorF &v, int step ) {
    int len = v.length();
    VectorF result;
    for (int i = 0; i < len; i += step)
        result.append( v[ i ] );
    return result;
}


/// create a vector of random elements sampled uniformly between min and max
VectorF randomVectorF( int len, float min, float max ) {
    assertDebug( len );
    VectorF result( len );
    for (int i = 0; i < len; i++) 
        result[ i ] = randomFloat( min, max );
    return result;
}


/// create a random permutation of integers in [0, len - 1]
VectorI randomPermutation( int len ) {
    VectorI result( len );
    for (int i = 0; i < len; i++)
        result[ i ] = i;
    for (int i = 0; i < len - 1; i++) {
        int index = randomInt( i, len - 1 );

        // swap result[ i ] and result[ index ]
        int val = result[ i ];
        result[ i ] = result[ index ];
        result[ index ] = val;
    }
    return result;
}


/// apply a gaussan smoothing kernel to the vector values
template <typename T> Vector<T> gaussFilter( const Vector<T> &v, float sigma ) {
    int len = v.length();
    Vector<T> result( len );

    // generate Gaussian table
    int tableRadius = (int) sigma * 3;
    double gFactor = gaussFactor( sigma );
    int tableSize = tableRadius * 2 + 1;
    double *table = new double[ tableSize ];
    for (int i = 0; i < tableSize; i++) {
        double x = (double) (i - tableRadius);
        table[ i ] = gauss( x * x, gFactor );
    }

    // apply filter
    // fix(faster): make faster
    for (int i = 0; i < len; i++) {
        double sum = 0, sumWt = 0;
        for (int j = -tableRadius; j <= tableRadius; j++) {
            int srcInd = i + j;
            if (srcInd >= 0 && srcInd < len) {
                double wt = table[ j + tableRadius ];
                sum += wt * v[ srcInd ];
                sumWt += wt;
            }
        }
        result[ i ] = (T) (sum / sumWt);
    }

    // clean up
    delete [] table;
    return result;
}
template Vector<float> gaussFilter( const Vector<float> &v, float sigma );
template Vector<double> gaussFilter( const Vector<double> &v, float sigma );


/// appply a bilateral filter to the data
VectorF bilateralFilter( const VectorF &v, float timeSigma, float valueSigma ) {
    int len = v.length();
    VectorF result( len );
    double timeFactor = gaussFactor( timeSigma );
    double valueFactor = gaussFactor( valueSigma );

    // generate Gaussian table
    int tableRadius = (int) timeSigma * 3;
    int tableSize = tableRadius * 2 + 1;
    double *table = new double[ tableSize ];
    for (int i = 0; i < tableSize; i++) {
        double x = (double) (i - tableRadius);
        table[ i ] = gauss( x * x, timeFactor );
    }

    // apply filter
    for (int i = 0; i < len; i++) {
        double sum = 0, sumWt = 0;
        double vCent = v[ i ];
        for (int j = -tableRadius; j <= tableRadius; j++) {
            int srcInd = i + j;
            if (srcInd >= 0 && srcInd < len) {
                double wt = table[ j + tableRadius ];
                double vDiff = v[ srcInd ] - vCent; // compute diff vs. center value
                wt *= gauss( vDiff * vDiff, valueFactor ); // add value weight
                sum += wt * v[ srcInd ];
                sumWt += wt;
            }
        }
        result[ i ] = (float) (sum / sumWt);
    }

    // clean up
    delete [] table;
    return result;
}


/// compute a gaussian smoothing kernel
VectorF gaussKernel( int tableRadius, float sigma ) {
    int tableSize = tableRadius * 2 + 1;
    float factor = gaussFactor( sigma );
    VectorF kernel( tableSize );
    for (int i = 0; i < tableSize; i++) {
        float x = (float) (i - tableRadius);
        kernel[ i ] = gauss( x * x, factor );
    }
    return kernel;
}


/// split a comma-delimited string of numbers
VectorI splitNumbersI( const String &s ) {
    Array<String> split = s.split( "," );
    VectorI result( split.count() );
    for (int i = 0; i < split.count(); i++) 
        result[ i ] = split[ i ].strip().toInt(); // fix(clean): make a to<type> template
    return result;
}


/// split a comma-delimited string of numbers
VectorF splitNumbersF( const String &s ) {
    Array<String> split = s.split( "," );
    VectorF result( split.count() );
    for (int i = 0; i < split.count(); i++) 
        result[ i ] = split[ i ].strip().toFloat(); // fix(clean): make a to<type> template
    return result;
}


/// create a comma-delimited string of numbers suitable for display
String toString( const VectorF &v, int minIndex, int maxIndex, const String &format ) {
    Array<String> strArr;
    minIndex = bound( minIndex, 0, v.length() - 1 );
    maxIndex = bound( maxIndex, 0, v.length() - 1 );
    for (int i = minIndex; i <= maxIndex; i++) 
        strArr.append( new String( sprintF( format.c_str(), v[ i ] )));
    return join( strArr, "," );
}


/// create a comma-delimited string of numbers suitable for display
String toString( const VectorI &v, int minIndex, int maxIndex ) {
    Array<String> strArr;
    minIndex = bound( minIndex, 0, v.length() - 1 );
    maxIndex = bound( maxIndex, 0, v.length() - 1 );
    for (int i = minIndex; i <= maxIndex; i++) 
        strArr.append( new String( sprintF( "%d", v[ i ] )));
    return join( strArr, "," );
}


//-------------------------------------------
// STATISTICS
//-------------------------------------------


/// compute standard deviation
template <typename T> T stDev( const Vector<T> &v, T mean ) {
    int len = v.length();
    T sumSqDiff = 0;
    for (int i = 0; i < len; i++) {
        T diff = v[ i ] - mean;
        sumSqDiff += diff * diff;
    }
    return (T) sqrt( (double) sumSqDiff / (double) (len - 1.0));
}
template float stDev( const Vector<float> &v, float mean );
template double stDev( const Vector<double> &v, double mean );


/// correlation between vectors
float correlationCoef( const VectorF &v1, const VectorF &v2 ) { 
    assertDebug( v1.length() == v2.length() );
    float len = (float) v1.length();
    float mean1 = v1.mean();
    float mean2 = v2.mean();
    float stDev1 = stDev( v1, mean1 );
    float stDev2 = stDev( v2, mean2 );
    if (stDev1 < 0.0001 || stDev2 < 0.0001)
        return 0;
    return (dot( v1, v2 ) - len * mean1 * mean2) / ((len - 1.0f) * stDev1 * stDev2);
}


/// correletion between vector and self, for each offset up to the max offset
// see: http://en.wikipedia.org/wiki/Correlogram#Estimation_of_autocorrelations
VectorF autoCorrelation( const VectorF &v, int maxLag ) {
    int len = v.length();
    float mean = v.mean();
    VectorF autoCor( maxLag + 1 );
    autoCor[ 0 ] = 1.0f;

    // compute sum variance
    float sumVar = 0;
    for (int i = 0; i < len; i++) {
        float diff = v[ i ] - mean;
        sumVar += diff * diff;
    }

    // compute auto-correlation for each lag
    for (int lag = 1; lag <= maxLag; lag++) {
        float sumAutoCov = 0;
        for (int i = 0; i < len - lag; i++) 
            sumAutoCov += (v[ i ] - mean) * (v[ i + lag ] - mean);
        autoCor[ lag ] = sumAutoCov / sumVar;
    }
    return autoCor;
}


/// compute entropy of vector values (assumes sum of v is 1)
double entropy( const VectorD &v ) {
    double sum = 0;
    for (int i = 0; i < v.length(); i++) {
        double val = v[ i ];
        if (val > 1e-10)
            sum -= val * log( val );
    }
    return sum;
}


/// compute median value
template <typename T> T median( const Vector<T> &v ) {
    assertDebug( v.length() );
    VectorI sortInd = sortIndex( v );
    return v[ sortInd[ v.length() / 2 ] ];
}
template float median( const Vector<float> &v );
template double median( const Vector<double> &v );


/// compute median absolute difference from the median
float medianDeviation( const VectorF &v ) {
    int len = v.length();
    float m = median( v );
    VectorF absDiff( len );
    for (int i = 0; i < len; i++)
        absDiff[ i ] = fAbs( v[ i ] - m );
    return median( absDiff );
}


/// computes min, 1st quantile, median, 3rd quantile, and max of the vector v
// fix(later): off-by-one errors?
template <typename T> void quantiles( const Vector<T> &v, T &min, T &per25, T &median, T &per75, T &max ) {
    assertDebug( v.length() );
    VectorI sortInd = sortIndex( v );
    int len = v.length();
    min = v[ sortInd[ 0 ] ];
    per25 = v[ sortInd[ len / 4 ] ];
    median = v[ sortInd[ len / 2 ] ];
    per75 = v[ sortInd[ len * 3 / 4 ] ];
    max = v[ sortInd[ len - 1 ] ];
}
template void quantiles<int>( const Vector<int> &v, int &min, int &per25, int &median, int &per75, int &max );
template void quantiles<float>( const Vector<float> &v, float &min, float &per25, float &median, float &per75, float &max );


/// assign all values above median to 1, below (or equal) to -1 
void medianSplit( VectorF &v ) {
    assertDebug( v.length() );
    int len = v.length();
    VectorI sortInd = sortIndex( v );
    float thresh = v[ sortInd[ len / 2 ] ];
    disp( 2, "median split thresh: %f, min: %f, mean: %f, max: %f", thresh, v.min(), v.mean(), v.max() );
    for (int i = 0; i < len; i++) {
        if (v[ i ] > thresh)
            v[ i ] = 1;
        else
            v[ i ] = -1;
    }
}


/// assign all values above given percentile to 1, below (or equal) to -1 
void percentileThresh( VectorF &v, float frac ) {
    assertDebug( v.length() );
    int len = v.length();
    VectorI sortInd = sortIndex( v );
    int splitIndex = round( (float) len * frac );
    splitIndex = bound( splitIndex, 0, len - 1 );
    float thresh = v[ sortInd[ splitIndex ] ];
    for (int i = 0; i < len; i++) {
        if (v[ i ] > thresh)
            v[ i ] = 1;
        else
            v[ i ] = -1;
    }
}


// returns number of non-zero entries
template <typename T> int nonZeroCount( const Vector<T> &v ) {
    int nonZeroCount = 0, len = v.length();
    for (int i = 0; i < len; i++)
        if (v[ i ])
            nonZeroCount++;
    return nonZeroCount;
}
template int nonZeroCount( const Vector<int> &v );
template int nonZeroCount( const Vector<float> &v );
template int nonZeroCount( const Vector<double> &v );


/// returns number of entries within specified bounds 
int countInBounds( const VectorF &v, float min, float max ) {
    int len = v.length();
    int count = 0;
    for (int i = 0; i < len; i++)
        if (v[ i ] >= min && v[ i ] <= max)
            count++;
    return count;
}


/// returns index of minimum element
int argMin( VectorD &v ) {
    assertDebug( v.length() );
    int minIndex = 0;
    double minValue = v[ 0 ];
    for (int i = 1; i < v.length(); i++) {
        double val = v[ i ];
        if (val < minValue) {
            minValue = val;
            minIndex = i;
        }
    }
    return minIndex;
}



/// returns index of minimum element
int argMax( VectorD &v ) {
    assertDebug( v.length() );
    int maxIndex = 0;
    double maxValue = v[ 0 ];
    for (int i = 1; i < v.length(); i++) {
        double val = v[ i ];
        if (val > maxValue) {
            maxValue = val;
            maxIndex = i;
        }
    }
    return maxIndex;
}


/// returns index of nearest local minimum
int nearestMinIndex( const VectorF &v, int startIndex ) {
    assertDebug( v.length() );
    assertDebug( startIndex >= 0 && startIndex < v.length() );
    int index = startIndex;
    float val = v[ index ];
    int len = v.length();
    while (true) {
        if (index > 0 && v[ index - 1 ] < val) {
            index--;
        } else if (index + 1 < len && v[ index + 1 ] < val) {
            index++;
        } else {
            break;
        }
        val = v[ index ];
    }
    return index;
}


/// returns index of nearest local maximum
int nearestMaxIndex( const VectorF &v, int startIndex ) {
    assertDebug( v.length() );
    assertDebug( startIndex >= 0 && startIndex < v.length() );
    int index = startIndex;
    float val = v[ index ];
    int len = v.length();
    while (true) {
        if (index > 0 && v[ index - 1 ] > val) {
            index--;
        } else if (index + 1 < len && v[ index + 1 ] > val) {
            index++;
        } else {
            break;
        }
        val = v[ index ];
    }
    return index;
}


//-------------------------------------------
// VECTOR STATS CLASS
//-------------------------------------------


/// the VectorStats class represents a set of stats about a particular vector
VectorStats::VectorStats( const VectorF &v ) {
    m_totalCount = v.length();
    m_posCount = 0;
    m_negCount = 0;
    m_zeroCount = 0;
    for (int i = 0; i < v.length(); i++) {
        float val = v[ i ];
        if (val > 0)
            m_posCount++;
        else if (val < 0)
            m_negCount++;
        else
            m_zeroCount++;
    }
    m_mean = v.mean();
    m_stDev = sbl::stDev( v, (float) m_mean );
    float vMin = 0, vPer25 = 0, vMedian = 0, vPer75 = 0, vMax = 0;
    quantiles( v, vMin, vPer25, vMedian, vPer75, vMax );
    m_min = vMin;
    m_per25 = vPer25;
    m_median = vMedian;
    m_per75 = vPer75;
    m_max = vMax;
}


/// the VectorStats class represents a set of stats about a particular vector
VectorStats::VectorStats( const VectorD &v ) {
    m_totalCount = v.length();
    m_posCount = 0;
    m_negCount = 0;
    m_zeroCount = 0;
    for (int i = 0; i < v.length(); i++) {
        double val = v[ i ];
        if (val > 0)
            m_posCount++;
        else if (val < 0)
            m_negCount++;
        else
            m_zeroCount++;
    }
    m_mean = v.mean();
    m_stDev = sbl::stDev( v, m_mean );
    double vMin = 0, vPer25 = 0, vMedian = 0, vPer75 = 0, vMax = 0;
    quantiles( v, vMin, vPer25, vMedian, vPer75, vMax );
    m_min = vMin;
    m_per25 = vPer25;
    m_median = vMedian;
    m_per75 = vPer75;
    m_max = vMax;
}


/// the VectorStats class represents a set of stats about a particular vector
VectorStats::VectorStats( const VectorI &v ) {
    m_totalCount = v.length();
    m_posCount = 0;
    m_negCount = 0;
    m_zeroCount = 0;
    for (int i = 0; i < v.length(); i++) {
        int val = v[ i ];
        if (val > 0)
            m_posCount++;
        else if (val < 0)
            m_negCount++;
        else
            m_zeroCount++;
    }
    m_mean = v.mean();
    m_stDev = sbl::stDev( toDouble( v ), m_mean );
    int vMin = 0, vPer25 = 0, vMedian = 0, vPer75 = 0, vMax = 0;
    quantiles( v, vMin, vPer25, vMedian, vPer75, vMax );
    m_min = vMin;
    m_per25 = vPer25;
    m_median = vMedian;
    m_per75 = vPer75;
    m_max = vMax;
}


//-------------------------------------------
// TEST CODE
//-------------------------------------------


// test the normalize function
bool testNormalize() {
    VectorF v = randomVectorF( 10, 1, 2 );
    normalize( v );
    float len1 = sqrtf( dot( v, v ) );
    unitAssert( fAbs( len1 - 1.0f ) < 1e-6 );
    normalize( v, 2.0 );
    float len2 = sqrtf( dot( v, v ) );
    unitAssert( fAbs( len2 - 2.0f ) < 1e-6 );
    v.clear( 0 );
    normalize( v );
    float len3 = sqrtf( dot( v, v ) );
    unitAssert( fAbs( len3 ) < 1e-6 );
    return true;
}


// test the sortIndex function
bool testSortIndex() {
    VectorF v1( 1.0f, 2.0f, 3.0f, 4.0f, 5.0f );
    VectorI sortIndex1 = sortIndex( v1 );
    unitAssert( sortIndex1 == VectorI( 0, 1, 2, 3, 4 ) );
    VectorF v2( 10.0f, 8.0f, 0.0f, -4.0f, -5.0f );
    VectorI sortIndex2 = sortIndex( v2 );
    unitAssert( sortIndex2 == VectorI( 4, 3, 2, 1, 0 ) );
    VectorF v3( 10.0f, 28.0f, 0.0f, 4.0f, 5.0f );
    VectorI sortIndex3 = sortIndex( v3 );
    unitAssert( sortIndex3 == VectorI( 2, 3, 4, 0, 1 ) );
    return true;
}


// test the sequenceI function
bool testSequenceI() {
    VectorI v1( 1, 2, 3, 4, 5 );
    unitAssert( v1 == sequenceI( 1, 5 ) );
    VectorI v2( 0, 1, 2, 3, 4 );
    unitAssert( v2 == sequenceI( 0, 4 ) );
    VectorI v3( -2, -1, 0, 1, 2 );
    unitAssert( v3 == sequenceI( -2, 2 ) );
    VectorI v4( 10, 11, 12, 13, 14 );
    unitAssert( v4 == sequenceI( 10, 14 ) );
    VectorI seq = sequenceI( 0, 99 );
    unitAssert( seq.length() == 100 );
    for (int i = 0; i < 100; i++)
        unitAssert( seq[ i ] == i );
    return true;
}


// register commands, etc. defined in this module
void initVectorUtil() {
    registerUnitTest( testNormalize );
    registerUnitTest( testSortIndex );
    registerUnitTest( testSequenceI );
}


} // end namespace sbl

