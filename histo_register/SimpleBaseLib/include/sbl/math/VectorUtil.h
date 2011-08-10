#ifndef _SBL_VECTOR_UTIL_H_
#define _SBL_VECTOR_UTIL_H_
#include <sbl/core/String.h>
#include <sbl/math/Vector.h>
namespace sbl {


/*! \file VectorUtil.h
    \brief The VectorUtil module provides various functions for creating, manipulating,
    and analyzing vectors (represented by Vector class instances).
*/


// register commands, etc. defined in this module
void initVectorUtil();


//-------------------------------------------
// VECTOR MATH OPS
//-------------------------------------------


/// multiply vector by scalar
void multiply( const VectorF &src, float val, VectorF &dest );


/// add scalar to vector
void add( const VectorF &src, float val, VectorF &dest );
void add( const VectorI &src, int val, VectorI &dest );


/// add vector to vector
void add( const VectorF &v1, const VectorF &v2, VectorF &dest );


/// subtract vector from vector
VectorF subtract( const VectorF &v1, const VectorF &v2 );


/// vector dot product
float dot( const VectorF &v1, const VectorF &v2 );


/// make vector have unit length
void normalize( VectorF &v, float scale = 1.0f );


/// scale and shift vector so that min = 0, max = 1 
// (assuming not all elements are equal)
void normalizeBounds( VectorF &v );


/// clamp all values to be within given range
void clamp( VectorF &v, float min, float max );


//-------------------------------------------
// VECTOR COMPARISON
//-------------------------------------------


/// squared euclidean distance between vectors (sum of squared differences)
int distSqd( const VectorU &v1, const VectorU &v2 );
float distSqd( const VectorF &v1, const VectorF &v2 );
float distSqd( const float *p1, const float *p2, int len );


/// squared euclidean distance between vectors (sum of squared differences);
/// stops calculating if reaches maxDistSqd
float distSqd( const VectorF &v1, const VectorF &v2, float maxDistSqd );


/// sum of absolute value of differences (element-wise)
float sumAbsDiff( const VectorF &v1, const VectorF &v2 );


/// returns 1 if vectors are identical (or parallel), 0 if orthogonal 
float cosineComparison( const VectorF &v1, const VectorF &v2 );


//-------------------------------------------
// VECTOR SORTING
//-------------------------------------------


/// returns indices of sorted elements
// fix(later): should use templates
VectorI sortIndex( const VectorD &v );
VectorI sortIndex( const VectorF &v );
VectorI sortIndex( const VectorI &v );


/// returns indices of elements, reverse sorted
// fix(later): should use templates
VectorI reverseSortIndex( const VectorD &v );
VectorI reverseSortIndex( const VectorF &v );
VectorI reverseSortIndex( const VectorI &v );


/// reverse sorts a vector (from high to low)
void reverseSort( VectorF &v);


//-------------------------------------------
// VECTOR CREATION / CONVERSION
//-------------------------------------------


/// convert vector to float vector
template <typename T> VectorF toFloat( const Vector<T> &v );


/// convert vector to double vector
template <typename T> VectorD toDouble( const Vector<T> &v );


/// create sequence of integers in [minItem, maxItem]
VectorI sequenceI( int minItem, int maxItem );


/// reverse the order of the elements in the vector
VectorF reverse( const VectorF &v );


/// re-order the elements using the given permutation vector
template <typename T> Vector<T> reorder( const Vector<T> &v, const VectorI &order );


/// extract a contiguous subset of a vector
template <typename T> Vector<T> subset( const Vector<T> &v, int minIndex, int maxIndex );


/// extract a subset consisting of each item where mask is non-zero
VectorF subset( const VectorF &v, const VectorI &mask );


/// extract a subset by stepping over the vector
VectorF stepSubset( const VectorF &v, int step );


/// create a vector of random elements sampled uniformly between min and max
VectorF randomVectorF( int len, float min, float max );


/// create a random permutation of integers in [0, len - 1]
VectorI randomPermutation( int len );


/// apply a gaussan smoothing kernel to the vector values
template <typename T> Vector<T> gaussFilter( const Vector<T> &v, float sigma );


/// appply a bilateral filter to the data
VectorF bilateralFilter( const VectorF &v, float timeSigma, float valueSigma );


/// compute a gaussian smoothing kernel
VectorF gaussKernel( int tableRadius, float sigma );


/// split a comma-delimited string of numbers
VectorI splitNumbersI( const String &s );
VectorF splitNumbersF( const String &s );


/// create a comma-delimited string of numbers suitable for display
String toString( const VectorF &v, int minIndex = 0, int maxIndex = 100, const String &format = "%4.2f" );
String toString( const VectorI &v, int minIndex = 0, int maxIndex = 100 );


//-------------------------------------------
// VECTOR STATISTICS
//-------------------------------------------


/// returns number of non-zero entries
template <typename T> int nonZeroCount( const Vector<T> &v );


/// returns number of entries within specified bounds 
int countInBounds( const VectorF &v, float min, float max );


/// compute median value
template <typename T> T median( const Vector<T> &v );


/// compute standard deviation
template <typename T> T stDev( const Vector<T> &v, T mean );


/// compute mean of absolute values of elements
float meanAbs( const VectorF &v1 );


/// correlation between vectors
float correlationCoef( const VectorF &v1, const VectorF &v2 );


/// correletion between vector and self, for each offset up to the max offset
VectorF autoCorrelation( const VectorF &v, int maxOffset );


/// compute entropy of vector values (assumes sum of v is 1)
double entropy( const VectorD &v );


/// compute median absolute difference from the median
float medianDeviation( const VectorF &v );


/// computes min, 1st quantile, median, 3rd quantile, and max of the vector v
// fix(later): off-by-one errors?
template <typename T> void quantiles( const Vector<T> &v, T &min, T &per25, T &median, T &per75, T &max );


/// assign all values above median to 1, below (or equal) to -1 
void medianThresh( VectorF &v );


/// assign all values above given percentile to 1, below (or equal) to -1 
void percentileThresh( VectorF &v, float frac );


/// returns index of minimum element
int argMin( VectorD &v );


/// returns index of maximum element
int argMax( VectorD &v );


/// returns index of nearest local minimum
int nearestMinIndex( const VectorF &v, int startIndex );


/// returns index of nearest local maximum
int nearestMaxIndex( const VectorF &v, int startIndex );


// don't want min/max macros
#undef min
#undef max


//-------------------------------------------
// VECTOR STATS CLASS
//-------------------------------------------


/// The VectorStats class represents a set of stats about a particular vector.
// fix(clean): make template?
class VectorStats {
public:

    /// compute stats
    VectorStats( const VectorF &v );
    VectorStats( const VectorD &v );
    VectorStats( const VectorI &v );
    
    /// access the stats
    inline double min() const { return m_min; }
    inline double per25() const { return m_per25; }
    inline double median() const { return m_median; }
    inline double per75() const { return m_per75; }
    inline double max() const { return m_max; }
    inline double mean() const { return m_mean; }
    inline double stDev() const { return m_stDev; }
    inline int zeroCount() const { return m_zeroCount; }
    inline int posCount() const { return m_posCount; }
    inline int negCount() const { return m_negCount; }
    inline int totalCount() const { return m_totalCount; }

private:

    /// the statistics
    double m_min;
    double m_per25;
    double m_median;
    double m_per75;
    double m_max;
    double m_mean;
    double m_stDev;
    int m_zeroCount;
    int m_posCount;
    int m_negCount;
    int m_totalCount;

    // disable copy constructor and assignment operator
    VectorStats( const VectorStats &x );
    VectorStats &operator=( const VectorStats &x );
};


} // end namespace sbl
#endif // _SBL_VECTOR_UTIL_H_

