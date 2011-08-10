#ifndef _SBL_MATH_UTIL_H_
#define _SBL_MATH_UTIL_H_
#include <math.h>
namespace sbl {


/*! \file MathUtil.h
    \brief The MathUtil module provides a set of basic numeric utility functions.
*/


/// basic min/max functions 
// fix(clean): use templates; how do we name to avoid conflicts with c libraries
#undef min
#undef max
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }
inline float min( float a, float b ) { return a < b ? a : b; }
inline float max( float a, float b ) { return a > b ? a : b; }
inline double min( double a, double b ) { return a < b ? a : b; }
inline double max( double a, double b ) { return a > b ? a : b; }
inline float min( float a, float b, float c ) { return min( a, min( b, c )); }
inline float max( float a, float b, float c ) { return max( a, max( b, c )); }
inline int min( int a, int b, int c ) { return min( a, min( b, c )); }
inline int max( int a, int b, int c ) { return max( a, max( b, c )); }


/// bound the input to lie within the given range [min, max]
inline int bound( int v, int min, int max ) { if (v < min) v = min; if (v > max) v = max; return v; }
inline float bound( float v, float min, float max ) { if (v < min) v = min; if (v > max) v = max; return v; }
inline double bound( double v, double min, double max ) { if (v < min) v = min; if (v > max) v = max; return v; }


/// rounds a float or double to the nearest int
inline int round( float val ) { return val > 0 ? (int) (val + 0.5f) : (int) (val - 0.5f); }
inline int round( double val ) { return val > 0 ? (int) (val + 0.5) : (int) (val - 0.5); }


/// absolute value
inline int iAbs( int a ) { return a > 0 ? a : -a; }
inline float fAbs( float a ) { return a > 0 ? a : -a; }
inline double dAbs( double a ) { return a > 0 ? a : -a; }


/// convert integer to bool
inline bool makeBool( int v ) { return v ? true : false; }


/// convert degrees to radians
inline float degToRad( float v ) { return v * 3.1415926535f / 180.0f; }


/// compute non-normalized value of Gaussian density function
inline float gauss( float distSqd, float gaussFactor ) { return expf( distSqd * gaussFactor ); }
inline double gauss( double distSqd, double gaussFactor ) { return exp( distSqd * gaussFactor ); }


/// compute factor used by gauss()
inline float gaussFactor( float sigma ) { return -0.5f / (sigma * sigma); }
inline double gaussFactor( double sigma ) { return -0.5 / (sigma * sigma); }


/// swap two values using an intermediate value
template <typename T> inline void swap( T &v1, T &v2 ) { T t = v2; v2 = v1; v1 = t; }


/// returns a random float between 0.0 and 1.0;
/// note: this isn't a very good random number
float randomFloat();


/// returns a random float between the specified bounds;
/// note: this isn't a very good random number
float randomFloat( float min, float max );


/// returns a random integer within specified bounds;
/// note: this isn't a very good random number
int randomInt( int min, int max );


/// returns 1 or -1;
/// note: this isn't a very good random number
inline int randomSign() { return randomInt( 0, 1 ) ? 1 : -1; }


/// set the random number generator seed
void randomSeed( int seed );


} // end namespace sbl
#endif // _SBL_MATH_UTIL_H_

