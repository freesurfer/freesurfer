
#ifndef H_MISC_MATHS_H
#define H_MISC_MATHS_H

#include <cmath>
#include <vector>

// defines a Gaussian kernel with given variance
// for now, only automatic trimming
std::vector<double> gauss_kernel(double, int max_len=-1);

// applies a kernel (odd len vector) to the 1D signal vin
std::vector<double> smooth_kernel(std::vector<double> vin,
                                  std::vector<double> vker);

template<class T>
T average(const std::vector<T>& data);

template<class T>
T stdev(const std::vector<T>& data);

//---------------------------------

template<class T>
T average(const std::vector<T>& data)
{
  T sum = (T)0;
  int count = 0;
  for ( typename std::vector<T>::const_iterator cit = data.begin();
        cit != data.end(); cit++ )
  {
    sum += *cit;
    count++;
  }
  if ( count > 0 )
    return sum/count;

  return sum;
}

template<class T>
T stdev(const std::vector<T>& data)
{
  T sum = average(data);
  T sum2 = (T)0;
  int count = 0;
  for ( typename std::vector<T>::const_iterator cit = data.begin();
        cit != data.end(); cit++ )
  {
    T tbuf = *cit - sum;
    sum2 += tbuf * tbuf;
    count++;
  }
  if ( count > 0 )
    return sqrt( sum2/count );

  return sum2;
}

//---------------------------------------------
//
// Functor returns True if condition is yes
//          false if condition is false
//
// the returned value is the maximum for which
//     the condition is the same as on the left side
template<class Functor>
double FindLeftZeroCrossing(const Functor& f,
                            double dmin,
                            double dmax,
                            double deps = 1.0e-3)
{
  bool bConditionToKeep = f(dmin);
  if ( f(dmax) == bConditionToKeep )
    throw std::string("FindLeftZeroCrossing - incorrect boundaries");

  double dmid;
  unsigned int noIterations=0, maxIterations=20;

  while ( dmax - dmin > deps && ++noIterations < maxIterations )
  {
    dmid = .5 * (dmin + dmax);
    if ( f(dmid) == bConditionToKeep ) dmin = dmid;
    else dmax = dmid;
  } // next iteration

  if ( dmax - dmin > deps )
    throw std::string( "No convergence finding optimal point");

  return dmin;
}


#endif

