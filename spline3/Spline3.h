/**
 * @brief A simple class representing cubic splines
 *
 */

/*
 * Original Author: Martin Reuter
 * Nov 4, 2014
 *
 */
 
#ifndef Spline3_H
#define Spline3_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cassert>

/** \class Spline3
 * \brief A simple class representing cubic splines
 * This class is designed to perform fast evaluations. It is possible to
 * chache the knot vector (x) for spline creation and also a new knot
 * vector (xnew) for resampling. Then for a given y data, the spline can 
 * be created quickly and the evaluation at xnew is also fast to 
 * generate ynew via eval()
 *
 * Example:
 *
 *   Spline3 S;
 *   S.preCacheX(x);
 *   S.preCacheXnew(xnew);
 *   
 *   // in a loop with lots of different y data:
 *   S.interp(y); // using cached x
 *   ynew = S.eval(y);  // using cached xnew, pass reference to y again (not cached earlier)
 *
 */
template <class T>
class Spline3
{
public:
  Spline3();
      
  //! Construct a cubic spline from x and y data in a single step
  void interp( const std::vector < T > & x , const std::vector < T > & y );
  void interp( unsigned int N, const T x[] , const T y[] );

  //! Cache everything related to x knots
  void preCacheX( const std::vector < T > & x );
  void preCacheX( unsigned int N, const T x[] );
  
  //! Construct a cubic spline from y (x information needs to be cached before)
  void interp( const std::vector < T > & y );
  void interp( unsigned int N, const T y[] );
    
    
  //! Evaluate spline at locations in xnew in a single step
  const std::vector < T > & eval( const std::vector <T > & xnew , const std::vector <T > & yold );
  T* eval( unsigned int M, const T xnew[] , unsigned int N, const T yold[], T ynew[]  );

  //! Cache everything related to xnew knots
  void preCacheXnew( const std::vector <T > & xnew );
  void preCacheXnew( unsigned int M, const T xnew[] );

  //! Evaluate spline at locations in xnew (from cache)
  const std::vector < T > & eval( const std::vector <T > & yold );
  T* eval( unsigned int N, const T yold[], T ynew[] );
  
  const std::vector < T > & getYPP(){return ypp;};
  
private:

  //! Computes and caches interval lengths
  void computeHi( const std::vector <T > & x);
  
  //! Computes and caches tri-diagnoal matrix M
  void computeM( const std::vector <T > & x);
  
  //! Computes and caches values needed to solve M x = b
  void cacheMi();
  
  //! Computes 2nd derivatives ypp via M ypp = b (based on cached values M and b)
  void computeYpp();

  //! Number of data points for interpolation
  unsigned int n;
  
  //! Cache (computed from x in preCacheX and subfunctions)
  std::vector < T > xc;   // size n   (allocated and set in preCacheX as copy of input x)
  std::vector < T > h;    // size n-1 (allocated and set in computeHi)
  std::vector < T > hi;   // size n-1 (allocated and set in computeHi)
  std::vector < T > M0;   // size n   (allocated and set in computeM)
  std::vector < T > M1;   // size n   (allocated and set in computeM)
  std::vector < T > M2;   // size n   (allocated and set in computeM)
  //std::vector < T > Mi1;  // size n   (allocated and set in cacheMi)
  std::vector < T > xm;   // size n   (allocated and set in cacheMi)
    
  //! Cache (computed from y, allocated in preCacheX, set in interp)
  std::vector < T > yd;   // size n-1 (allocated empty in preCacheX)
  std::vector < T > ydhi; // size n-1 (allocated empty in preCacheX)
  std::vector < T > b;    // size n   (allocated empty in preCacheX)
  std::vector < T > ypp;  // size n   (allocated empty in preCacheX)
  
  //! Cache (computed from xnew in preCacheXnew)
  unsigned int m;
  std::vector < unsigned int > intervals;  // size m=length(xnew) (allocated and set in preCacheXnew)
  std::vector < T > d;    // size m  (allocated and set in preCacheXnew)
  std::vector < T > ynew; // size m  (allocated empty in preCacheXnew)

};

template <class T>
Spline3<T>::Spline3():n(0)
{}

/** Interpolates the cubic spline into x and y (by computing the Ypp). 
    The boundary condition used: spline is quadratic in the first and 
    last interval. Does not create a copy of y (although it is needed
    later for eval) as this function is time sensitive. */
template <class T>
void Spline3<T>::interp(const std::vector <T > & x , const std::vector <T > & y )
{
  preCacheX(x);
  interp(y);
}
template <class T>
void Spline3<T>::interp( unsigned int N, const T x[] , const T y[] )
{
  preCacheX(N,x);
  interp(N,y);
}

/** Function to Cache everything that depends on x. This also
    creates a copy of x (as it is needed later and we have time here)
    and also allocates the memory for everything y related. */
template <class T>
void Spline3<T>::preCacheX(const std::vector <T > & x)
{
  n = x.size();
  if (n<3)
  {
    std::cerr << "Spline3 ERROR: we need at least 3 knots" << std::endl;
    exit(1); 
  }
  xc = x;
  computeHi(x);
  computeM(x);
  cacheMi();

  // allocate space for b (right hand side) and ypp (second derivatives)
  b.resize(n);
  ypp.resize(n);
  yd.resize(n-1);
  ydhi.resize(n-1);
}

template <class T>
void Spline3<T>::preCacheX( unsigned int N, const T x[] )
{
  n = N;
  if (n<3)
  {
    std::cerr << "Spline3 ERROR: we need at least 3 knots" << std::endl;
    exit(1); 
  }
  // copy x
  xc.resize(N);
  memcpy( &xc[0], x, N*sizeof(T) );
  computeHi(xc);
  computeM(xc);
  cacheMi();

  // allocate space for b (right hand side) and ypp (second derivatives)
  b.resize(n);
  ypp.resize(n);
  yd.resize(n-1);
  ydhi.resize(n-1);
} 

/** Interpolates the cubic spline into x and y (by computing the Ypp). 
    preCacheX(x) needs to be called beforehand!
    The boundary condition used: spline is quadratic in the first and 
    last interval. Does not create a copy of y (although it is needed
    later for eval) as this function is time sensitive. */
template <class T>
void Spline3<T>::interp( const std::vector <T > & y )
{
  assert(y.size() == n);
  //bcond: spline quadratic in first and last interval
  b[0]   = 0.0;
  b[n-1] = 0.0;
  // avoid re-computation of stuff
  yd[0] = y[1] - y[0];
  ydhi[0] = yd[0] * hi[0];
  for (unsigned int i = 1; i<n-1;i++)
  {
    yd[i]   = y[i+1] - y[i];
    ydhi[i] = yd[i] * hi[i]; 
    b[i]    = ydhi[i] - ydhi[i-1];
  }
  
  // now we should have Mi (from cache) and b, time to compute ypp:
  computeYpp();
}
template <class T>
void Spline3<T>::interp( unsigned int N, const T y[] )
{
  assert(n==N);
  //bcond: spline quadratic in first and last interval
  b[0]   = 0.0;
  b[n-1] = 0.0;
  // avoid re-computation of stuff
  yd[0] = y[1] - y[0];
  ydhi[0] = yd[0] * hi[0];
  for (unsigned int i = 1; i<n-1;i++)
  {
    yd[i]   = y[i+1] - y[i];
    ydhi[i] = yd[i] * hi[i]; 
    b[i]    = ydhi[i] - ydhi[i-1];
  }
  
  // now we should have Mi (from cache) and b, time to compute ypp:
  computeYpp();
}

/** Evaluates the spline at xnew locations (interp needs to be called
    before to fit the spline). This function requires yold to be passed again,
    as it is needed for eval and we did not create a copy earlier to save time. */
template <class T>
const std::vector <T> & Spline3<T>::eval( const std::vector < T > & xnew , const std::vector < T > & yold)
{
  preCacheXnew(xnew);
  return eval(yold);
}

template <class T>
T* Spline3<T>::eval( unsigned int M, const T xnew[] , unsigned int N, const T yold[], T ynew[] )
{
  preCacheXnew(M,xnew);
  return eval(N,yold,ynew);
}


/** Caches everything related to xnew (the values used for the
    spline evaluation). It specifically sets m (size of xnew)
    and locates intervals of each xnew value in the old x data,
   as well as the offsets within that interval (d). */
template <class T>
void Spline3<T>::preCacheXnew( const std::vector < T > & xnew )
{
  m = xnew.size();
  intervals.resize(m);
  d.resize(m); 
  unsigned int i,j;
  for ( j = 0; j< m; j++)
  {
    // determine intervals
    intervals[j] = n-2;
    for ( i = 0; i < n-1; i++ )
    {
      if ( xnew[j] < xc[i+1] )
      {
        intervals[j] = i;
        break;
      }
    }    
    d[j] = xnew[j] - xc[intervals[j]];
  }  
  ynew.resize(m);
}
template <class T>
void Spline3<T>::preCacheXnew( unsigned int M, const T xnew[] )
{
  m = M;
  intervals.resize(m);
  d.resize(m); 
  unsigned int i,j;
  for ( j = 0; j< m; j++)
  {
    // determine intervals
    intervals[j] = n-2;
    for ( i = 0; i < n-1; i++ )
    {
      if ( xnew[j] < xc[i+1] )
      {
        intervals[j] = i;
        break;
      }
    }    
    d[j] = xnew[j] - xc[intervals[j]];
  }  
  ynew.resize(m);
}

/** Evaluates the spline at the xnew locations (preCacheXnew has
    to be called before and the spline needs to be constructed 
    before that via preChaceX and interp). Reference to yold needs
    to be passed again (it was not copyied earlier to save time).
    Derivatives can easily be computed here (but are not!). */
template <class T>
const std::vector <T> & Spline3<T>::eval( const std::vector < T > & yold )
{
  for (unsigned int j = 0; j< m; j++)
  {
    unsigned int & ival = intervals[j];
    T & dval = d[j];
    T & hval = h[ival];
    //T & hival = hi[ival];
    T & ydhiival = ydhi[ival];
    
    ynew[j] = yold[ival]
      + dval * ( ydhiival
             - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * hval
      + dval * ( 0.5 * ypp[ival]
      + dval * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * hval ) ) ) );
  
  // Derivatives could be computed like this:
//  yx = ( ydhiival
//    - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * hval
//    + dval * ( ypp[ival]
//    + dval * ( 0.5 * ( ypp[ival+1] - ypp[ival] ) * hival ) );
//
//  yxx = ypp[ival] + dval * ( ypp[ival+1] - ypp[ival] ) * hival;
  }
  return ynew;
}

template <class T>
T*  Spline3<T>::eval( unsigned int N, const T yold[], T ynew[] )
{
  if (ynew == NULL) ynew = (T*)malloc(m*sizeof(T));
  for (unsigned int j = 0; j< m; j++)
  {
    unsigned int & ival = intervals[j];
    T & dval = d[j];
    T & hval = h[ival];
    //T & hival = hi[ival];
    T & ydhiival = ydhi[ival];
    
    ynew[j] = yold[ival]
      + dval * ( ydhiival
             - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * hval
      + dval * ( 0.5 * ypp[ival]
      + dval * ( ( ypp[ival+1] - ypp[ival] ) / ( 6.0 * hval ) ) ) );
  
  // Derivatives could be computed like this:
//  yx = ( ydhiival
//    - ( ypp[ival+1] / 6.0 + ypp[ival] / 3.0 ) * hval
//    + dval * ( ypp[ival]
//    + dval * ( 0.5 * ( ypp[ival+1] - ypp[ival] ) * hival ) );
//
//  yxx = ypp[ival] + dval * ( ypp[ival+1] - ypp[ival] ) * hival;
  }
  return ynew;
}

/** Compute interval length (reciprocal) based on x.
    Interval h[i] is between x[i] and x[i+1] */
template <class T>
void Spline3<T>::computeHi(const std::vector <T > & x)
{
  unsigned int nm1 = n - 1;
  hi.resize(nm1);
  h.resize(nm1);
  for (unsigned int i = 0; i< nm1 ; i++)
  {
    h[i]  = x[i+1] - x[i];
    hi[i] = 1.0 / h[i];
  }
}

/** Construct tri-diagonal matrix M (M0,M1,M2) based on x.
    bcond: assuming cubic spline is qadratic over the first and last interval */
template <class T>
void Spline3<T>::computeM(const std::vector <T > & x)
{
  unsigned int nm1 = n - 1;
  M0.resize(n);
  M1.resize(n);
  M2.resize(n);
  
  // boundary conditions (quadratic over first and last interval)
  M0[0]   =  0.0;
  M1[0]   =  1.0;
  M2[0]   = -1.0;
  M0[nm1] = -1.0;
  M1[nm1] =  1.0;
  M2[nm1] =  0.0;
  
  for (unsigned int i = 1; i<nm1; i++)
  {
    M0[i] = ( x[i+1] - x[i]   ) / 6.0;
    M1[i] = ( x[i+1] - x[i-1] ) / 3.0;
    M2[i] = ( x[i]   - x[i-1] ) / 6.0;  
  }
  
}

/** Chache xm based on M (coefficients necessary for solving the system M x = b).
    also M1 is being updated. */
template <class T>
void Spline3<T>::cacheMi()
{

  xm.resize(n-1);
  //Mi1.resize(n);
  //Mi1[0] = 0.0;

  for (unsigned int i = 1; i < n-1 ; i++ )
  {
    xm[i-1] = M0[i] / M1[i-1];
    M1[i]   = M1[i] - xm[i-1] * M2[i-1];
    //Mi1[i] = 1.0 / (M1[i] - xm[i-1] * M2[i-1]);
  }
}


/** Compute second derivatives ypp at the knot points.
    This uses xm, M0, M1, M2 and b. */
template <class T>
void Spline3<T>::computeYpp()
{
  int i; // unsigned does not work, for second loop below!
  for ( i = 1; i < (int)n - 1; i++ )
  {
    b[i] = b[i] - xm[i-1] * b[i-1];
  }
  
  T xmult = M0[n-1] / M1[n-2];
  M1[n-1] = M1[n-1] - xmult * M2[n-2];
  ypp[n-1] = ( b[n-1] - xmult * b[n-2] ) / M1[n-1];
  ypp[n-2] = ( b[n-2] - M2[n-2] * ypp[n-1] ) / M1[n-2];
  for ( i = (int)n - 3; 0 <= i; i-- )
  {
    ypp[i] = ( b[i] - M2[i] * ypp[i+1] ) / M1[i];
  }
}

#endif
