/**
 * @file Spline3.h
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
class Spline3
{
public:
  Spline3();
      
  //! Construct a cubic spline from x and y data in a single step
  void interp( const std::vector <double > & x , const std::vector <double > & y );

  //! Cache everything related to x knots
  void preCacheX( const std::vector <double > & x );
  
  //! Construct a cubic spline from y (x information needs to be cached before)
  void interp( const std::vector <double > & y );
    
    
  //! Evaluate spline at locations in xnew in a single step
  const std::vector < double > & eval( const std::vector <double > & xnew , const std::vector <double > & yold );

  //! Cache everything related to xnew knots
  void preCacheXnew( const std::vector <double > & xnew );

  //! Evaluate spline at locations in xnew (from cache)
  const std::vector < double > & eval( const std::vector <double > & yold );
  
  const std::vector < double > & getYPP(){return ypp;};
  
private:

  //! Computes and caches interval lengths
  void computeHi( const std::vector <double > & x);
  
  //! Computes and caches tri-diagnoal matrix M
  void computeM( const std::vector <double > & x);
  
  //! Computes and caches values needed to solve M x = b
  void cacheMi();
  
  //! Computes 2nd derivatives ypp via M ypp = b (based on cached values M and b)
  void computeYpp();

  //! Number of data points for interpolation
  unsigned int n;
  
  //! Cache (computed from x in preCacheX and subfunctions)
  std::vector < double > xc;   // size n   (allocated and set in preCacheX as copy of input x)
  std::vector < double > h;    // size n-1 (allocated and set in computeHi)
  std::vector < double > hi;   // size n-1 (allocated and set in computeHi)
  std::vector < double > M0;   // size n   (allocated and set in computeM)
  std::vector < double > M1;   // size n   (allocated and set in computeM)
  std::vector < double > M2;   // size n   (allocated and set in computeM)
  //std::vector < double > Mi1;  // size n   (allocated and set in cacheMi)
  std::vector < double > xm;   // size n   (allocated and set in cacheMi)
    
  //! Cache (computed from y, allocated in preCacheX, set in interp)
  std::vector < double > yd;   // size n-1 (allocated empty in preCacheX)
  std::vector < double > ydhi; // size n-1 (allocated empty in preCacheX)
  std::vector < double > b;    // size n   (allocated empty in preCacheX)
  std::vector < double > ypp;  // size n   (allocated empty in preCacheX)
  
  //! Cache (computed from xnew in preCacheXnew)
  unsigned int m;
  std::vector < unsigned int > intervals;  // size m=length(xnew) (allocated and set in preCacheXnew)
  std::vector < double > d;    // size m  (allocated and set in preCacheXnew)
  std::vector < double > ynew; // size m  (allocated empty in preCacheXnew)

};

Spline3::Spline3():n(0)
{}

/** Interpolates the cubic spline into x and y (by computing the Ypp). 
    The boundary condition used: spline is quadratic in the first and 
    last interval. Does not create a copy of y (although it is needed
    later for eval) as this function is time sensitive. */
void Spline3::interp(const std::vector <double > & x , const std::vector <double > & y )
{
  preCacheX(x);
  interp(y);
}

/** Function to Cache everything that depends on x. This also
    creates a copy of x (as it is needed later and we have time here)
    and also allocates the memory for everything y related. */
void Spline3::preCacheX(const std::vector <double > & x)
{
  n = x.size();
  xc = x;
  if (n<3)
  {
    std::cerr << "Spline3 ERROR: we need at least 3 knots" << std::endl;
    exit(1); 
  }
  computeHi(x);
  computeM(x);
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
void Spline3::interp( const std::vector <double > & y )
{
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
const std::vector <double> & Spline3::eval( const std::vector < double > & xnew , const std::vector < double > & yold)
{
  preCacheXnew(xnew);
  return eval(yold);
}

/** Caches everything related to xnew (the values used for the
    spline evaluation). It specifically sets m (size of xnew)
    and locates intervals of each xnew value in the old x data,
   as well as the offsets within that interval (d). */
void Spline3::preCacheXnew( const std::vector < double > & xnew )
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

/** Evaluates the spline at the xnew locations (preCacheXnew has
    to be called before and the spline needs to be constructed 
    before that via preChaceX and interp). Reference to yold needs
    to be passed again (it was not copyied earlier to save time).
    Derivatives can easily be computed here (but are not!). */
const std::vector <double> & Spline3::eval( const std::vector < double > & yold )
{
  for (unsigned int j = 0; j< m; j++)
  {
    unsigned int & ival = intervals[j];
    double & dval = d[j];
    double & hval = h[ival];
    //double & hival = hi[ival];
    double & ydhiival = ydhi[ival];
    
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
void Spline3::computeHi(const std::vector <double > & x)
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
void Spline3::computeM(const std::vector <double > & x)
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
void Spline3::cacheMi()
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
void Spline3::computeYpp()
{
  int i; // unsigned does not work, for second loop below!
  for ( i = 1; i < (int)n - 1; i++ )
  {
    b[i] = b[i] - xm[i-1] * b[i-1];
  }
  
  double xmult = M0[n-1] / M1[n-2];
  M1[n-1] = M1[n-1] - xmult * M2[n-2];
  ypp[n-1] = ( b[n-1] - xmult * b[n-2] ) / M1[n-1];
  ypp[n-2] = ( b[n-2] - M2[n-2] * ypp[n-1] ) / M1[n-2];
  for ( i = (int)n - 3; 0 <= i; i-- )
  {
    ypp[i] = ( b[i] - M2[i] * ypp[i+1] ) / M1[i];
  }
}

#endif
