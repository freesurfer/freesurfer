/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef TOPOLOGY_GLOBALS_H
#define TOPOLOGY_GLOBALS_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cassert>
using namespace std;

#define PRINT_ERROR_MODE 0
#define PRINT_MODE 0

//#define ASSERT(exp)  assert(exp)
#define ASSERT(exp) check(exp)

//////////////////////////////////////////////////////////////////
// macros

#define __DOT(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define __CROSS(a,b,d) (d[0]=a[1]*b[2]-b[1]*a[2], \
   d[1]=a[2]*b[0]-b[2]*a[0], \
   d[2]=a[0]*b[1]-b[0]*a[1])

#ifndef __SIGN
#define __SIGN(x) (((x)>0)? 1.0 : -1.0 )
#endif
#ifndef __SQR3
#define __SQR3(x) (SQR(x[0])+SQR(x[1])+SQR(x[2]))
#endif
#ifndef __NORM3
#define __NORM3(x) (sqrt(SQR3(x)))
#endif
#ifndef __MAX
#define __MAX(a,b) (((a)<(b))?(b):(a))
#endif
#ifndef __MIN
#define __MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef __MAX3
#define __MAX3(a,b,c) (MAX(a,MAX(b,c)))
#endif
#ifndef __MIN3
#define __MIN3(a,b,c) (MIN(a,MIN(b,c)))
#endif

inline double __SQR(double x)
{
  return x*x;
}

inline double __norm(double x,double y, double z)
{
  return sqrt(x*x+y*y+z*z);
}

inline void __normalize(double &x,double &y, double &z)
{
  double n=__norm(x,y,z);
  x /= n;
  y /= n;
  z /= n;
}

inline double __dot(double x[3],double y[3])
{
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}

inline void __cross(double x[3], double y[3], double *z)
{
  z[0] = x[1]*y[2]-x[2]*y[1];
  z[1] = x[2]*y[0]-x[0]*y[2];
  z[2] = x[0]*y[1]-x[1]*y[0];
}

// a random number in the range 0 to nmax
int Random(int nmax);

void check(bool exp);

void ErrorExit(string s);

#endif
