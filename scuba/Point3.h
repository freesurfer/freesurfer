#ifndef Point3_h
#define Point3_h

#include <iostream>


template <typename T>
class Point3 { 
 public:
  Point3 () {}
  Point3 ( T iX, T iY, T iZ ) { m[0] = iX; m[1] = iY; m[2] = iZ;}
  Point3 ( T iXYZ[3] ) { m[0] = iXYZ[0]; m[1] = iXYZ[1]; m[2] = iXYZ[2]; }
  void Set ( T iX, T iY, T iZ ) { m[0] = iX; m[1] = iY; m[2] = iZ;}
  T x() { return m[0]; }
  T y() { return m[1]; }
  T z() { return m[2]; }
  T m[3];
};

template <typename T>
std::ostream& operator << ( std::ostream&, Point3<T> iPoint  );

#endif
