#ifndef Point2_h
#define Point2_h

#include <iostream>


template <typename T>
class Point2 { 
 public:
  Point2 () {}
  Point2 ( T iX, T iY ) { m[0] = iX; m[1] = iY; }
  Point2 ( T const iXY[2] ) { m[0] = iXY[0]; m[1] = iXY[1]; }
  void Set ( T iX, T iY ) { m[0] = iX; m[1] = iY; }
  void Set ( T iXY[2] ) { m[0] = iXY[0]; m[1] = iXY[1]; }
  void Set ( Point2<T>& i ) { m[0] = i.x(); m[1] = i.y(); }
  T& operator[]( int i ) { return m[i]; }
  T* xy() { return m; }
  T x () { return m[0]; }
  T y () { return m[1]; }
  T m[2];
};


template <typename T>
std::ostream& operator << ( std::ostream&, Point2<T> iPoint  );

#endif
