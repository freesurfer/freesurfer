#ifndef Utilities_h
#define Utilities_h

#include <list>
#include "Point2.h"
#include "Point3.h"


class Utilities {

 public:
  
  static void FindPointsOnLine2d ( int iPointA[2], int iPointB[2],
				   int iThickness,
				   std::list<Point2<int> >& oPoints );

  static float DistanceFromLineToPoint3f ( Point3<float>& iLineA, 
					   Point3<float>& iLineB,
					   Point3<float>& iPoint );
};


#endif
