#include "string_fixed.h"
#include <stdexcept>
#include "Utilities.h"
#include "VectorOps.h"

using namespace std;

void
Utilities::FindPointsOnLine2d ( int iPointA[2], int iPointB[2],
				int iThickness,
				list<Point2<int> >& oPoints ) {

  int dx = iPointB[0] - iPointA[0];
  int ax = abs(dx) * 2;
  int sx = dx > 0 ? 1 : dx < 0 ? -1 : 0; // sign of dx

  int dy = iPointB[1] - iPointA[1];
  int ay = abs(dy) * 2;
  int sy = dy > 0 ? 1 : dy < 0 ? -1 : 0; // sign of dy

  int cur[2];
  cur[0] = iPointA[0];
  cur[1] = iPointA[1];

  if( ax > ay ) {

    int d = ay - (ax / 2);
    while( cur[0] != iPointB[0] ) {
      
      Point2<int> point( cur[0], cur[1] );
      oPoints.push_back( point );

      if( d >= 0 ) {
	cur[1] += sy;
	d -= ax;
      } 
      cur[0] += sx;
      d += ay;
    }

  } else {

    int d = ax - (ay / 2);
    while( cur[1] != iPointB[1] ) {

      Point2<int> point( cur[0], cur[1] );
      oPoints.push_back( point );

      if( d >= 0 ) {
	cur[0] += sx;
	d -= ay;
      } 
      cur[1] += sy;
      d += ax;
    }

  }
}

float
Utilities::DistanceFromLineToPoint3f ( Point3<float>& iLineA, 
				       Point3<float>& iLineB,
				       Point3<float>& iPoint ) {
  
  // iLineA is A, iLineB is B, iPoint is P. Find line segment lengths.
  float AB = Distance( iLineA, iLineB );
  float AP = Distance( iLineA, iPoint );
  float BP = Distance( iLineB, iPoint );

  // These points form a triangle. Find the semiperimeter.
  float semiperimeter = ( AB + AP + BP ) / 2.0;

  // Area of triangle.
  float area = semiperimeter * (semiperimeter - AB) *
    (semiperimeter - AP) * (semiperimeter - BP);

  // h is height of this triangle, which the distance from P to the
  // base or the line AB.
  float h = 2 * (area / AB);

  return h;
}
