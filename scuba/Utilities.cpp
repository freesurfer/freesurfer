#include <stdexcept>
#include "Utilities.h"

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

