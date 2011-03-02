/**
 * @file  Utilities.cpp
 * @brief Collection of utilities
 *
 * Some utilities that don't fit in anywhere else.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.12 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#include "string_fixed.h"
#include <stdexcept>
#ifdef __cplusplus
extern "C"
{
#endif
#include <math.h>
#include <stdlib.h> // abs
#ifdef __cplusplus
}
#endif
#include "Utilities.h"
#include "VectorOps.h"

using namespace std;

void
Utilities::FindPointsOnLine2d ( int const iPointA[2],
				int const iPointB[2],
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

  if ( ax > ay ) {

    int d = ay - (ax / 2);
    while ( cur[0] != iPointB[0] ) {

      Point2<int> point( cur[0], cur[1] );
      oPoints.push_back( point );

      if ( d >= 0 ) {
        cur[1] += sy;
        d -= ax;
      }
      cur[0] += sx;
      d += ay;
    }

  } else {

    int d = ax - (ay / 2);
    while ( cur[1] != iPointB[1] ) {

      Point2<int> point( cur[0], cur[1] );
      oPoints.push_back( point );

      if ( d >= 0 ) {
        cur[0] += sx;
        d -= ay;
      }
      cur[1] += sy;
      d += ax;
    }

  }
}

float
Utilities::DistanceFromLineToPoint3f ( Point3<float> const& iLineA,
                                       Point3<float> const& iLineB,
                                       Point3<float> const& iPoint ) {

  // If this line segment is really short, just use the point to point
  // distance function.
  if ( fabs(iLineA[0] - iLineB[0]) < 0.001 &&
       fabs(iLineA[1] - iLineB[1]) < 0.001 &&
       fabs(iLineA[2] - iLineB[2]) < 0.001 ) {
    return VectorOps::Distance( iLineB, iPoint );
  }

  // iLineA is A, iLineB is B, iPoint is P. Find line segment lengths.
  float AB = VectorOps::Distance( iLineA, iLineB );
  float AP = VectorOps::Distance( iLineA, iPoint );
  float BP = VectorOps::Distance( iLineB, iPoint );

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

float
Utilities::DistanceFromSegmentToPoint3f ( Point3<float> const& iLineA,
					  Point3<float> const& iLineB,
					  Point3<float> const& iPoint ) {

  // Find the distance from the point to the line formed by the two
  // segments. Then find the distance from the point to the segment
  // endpoints. If the distance to the line is < the distance to the
  // closest segment, return the distance to the segment point
  // instead.
  float distanceToLine = DistanceFromLineToPoint3f( iLineA, iLineB, iPoint );
  float distanceToA = VectorOps::Distance( iLineA, iPoint );
  float distanceToB = VectorOps::Distance( iLineB, iPoint );

  if ( distanceToLine < distanceToA && distanceToLine < distanceToB ) {
    return (distanceToA < distanceToB ? distanceToA : distanceToB);
  } else {
    return distanceToLine;
  }
}


// Following is written by
// Paul J. Weiss, http://www.codeproject.com/string/stringsplit.asp
int
Utilities::SplitString( string const& input,
                        string const& delimiter,
                        vector<string>& results ) {

  int iPos = 0;
  int newPos = -1;
  int sizeS2 = delimiter.size();
  int isize = input.size();

  vector<int> positions;

  newPos = input.find (delimiter, 0);

  if ( newPos < 0 ) {
    return 0;
  }

  int numFound = 0;

  while ( newPos > iPos ) {
    numFound++;
    positions.push_back(newPos);
    iPos = newPos;
    newPos = input.find (delimiter, iPos+sizeS2+1);
  }

  for ( int i=0; i <= (int)positions.size(); i++ ) {
    string s;
    if ( i == 0 ) {
      s = input.substr( i, positions[i] );
    }
    int offset = positions[i-1] + sizeS2;
    if ( offset < isize ) {
      if ( i == (int)positions.size() ) {
        s = input.substr(offset);
      } else if ( i > 0 ) {
        s = input.substr( positions[i-1] + sizeS2,
                          positions[i] - positions[i-1] - sizeS2 );
      }
    }
    if ( s.size() > 0 ) {
      results.push_back(s);
    }
  }
  return numFound;
}

