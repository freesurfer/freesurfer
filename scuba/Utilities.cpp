/**
 * @file  Utilities.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.9 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include "string_fixed.h"
#include <stdexcept>
#include <math.h>
#include "Utilities.h"
#include "VectorOps.h"

using namespace std;

void
Utilities::FindPointsOnLine2d ( int iPointA[2], int iPointB[2],
                                int,
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

void
Utilities::FindPointsOnLine3d ( int iPointA[3], int iPointB[3],
                                list<Point3<int> >& oPoints ) {

  float dx = iPointB[0] - iPointA[0];
  float dy = iPointB[1] - iPointA[1];
  float dz = iPointB[2] - iPointA[2];

  float l = sqrt( dx * dx + dy * dy + dz * dz );

  float sx = dx / l;
  float sy = dy / l;
  float sz = dz / l;

  float cur[3];
  cur[0] = iPointA[0];
  cur[1] = iPointA[1];
  cur[2] = iPointA[2];

  while ( (int)rint(cur[0]) != iPointB[0] ||
          (int)rint(cur[1]) != iPointB[1] ||
          (int)rint(cur[2]) != iPointB[2] ) {

    Point3<int> p( (int)rint(cur[0]), (int)rint(cur[1]), (int)rint(cur[2]) );
    oPoints.push_back( p );

    cur[0] += sx;
    cur[1] += sy;
    cur[2] += sz;
  }

  // Add last point.
  Point3<int> p( iPointB[0], iPointB[1], iPointB[2] );
  oPoints.push_back( p );
}

void
Utilities::FindPointsOnLine3f ( float iPointA[3], float iPointB[3],
                                list<Point3<float> >& oPoints ) {

  float dx = iPointB[0] - iPointA[0];
  float dy = iPointB[1] - iPointA[1];
  float dz = iPointB[2] - iPointA[2];

  float max = fabs(dx) > fabs(dy) ? fabs(dx) : fabs(dy);
  max = fabs(max) > fabs(dz) ? fabs(max) : fabs(dz);
  float sx = dx / max;
  float sy = dy / max;
  float sz = dz / max;


  float cur[3];
  cur[0] = iPointA[0];
  cur[1] = iPointA[1];
  cur[2] = iPointA[2];

  for ( int n = 0; n < max; n++ ) {

    Point3<float> p( cur[0], cur[1], cur[2] );
    oPoints.push_back( p );

    cur[0] += sx;
    cur[1] += sy;
    cur[2] += sz;
  }

  // Add last point.
  Point3<float> p( iPointB[0], iPointB[1], iPointB[2] );
  oPoints.push_back( p );
}

float
Utilities::DistanceFromLineToPoint3f ( Point3<float>& iLineA,
                                       Point3<float>& iLineB,
                                       Point3<float>& iPoint ) {

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
Utilities::DistanceFromSegmentToPoint3f ( Point3<float>& iLineA,
    Point3<float>& iLineB,
    Point3<float>& iPoint ) {

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
Utilities::SplitString( const string& input,
                        const string& delimiter,
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

