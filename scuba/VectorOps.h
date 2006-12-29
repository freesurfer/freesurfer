/**
 * @file  VectorOps.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.7 $
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


#ifndef VectorOps_h
#define VectorOps_h

#include <string>
#include "Point3.h"

Point3<float> operator*( float s, Point3<float>& v );
//Point3<float> operator*( float s, Point3<float> v );
Point3<float> operator*( Point3<float>& v, float s );
Point3<float> operator/( Point3<float>& v, float s);
Point3<float> operator+( Point3<float>& v, Point3<float>& u );
Point3<float> operator+( Point3<float> v, Point3<float> u );
Point3<float> operator-( Point3<float>& v, Point3<float>& u );
bool operator==( Point3<float>& v, Point3<float>& u );
bool operator!=( Point3<float>& v, Point3<float>& u );

class VectorOps {

public:

  static float Dot( Point3<float>& v1, Point3<float>& v2 );
  static float Length ( Point3<float>& v );
  static float Distance ( Point3<float>& v, Point3<float>& w );
  static Point3<float> Cross( Point3<float>& u, Point3<float>& v);
  static float TripleScale ( Point3<float>& a, Point3<float>& b,
                             Point3<float>& c );
  static int AreVectorsParallel ( Point3<float>& u, Point3<float>& v );
  static double RadsBetweenVectors ( Point3<float>& u, Point3<float>& v );
  static Point3<float> Normalize ( Point3<float>& u );

  static float PerpProduct ( Point3<float>& u, Point3<float>& v );

  enum  IntersectionResult {
    segmentInPlane, dontIntersect,
    segmentParallelToPlane, intersect
  } ;
  static std::string IntersectionResultToString ( IntersectionResult iR );

  // Tests a point and a collinear segment and returns if point is in
  // the segment. Note that the point must already be collinear with
  // the segment.
  static IntersectionResult
  PointInSegment ( Point3<float>& p1,
                   Point3<float>& q1, Point3<float> q2 );

  static IntersectionResult
  SegmentIntersectsPlane ( Point3<float>& p1, Point3<float>& p2,
                           Point3<float>& plane, Point3<float>& n,
                           Point3<float>& oIntersection );
  static IntersectionResult
  SegmentIntersectsSegment ( Point3<float>& p1, Point3<float>& p2,
                             Point3<float>& q1, Point3<float>& q2,
                             Point3<float>& oIntersection );
  static const float Epsilon;
};

#endif
