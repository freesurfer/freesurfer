/**
 * @file  VectorOps.h
 * @brief Math operations on 3D vectors
 *
 * A collection of operations that can be done on 3D vectors
 * represented by form of Point3<floats>. Operators are global scope
 * but the named functions are int he VectorOps namepsace.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.9 $
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


#ifndef VectorOps_h
#define VectorOps_h

#include <string>
#include "Point3.h"


// Global operators
Point3<float> operator*( float s, Point3<float> const& v );
Point3<float> operator*( Point3<float> const& v, float s );
Point3<float> operator/( Point3<float> const& v, float s);
Point3<float> operator+( Point3<float> const& v, Point3<float> const& u );
Point3<float> operator-( Point3<float> const& v, Point3<float> const& u );
bool operator==( Point3<float> const& v, Point3<float> const& u );
bool operator!=( Point3<float> const& v, Point3<float> const& u );

// These are all static functions, we just use the class for the
// namespace.
class VectorOps {

public:
  
  // Length of a vector.
  static float Length ( Point3<float> const& v );

  // Distance between two points. (Not really a vector op but hey.)
  static float Distance ( Point3<float> const& p1,
			  Point3<float> const& p2 );

  // Normalize a vector and return the result. Does not modify the
  // argument.
  static Point3<float> Normalize ( Point3<float> const& u );

  // Dot product.
  static float Dot( Point3<float> const& u,
		    Point3<float> const& v );

  // Cross vector.
  static Point3<float> Cross( Point3<float> const& u,
			      Point3<float> const& v );

  // Triple scalar prodcuct of three vectors. 
  static float TripleScalar ( Point3<float> const& u,
			      Point3<float> const& v,
			      Point3<float> const& w );

  // Returns whether the two vectors are parallel.
  static bool AreVectorsParallel ( Point3<float> const& u, 
				   Point3<float> const& v );

  // Return the radians between two vectors.
  static double RadsBetweenVectors ( Point3<float> const& u,
				     Point3<float> const& v );

  // Perp dot product.
  static float PerpDotProduct ( Point3<float> const& u,
				Point3<float> const& v );

  // Constants used in our intersection functions.
  enum IntersectionResult {
    segmentInPlane, dontIntersect, segmentParallelToPlane, intersect
  };

  // Convert an intersection result to a string.
  static std::string IntersectionResultToString ( IntersectionResult iR );

  // Tests a point and a collinear segment and returns if point is in
  // the segment. Note that the point must already be collinear with
  // the segment. p1 is the point and the segment is formed from q1
  // and q2.
  static IntersectionResult
    PointInSegment ( Point3<float> const& p1,
		     Point3<float> const& q1, Point3<float> const& q2 );
  
  // Tests for intersection between a segment and a plane. The segment
  // is formed be q1 and q2, the plane is formed by the point p1 and
  // the normal n, and the intersection point is returned in
  // oIntersection if the return value is "intersect".
  static IntersectionResult
    SegmentIntersectsPlane( Point3<float> const& q1, Point3<float> const& q2,
			    Point3<float> const& p1, Point3<float> const& n,
			    Point3<float>& oIntersection );

  // Tests for intersection between two segments. The first segment is
  // formed by p1 and p2 and the second is from q1 and q2. The
  // intersection point is returned in oIntersection if the return
  // value is "intersect".
  static IntersectionResult
    SegmentIntersectsSegment( Point3<float> const& p1, Point3<float> const& p2,
			      Point3<float> const& q1, Point3<float> const& q2,
			      Point3<float>& oIntersection );
  static const float Epsilon;
};

#endif
