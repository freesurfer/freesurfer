/**
 * @file  VectorOps.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:15 $
 *    $Revision: 1.8 $
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


#include <math.h>
#include "VectorOps.h"

const float VectorOps::Epsilon = 0.00001;

Point3<float> operator* ( float s, Point3<float>& v ) {
  return Point3<float>( v[0]*s, v[1]*s, v[2]*s );
}

#if 0
Point3<float> operator* ( float s, Point3<float> v ) {
  return Point3<float>( v[0]*s, v[1]*s, v[2]*s );
}
#endif

Point3<float> operator* ( Point3<float>& v, float s ) {
  return Point3<float>( v[0]*s, v[1]*s, v[2]*s );
}

Point3<float> operator/ ( Point3<float>& v, float s) {
  return Point3<float>( v[0]/s, v[1]/s, v[2]/s );
}

Point3<float> operator+ ( Point3<float>& v, Point3<float>& u ) {
  return Point3<float>( v[0]+u[0], v[1]+u[1], v[2]+u[2] );
}

Point3<float> operator+ ( Point3<float> v, Point3<float> u ) {
  return Point3<float>( v[0]+u[0], v[1]+u[1], v[2]+u[2] );
}

Point3<float> operator- ( Point3<float>& v, Point3<float>& u ) {
  return Point3<float>( v[0]-u[0], v[1]-u[1], v[2]-u[2] );
}

bool operator== ( Point3<float>& v, Point3<float>& u ) {
  return ( fabs(v[0] - u[0]) < VectorOps::Epsilon &&
           fabs(v[1] - u[1]) < VectorOps::Epsilon &&
           fabs(v[2] - u[2]) < VectorOps::Epsilon );
}

bool operator!= ( Point3<float>& v, Point3<float>& u ) {
  return !(v == u);
}

#if 0
Point3<float> operator- ( Point3<float> v, Point3<float> u ) {
  return Point3<float>( v[0]-u[0], v[1]-u[1], v[2]-u[2] );
}
#endif

float
VectorOps::Dot ( Point3<float>& v1, Point3<float>& v2 ) {
  return( v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] );
}

float
VectorOps::Length ( Point3<float>& v ) {
  return( sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ) );
}

float
VectorOps::Distance ( Point3<float>& v, Point3<float>& w ) {
  return( sqrt( (v[0] - w[0]) * (v[0] - w[0]) +
                (v[1] - w[1]) * (v[1] - w[1]) +
                (v[2] - w[2]) * (v[2] - w[2]) ) );
}

Point3<float>
VectorOps::Cross ( Point3<float>& u, Point3<float>& v) {
  return Point3<float>( u[1]*v[2] - u[2]*v[1],
                        u[2]*v[0] - u[0]*v[2],
                        u[0]*v[1] - u[1]*v[0] );
}

float
VectorOps::TripleScale ( Point3<float>& a, Point3<float>& b, Point3<float>& c ) {
  return( c[0]*a[1]*b[2] - c[0]*a[2]*b[1] +
          c[1]*a[2]*b[0] - c[1]*a[0]*b[2] +
          c[2]*a[0]*b[1] - c[2]*a[1]*b[0] );
}

int
VectorOps::AreVectorsParallel ( Point3<float>& u, Point3<float>& v ) {
  Point3<float> theCross = Cross( u, v );
  return( theCross[0] == 0 && theCross[1] == 0 && theCross[2] == 0 );
}

double
VectorOps::RadsBetweenVectors ( Point3<float>& u, Point3<float>& v ) {
  double dot = Dot(u,v);
  double lu  = Length(u);
  double lv  = Length(v);
  double costheta = dot / (lu*lv);
  if ( fabs( costheta - -1.0 ) < Epsilon ) return M_PI;
  return acos( costheta );
}

Point3<float>
VectorOps::Normalize ( Point3<float>& u ) {
  float length = Length(u);
  if ( 0 == length ) return u;
  return Point3<float>( u[0]/length, u[1]/length, u[2]/length );
}

VectorOps::IntersectionResult
VectorOps::PointInSegment ( Point3<float>& p,
                            Point3<float>& q1, Point3<float> q2 ) {

  if ( q1[0] != q2[0] ) {
    if ( q1[0] <= p[0] && p[0] <= q2[0] ) return intersect;
    if ( q1[0] >= p[0] && p[0] >= q2[0] ) return intersect;
  } else if ( q2[1] != q2[1] ) {
    if ( q1[1] <= p[1] && p[1] <= q2[1] ) return intersect;
    if ( q1[1] >= p[1] && p[1] >= q2[1] ) return intersect;
  } else {
    if ( q1[2] <= p[2] && p[2] <= q2[2] ) return intersect;
    if ( q1[2] >= p[2] && p[2] >= q2[2] ) return intersect;
  }
  return dontIntersect;
}


VectorOps::IntersectionResult
VectorOps::SegmentIntersectsPlane ( Point3<float>& p1, Point3<float>& p2,
                                    Point3<float>& plane, Point3<float>& n,
                                    Point3<float>& oIntersection ) {
  Point3<float> u;
  u = p2 - p1;
  Point3<float> w = p1 - plane;

  float D = Dot( n, u );
  float N = -Dot( n, w );

  if (fabs(D) < Epsilon) {          // segment is parallel to plane
    if (fabs(N) < Epsilon)          // segment lies in plane
      return segmentInPlane;
    else
      return segmentParallelToPlane; // no intersection
  }
  // they are not parallel
  // compute intersect param
  float sI = N / D;
  if (sI < 0 || sI > 1)
    return dontIntersect;               // no intersection

  oIntersection = p1 + sI * u;    // compute segment intersect point

  return intersect;
}

float
VectorOps::PerpProduct ( Point3<float>& u, Point3<float>& v ) {

  return( u[1]*v[2] - u[2]*v[1] +
          u[2]*v[0] - u[0]*v[2] +
          u[0]*v[1] - u[1]*v[0] );
}

VectorOps::IntersectionResult
VectorOps::SegmentIntersectsSegment ( Point3<float>& p1, Point3<float>& p2,
                                      Point3<float>& q1, Point3<float>& q2,
                                      Point3<float>& oIntersection ) {

  Point3<float> u = p2 - p1;
  Point3<float> v = q2 - q1;
  Point3<float> w = p1 - q1;
  float D = PerpProduct( u, v );

  // This is true if they are parallel or either one is a point.
  if ( fabs(D) < Epsilon ) {

    // This is true if they are not collinear, so no intersection.
    if ( PerpProduct(u,w) != 0 ||
         PerpProduct(v,w) != 0 ) {
      return dontIntersect;
    }

    // Check if they are points.
    float du = Dot( u, u );
    float dv = Dot( v, v );
    if ( 0 == du && 0 == dv ) {
      // Check if they are the same point. If not, no intersection, if
      // so, they intersect at the point.
      if ( p1 != q1 ) {
        return dontIntersect;
      } else {
        oIntersection = p1;
        return intersect;
      }
    }

    // p is the single point. Check if it's in q. If so, return p.
    if ( 0 == du ) {
      if ( PointInSegment( p1, q1, q2 ) == intersect ) {
        oIntersection = p1;
        return intersect;
      } else {
        return dontIntersect;
      }
    }
    // q is the single point. Check if it's in p. If so, return q.
    if ( 0 == dv ) {
      if ( PointInSegment( q1, p1, p2 ) == intersect ) {
        oIntersection = q1;
        return intersect;
      } else {
        return dontIntersect;
      }
    }

    // They are collinear. They may overlap or not.
    // Find the endpoints of p in the equation of q.
    float t0, t1;
    Point3<float> w2 = p2 - q1;
    if ( 0 != v[0] ) {
      t0 = w[0] / v[0];
      t1 = w2[0] / v[0];
    } else if ( 0 != v[1] ) {
      t0 = w[1] / v[1];
      t1 = w2[1] / v[1];
    } else {
      t0 = w[2] / v[2];
      t1 = w2[2] / v[2];
    }
    if ( t0 > t1 ) { // Make sure t0 < t1.
      float t = t0;
      t0 = t1;
      t1 = t;
    }
    if ( t0 > 1 || t1 < 0 ) {
      return dontIntersect;  // No overlap
    }
    t0 = t0 < 0 ? 0 : t0; // Clip t0 to 0.
    t1 = t1 > 1 ? 1 : t1; // Clip t1 to 1.
    if ( t0 == t1 ) { // Intersection is a point.
      oIntersection = q1 + (t0 * v);
      return intersect;
    }

    // Intersection is a segment.
    oIntersection = q1 + (t0 * v);
    // oIntersection2 = q[0] + t1 * v;
    return intersect;
  }

  // Segments are skew and may intersect. Get intersect parameter for
  // p and q and see if they intersect.
  float pI = PerpProduct( v, w ) / D;
  if ( pI < 0 || pI > 1 )
    return dontIntersect;

  float qI = PerpProduct( u, w ) / D;
  if ( qI < 0 || qI > 1 )
    return dontIntersect;

  oIntersection = p1 + (pI * u);
  return intersect;
}

std::string
VectorOps::IntersectionResultToString ( IntersectionResult iR ) {
  switch ( iR ) {
  case segmentInPlane:
    return "segmentInPlane";
    break;
  case dontIntersect:
    return "dontIntersect";
    break;
  case segmentParallelToPlane:
    return "segmentParallelToPlane";
    break;
  case intersect:
    return "intersect";
    break;
  }
  return "unknown";
}
