#include <math.h>
#include "VectorOps.h"

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
  return acos(  Dot(u,v) / ( Length(u) * Length(v) )  );
}

Point3<float> 
VectorOps::Normalize ( Point3<float>& u ) {
  float length = Length(u);
  return Point3<float>( u[0]/length, u[1]/length, u[2]/length );
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
  
  if (fabs(D) < 0.00001) {          // segment is parallel to plane
    if (fabs(N) < 0.00001)          // segment lies in plane
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
