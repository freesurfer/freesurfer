#include <math.h>
#include "VectorOps.h"

float DotVectors ( Point3<float>& v1, Point3<float>& v2 ) {
  return( v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] );
}

float Length ( Point3<float>& v ) {
  return( sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ) );
}
  
Point3<float> operator* ( float s, Point3<float>& v ) {
  return Point3<float>( v[0]*s, v[1]*s, v[2]*s );
}

Point3<float> operator* ( float s, Point3<float> v ) {
  return Point3<float>( v[0]*s, v[1]*s, v[2]*s );
}

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

Point3<float> operator- ( Point3<float> v, Point3<float> u ) {
  return Point3<float>( v[0]-u[0], v[1]-u[1], v[2]-u[2] );
}

Point3<float> CrossVectors ( Point3<float>& u, Point3<float>& v) {
  return Point3<float>( u[1]*v[2] - u[2]*v[1],
		  u[2]*v[0] - u[0]*v[2],
		  u[0]*v[1] - u[1]*v[0] );
}

float TripleScaleVectors ( Point3<float>& a, Point3<float>& b, Point3<float>& c ) {
  return( c[0]*a[1]*b[2] - c[0]*a[2]*b[1] + 
	  c[1]*a[2]*b[0] - c[1]*a[0]*b[2] +
	  c[2]*a[0]*b[1] - c[2]*a[1]*b[0] );
}

int AreVectorsParallel ( Point3<float>& u, Point3<float>& v ) {
  Point3<float> theCross = CrossVectors( u, v );
  return( theCross[0] == 0 && theCross[1] == 0 && theCross[2] == 0 );
}

double RadsBetweenVectors ( Point3<float>& u, Point3<float>& v ) {
  return acos(  DotVectors(u,v) / ( Length(u) * Length(v) )  );
}

Point3<float> NormalizeVector ( Point3<float>& u ) {
  float length = Length(u);
  return Point3<float>( u[0]/length, u[1]/length, u[2]/length );
}
