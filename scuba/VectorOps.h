#include "Point3.h"

float DotVectors( Point3<float>& v1, Point3<float>& v2 );
float Length ( Point3<float>& v );
float Distance ( Point3<float>& v, Point3<float>& w );
Point3<float> operator*( float s, Point3<float>& v );
Point3<float> operator*( float s, Point3<float> v );
Point3<float> operator*( Point3<float>& v, float s );
Point3<float> operator/( Point3<float>& v, float s);
Point3<float> operator+( Point3<float>& v, Point3<float>& u );
Point3<float> operator+( Point3<float> v, Point3<float> u );
Point3<float> operator-( Point3<float>& v, Point3<float>& u );
Point3<float> operator-( Point3<float> v, Point3<float> u );
Point3<float> CrossVectors( Point3<float>& u, Point3<float>& v);
float TripleScaleVectors ( Point3<float>& a, Point3<float>& b, Point3<float>& c );
int AreVectorsParallel ( Point3<float>& u, Point3<float>& v );
double RadsBetweenVectors ( Point3<float>& u, Point3<float>& v );
Point3<float> NormalizeVector ( Point3<float>& u );
