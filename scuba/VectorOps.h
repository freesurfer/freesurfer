#include "Point3.h"

Point3<float> operator*( float s, Point3<float>& v );
//Point3<float> operator*( float s, Point3<float> v );
Point3<float> operator*( Point3<float>& v, float s );
Point3<float> operator/( Point3<float>& v, float s);
Point3<float> operator+( Point3<float>& v, Point3<float>& u );
Point3<float> operator+( Point3<float> v, Point3<float> u );
Point3<float> operator-( Point3<float>& v, Point3<float>& u );
//Point3<float> operator-( Point3<float> v, Point3<float> u );

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

  enum  IntersectionResult { segmentInPlane, dontIntersect,
			     segmentParallelToPlane, intersect } ;
  static IntersectionResult 
    SegmentIntersectsPlane ( Point3<float>& p1, Point3<float>& p2,
			     Point3<float>& plane, Point3<float>& n,
			     Point3<float>& oIntersection );

};
