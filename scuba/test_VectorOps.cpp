/**
 * @file  test_VectorOps.cpp
 * @brief test_VectorOps.cpp
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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


#include <sstream>
#include "VectorOps.h"
#include "Scuba-impl.h"
#include "Listener.h"

const char* Progname = "test_VectorOps";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


class VectorOpsTester {
public:
  void Test();
};

void
VectorOpsTester::Test () {

  stringstream ssError;

  try {

    Point3<float> p;
    Point3<float> q;

    p.Set( 1, 2, 3 );
    q.Set( 4, 5, 6 );
    Assert( !(p == q), "== didn't work." );
    Assert( (p != q), "!= didn't work." );
    p.Set( 1, 2, 3 );
    q.Set( 1, 2, 3 );
    Assert( (p == q), "== didn't work." );
    Assert( !(p != q), "!= didn't work." );


    const int cVals = 10;
    float vals[cVals];
    for ( int n = 0; n < cVals; n++ ) {
      vals[n] = (float)random() / (float)random();
    }

    p.Set( vals[0], vals[1], vals[2] );
    q.Set( vals[3], vals[4], vals[5] );

    float dot = VectorOps::Dot( p, q );
    Assert( (fabs( dot - (p[0] * q[0] + p[1] * q[1] + p[2] * q[2])) < 0.0001),
            "dot didn't work." );



    p.Set( 0, 0, 0 );
    q.Set( 0, 0, 3 );
    Point3<float> plane( 0, 0, 2 );
    Point3<float> n( 0, 0, 1 );
    Point3<float> x;

    VectorOps::IntersectionResult rIntersect =
      VectorOps::SegmentIntersectsPlane( p, q, plane, n, x );
    Assert( (VectorOps::intersect == rIntersect),
            "SegmentIntersectsPlane failed, incorrect result" );
    Assert( (x[0] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[1] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[2] == 2),
            "SegmentIntersectsPlane failed, x was wrong" );

    p.Set( 0, 0, 0 );
    q.Set( 0, 0, 3 );
    plane.Set( 0, 0, 4 );
    n.Set( 0, 0, 1 );
    rIntersect = VectorOps::SegmentIntersectsPlane( p, q, plane, n, x );
    Assert( (VectorOps::dontIntersect == rIntersect),
            "SegmentIntersectsPlane failed, incorrect result" );
    Assert( (x[0] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[1] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[2] == 2),
            "SegmentIntersectsPlane failed, x was wrong" );

    p.Set( 0, 0, 0 );
    q.Set( 0, 0, 1 );
    plane.Set( 0, 0, 1 );
    n.Set( 0, 1, 0 );
    rIntersect = VectorOps::SegmentIntersectsPlane( p, q, plane, n, x );
    Assert( (VectorOps::segmentInPlane == rIntersect),
            "SegmentIntersectsPlane failed, incorrect result" );
    Assert( (x[0] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[1] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[2] == 2),
            "SegmentIntersectsPlane failed, x was wrong" );

    p.Set( 0, 0, 0 );
    q.Set( 0, 0, 1 );
    plane.Set( 0, 1, 0 );
    n.Set( 0, 1, 0 );
    rIntersect = VectorOps::SegmentIntersectsPlane( p, q, plane, n, x );
    Assert( (VectorOps::segmentParallelToPlane == rIntersect),
            "SegmentIntersectsPlane failed, incorrect result" );
    Assert( (x[0] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[1] == 0),
            "SegmentIntersectsPlane failed, x was wrong" );
    Assert( (x[2] == 2),
            "SegmentIntersectsPlane failed, x was wrong" );

    p.Set( 1, 0, 0 );
    q.Set( 0, 1, 0 );
    double rads = VectorOps::RadsBetweenVectors( p, q );
    Assert( fabs(rads - M_PI/2.0) < 0.0001,
            "RadsBetweenVectors was wrong" );

    Point3<float> sq[4];
    sq[0].Set( 0, 0, 0 );
    sq[1].Set( 4, 0, 0 );
    sq[2].Set( 4, 4, 0 );
    sq[3].Set( 0, 4, 0 );

    p.Set( 1.402, 3.2, 0 );
    Point3<float> v[4];
    v[0] = sq[0] - p;
    v[1] = sq[1] - p;
    v[2] = sq[2] - p;
    v[3] = sq[3] - p;

    double angle[4];
    angle[0] = VectorOps::RadsBetweenVectors( v[0], v[1] );
    angle[1] = VectorOps::RadsBetweenVectors( v[1], v[2] );
    angle[2] = VectorOps::RadsBetweenVectors( v[2], v[3] );
    angle[3] = VectorOps::RadsBetweenVectors( v[3], v[0] );
    double sum = angle[0] + angle[1] + angle[2] + angle[3];

    Assert( fabs(sum - M_PI*2.0) < 0.0001,
            "angle sum was wrong" );


    Point3<float> p1, p2, q1, q2;
    p1.Set( 0, 0, 0 );
    p2.Set( 0, 0, 2 );
    q1.Set( -1, 0, 1 );
    q2.Set( 1, 0, 1 );
    rIntersect = VectorOps::SegmentIntersectsSegment( p1, p2, q1, q2, x );
    Assert( (rIntersect == VectorOps::intersect),
            "SegmentIntersectsSegment didn't intersect" );
    {
      stringstream ssError;
      ssError << "SegmentIntersectsSegment returned incorrect intersect: "
      << x << endl;
      Assert( (x[0] == 0 && x[1] == 0 && x[2] == 1), ssError.str() );
    }

    p1.Set( 0, 0, 0 );
    p2.Set( 0, 0, 2 );
    q1.Set( 0, 0, 2 );
    q2.Set( 0, 0, 4 );
    rIntersect = VectorOps::SegmentIntersectsSegment( p1, p2, q1, q2, x );
    Assert( (rIntersect == VectorOps::intersect),
            "SegmentIntersectsSegment didn't intersect" );
    {
      stringstream ssError;
      ssError << "SegmentIntersectsSegment returned incorrect intersect: "
      << x << endl;
      Assert( (x[0] == 0 && x[1] == 0 && x[2] == 2), ssError.str() );
    }

    p1.Set( 0, 0, 0 );
    p2.Set( 0, 0, 2 );
    q1.Set( 0, 0, 2 );
    q2.Set( 0, 0, 2 );
    rIntersect = VectorOps::SegmentIntersectsSegment( p1, p2, q1, q2, x );
    Assert( (rIntersect == VectorOps::intersect),
            "SegmentIntersectsSegment didn't intersect" );
    {
      stringstream ssError;
      ssError << "SegmentIntersectsSegment returned incorrect intersect: "
      << x << endl;
      Assert( (x[0] == 0 && x[1] == 0 && x[2] == 2), ssError.str() );
    }

    p1.Set( 0, 0, 2 );
    p2.Set( 0, 0, 2 );
    q1.Set( 0, 0, 2 );
    q2.Set( 0, 0, 2 );
    rIntersect = VectorOps::SegmentIntersectsSegment( p1, p2, q1, q2, x );
    Assert( (rIntersect == VectorOps::intersect),
            "SegmentIntersectsSegment didn't intersect" );
    {
      stringstream ssError;
      ssError << "SegmentIntersectsSegment returned incorrect intersect: "
      << x << endl;
      Assert( (x[0] == 0 && x[1] == 0 && x[2] == 2), ssError.str() );
    }


  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    VectorOpsTester tester0;
    for ( int i = 0; i < 1; i++ ) {
      tester0.Test();
    }

  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

