/**
 * @file  test_Matrix44.cpp
 * @brief test Matrix 4x4 routines
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
#include "Matrix44.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_Matrix44";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


#define VFEQUAL(v,a,b,c) \
   (FEQUAL(((v)[0]),a) && FEQUAL(((v)[1]),b) && FEQUAL(((v)[2]),c))

#define FEQUAL2(a,b) ((a) - (b) < 0.001)

class Matrix44Tester {
public:
  void Test();
};

void
Matrix44Tester::Test () {

  stringstream ssError;

  try {

    Matrix44 transform;

    // Make sure starts out as identity.
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        if ( r == c ) {
          Assert((transform.GetCR(c,r) == 1),
                 "Didn't start out as identity");
        } else {
          Assert((transform.GetCR(c,r) == 0),
                 "Didn't start out as identity");
        }
      }
    }

    // Set a bunch of values, make sure they got set.
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        transform.SetCR( c, r, (r * 4) + c);
      }
    }
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        float value = transform.GetCR( c, r );
        ssError << "Check failed for " << c << ", " << r;
        Assert((value - (r*4) == c &&
                (value - c) / 4 == r),
               ssError.str() );
      }
    }

    // Set to identity, make sure multing a vector returns same vector.
    transform.MakeIdentity();
    float in[3], out[3];
    in[0] = 5;
    in[1] = 35.67;
    in[2] = 1000;
    transform.MultiplyVector3( in, out );
    Assert((in[0] == out[0] && in[1] == out[1] && in[2] == out[2]),
           "Identity mult check failed");


    // Set to a scale matrix, make sure multing a vector returns
    // correct response.
    transform.SetMatrix( 5, 0, 0, 0,
                         0, 5, 0, 0,
                         0, 0, 5, 0,
                         0, 0, 0, 1 );

    in[0] = 5;
    in[1] = 6;
    in[2] = 7;
    transform.MultiplyVector3( in, out );
    Assert((in[0]*5 == out[0] && in[1]*5 == out[1] && in[2]*5 == out[2]),
           "Scale mult check failed");



    // Test rotation by making rotation matrices and transforming
    // points.
    Point3<float> p( 1, 0, 0 );
    Point3<float> q;

    transform.MakeXRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,1,0,0), "X Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,-1,0), "X Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,-1), "X Rotation failed");

    transform.MakeYRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,-1,0,0), "Y Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,1,0), "Y Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,-1), "Y Rotation failed");

    transform.MakeZRotation( M_PI );
    p.Set( 1, 0, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,-1,0,0), "Z Rotation failed");
    p.Set( 0, 1, 0 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,-1,0), "Z Rotation failed");
    p.Set( 0, 0, 1 );
    transform.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,1), "Z Rotation failed");

    // Test an identity transform.
    Matrix44 m;
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        m.SetCR( c, r, (r * 4) + c);
      }
    }
    Matrix44 id;
    id.MakeIdentity();
    Matrix44 n = m * id;
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        Assert( (FEQUAL(n(c,r),m(c,r))), "Operator* failed" );
      }
    }

    // Test rotation around a vector.
    p.Set( 0, 1, 0 );
    Point3<float> v( 1, 0, 0 );
    m.MakeRotation( p.xyz(), v.xyz(), M_PI );
    p.Set( 1, 0, 0 );
    m.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,1,2,0), "Rotation failed p " << p << " q " << q
           << endl << "m " << m);

    p.Set( 0, 0, 1 );
    v.Set( 0, 1, 0 );
    m.MakeRotation( p.xyz(), v.xyz(), M_PI );
    p.Set( 0, 0, 0 );
    m.MultiplyVector3( p.xyz(), q.xyz() );
    Assert(VFEQUAL(q,0,0,2), "Rotation failed");


    // Assignment operator.
    Matrix44 a1;
    a1.SetMatrix(  1 , 2,  3,  4,
                   5,  6,  7,  8,
                   9, 10, 11, 12,
                   13, 14, 15, 16 );
    Matrix44 a2;
    a2 = a1;
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        Assert( (a2(c,r) == a1(c,r)), "Operator= failed" );
      }
    }

    // Multiply operator.
    a1.SetMatrix( 0, 1, 0, 1,
                  1, 0, 1, 0,
                  0, 1, 0, 1,
                  1, 0, 1, 0 );
    Matrix44 a3 = a2 * a1;
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        float sum = 0;
        stringstream ssSum;
        for ( int i = 0; i < 4; i++ ) {
          ssSum << "(a1("<<c<<","<<i<<") "<<a1(c,i)
          << " * a2("<<i<<","<<c<<") "<<a2(i,c)
          << " = "<<a1(c,i) * a2(i,c) << ") + ";
          sum += a1(c,i) * a2(i,r);
        }
        ssSum << " = " << sum << endl;
        if ( !(FEQUAL(a3(c,r),sum)) ) {

          stringstream ss;
          ss << "Line " << __LINE__ << ": Operator* failed"<< endl;
          ss << a1;
          ss << a2;
          ss << a3;
          ss << "mult failed on a3(" << c << "," << r << "), sum should be "
          << sum << endl;
          ss << ssSum.str();
          throw runtime_error( ss.str() );
        }
      }
    }


    // Extract translate.
    m.SetMatrix(  1,  2,  3,  4,
                  5,  6,  7,  8,
                  9, 10, 11, 12,
                  13, 14, 15, 16 );
    Matrix44 t = m.ExtractTranslation();
    if ( !(t(0,0) == 1  &&  t(1,0) == 0  &&  t(2,0) == 0  &&  t(3,0) == 4  &&
           t(0,1) == 0  &&  t(1,1) == 1  &&  t(2,1) == 0  &&  t(3,1) == 8  &&
           t(0,2) == 0  &&  t(1,2) == 0  &&  t(2,2) == 1  &&  t(3,2) == 12  &&
           t(0,3) == 0  &&  t(1,3) == 0  &&  t(2,3) == 0  &&  t(3,3) == 1) ) {

      stringstream ss;
      ss << "Line " << __LINE__ << ": ExtractTranslation failed"<< endl;
      ss << m;
      ss << t;
      throw runtime_error( ss.str() );
    }

    // Extract scale.
    m.SetMatrix( 5, 0, 0, 0,
                 0, 5, 0, 0,
                 0, 0, 5, 0,
                 0, 0, 0, 1 );
    Matrix44 s = m.ExtractScale();
    if ( !(s(0,0) == 5  &&  s(1,0) == 0  &&  s(2,0) == 0  &&  s(3,0) == 0  &&
           s(0,1) == 0  &&  s(1,1) == 5  &&  s(2,1) == 0  &&  s(3,1) == 0  &&
           s(0,2) == 0  &&  s(1,2) == 0  &&  s(2,2) == 5  &&  s(3,2) == 0  &&
           s(0,3) == 0  &&  s(1,3) == 0  &&  s(2,3) == 0  &&  s(3,3) == 1) ) {

      stringstream ss;
      ss << "Line " << __LINE__ << ": ExtractScale failed"<< endl;
      ss << m;
      ss << s;
      throw runtime_error( ss.str() );
    }


    // Inverse.
    m.SetMatrix( -0.263203142774843, -0.92207253201611, -0.283736411718422, -3.35437752965114,
                 0.960353372502555, -0.222402361044052, -0.168102914088361, 19.7297470143766,
                 0.0918994317523974, -0.316732435193226, 0.944052466200981, -17.5292922485199,
                 0, 0, 0, 1 );
    Matrix44 inv = m.Inverse();

    Matrix44 i = m * inv;
    if ( !(FEQUAL2((i(0,0)),1)  &&  FEQUAL2((i(1,0)),1)  &&  FEQUAL2((i(2,0)),1)  &&  FEQUAL2((i(3,0)),1)  &&
           FEQUAL2((i(0,1)),1)  &&  FEQUAL2((i(1,1)),1)  &&  FEQUAL2((i(2,1)),1)  &&  FEQUAL2((i(3,1)),1)  &&
           FEQUAL2((i(0,2)),1)  &&  FEQUAL2((i(1,2)),1)  &&  FEQUAL2((i(2,2)),1)  &&  FEQUAL2((i(3,2)),1)  &&
           FEQUAL2((i(0,3)),1)  &&  FEQUAL2((i(1,3)),1)  &&  FEQUAL2((i(2,3)),1)  &&  FEQUAL2((i(3,3)),1))) {

      stringstream ss;
      ss << "Line " << __LINE__ << ": Inverse failed"<< endl;
      ss << inv;
      ss << i;
      throw runtime_error( ss.str() );
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

    Matrix44Tester tester0;
    tester0.Test();
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

