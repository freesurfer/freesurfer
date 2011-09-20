/**
 * @file  test_Transform44.cpp
 * @brief test 4x4 Transform routines
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/09/20 19:06:38 $
 *    $Revision: 1.11 $
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
#include <fstream>
#include "Transform44.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_Transform44";

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

class Transform44Tester {
public:
  void Test();
};

void
Transform44Tester::Test () {

  stringstream ssError;

  try {

    Transform44 transform;

    // Set to identity, make sure multing a vector returns same vector.
    float in[3], out[3];
    in[0] = 5;
    in[1] = 35.67;
    in[2] = 1000;
    transform.MultiplyVector3( in, out );
    Assert((in[0] == out[0] && in[1] == out[1] && in[2] == out[2]),
           "Identity mult check failed");


    // Set to a scale matrix, make sure multing a vector returns
    // correct response.
    transform.SetMainTransform( 5, 0, 0, 0,
                                0, 5, 0, 0,
                                0, 0, 5, 0,
                                0, 0, 0, 1 );

    in[0] = 5;
    in[1] = 6;
    in[2] = 7;
    transform.MultiplyVector3( in, out );
    if ( !(in[0]*5 == out[0] && in[1]*5 == out[1] && in[2]*5 == out[2]) ) {
      ssError << "Scale mult check failed" << endl
      << transform << endl
      << "out " << Point3<float>(out) << endl;
      throw(runtime_error(ssError.str()));
    }

    // Check the inverse.
    transform.InvMultiplyVector3( out, in );
    Assert((FEQUAL( in[0], 5.0 ) &&
            FEQUAL( in[1], 6.0 ) && FEQUAL( in[2], 7.0) ),
           "Inv scale mult check failed");

    // Try loading an LTA.
    string fnLTA = "test_data/testTransform44.lta";
    ifstream fLTA( fnLTA.c_str(), ios::in );
    if ( !fLTA ) {
      cerr << "WARNING: File " + fnLTA + " not found, test skipped." << endl;
    } else {

      Transform44 l;
      l.LoadFromLTAFile( fnLTA );
      float l00,l10,l20,l30,l01,l11,l21,l31,l02,l12,l22,l32,l03,l13,l23,l33;
      l00=l(0,0) -  0.999194562;
      l10=l(1,0) -  0.0126949335;
      l20=l(2,0) -  0.0380673669;
      l30=l(3,0) - -7.12072897;
      l01=l(0,1) - -0.0120223239;
      l11=l(1,1) -  0.999768436;
      l21=l(2,1) - -0.0178460851;
      l31=l(3,1) -  2.41395235;
      l02=l(0,2) - -0.0382851064;
      l12=l(1,2) -  0.0173740536;
      l22=l(2,2) -  0.99911582;
      l32=l(3,2) -  1.43621171;
      l03=l(0,3) -  0;
      l13=l(1,3) -  0;
      l23=l(2,3) -  0;
      l33=l(3,3) -  1;
      float tolerance = 0.00000001;
      Assert((l00  < tolerance && 
              l10  < tolerance && 
              l20  < tolerance && 
              l30  < tolerance &&
              l01  < tolerance && 
              l11  < tolerance && 
              l21  < tolerance && 
              l31  < tolerance &&
              l02  < tolerance && 
              l12  < tolerance && 
              l22  < tolerance && 
              l32  < tolerance &&
              l03  < tolerance && 
              l13  < tolerance && 
              l23  < tolerance && 
              l33  < tolerance),
             "LTA didn't load properly.");
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

    Transform44Tester tester0;
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

