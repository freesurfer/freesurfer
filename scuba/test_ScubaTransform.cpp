/**
 * @file  test_ScubaTransform.cpp
 * @brief test ScubaTransform class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.13 $
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


#include <fstream>
#include "ScubaTransform.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_ScubaTransform";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      stringstream ssError; \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
      << "Command: " << sCommand << endl \
      << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


#define VFEQUAL(v,a,b,c) \
   (FEQUAL(((v)[0]),a) && FEQUAL(((v)[1]),b) && FEQUAL(((v)[2]),c))

class ScubaTransformTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
ScubaTransformTester::Test ( Tcl_Interp* iInterp ) {

  try {

    ScubaTransform transform;

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
      stringstream ssError;
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


    // Try the test case from VolumeCollection.
    ScubaTransform t1;
    ScubaTransform* mDataToWorldTransform =
      &ScubaTransform::FindByID( t1.GetID() );
    ScubaTransform mDataToIndexTransform;
    Transform44 mWorldToIndexTransform;

    mDataToIndexTransform.SetMainTransform( -1, 0,  0, 128,
                                            0, 0, -1, 128,
                                            0, 1,  0, 128,
                                            0, 0,  0,   1 );
    mDataToWorldTransform->SetMainTransform( 2, 0,  0,   0,
        0, 2,  0,   0,
        0, 0,  2,   0,
        0, 0,  0,   1 );

    Transform44 worldToData = mDataToWorldTransform->Inverse();
    mWorldToIndexTransform = worldToData;
    mWorldToIndexTransform.ApplyTransform( mDataToIndexTransform );

    Transform44& t = mWorldToIndexTransform;
    if ( !(t(0,0) == -0.5 && t(1,0) == 0 && t(2,0) == 0 && t(3,0) == 128 &&
           t(0,1) == 0 && t(1,1) == 0 && t(2,1) == -0.5 && t(3,1) == 128 &&
           t(0,2) == 0 && t(1,2) == 0.5 && t(2,2) == 0 && t(3,2) == 128 &&
           t(0,3) == 0 && t(1,3) == 0 && t(2,3) == 0 && t(3,3) == 1) ) {

      cerr << t << endl;
      Assert( 0, "ApplyTransform case didn't work." );
    }



    // Try the tcl commands.
    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "SetTransformValues %d {0 1 2 3 4 5 6 7 8 9 "
             "10 11 12 13 14 15}", transform.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        float value = transform.m( c, r );
        stringstream ssError;
        ssError << "TCL set check failed for " << c << ", " << r;
        Assert((value - (r*4) == c &&
                (value - c) / 4 == r),
               ssError.str() );
      }
    }

    sprintf( sCommand, "GetTransformValues %d", transform.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    const char* sTclResult = Tcl_GetStringResult( iInterp );
    stringstream ssResult( sTclResult );
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++ ) {
        float value;
        ssResult >> value;
        stringstream ssError;
        ssError << "TCL set check failed for " << c << ", " << r;
        Assert((value - (r*4) == c &&
                (value - c) / 4 == r),
               ssError.str() );
      }
    }

    string fnLTA = "test_data/testScubaTransform.lta";
    ifstream fLTA( fnLTA.c_str(), ios::in );
    if ( !fLTA ) {
      cerr << "WARNING: File " + fnLTA + " not found, test skipped." << endl;
    } else {

      fLTA.close();
      ScubaTransform l;
      sprintf( sCommand, "LoadTransformFromLTAFile %d %s",
               l.GetID(), fnLTA.c_str() );
      rTcl = Tcl_Eval( iInterp, sCommand );
      AssertTclOK( rTcl );
      Assert((l(0,0) == 1 && l(1,0) == 2 && l(2,0) == 3 && l(3,0) == 4 &&
              l(0,1) == 5 && l(1,1) == 6 && l(2,1) == 7 && l(3,1) == 8 &&
              l(0,2) == 9 && l(1,2) == 10 && l(2,2) == 11 && l(3,2) == 12 &&
              l(0,3) == 13 && l(1,3) == 14 && l(2,3) == 15 && l(3,3) == 16),
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

    Tcl_Interp* interp = Tcl_CreateInterp();
    Assert( interp, "Tcl_CreateInterp returned null" );

    int rTcl = Tcl_Init( interp );
    Assert( TCL_OK == rTcl, "Tcl_Init returned not TCL_OK" );

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.SetOutputStreamToCerr();
    commandMgr.Start( interp );


    ScubaTransformTester tester0;
    tester0.Test( interp );


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

