/**
 * @file  test_ScubaTransform.cpp
 * @brief test ScubaTransform class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/09/20 19:06:38 $
 *    $Revision: 1.14 $
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

    const float answer[16] = {0.999194562, 0.0126949335, 0.0380673669, -7.12072897, -0.0120223239, 0.999768436, -0.0178460851, 2.41395235, -0.0382851064, 0.0173740536, 0.99911582, 1.43621171, 0, 0, 0, 1};
    sprintf( sCommand, "SetTransformValues %d {0.999194562 0.0126949335 0.0380673669 -7.12072897 -0.0120223239 0.999768436 -0.0178460851 2.41395235 -0.0382851064 0.0173740536 0.99911582 1.43621171 0 0 0 1}", transform.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    int i=0;
    for ( int r = 0; r < 4; r++ ) {
      for ( int c = 0; c < 4; c++, i++ ) {
        float value = transform.m( c, r );
        stringstream ssError;
        ssError << "TCL set check failed for " << c << ", " << r;
        Assert((value - answer[i]) < 0.00000001 , ssError.str() );
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
//        Assert((value - answer[i]) < 0.00000001 , ssError.str() );
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

