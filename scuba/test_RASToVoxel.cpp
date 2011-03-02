/**
 * @file  test_RASToVoxel.cpp
 * @brief test RASToVoxel class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
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


#include <stdlib.h>
#include <iomanip>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "Transform44.h"
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

using namespace std;

const char* Progname = "test_RASToVoxel";

#define FLOAT_EQUAL(x,y) (fabs((x)-(y))<0.0001)

#define USEFLOORTOROUND 1

void TestCoord ( Transform44& rasToVoxel, Point3<float>& ras, bool ibPrint ) {


  Point3<int> idxi;
  Point3<float> idxf;
  Point3<float> rasFromIdxICheck;
  Point3<float> rasFromIdxFCheck;
  Point3<int> idxiFromRASFromIdxICheck;
  Point3<float> idxfFromRASFromIdxICheck;

  rasToVoxel.MultiplyVector3( ras.xyz(), idxi.xyz() );
  rasToVoxel.MultiplyVector3( ras.xyz(), idxf.xyz() );

  if ( ibPrint )
    cerr << setprecision(20) << "idxi " << idxi.x() << ", "
    << idxi.y() << ", " << idxi.z() << " "
    << "idxf " << idxf.x() << ", "
    << idxf.y() << ", " << idxf.z() << endl;

#if USEFLOORTOROUND
  if ( idxi.x() != (int) floor( idxf.x() + 0.5 ) ||
       idxi.y() != (int) floor( idxf.y() + 0.5 ) ||
       idxi.z() != (int) floor( idxf.z() + 0.5 ) ) {
#else
  if ( idxi.x() != (int) floor( idxf.x() ) ||
       idxi.y() != (int) floor( idxf.y() ) ||
       idxi.z() != (int) floor( idxf.z() ) ) {
#endif

    stringstream ssError;
    ssError << "RAS -> index float didn't match int: "
    << ras << " -> " << idxf << ", " << idxi;
    throw runtime_error( ssError.str() );
  }

  rasToVoxel.InvMultiplyVector3( idxf.xyz(), rasFromIdxFCheck.xyz() );

  if ( !FLOAT_EQUAL(rasFromIdxFCheck.x(),ras.x()) ||
       !FLOAT_EQUAL(rasFromIdxFCheck.y(),ras.y()) ||
       !FLOAT_EQUAL(rasFromIdxFCheck.z(),ras.z()) ) {

    stringstream ssError;
    ssError << "RAS -> index float -> RAS didn't match orig RAS: "
    << ras << " -> " << idxf << " -> " << rasFromIdxFCheck;
    throw runtime_error( ssError.str() );
  }

  rasToVoxel.InvMultiplyVector3( idxi.xyz(), rasFromIdxICheck.xyz() );

  rasToVoxel.MultiplyVector3( rasFromIdxICheck.xyz(),
                              idxiFromRASFromIdxICheck.xyz() );
  rasToVoxel.MultiplyVector3( rasFromIdxICheck.xyz(),
                              idxfFromRASFromIdxICheck.xyz() );

#if USEFLOORTOROUND
  if ( idxiFromRASFromIdxICheck.x() !=
       (int) floor( idxfFromRASFromIdxICheck.x() + 0.5 ) ||
       idxiFromRASFromIdxICheck.y() !=
       (int) floor( idxfFromRASFromIdxICheck.y() + 0.5 ) ||
       idxiFromRASFromIdxICheck.z() !=
       (int) floor( idxfFromRASFromIdxICheck.z() + 0.5 ) ) {
#else
  if ( idxiFromRASFromIdxICheck.x() !=
       (int) floor( idxfFromRASFromIdxICheck.x() ) ||
       idxiFromRASFromIdxICheck.y() !=
       (int) floor( idxfFromRASFromIdxICheck.y() ) ||
       idxiFromRASFromIdxICheck.z() !=
       (int) floor( idxfFromRASFromIdxICheck.z() ) ) {
#endif

    stringstream ssError;
    ssError << "RAS -> index int -> RAS -> index float "
    << "didn't match int: "
    << ras << " -> " << idxf << " -> " << rasFromIdxICheck
    << " -> " << idxiFromRASFromIdxICheck << ", "
    << idxfFromRASFromIdxICheck;
    throw runtime_error( ssError.str() );
  }

  if ( idxiFromRASFromIdxICheck.x() != idxi.x() ||
       idxiFromRASFromIdxICheck.y() != idxi.y() ||
       idxiFromRASFromIdxICheck.z() != idxi.z() ) {

    Point3<float> idxf2, ras2, idxf3;
    idxf2[0] = idxi[0];
    idxf2[1] = idxi[1];
    idxf2[2] = idxi[2];
    rasToVoxel.InvMultiplyVector3( idxf2.xyz(), ras2.xyz() );
    rasToVoxel.MultiplyVector3( ras2.xyz(), idxf3.xyz() );
    stringstream ssError;
    ssError << "index int -> RAS -> index, "
    << "index ints didn't match:\n"
    << setprecision(20)
    << "\t" << idxi << " -> " << rasFromIdxICheck
    << " -> " << idxiFromRASFromIdxICheck << "\n"
    << "\t" << idxf2 << " -> " << ras2 << " -> " << idxf3 << endl;
    throw runtime_error( ssError.str() );
  }

  if ( ibPrint )
    cerr << "RAS -> index int, float -> RAS -> index int, index float " << endl
    << ras << " -> " << idxi << ", " << idxf
    << " -> " << rasFromIdxICheck
    << " -> " << idxiFromRASFromIdxICheck << ", "
    << idxfFromRASFromIdxICheck << endl;
}

int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

#if 0
    Transform44 rasToVoxel;
    rasToVoxel.SetMainTransform( -9.78514, 0.0825722, 1.09133, 189.96,
                                 -3.04974e-08, -9.81809, 0.742855, 286.028,
                                 1.11155, 0.749787, 9.90971, 113.646,
                                 0,      0,      0,      1 );
#endif

    string fnMRI = "test_data/bertT1.mgz";
    if ( argc > 1 ) {
      fnMRI = argv[1];
    }

    VolumeCollection* vol = new VolumeCollection;
    try {
      vol->SetFileName( fnMRI );
      vol->LoadVolume();
    } catch ( runtime_error& e ) {
      cerr << e.what() << endl;
      exit( 1 );
    } catch ( ... ) {
      cerr << "Couldn't load volume" << endl;
      exit( 1 );
    }

    Transform44 rasToVoxel;
    rasToVoxel.SetMainTransform ( vol->GetWorldToIndexTransform() );
    delete vol;

    Point3<float> ras;
    if ( argc == 5 ) {
      //    ras.Set (0.364583,2.92578,18.4121);
      ras.Set(0.0595181,-52.65,12.19);
      ras.Set( atof(argv[2]), atof(argv[3]), atof(argv[4]) );
      TestCoord( rasToVoxel, ras, true );

    } else {

      for ( float z = -10.0; z < 10.0; z += 0.0490005 ) {
        for ( float y = -10.0; y < 10.0; y += 0.0490003 ) {
          for ( float x = -10.0; x < 10.0; x += 0.0490001 ) {
            ras.Set( x, y, z );
            TestCoord( rasToVoxel, ras, false );
          }
        }
      }

      for ( int nTrial = 0; nTrial < 1000; nTrial++ ) {

        do {
          ras.Set( (float)random()/(float)random(),
                   (float)random()/(float)random(),
                   (float)random()/(float)random() );
        } while ( ras.x() < -50.0 || ras.x() > 50.0 ||
                  ras.y() < -50.0 || ras.y() > 50.0 ||
                  ras.z() < -50.0 || ras.z() > 50.0 );

        TestCoord( rasToVoxel, ras, false );
      }
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
