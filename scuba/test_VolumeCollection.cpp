/**
 * @file  test_VolumeCollection.cpp
 * @brief test VolumeCollection class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.23 $
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
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <fstream>
#include "VolumeCollection.h"
#include "DataManager.h"
extern "C" {
#include "mri.h"
#include "macros.h" // FEQUAL
}
#include "Scuba-impl.h"

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


#define AssertTclOK(x) \
    if( TCL_OK != (x) ) { \
      ssError << "Tcl_Eval returned not TCL_OK: " << endl  \
      << "Command: " << sCommand << endl \
      << "Result: " << iInterp->result; \
      throw runtime_error( ssError.str() ); \
    } \


using namespace std;

const char* Progname = "test_VolumeCollection";



class VolumeCollectionTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
VolumeCollectionTester::Test ( Tcl_Interp* iInterp ) {

  stringstream ssError;

  try {

    string fnMRI = "test_data/bertT1.mgz";
    VolumeCollection* vol = new VolumeCollection();
    vol->SetFileName( fnMRI );
    MRI* mri = const_cast<MRI*>(vol->GetMRI());

    Assert( (vol->GetTypeDescription() == "Volume"),
            "GetTypeDescription didn't return Volume" );

    DataManager dataMgr = DataManager::GetManager();
    MRILoader mriLoader = dataMgr.GetMRILoader();
    Assert( 1 == mriLoader.CountLoaded(),
            "CountLoaded didn't return 1" );
    Assert( 1 == mriLoader.CountReferences(mri),
            "CountReferences didn't return 1" );

    char* fnMRIC = strdup( fnMRI.c_str() );
    MRI* mriComp = MRIread( fnMRIC );

    Assert( (MRImatch( mriComp, mri )), "MRImatch failed for load" );

    MRIfree( &mriComp );


    // Save it in /tmp, load it, and match it again.
    string fnSave( "/tmp/test.mgz" );
    vol->Save( fnSave );

    VolumeCollection* testVol = new VolumeCollection();
    testVol->SetFileName( fnSave );
    MRI* testMri = const_cast<MRI*>(testVol->GetMRI());
    Assert( (MRImatch( testMri, mri )), "MRImatch failed for load after save");



    // Make an ROI and make sure it's a volume ROI.
    try {
      int roiID = vol->NewROI();
      ScubaROIVolume* roi =
        dynamic_cast<ScubaROIVolume*>(&ScubaROI::FindByID( roiID ));
      roi = NULL;
    } catch (...) {
      throw( runtime_error("typecast failed for NewROI") );
    }


    // Try our conversions.
    Point3<float> world;
    Point3<float> data;
    Point3<int> index;
    world.Set( -50, 0, -80 );
    vol->RASToMRIIndex( world.xyz(), index.xyz() );
    {
      stringstream ssError;
      ssError << "RASToMRIIndex failed. world "
      << world << " index " << index;
      Assert( (index.x() == 178 && index.y() == 208 && index.z() == 128),
              ssError.str() );
    }

    // Set a transform that scales the volume up by 2x in the world.
    ScubaTransform dataTransform;
    dataTransform.SetMainTransform( 2, 0, 0, 0,
                                    0, 2, 0, 0,
                                    0, 0, 2, 0,
                                    0, 0, 0, 1 );
    vol->SetDataToWorldTransform( dataTransform.GetID() );

    world.Set( -50, 0, -80 );
    vol->RASToDataRAS( world.xyz(), data.xyz() );
    {
      stringstream ssError;
      ssError << "RASToDataRAS failed. world "
      << world << " data " << data;
      Assert( ((FEQUAL(data.x(),-25)) &&
               (FEQUAL(data.y(),0)) &&
               (FEQUAL(data.z(),-40))),
              ssError.str() );
    }

    vol->RASToMRIIndex( world.xyz(), index.xyz() );

    if ( index.x() != 153 || index.y() != 168 || index.z() != 128 ) {
      cerr << "RASToMRIIndex with data transform failed. world "
      << world << " index " << index << endl;
      throw( runtime_error( "failed" ) );
    }


    world.Set( -50, 0, -80 );
    VolumeLocation loc( vol->MakeVolumeLocationFromRAS( world.xyz() ) );
    if ( !vol->IsInBounds( loc ) ) {
      stringstream ssError;
      ssError << "IsInBounds failed. world " << world;
      throw( runtime_error( ssError.str() ) );
    }
    world.Set( -1000, 0, 0 );
    VolumeLocation loc2( vol->MakeVolumeLocationFromRAS( world.xyz() ) );
    if ( vol->IsInBounds( loc2 ) ) {
      stringstream ssError;
      ssError << "IsInBounds failed. world " << world;
      throw( runtime_error( ssError.str() ) );
    }


    dataTransform.SetMainTransform( 2, 0, 0, 0,
                                    0, 2, 0, 0,
                                    0, 0, 2, 0,
                                    0, 0, 0, 1 );
    vol->SetDataToWorldTransform( dataTransform.GetID() );
    world.Set( 0, -1000, -254 );
    VolumeLocation loc3( vol->MakeVolumeLocationFromRAS( world.xyz() ) );
    if ( vol->IsInBounds( loc3 ) ) {
      stringstream ssError;
      vol->RASToMRIIndex( world.xyz(), index.xyz() );
      ssError << "IsRASInMRIBounds failed. world " << world
      << " index " << index;
      throw( runtime_error( ssError.str() ) );
    }


    {
      // Create a new one from template.
      VolumeCollection* vol2 = new VolumeCollection();
      vol2->MakeUsingTemplate( vol->GetID(), -1 );
      Assert( (vol->mVoxelSize[0] == vol2->mVoxelSize[0] &&
	       vol->mVoxelSize[1] == vol2->mVoxelSize[1] &&
	       vol->mVoxelSize[2] == vol2->mVoxelSize[2]),
	      "NewUsingTemplate failed, vol2 didn't match vol's voxelsize" );
      Assert( (vol->GetDataType() == vol2->GetDataType()),
	    "NewUsingTemplate(-1) failed, vol2 didn't match vol's data type" );
      delete vol2;
    }
    
    {
      VolumeCollection* vol2 = new VolumeCollection();
      vol2->MakeUsingTemplate( vol->GetID(), MRI_FLOAT );
      Assert( (vol->mVoxelSize[0] == vol2->mVoxelSize[0] &&
	       vol->mVoxelSize[1] == vol2->mVoxelSize[1] &&
	       vol->mVoxelSize[2] == vol2->mVoxelSize[2]),
	      "NewUsingTemplate failed, vol2 didn't match vol's voxelsize" );
      Assert( (MRI_FLOAT == vol2->GetDataType()),
	      "NewUsingTemplate(float) failed, vol2 wasn't a float" );
      
      int idx[3] = { 0, 0, 0 };
      VolumeLocation loc( vol2->MakeVolumeLocationFromIndex( idx ) );
      vol2->SetMRIValue( loc, 0.6789 );
      float value = vol2->GetMRINearestValue(loc);

      stringstream ssError;
      ssError << "NewUsingTemplate(float) failed value comparison, "
	      << "was expecting 0.6789 but got " << value;
      Assert( (fabs(0.6789 - value) < 0.00001), ssError.str() );
      delete vol2;
    }

    {
      VolumeCollection* vol2 = NULL;
      try {
	vol2 = new VolumeCollection();
	vol2->MakeUsingTemplate( vol->GetID(), 10000 );
	throw runtime_error( "MakeUsingTemplate(10000) didn't throw an error");
      }
      catch(...) {}
      delete vol2;
    }

    dataTransform.SetMainTransform( 1, 0, 0, 0,
                                    0, 1, 0, 0,
                                    0, 0, 1, 0,
                                    0, 0, 0, 1 );
    vol->SetDataToWorldTransform( dataTransform.GetID() );

    // FindRASPointsInSquare
    {
      Point3<float> sqRAS[4], cRAS;
      sqRAS[0].Set( 0, -3, 71 );
      sqRAS[1].Set( 0, 1, 71 );
      sqRAS[2].Set( 0, 1, 67 );
      sqRAS[3].Set( 0, -3, 67 );
      cRAS.Set( 0, -1, 69 );
      list<Point3<float> > points;
      vol->FindRASPointsInSquare( cRAS.xyz(),
                                  sqRAS[0].xyz(), sqRAS[1].xyz(),
                                  sqRAS[2].xyz(), sqRAS[3].xyz(),
                                  0, points );
      list<Point3<float> >::iterator tPoint;
      for ( tPoint = points.begin(); tPoint != points.end(); ++tPoint ) {
        Point3<float> pRAS = *tPoint;

        stringstream ssMsg;
        ssMsg << "Failed: " << pRAS << " outside of square";

        Assert( ( pRAS[0] == 0 &&
                  (pRAS[1] >= -3 && pRAS[1] <= 1) &&
                  (pRAS[2] >= 67 && pRAS[2] <= 71) ),
                ssMsg.str() );
      }
    }


    // VoxelIntersectsSegment
    {
      Point3<int> idx;
      Point3<float> segIdxA, segIdxB;
      Point3<float> intIdx;

      idx.Set( 5, 5, 5 );

      float aSegments[6][3] = { {3, 5.5, 5.5}, {6, 5.5, 5.5},
                                {5.5, 3, 5.5}, {5.5, 6, 5.5},
                                {5.5, 5.5, 3}, {5.5, 5.5, 6} };
      for ( int nSegment = 0; nSegment < 6; nSegment += 2 ) {

        segIdxA.Set( aSegments[nSegment] );
        segIdxB.Set( aSegments[nSegment+1] );

        VectorOps::IntersectionResult rInt =
          vol->VoxelIntersectsSegment( idx, segIdxA, segIdxB, intIdx );

        if ( VectorOps::intersect != rInt ) {
          cerr << "Failed VoxelIntersectsSegment test: idx " << idx
          << ", seg " << segIdxA << ", " << segIdxB << endl
          << "\tDidn't intersect" << endl;
          throw runtime_error("failed");
        }
      }

      segIdxA.Set( 0, 5.5, 5.5 );
      segIdxB.Set( 4, 5.5, 5.5 );
      VectorOps::IntersectionResult rInt =
        vol->VoxelIntersectsSegment( idx, segIdxA, segIdxB, intIdx );

      if ( VectorOps::dontIntersect != rInt ) {
        cerr << "Failed VoxelIntersectsSegment test: idx " << idx
        << ", seg " << segIdxA << ", " << segIdxB << endl
        << "\tIntersected" << endl;
        throw runtime_error("failed");
      }
    }

    // FindRASPointsOnSegment
    {

    }


    // GetVoxelsWithValue
    {
      // This is a 5cubed volume whose values are set to the x
      // coordinate. So for x=3,y=0..4,z=0..4, value = 3.
      string fnVol = "test_data/testVolumeCollection-GetVoxelsWithValue.mgh";
      ifstream fVol( fnMRI.c_str(), ios::in );
      if ( !fVol ) {
        throw runtime_error("Couldn't find necessary test data.");
      }
      fVol.close();

      VolumeCollection* vol = new VolumeCollection();
      vol->SetFileName( fnVol );
      vol->LoadVolume();

      // Get the values 0-4 and make sure we got the right voxels.
      for ( int nStructure = 0; nStructure < 5; nStructure++ ) {

        list<VolumeLocation> lLocations;
        vol->GetVoxelsWithValue( nStructure, lLocations );

        Volume3<bool> bGot( 5, 5, 5, false );
        list<VolumeLocation>::iterator tLocation;
        for ( tLocation = lLocations.begin(); tLocation != lLocations.end();
              ++tLocation ) {

          VolumeLocation loc = *(tLocation);
          bGot.Set( loc.Index()[0], loc.Index()[1], loc.Index()[2], true );
        }

        for ( int nZ = 0; nZ < 5; nZ++ ) {
          for ( int nY = 0; nY < 5; nY++ ) {
            for ( int nX = 0; nX < 5; nX++ ) {
              if ( nX == nStructure && !bGot.Get( nX, nY, nZ ) ) {
                stringstream ssErr;
                ssErr << "Failed GetVoxelsWithValue test: "
                << " nStructure = " << nStructure
                << " index " << Point3<int>(nX,nY,nZ)
                << " - was supposed to get voxel but didn't";
                throw runtime_error(ssErr.str());
              }
              if ( nX != nStructure && bGot.Get( nX, nY, nZ ) ) {
                stringstream ssErr;
                ssErr << "Failed GetVoxelsWithValue test: "
                << " nStructure = " << nStructure
                << " index " << Point3<int>(nX,nY,nZ)
                << " - wasn't supposed to get voxel but did";
                throw runtime_error(ssErr.str());
              }
            }
          }
        }
      }

      delete vol;
    }

    // GetAverageValue
    {
      // We'll use the same volume as with GetVoxelsWithValue.
      string fnVol = "test_data/testVolumeCollection-GetVoxelsWithValue.mgh";
      ifstream fVol( fnMRI.c_str(), ios::in );
      if ( !fVol ) {
        throw runtime_error("Couldn't find necessary test data.");
      }
      fVol.close();

      VolumeCollection* vol = new VolumeCollection();
      vol->SetFileName( fnVol );
      vol->LoadVolume();

      // Get values of all voxels in plane x=3 and make sure it's 3.
      list<VolumeLocation> lVoxels;
      for ( int nZ = 0; nZ < 5; nZ++ ) {
        for ( int nY = 0; nY < 5; nY++ ) {
          int index[3] = { 3, nY, nZ };
          VolumeLocation loc( vol->MakeVolumeLocationFromIndex( index ) );
          lVoxels.push_back( loc );
        }
      }
      float average = vol->GetAverageValue( lVoxels );
      if ( average != 3.0 ) {
        stringstream ssErr;
        ssErr << "Failed GetAverageValue: Getting all voxels in x=3, "
        << "average should have been 3, but was " << average;
        throw runtime_error( ssErr.str() );
      }

      // Get values of 5 voxels, one from each x plane
      // (val=0,1,2,3,4), and make sure it's 2.
      lVoxels.clear();
      for ( int nX = 0; nX < 5; nX++ ) {
        int index[3] = { nX, 3, 3 };
        VolumeLocation loc( vol->MakeVolumeLocationFromIndex( index ) );
        lVoxels.push_back( loc );
      }
      average = vol->GetAverageValue( lVoxels );
      if ( average != 2.0 ) {
        stringstream ssErr;
        ssErr << "Failed GetAverageValue: Getting voxels in different planes, "
        << "average should have been 2, but was " << average;
        throw runtime_error( ssErr.str() );
      }

      // Make sure we get an error for no voxels.
      lVoxels.clear();
      bool bDidntThrow = false;
      try {
        average = vol->GetAverageValue( lVoxels );
        bDidntThrow = true;
      } catch ( exception& e ) {}
      if ( bDidntThrow ) {
        stringstream ssErr;
        ssErr << "Failed GetAverageValue: Didn't throw with empty list";
        throw runtime_error( ssErr.str() );
      }

      delete vol;
    }


    // Check the tcl commands.
    char sCommand[1024];
    int rTcl;

    int id = vol->GetID();
    string fnTest = "test-name";
    sprintf( sCommand, "SetVolumeCollectionFileName %d test-name", id );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    Assert( (vol->mfnMRI == fnTest),
            "Setting file name via tcl didn't work" );

    vol->SetDataToWorldTransform( 0 );

    delete vol;
    delete testVol;
  } catch ( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed." << endl;
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


    VolumeCollectionTester tester0;
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
