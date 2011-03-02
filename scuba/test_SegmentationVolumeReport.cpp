/**
 * @file  test_SegmentationVolumeReport.cpp
 * @brief test SegmentationVolumeReport class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.8 $
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
#include "SegmentationVolumeReport.h"
extern "C" {
#include "mri.h"
}
#include "Scuba-impl.h"

// Uses test files
//
// test_data/testSEgmentationVolumeReportData-Seg.mgz
//   This is a 25 cubed volume with 1s mm voxels. On slice s=-11, all
//   values are on. On s=-10, all are 2, etc, up to s=13 with all
//   values 25.

// test_data/testSEgmentationVolumeReportData-Int.mgz
//   Same size as the seg volume, except the intensity for each
//   corresponding segmentation area is 2x the segmentation value.

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


using namespace std;

const char* Progname = "test_SegmentationVolumeReport";

class SegmentationVolumeReportTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void
SegmentationVolumeReportTester::Test ( Tcl_Interp* iInterp ) {

  try {


    // Load our seg volume.
    string fnSegVolume = "test_data/testSegmentationVolumeReportData-Seg.mgz";
    VolumeCollection seg;
    seg.SetFileName( fnSegVolume );
    seg.LoadVolume();
    seg.SetLabel( "Seg" );

    // Load our intensity volume.
    string fnIntVolume = "test_data/testSegmentationVolumeReportData-Int.mgz";
    VolumeCollection vol;
    vol.SetFileName( fnIntVolume );
    vol.LoadVolume();
    vol.SetLabel( "Int" );


    // Load our LUT.
    string fnLUT = "test_data/TestLUT.txt";
    ScubaColorLUT lut;
    lut.UseFile( fnLUT );

    // Set up the report.
    SegmentationVolumeReport& report =
      SegmentationVolumeReport::GetReport();

    report.SetSegmentation( seg );
    if ( NULL == report.mSegVol ) {
      stringstream ssError;
      ssError << "Error on SetSegmentation, mSegVol was NULL";
      throw runtime_error( ssError.str() );
    }
    if ( report.mSegVol->GetID() != seg.GetID() ) {
      stringstream ssError;
      ssError << "Error on SetSegmentation, mSegVol was the wrong volume (should be ID " << seg.GetID() << " but was " << report.mSegVol->GetID();
      throw runtime_error( ssError.str() );
    }

    report.DontUseROI();
    if ( report.mbUseROI ) {
      stringstream ssError;
      ssError << "Error on DontUseROI, mbUseROI was true";
      throw runtime_error( ssError.str() );
    }


    report.SetColorLUT( lut );
    if ( NULL == report.mLUT ) {
      stringstream ssError;
      ssError << "Error on SetColorLUT, mROI was NULL";
      throw runtime_error( ssError.str() );
    }
    if ( report.mLUT->GetID() != lut.GetID() ) {
      stringstream ssError;
      ssError << "Error on SetColorLUT, id didn't match";
      throw runtime_error( ssError.str() );
    }

    // Add 1-5 but not 3.
    report.AddSegmentationStructure( 1 );
    report.AddSegmentationStructure( 2 );
    report.AddSegmentationStructure( 4 );
    report.AddSegmentationStructure( 5 );
    map<int,bool> structureMap;
    list<int>::iterator tStructure;
    for ( tStructure = report.mlStructures.begin();
          tStructure != report.mlStructures.end(); ++tStructure ) {
      int nStructure = *tStructure;
      if ( nStructure != 1 && nStructure != 2 &&
           nStructure != 4 && nStructure != 5 ) {
        stringstream ssError;
        ssError << "Error on AddSegmentationStructure, added an unknown structure " << nStructure;
        throw runtime_error( ssError.str() );
      }
      structureMap[nStructure] = true;
    }
    if ( !(structureMap[1] && structureMap[2] &&
           structureMap[4] && structureMap[5]) ) {
      stringstream ssError;
      ssError << "Error in AddSegmentationStructure, didn't add all structures";
      throw runtime_error( ssError.str() );
    }

    // Test handling of undefined structures.
    report.AddSegmentationStructure( 200 );

    // Add the intensity volume. Also add the seg vol as an additional
    // intensity volume.
    report.AddIntensityVolume( vol );
    report.AddIntensityVolume( seg );
    map<int,bool> volsLoadedMap;
    list<VolumeCollection*>::iterator tVolume;
    for ( tVolume = report.mlIntVols.begin();
          tVolume != report.mlIntVols.end(); ++tVolume ) {
      VolumeCollection* testVol = *tVolume;
      int volID = testVol->GetID();
      if ( volID != vol.GetID() && volID != seg.GetID() ) {
        stringstream ssError;
        ssError << "Error in AddIntensityVolume, added a volume with an unknown id " << volID;
        throw runtime_error( ssError.str() );
      }
      volsLoadedMap[volID] = true;
    }
    if ( !(volsLoadedMap[vol.GetID()] && volsLoadedMap[seg.GetID()]) ) {
      stringstream ssError;
      ssError << "Error in AddIntensityVolume, didn't add both volumes";
      throw runtime_error( ssError.str() );
    }

    report.MakeVolumeReport();


    for ( int nStructure = 1; nStructure <= 5; nStructure++ ) {

      // Check the number of voxels we got.
      float expectedSize;
      if ( nStructure == 3 ) {
        expectedSize = 0;
      } else {
        expectedSize = 25 * 25;
      }
      if ( report.mStructureToVolumeMap[nStructure] != expectedSize ) {
        stringstream ssError;
        ssError << "Error on report for seg " << nStructure
        << ": expectedSize was " << expectedSize << ", but got "
        << report.mStructureToVolumeMap[nStructure];
        throw runtime_error( ssError.str() );
      }

      // Check the intensity values (for the vol, should be 2x the
      // segmentation index, and for the seg, should be ==).
      float expectedIntensityAverage;
      if ( nStructure == 3 ) {
        expectedIntensityAverage = 0;
      } else {
        expectedIntensityAverage = nStructure * 2.0;
      }
      if ( report.mVolumeToIntensityAverageMap[&vol][nStructure] !=
           expectedIntensityAverage ) {
        stringstream ssError;
        ssError << "Error on report for seg " << nStructure
        << ": expectedIntensityAverage for intensity vol " << &vol
        << " was " << expectedIntensityAverage << ", but got "
        << report.mVolumeToIntensityAverageMap[&vol][nStructure];
        throw runtime_error( ssError.str() );
      }

      if ( nStructure == 3 ) {
        expectedIntensityAverage = 0;
      } else {
        expectedIntensityAverage = nStructure;
      }
      if ( report.mVolumeToIntensityAverageMap[&seg][nStructure] !=
           expectedIntensityAverage ) {
        stringstream ssError;
        ssError << "Error on report for seg " << nStructure
        << ": expectedIntensityAverage for seg vol " << &seg
        << " was " << expectedIntensityAverage << ", but got "
        << report.mVolumeToIntensityAverageMap[&seg][nStructure];
        throw runtime_error( ssError.str() );
      }
    }

    // Make the report.
    report.MakeVolumeReport( "/tmp/testReport.txt" );

    // Compare it with the correct one.
    int rDiff =
      system( "diff /tmp/testReport.txt testSegmentationVolumeReport-correct.txt" );
    if ( rDiff != 0 ) {
      throw runtime_error( "diff failed on testReport" );
    }

    // Make the report.
    report.MakeIntensityReport( "/tmp/testIntReport.txt" );

    // Compare it with the correct one.
    rDiff =
      system( "diff /tmp/testIntReport.txt testSegmentationVolumeReportIntReport-correct.txt" );
    if ( rDiff != 0 ) {
      throw runtime_error( "diff failed on testIntReport" );
    }

    // Test the tcl functions now. Just call them all and check the
    // results. Need to clear the report first.
    report.Clear();
    if ( NULL != report.mSegVol ||
         report.mlIntVols.size() != 0 ||
         report.mbUseROI ||
         NULL != report.mROIVol ||
         NULL != report.mROI ||
         NULL != report.mLUT ||
         report.mlStructures.size() != 0 ||
         !report.mbReportDirty ||
         report.mVolumeToIntensityAverageMap.size() != 0 ||
         report.mStructureToVolumeVoxelListMap.size() != 0 ) {
      stringstream ssError;
      ssError << "Error on Clear: not cleared";
      throw runtime_error( ssError.str() );
    }

    char sCommand[1024];
    int rTcl;

    sprintf( sCommand, "SetSegVolReportSegmentation %d", seg.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

    if ( NULL == report.mSegVol ) {
      stringstream ssError;
      ssError << "Error on Tcl SetSegVolReportSegmentation, mSegVol was NULL";
      throw runtime_error( ssError.str() );
    }
    if ( report.mSegVol->GetID() != seg.GetID() ) {
      stringstream ssError;
      ssError << "Error on Tcl SetSegVolReportSegmentation, mSegVol was the wrong volume (should be ID " << seg.GetID() << " but was " << report.mSegVol->GetID();
      throw runtime_error( ssError.str() );
    }

    sprintf( sCommand, "AddSegVolReportIntensityVolume %d", vol.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

    sprintf( sCommand, "AddSegVolReportIntensityVolume %d", seg.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );

    volsLoadedMap.clear();
    for ( tVolume = report.mlIntVols.begin();
          tVolume != report.mlIntVols.end(); ++tVolume ) {
      VolumeCollection* testVol = *tVolume;
      int volID = testVol->GetID();
      if ( volID != vol.GetID() && volID != seg.GetID() ) {
        stringstream ssError;
        ssError << "Error in Tcl AddSegVolReportIntensityVolume, added a volume with an unknown id " << volID;
        throw runtime_error( ssError.str() );
      }
      volsLoadedMap[volID] = true;
    }
    if ( !(volsLoadedMap[vol.GetID()] && volsLoadedMap[seg.GetID()]) ) {
      stringstream ssError;
      ssError << "Error in Tcl AddSegVolReportIntensityVolume, didn't add both volumes";
      throw runtime_error( ssError.str() );
    }

    sprintf( sCommand, "AddSegmentationToSegVolReport 1" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sprintf( sCommand, "AddSegmentationToSegVolReport 2" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sprintf( sCommand, "AddSegmentationToSegVolReport 4" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    sprintf( sCommand, "AddSegmentationToSegVolReport 5" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    structureMap.clear();
    for ( tStructure = report.mlStructures.begin();
          tStructure != report.mlStructures.end(); ++tStructure ) {
      int nStructure = *tStructure;
      if ( nStructure != 1 && nStructure != 2 &&
           nStructure != 4 && nStructure != 5 ) {
        stringstream ssError;
        ssError << "Error on Tcl AddSegmentationToSegVolReport, added an unknown structure " << nStructure;
        throw runtime_error( ssError.str() );
      }
      structureMap[nStructure] = true;
    }
    if ( !(structureMap[1] && structureMap[2] &&
           structureMap[4] && structureMap[5]) ) {
      stringstream ssError;
      ssError << "Error in Tcl AddSegmentationToSegVolReport, didn't add all structures";
      throw runtime_error( ssError.str() );
    }

    report.AddSegmentationStructure( 200 );

    sprintf( sCommand, "SetSegVolReportLUT %d", lut.GetID() );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    if ( NULL == report.mLUT ) {
      stringstream ssError;
      ssError << "Error on Tcl SetSegVolReportLUT, mLUT was NULL";
      throw runtime_error( ssError.str() );
    }
    if ( report.mLUT->GetID() != lut.GetID() ) {
      stringstream ssError;
      ssError << "Error on Tcl SetSegVolReportLUT, id didn't match";
      throw runtime_error( ssError.str() );
    }

    sprintf( sCommand, "MakeSegVolReport /tmp/testReport2.txt" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    rDiff =
      system( "diff /tmp/testReport2.txt testSegmentationVolumeReport-correct.txt" );
    if ( rDiff != 0 ) {
      throw runtime_error( "diff failed on testReport in Tcl MakeSegVolReport" );
    }

    sprintf( sCommand, "MakeSegVolIntensityReport /tmp/testIntReport2.txt" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    rDiff =
      system( "diff /tmp/testIntReport2.txt testSegmentationVolumeReportIntReport-correct.txt" );
    if ( rDiff != 0 ) {
      throw runtime_error( "diff failed on testIntReport in Tcl MakeSegVolIntensityReport" );
    }

    sprintf( sCommand, "ClearSegVolReport" );
    rTcl = Tcl_Eval( iInterp, sCommand );
    AssertTclOK( rTcl );
    if ( NULL != report.mSegVol ||
         report.mlIntVols.size() != 0 ||
         report.mbUseROI ||
         NULL != report.mROIVol ||
         NULL != report.mROI ||
         NULL != report.mLUT ||
         report.mlStructures.size() != 0 ||
         !report.mbReportDirty ||
         report.mVolumeToIntensityAverageMap.size() != 0 ||
         report.mStructureToVolumeVoxelListMap.size() != 0 ) {
      stringstream ssError;
      ssError << "Error on Tcl ClearSegVolReport: not cleared";
      throw runtime_error( ssError.str() );
    }

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


    SegmentationVolumeReportTester tester0;
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
