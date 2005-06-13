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

char* Progname = "test_SegmentationVolumeReport";



class SegmentationVolumeReportTester {
public:
  void Test( Tcl_Interp* iInterp );
};

void 
SegmentationVolumeReportTester::Test ( Tcl_Interp* iInterp ) {

  try {


    char* testDataPath = getenv("FSDEV_TEST_DATA");
    if( NULL == testDataPath ) {
      throw
	runtime_error("Couldn't load test data: FSDEV_TEST_DATA not defined" );
    }
    string fnTestDataPath( testDataPath );

    string fnSegVolume = fnTestDataPath + "/anatomical/" + 
      "testSegmentationVolumeReportData-Seg.mgz";

    VolumeCollection seg;
    seg.SetFileName( fnSegVolume );
    seg.LoadVolume();

    string fnLUT = fnTestDataPath + "/lut/TestLUT.txt";

    ScubaColorLUT lut;
    lut.UseFile( fnLUT );

    SegmentationVolumeReport& report = 
      SegmentationVolumeReport::GetReport();
    
    report.SetSegmentation( seg );
    report.DontUseROI();
    
    report.SetColorLUT( lut );
    report.AddSegmentationStructure( 1 );

    report.MakeVolumeReport();

    float expectedSize = 25 * 25;
    if( report.mStructureToVolumeMap[1] != expectedSize ) {
      stringstream ssError;
      ssError << "Error on report for seg 1: expectedSize was " 
	      << expectedSize << ", but got " 
	      << report.mStructureToVolumeMap[1];
      throw runtime_error( ssError.str() );
    }

  }
  catch( exception& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
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

 
  }
  catch( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}
