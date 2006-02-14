#include <stdexcept>
#include <sstream>
#include <iostream>

extern "C" {
#include "mri.h"
#include "NrrdIO.h"
}

using namespace std;

char *Progname = NULL; // to satisfy mriio.c's extern Progname declaration

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  throw runtime_error( ss.str() ); \
  }

class MriioTester {
public:
  void TestNrrdIO ();
};


void 
MriioTester::TestNrrdIO () {
  cerr << "Check that mriNrrdRead reads foolc.nrrd...";
  Assert(MRIreadType("test_mriio_data/foolc.nrrd", NRRD_FILE) != NULL,
	 "failed.");
  cerr << "passed." << endl;
}


int main ( int argc, char** argv ) {
  //Progname = argv[0];

  cerr << "Beginning tests..." << endl;

  try {

    MriioTester tester;
    tester.TestNrrdIO();
 
  } catch( runtime_error& e ) {
    cerr << "failed " << endl << "exception: " << e.what() << endl;
    exit( 1 );
  } catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

