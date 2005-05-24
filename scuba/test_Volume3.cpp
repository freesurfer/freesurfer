#include <sstream>
#include "Volume3.h"
#include "Scuba-impl.h"

char* Progname = "test_Volume3";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


class VolumeTester {
public:
  void Test();
};

void 
VolumeTester::Test () {

  stringstream ssError;

  try {

    Volume3<bool> b(4, 4, 4, false );
    


  }
  catch( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    VolumeTester tester0;
    for( int i = 0; i < 1; i++ ) {
      tester0.Test();
    }

 
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

