#include <sstream>
#include "ShortestPathFinder.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

char* Progname = "test_ScubaTransform";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


class ShortestPathFinderTester {
public:
  void Test();
};

float const maEdgeCosts[10][10] = {
  {  9,   9,   9,   9,   9,   9,   9,   9,   9,   9},
  {  9,   1,   1,   1,   1,   1,   1,   1,   1,   9},
  {  9,   7,   1,   9,   1,   1,   1,   9,   1,   9},
  {  9,   1,   7,   1,   9,   1,   1,   1,   1,   9},
  {  9,   5,   3,   7,   1,   9,   1,   9,   1,   9},
  {  9,   5,   3,   2,   7,   1,   9,   9,   1,   9},
  {  9,   3,   3,   5,   5,   7,   1,   9,   1,   9},
  {  9,   3,   3,   3,   3,   5,   7,   1,   1,   9},
  {  9,   5,   5,   5,   3,   3,   1,   7,   1,   9},
  {  9,   9,   9,   9,   9,   9,   9,   9,   9,   9} };

class TestShortestPathFinder : public ShortestPathFinder {
public:
  TestShortestPathFinder () {
    SetDimensions( 10, 10, 9 );
  }

  float GetEdgeCost( Point2<int> &iPoint ) {
    return maEdgeCosts[iPoint.y()][iPoint.x()];
  }

};


void 
ShortestPathFinderTester::Test () {

  stringstream ssError;

  try {

    TestShortestPathFinder finder;
    finder.SetDebug( true );

    Point2<int> begin( 1, 1 );
    Point2<int> end( 8, 1 );
    list<Point2<int> > points;
    finder.FindPath( begin, end, points );
    
  }
  catch( runtime_error e ) {
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

    ShortestPathFinderTester tester0;
    tester0.Test();

 
  }
  catch( runtime_error e ) {
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

