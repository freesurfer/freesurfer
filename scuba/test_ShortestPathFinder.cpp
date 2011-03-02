/**
 * @file  test_ShortestPathFinder.cpp
 * @brief test ShortestPathFinder class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.7 $
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
#include "ShortestPathFinder.h"
extern "C" {
#include "macros.h"
}
#include "Scuba-impl.h"

const char* Progname = "test_ShortestPathFinder";

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

// These edge costs set us up so that we should be straight from 2,2 to 8,8
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
  {  9,   9,   9,   9,   9,   9,   9,   9,   9,   9}
};

class TestShortestPathFinder : public ShortestPathFinder {
public:
  TestShortestPathFinder () {
    SetDimensions( 10, 10 );
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

    Point2<int> begin( 1, 1 );
    Point2<int> end( 8, 8 );
    list<Point2<int> > points;

    finder.SetStartPoint( begin );
    finder.FindPathToEndPoint( end, points );

    // Our path should be straight from 2,2 to 8,8.
    list<Point2<int> >::iterator tPoint;
    int nPoint = 0;
    for ( tPoint = points.begin(); tPoint != points.end(); ++tPoint, nPoint++ ) {
      Point2<int> p = *tPoint;
      stringstream ssError;
      bool bGood = false;
      switch ( nPoint ) {
      case 0:
        bGood = p[0] == 2 && p[1] == 2;
        break;
      case 1:
        bGood = p[0] == 3 && p[1] == 3;
        break;
      case 2:
        bGood = p[0] == 4 && p[1] == 4;
        break;
      case 3:
        bGood = p[0] == 5 && p[1] == 5;
        break;
      case 4:
        bGood = p[0] == 6 && p[1] == 6;
        break;
      case 5:
        bGood = p[0] == 7 && p[1] == 7;
        break;
      case 6:
        bGood = p[0] == 8 && p[1] == 8;
        break;
      default:
        ssError << "Too many points!" << endl;
        break;
      }
      ssError << "Wrong point at index " << nPoint << ", was " << p;
      Assert( (bGood), ssError.str() );
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

    ShortestPathFinderTester tester0;
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

