/**
 * @file  test_Utilities.cpp
 * @brief test Utilities class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/11 00:06:04 $
 *    $Revision: 1.6.2.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdexcept>
#include <sstream>

#include "Utilities.h"
#include "Volume3.h"
#include "Point2.h"
#include "Point3.h"
#include "Scuba-impl.h"

extern "C" {
#include "macros.h"
}


const char* Progname = "test_Utilities";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }

class UtilitiesTester {
public:
  void Test ();
};


void
UtilitiesTester::Test () {

  try {

    Volume3<bool> bVol( 10, 15, 20, false );
    for ( int nZ = 0; nZ < 20; nZ++ ) {
      for ( int nY = 0; nY < 15; nY++ ) {
        for ( int nX = 0; nX < 10; nX++ ) {

          {
            stringstream ssError;
            bool bVal = bVol.Get( nX, nY, nZ );
            ssError << "initial value was not false  for " << nX
            << ", " << nY << ", " << nZ;
            Assert( (false == bVal), ssError.str() );
          }

          {
            stringstream ssError;
            ssError << "set failed for " << nX << ", " << nY << ", " << nZ;
            bVol.Set( nX, nY, nZ, true );
          }

          {
            stringstream ssError;
            bool bVal = bVol.Get( nX, nY, nZ );
            ssError << "set value was not true for " << nX
            << ", " << nY << ", " << nZ;
            Assert( (true == bVal), ssError.str() );
          }
        }
      }
    }

    try {
      bVol.Get( 11, 0, 0 );
      throw runtime_error( "didn't throw for illegal x access" );
    } catch (...) {}
    try {
      bVol.Get( 0, 16, 0 );
      throw runtime_error( "didn't throw for illegal x access" );
    } catch (...) {}
    try {
      bVol.Get( 0, 0, 21 );
      throw runtime_error( "didn't throw for illegal x access" );
    } catch (...) {}



    int begin[2] = {0, 0};
    int end[2] = {1, 1};
    list<Point2<int> > points;
    Utilities::FindPointsOnLine2d( begin, end, 1, points );
    {
      Point3<int> a( 0, 0, 0 );
      Point3<int> b( 5, 6, 7 );
      list<Point3<int> > l;
      Utilities::FindPointsOnLine3d( a.xyz(), b.xyz(), l );

      // Our list should be: (0,0,0) (0,1,1) (1,1,1) (1,2,2) (2,2,3)
      // (2,3,3) (3,3,4) (3,4,5) (4,5,5) (4,5,6) (5,6,7)
      int nPoint = 0;
      list<Point3<int> >::iterator t;
      for ( t = l.begin(); t != l.end(); ++t, nPoint++ ) {
        stringstream ssError;
        bool bGood = false;
        Point3<int> p = *t;
        switch ( nPoint ) {
        case 0:
          bGood = p[0] == 0 && p[1] == 0 && p[2] == 0;
          break;
        case 1:
          bGood = p[0] == 0 && p[1] == 1 && p[2] == 1;
          break;
        case 2:
          bGood = p[0] == 1 && p[1] == 1 && p[2] == 1;
          break;
        case 3:
          bGood = p[0] == 1 && p[1] == 2 && p[2] == 2;
          break;
        case 4:
          bGood = p[0] == 2 && p[1] == 2 && p[2] == 3;
          break;
        case 5:
          bGood = p[0] == 2 && p[1] == 3 && p[2] == 3;
          break;
        case 6:
          bGood = p[0] == 3 && p[1] == 3 && p[2] == 4;
          break;
        case 7:
          bGood = p[0] == 3 && p[1] == 4 && p[2] == 5;
          break;
        case 8:
          bGood = p[0] == 4 && p[1] == 5 && p[2] == 5;
          break;
        case 9:
          bGood = p[0] == 4 && p[1] == 5 && p[2] == 6;
          break;
        case 10:
          bGood = p[0] == 5 && p[1] == 6 && p[2] == 7;
          break;
        default:
          ssError << "Too many points!";
          break;
        }
        ssError << "Incorrect point at position " << nPoint << endl;
        Assert( (bGood), ssError.str() );
      }


    }


    {
      string sToSplit = "option1:option2=value2:option3";
      string sDelimiter = ":";
      vector<string> lsResults;

      int cFound = Utilities::SplitString( sToSplit, sDelimiter, lsResults );

      {
        stringstream ssError;
        ssError << "cFound is incorrect: should be 2, is " << cFound;
        Assert( (2 == cFound), ssError.str() );
      }
      {
        stringstream ssError;
        ssError << "size of lsResults is incorrect: should be 3, is "
        << lsResults.size();
        Assert( (lsResults.size() == 3), ssError.str() );
      }
      string sResult = lsResults[0];
      {
        stringstream ssError;
        ssError << "Result 0 is incorrect, should be option1, is " << sResult;
        Assert( (sResult == "option1"), ssError.str() );
      }
      {
        stringstream ssError;
        sResult = lsResults[1];
        ssError << "Result 1 is incorrect, should be option1, is " << sResult;
        Assert( (sResult == "option2=value2"), ssError.str() );
      }
      {
        stringstream ssError;
        sResult = lsResults[2];
        ssError << "Result 2 is incorrect, should be option1, is " << sResult;
        Assert( (sResult == "option3"), ssError.str() );
      }
    }

    {
      string sToSplit = "option1";
      string sDelimiter = ":";
      vector<string> lsResults;

      int cFound = Utilities::SplitString( sToSplit, sDelimiter, lsResults );

      {
        stringstream ssError;
        ssError << "cFound is incorrect: should be 0, is " << cFound;
        Assert( (0 == cFound), ssError.str() );
      }
      {
        stringstream ssError;
        ssError << "size of lsResults is incorrect: should be 0, is "
        << lsResults.size();
        Assert( (lsResults.size() == 0), ssError.str() );
      }
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

    UtilitiesTester tester0;
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

