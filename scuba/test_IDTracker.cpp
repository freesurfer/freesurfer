/**
 * @file  test_IDTracker.cpp
 * @brief Test file for IDTracker
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


#include <iostream>
#include <stdexcept>
#include "IDTracker.h"
#include "Scuba-impl.h"

using namespace std;

const char* Progname = "test_IDTracker";

class A : public IDTracker<A> {
public:
  A() {}
};

class B : public IDTracker<B> {
public:
  B() {}
};

int main ( int argc, char** argv ) {

  try {

    {
      A a0;
      A a1;
      A a2;
      A a3;
      A a4;

      if ( !(a0.GetID() == 0) ||
           !(a1.GetID() == 1) ||
           !(a2.GetID() == 2) ||
           !(a3.GetID() == 3) ||
           !(a4.GetID() == 4) ) {
	stringstream ssError;
	ssError << "IDs didn't match. ID of a0=" << a0.GetID()
		<< " a1=" << a1.GetID()
		<< " a2=" << a2.GetID()
		<< " a3=" << a3.GetID()
		<< " a4=" << a4.GetID();
        throw logic_error( ssError.str() );
      }

      A& findA0 = A::FindByID( 0 );
      A& findA1 = A::FindByID( 1 );
      A& findA2 = A::FindByID( 2 );
      A& findA3 = A::FindByID( 3 );
      A& findA4 = A::FindByID( 4 );
      if ( !(a0.GetID() == findA0.GetID()) ||
           !(a1.GetID() == findA1.GetID()) ||
           !(a2.GetID() == findA2.GetID()) ||
           !(a3.GetID() == findA3.GetID()) ||
           !(a4.GetID() == findA4.GetID()) ) {
        throw logic_error("IDs from FindByID didn't match");
      }


      bool bGood = false;
      try {
        A findNonExistent = A::FindByID( 10 );
        bGood = false;
      } catch ( exception& e) {
        bGood = true;
      }
      if ( !bGood ) throw logic_error( "found a nonexistant object" );


      list<int> idList;
      A::GetIDList( idList );
      bool abFound[5];
      for ( int nObj = 0; nObj < 5; nObj++ ) {
        abFound[nObj] = false;
      }
      list<int>::iterator tID;
      for ( tID = idList.begin(); tID != idList.end(); ++tID ) {
        int id = *tID;
        if ( id >= 0 && id < 5 ) {
          abFound[id] = true;
        } else {
          throw logic_error( "found invalid id in id list" );
        }
      }
      for ( int nObj = 0; nObj < 5; nObj++ ) {
        if ( !abFound[nObj] ) {
          throw logic_error( "valid id not found in id list" );
        }
      }

    }

    bool bGood = false;
    try {
      A findNonExistent = A::FindByID( 0 );
      bGood = false;
    } catch ( exception& e) {
      bGood = true;
    }
    if ( !bGood ) throw logic_error( "found object after its deletion" );


    B* b0 = new B();
    list<int> bList;
    //    list<int>::iterator tID;
    bList.clear();
    B::GetIDList( bList );
    if ( bList.size() != 1 ) {
      throw runtime_error( "bList wasn't correct size after 1 init" );
    }

    B* b1 = new B;
    bList.clear();
    B::GetIDList( bList );
    if ( bList.size() != 2 ) {
      throw runtime_error( "bList wasn't correct size after 2 init" );
    }

    B* b2 = new B;
    bList.clear();
    B::GetIDList( bList );
    if ( bList.size() != 3 ) {
      throw runtime_error( "bList wasn't correct size after 3 init" );
    }

    delete b0;
    delete b1;
    delete b2;
    bList.clear();
    B::GetIDList( bList );
    if ( bList.size() != 0 ) {
      stringstream ssError;
      ssError << "ID list not empty after deleting all tracked objects:"
      << endl << "\t";
      list<int>::iterator tID;
      for ( tID = bList.begin(); tID != bList.end(); ++tID ) {
        int id = *tID;
        ssError << id << " ";
      }
      ssError << endl;
      throw runtime_error( ssError.str() );
    }

  } catch ( exception& e) {
    cerr << "failed: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "success" << endl;
  exit( 0 );
}
