#include <iostream>
#include <stdexcept>
#include "IDTracker.h"

using namespace std;


class A : public IDTracker<A> {
  
public:
  A() {}
};

template IDTracker<A>;
int IDTracker<A>::mNextID = 0;
std::map<int,A*> IDTracker<A>::mIDMap;

int main ( int argc, char** argv ) {
  
  try { 

    {
      A a0;
      A a1;
      A a2;
      A a3;
      A a4;
      
      if( !(a0.GetID() == 0) ||
	  !(a1.GetID() == 1) ||
	  !(a2.GetID() == 2) ||
	  !(a3.GetID() == 3) ||
	  !(a4.GetID() == 4) ) {
	throw logic_error("IDs didn't match");
      }
      
      A findA0 = A::FindByID( 0 );
      A findA1 = A::FindByID( 1 );
      A findA2 = A::FindByID( 2 );
      A findA3 = A::FindByID( 3 );
      A findA4 = A::FindByID( 4 );
      if( !(a0.GetID() == findA0.GetID()) ||
	  !(a1.GetID() == findA1.GetID()) ||
	  !(a2.GetID() == findA2.GetID()) ||
	  !(a3.GetID() == findA3.GetID()) ||
	  !(a4.GetID() == findA4.GetID()) ) {
	throw logic_error("IDs didn't match");
      }

      
      bool bGood = false;
      try {
	A findNonExistent = A::FindByID( 10 );
	bGood = false;
      }
      catch(exception e) {
	bGood = true;
      }
      if( !bGood ) throw logic_error( "found a nonexistant object" );


      list<int> idList;
      A::GetIDList( idList );
      bool abFound[5];
      for( int nObj = 0; nObj < 5; nObj++ ) {
	abFound[nObj] = false;
      }
      list<int>::iterator tID;
      for( tID = idList.begin(); tID != idList.end(); ++tID ) {
	int id = *tID;
	if( id >= 0 && id < 5 ) { abFound[id] = true; }
	else { throw logic_error( "found invalid id in id list" ); }
      }
      for( int nObj = 0; nObj < 5; nObj++ ) {
	if( !abFound[nObj] ) {
	  throw logic_error( "valid id not found in id list" ); 
	}
      }

    }

    bool bGood = false;
    try {
      A findNonExistent = A::FindByID( 0 );
      bGood = false;
    }
    catch(exception e) {
      bGood = true;
    }
    if( !bGood ) throw logic_error( "found object after its deletion" );

  }
  catch(exception e) {
    cerr << "failed: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "success" << endl;
  exit( 0 );
}
