
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "DataCollection.h"

#define Assert(x,e)   if(!(x)) { throw e; }

using namespace std;

char* Progname = "test_DataCollection";

class TestCollection : public DataCollection {

public:
  TestCollection( string isLabel ) :
    DataCollection( isLabel ) { 
  }
  
  virtual void GetInfoAtRAS( float const iX, float const iY, float const iZ,
			     std::list<string> olLabels,
			     std::list<string> olValues ) const {
    return;
  }
};


int main ( int argc, char** argv ) {

  try {

    TestCollection col1( "col1" );
    Assert( (col1.GetLabel() == "col1"), logic_error("col1 label incorrect") );
    TestCollection col2( "col2" );
    Assert( (col2.GetLabel() == "col2"), logic_error("col2 label incorrect") );
    TestCollection col3( "col3" );
    Assert( (col3.GetLabel() == "col3"), logic_error("col3 label incorrect") );
    Assert( (col1.GetID() != col2.GetID() != col3.GetID()),
	    logic_error("not unique IDs") );


    TestCollection* col4 = new TestCollection( "col4" );
    DataCollection::ID col4ID = col4->GetID();
    DataCollection col4comp = DataCollection::GetDataCollection( col4ID );
    Assert( (col4ID == col4comp.GetID()), 
	    logic_error("Didn't get correct collection") );
    Assert( (col4->GetLabel() == col4comp.GetLabel()), 
	    logic_error("Didn't get correct label") );

    delete col4;
    try {
      DataCollection col4comp = DataCollection::GetDataCollection( col4ID );
      throw( logic_error("Didn't throw on deleted collection") );
    }
    catch( ... ) {}
    

    
  }
  catch( exception e ) {
    cerr << "failed: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
