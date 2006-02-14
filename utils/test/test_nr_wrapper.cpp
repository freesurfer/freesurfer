#include <stdexcept>
#include <sstream>
#include <iostream>

// testing
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

extern "C" {
  #include "nr_wrapper.h"
  #include "matrix.h"
}

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }
  
class NRWrapperTest : public CppUnit::TestFixture {

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( NRWrapperTest );
  CPPUNIT_TEST( TestLUDecomp );
  CPPUNIT_TEST_SUITE_END();

private:
  MATRIX* mIdentityMatrix;

public:  
  // setUp is called automatically before each test
  void setUp();
  
  // tearDown is called automatically after each test
  void tearDown();

  void PrintMatrix( MATRIX* matrix );
  bool AreMatricesEqual( MATRIX *m1, MATRIX *m2 );
  
  void TestLUDecomp();
};

void
NRWrapperTest::setUp() {
  MATRIX *identityMatrix = NULL;
  identityMatrix = MatrixIdentity( 3, identityMatrix );
  mIdentityMatrix = identityMatrix;
}

void
NRWrapperTest::tearDown() {
  delete mIdentityMatrix;
}

bool 
NRWrapperTest::AreMatricesEqual( MATRIX *m1, MATRIX *m2 ) {
  bool areEqual = true;

  if( m1->rows == m2->rows &&
      m2->cols == m2->cols ) {  

    int numberOfRows =  m1->rows;
    int numberOfColumns = m1->cols;

    for( int row=0; row<numberOfRows; row++ ) {
      for( int column=0; column<numberOfColumns; column++ ) {    

        int index = column + row * numberOfColumns;
        if( m1->data[ index ] != m2->data[ index ] ) {
          areEqual = false;
        }        

      }
      cerr << endl;
    }
            
  } else {
    areEqual = false;    
  }
  return areEqual;
}

void 
NRWrapperTest::PrintMatrix( MATRIX *matrix ) {
  int numberOfRows =  matrix->rows;
  int numberOfColumns = matrix->cols;
  
  cerr << "number of rows:    " <<  numberOfRows << endl;  
  cerr << "number of columns: " <<  numberOfColumns << endl;  
  
  cerr << "matrix:" << endl;
  for( int row=0; row<numberOfRows; row++ ) {
    for( int column=0; column<numberOfColumns; column++ ) {
      cerr << "  " << matrix->data[ column + row*numberOfColumns ] << " ";
    }
    cerr << endl;
  }
}

void 
NRWrapperTest::TestLUDecomp() {
  MATRIX *identityMatrix = NULL;
  identityMatrix = MatrixIdentity( 3, identityMatrix );

  MATRIX *inverseMatrix = NULL;
  
  PrintMatrix( identityMatrix );

  inverseMatrix = MatrixInverse( identityMatrix, inverseMatrix );
  PrintMatrix( inverseMatrix );

  Assert( AreMatricesEqual( identityMatrix, inverseMatrix ),
    "matrices not equal...");
      
  Assert(true, "honky tonk...");
}

int main ( int argc, char** argv ) {
  
  int SUCCESS = 0;
  int FAIL = 1;
  
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( NRWrapperTest::suite() );
  
  if( runner.run() ) {
    exit ( SUCCESS );
  }
  
  exit( FAIL );
  
/*
  cerr << "Beginning tests..." << endl;

  try {

    NRWrapperTester tester0;
    tester0.TestLUDecomp();
 
  } catch( runtime_error& e ) {

    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );

  } catch( ... ) {

    cerr << "failed" << endl;
    exit( 1 );

  }

  cerr << "Success" << endl;

  exit( 0 );
*/  
}

