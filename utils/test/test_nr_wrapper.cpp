#include <stdexcept>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>

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
  
class NRWrapperTest : public CppUnit::TestFixture {

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( NRWrapperTest );
  
    CPPUNIT_TEST( TestLUDecomp );
    CPPUNIT_TEST( TestDFPMin );
    CPPUNIT_TEST( TestPowell );
    CPPUNIT_TEST( TestSVDCmp );
    CPPUNIT_TEST( TestTred2 );
    CPPUNIT_TEST( TestTQLI );
    CPPUNIT_TEST( TestRan1 );
    CPPUNIT_TEST( TestSplineAndSplint );
    CPPUNIT_TEST( TestLUMatrixInverse );

  CPPUNIT_TEST_SUITE_END();

private:
  MATRIX* mIdentityMatrix;
  MATRIX* mZeroMatrix;
  MATRIX* mWellFormedMatrix;
  MATRIX* mPascalMatrix;

public:
  static const float EPSILON = 1e-5;
  
  static const string TESTING_DIR;

/*
  static const string WELL_FORMED_MATRIX;
  static const string WELL_FORMED_LUDCMP;

  static const string PASCAL_MATRIX;
  static const string PASCAL_LUDCMP;
  static const string PASCAL_INVERSE;
*/  
  
  static const string SINE_X;
  static const string SINE_Y;
  
  static const string SINE_XX_COURSE;
  static const string SINE_SPLINE_YY_0_0_COURSE;
  static const string SINE_SPLINE_YY_5_5_COURSE;
  static const string SINE_SPLINE_YY_0_17_COURSE;
  static const string SINE_SPLINE_YY_19_2_COURSE;
  
  static const string SINE_XX_FINE;
  static const string SINE_SPLINE_YY_0_0_FINE;
  static const string SINE_SPLINE_YY_5_5_FINE;
  static const string SINE_SPLINE_YY_0_17_FINE;
  static const string SINE_SPLINE_YY_19_2_FINE;


  // setUp is called automatically before each test
  void setUp();
  
  // tearDown is called automatically after each test
  void tearDown();

  bool AreMatricesEqual( MATRIX *m1, MATRIX *m2 );
  bool AreSplinesEqual(MATRIX* x, MATRIX* y, MATRIX* xx, 
    string yyFile, float derivative1, float derivativeN);
  
  void TestLUDecomp();
  void TestDFPMin();
  void TestPowell();
  void TestSVDCmp();
  void TestTred2();
  void TestTQLI();
  void TestRan1();
  void TestSplineAndSplint();
  void TestLUMatrixInverse();

};

const string NRWrapperTest::TESTING_DIR = "test_nr_wrapper_data/";

/*
const string NRWrapperTest::WELL_FORMED_MATRIX = TESTING_DIR + "WellFormed.mat";
const string NRWrapperTest::WELL_FORMED_LUDCMP = TESTING_DIR + "WellFormedLUDCMP.mat";

const string NRWrapperTest::PASCAL_MATRIX = TESTING_DIR + "Pascal.mat";
const string NRWrapperTest::PASCAL_LUDCMP = TESTING_DIR + "PascalLUDCMP.mat";
const string NRWrapperTest::PASCAL_INVERSE = TESTING_DIR + "PascalInverse.mat";
*/

const string NRWrapperTest::SINE_X = TESTING_DIR + "SineX.mat";
const string NRWrapperTest::SINE_Y = TESTING_DIR + "SineY.mat";

const string NRWrapperTest::SINE_XX_COURSE = TESTING_DIR + "SineXXCourse.mat";
const string NRWrapperTest::SINE_SPLINE_YY_0_0_COURSE = TESTING_DIR + "SineSplineYY_0_0_Course.mat";
const string NRWrapperTest::SINE_SPLINE_YY_5_5_COURSE = TESTING_DIR + "SineSplineYY_5_5_Course.mat";
const string NRWrapperTest::SINE_SPLINE_YY_0_17_COURSE = TESTING_DIR + "SineSplineYY_0_17_Course.mat";
const string NRWrapperTest::SINE_SPLINE_YY_19_2_COURSE = TESTING_DIR + "SineSplineYY_19_2_Course.mat";

const string NRWrapperTest::SINE_XX_FINE = TESTING_DIR + "SineXXFine.mat";
const string NRWrapperTest::SINE_SPLINE_YY_0_0_FINE = TESTING_DIR + "SineSplineYY_0_0_Fine.mat";
const string NRWrapperTest::SINE_SPLINE_YY_5_5_FINE = TESTING_DIR + "SineSplineYY_5_5_Fine.mat";
const string NRWrapperTest::SINE_SPLINE_YY_0_17_FINE = TESTING_DIR + "SineSplineYY_0_17_Fine.mat";
const string NRWrapperTest::SINE_SPLINE_YY_19_2_FINE = TESTING_DIR + "SineSplineYY_19_2_Fine.mat";

void
NRWrapperTest::setUp() {
/*  
  MATRIX *identityMatrix = NULL;
  identityMatrix = MatrixIdentity( 3, identityMatrix );
  mIdentityMatrix = identityMatrix;
  
  MATRIX *zeroMatrix = MatrixAlloc( 3, 3, MATRIX_REAL );
  MatrixClear(zeroMatrix);
  mZeroMatrix = zeroMatrix;
  
  mWellFormedMatrix = MatrixRead( (char*)( WELL_FORMED_MATRIX.c_str() ) );
  
  mPascalMatrix = MatrixRead( (char*)( PASCAL_MATRIX.c_str() ) );
*/  
}

void
NRWrapperTest::tearDown() {
/*  
  delete mIdentityMatrix;
  delete mZeroMatrix;
*/  
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
    }
            
  } else {
    areEqual = false;    
  }
  return areEqual;
}

void 
NRWrapperTest::TestLUDecomp() {
/*
  int ERROR = -1;
    
//  MATRIX *inverseMatrix = NULL;
  inverseMatrix = MatrixInverse( mIdentityMatrix, inverseMatrix );  
  MatrixPrint(stdout, inverseMatrix);    
  CPPUNIT_ASSERT( AreMatricesEqual( mIdentityMatrix, inverseMatrix ) );
  
  MATRIX *matrix = mPascalMatrix;
  
  // LU Decomposition, takes
  int numberOfRows = matrix->rows;
  MATRIX *decomposedMatrix = MatrixCopy( matrix, NULL );
  float **decomposedMatrixData = decomposedMatrix->rptr;
  int *index = (int *)calloc( numberOfRows+1, sizeof( int ) );
  float d;

  cerr << "--- original" << endl;
  MatrixPrint( stdout, matrix );

  cerr << "+++ ludcmp calculated" << endl;

  int isError = ludcmp( decomposedMatrixData, numberOfRows, index, &d );
  if( isError == ERROR ) {
    cerr << "--- Error ----" << endl;
    MatrixPrint( stdout, decomposedMatrix );
  } else {
    MatrixPrint( stdout, decomposedMatrix );
  }
  
  cerr << "+++ ludcmp index" << endl;
  for(int i=0; i<numberOfRows; i++) {
    cerr << index[0] << endl;
  }
  
  cerr << "*** ludcmp matlab" << endl;  
  MATRIX* pascalLUDCMP = MatrixRead( (char*)( PASCAL_LUDCMP.c_str() ));
  MatrixPrint( stdout, pascalLUDCMP );
*/  
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestDFPMin() {
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestPowell() {
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestSVDCmp() {
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestTred2() {
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestTQLI() {
  CPPUNIT_FAIL("not implemented");
}

void 
NRWrapperTest::TestRan1() {

/*    
  long seed = -1L * (long)(abs((int)time(NULL))) ;

  for(int i=0; i<25; i++) {
    float val = ran1(&seed);  
    cerr << " " << val << " ";
  }
*/  
  
  CPPUNIT_FAIL("not implemented");
}

bool
NRWrapperTest::AreSplinesEqual(MATRIX* x, MATRIX* y, MATRIX* xx,
    string yyFile, float derivative1, float derivativeN) {    

  bool areEqual = true;
  
  float* inputX = x->data;
  float* inputY = y->data;
  
  int n = x->cols;      
  float output[300];
    
  MATRIX* yy = MatrixRead( (char*) ( yyFile.c_str() ) );
  spline(inputX, inputY, n, derivative1, derivativeN, output);

  for(int i=0; i<xx->cols; i++) {    
    float pointX = xx->data[i];
    float result;

    splint(inputX, inputY, output, n, pointX, &result);
    
    float expected = yy->data[i];
        
    if( fabs( result - expected ) > EPSILON) {
      areEqual = false;
    }
  }
  
  MatrixFree( &yy );
    
  return areEqual;
}

void 
NRWrapperTest::TestSplineAndSplint() {
  MATRIX* x = MatrixRead( (char*) ( SINE_X.c_str() ) );
  MATRIX* y = MatrixRead( (char*) ( SINE_Y.c_str() ) );

  MATRIX* xxCourse = MatrixRead( (char*) ( SINE_XX_COURSE.c_str() ) );
  MATRIX* xxFine = MatrixRead( (char*) ( SINE_XX_FINE.c_str() ) );

  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxCourse, SINE_SPLINE_YY_0_0_COURSE,
    0.0, 0.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxCourse, SINE_SPLINE_YY_5_5_COURSE,
    5.0, 5.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxCourse, SINE_SPLINE_YY_0_17_COURSE,
    0.0, 17.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxCourse, SINE_SPLINE_YY_19_2_COURSE,
    19.0, 2.0) );  
    
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxFine, SINE_SPLINE_YY_0_0_FINE,
    0.0, 0.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxFine, SINE_SPLINE_YY_5_5_FINE,
    5.0, 5.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxFine, SINE_SPLINE_YY_0_17_FINE,
    0.0, 17.0) );  
  CPPUNIT_ASSERT( AreSplinesEqual(x, y, xxFine, SINE_SPLINE_YY_19_2_FINE,
    19.0, 2.0) );  
    
  MatrixFree( &x );
  MatrixFree( &y ); 
  MatrixFree( &xxCourse );
  MatrixFree( &xxFine );
}

void
NRWrapperTest::TestLUMatrixInverse() {
  CPPUNIT_FAIL("not implemented");
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
}

