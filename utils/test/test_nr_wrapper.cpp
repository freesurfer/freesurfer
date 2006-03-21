#include <stdexcept>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>

// testing
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "nr_wrapper_open_source.h"

extern "C" {
  #include "nr_wrapper.h"
  #include "matrix.h"
  #include "error.h"
}
  
class NRWrapperTest : public CppUnit::TestFixture {

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( NRWrapperTest );
  
// tested by test_matrix::inverse
//   CPPUNIT_TEST( TestLUDecomp );
//   CPPUNIT_TEST( TestLUMatrixInverse );

    CPPUNIT_TEST( TestDFPMin );
    CPPUNIT_TEST( TestPowell );
    CPPUNIT_TEST( TestSVDcmp );

// tested by test_matrix::eigensystem
//    CPPUNIT_TEST( TestTred2 );
//    CPPUNIT_TEST( TestTQLI );

    CPPUNIT_TEST( TestRan1 );
    CPPUNIT_TEST( TestSplineAndSplInt );

  CPPUNIT_TEST_SUITE_END();

private:
  float mExpectedFret1;
  float *mExpectedP1;
  int mFunction1Size;

  float mExpectedFret2;
  float *mExpectedP2;
  int mFunction2Size;

  void SetUpMinimizationTestResults();
  void TearDownMinimizationTestResults();

public:
  static const bool IS_USING_VNL;

  static const int MATRICES_NOT_EQUAL;
  static const int MATRICES_EQUAL;
  static const int MATRICES_ERROR;
  
  static const float EPSILON = 5e-5;
  
  static const std::string TESTING_DIR;

  static const std::string PASCAL_MATRIX;
  static const std::string ZEROES_MATRIX;
  static const std::string IDENTITY_MATRIX;
  static const std::string SINGULAR_MATRIX;
  static const std::string ONE_MATRIX;
  static const std::string ONE_SMALL_MATRIX;
  static const std::string RANDOM_61_2;
  static const std::string RANDOM_5_11;
  
  static const std::string SINE_X;
  static const std::string SINE_Y;
  
  static const std::string SINE_XX_COURSE;
  static const std::string SINE_SPLINE_YY_0_0_COURSE;
  static const std::string SINE_SPLINE_YY_5_5_COURSE;
  static const std::string SINE_SPLINE_YY_0_17_COURSE;
  static const std::string SINE_SPLINE_YY_19_2_COURSE;
  
  static const std::string SINE_XX_FINE;
  static const std::string SINE_SPLINE_YY_0_0_FINE;
  static const std::string SINE_SPLINE_YY_5_5_FINE;
  static const std::string SINE_SPLINE_YY_0_17_FINE;
  static const std::string SINE_SPLINE_YY_19_2_FINE;

  static float TestFunction1( float *p );
  static void TestFunction1Gradient( float *p, float *g );
  
  static float TestFunction2( float *p );
  static void TestFunction2Gradient( float *p, float *g );

  static void StepFunction( int itno, float sse, void *parms, float *p );

  // setUp is called automatically before each test
  void setUp();
  
  // tearDown is called automatically after each test
  void tearDown();

  bool AreMatricesEqual( MATRIX *m1, MATRIX *m2, float tolerance,
    int numberOfRows, int numberOfColumns );

  bool AreSplinesEqual( MATRIX* x, MATRIX* y, MATRIX* xx, 
    std::string yyFile, float derivative1, float derivativeN );

  bool AreEqualWithinTolerance( const float expected, 
                                const float actual,
                                const float tolerance );

  bool AreNByNArraysEqual( float *expected, float *actual, int n, 
    float tolerance );

  bool IsBetween0And1( float x );

  MATRIX* ReconstructFromSVD( MATRIX* u, VECTOR* wVector, MATRIX* v );

  int TestSVDcmpHelper( std::string matrixFile, 
    int numberOfRows, int numberOfColumns );
      
  bool TestPowellHelper( const int numberOfParameters, float expectedFret,
                         float * expectedP, float (*function)(float []) );
                         
  bool TestDFPMinHelper( int numberOfParameters, 
    float expectedFret, float * expectedP, 
    float ( *function )(float []),
    void ( *functionGradient )(float [], float []),
    void ( *stepFunction )( int itno, float sse, void *parms, float *p ), 
    void *params );
          
  
  void TestDFPMin();
  void TestPowell();
  void TestSVDcmp();
  void TestRan1();
  void TestSplineAndSplInt();

};

const bool NRWrapperTest::IS_USING_VNL = true;

const int NRWrapperTest::MATRICES_NOT_EQUAL = 0;
const int NRWrapperTest::MATRICES_EQUAL = 1;
const int NRWrapperTest::MATRICES_ERROR = 2;

const std::string NRWrapperTest::TESTING_DIR = "test_nr_wrapper_data/";
 
const std::string NRWrapperTest::PASCAL_MATRIX = TESTING_DIR + "Pascal.mat";
const std::string NRWrapperTest::ZEROES_MATRIX = TESTING_DIR + "Zeroes.mat";
const std::string NRWrapperTest::IDENTITY_MATRIX = TESTING_DIR + "Identity.mat";
const std::string NRWrapperTest::SINGULAR_MATRIX = TESTING_DIR + "Singular.mat";
const std::string NRWrapperTest::ONE_MATRIX = TESTING_DIR + "One.mat";
const std::string NRWrapperTest::ONE_SMALL_MATRIX = 
  TESTING_DIR + "OneSmall.mat";

/**
 * This is actually a 61 x 61 matrix that has zeros outside of the 61 x 2 
 * bounds. 
 **/
const std::string NRWrapperTest::RANDOM_61_2 = TESTING_DIR + "Random_61_2.mat";
const std::string NRWrapperTest::RANDOM_5_11 = TESTING_DIR + "Random_5_11.mat";

const std::string NRWrapperTest::SINE_X = TESTING_DIR + "SineX.mat";
const std::string NRWrapperTest::SINE_Y = TESTING_DIR + "SineY.mat";

const std::string NRWrapperTest::SINE_XX_COURSE = 
  TESTING_DIR + "SineXXCourse.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_0_COURSE = 
  TESTING_DIR + "SineSplineYY_0_0_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_5_5_COURSE = 
  TESTING_DIR + "SineSplineYY_5_5_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_17_COURSE = 
  TESTING_DIR + "SineSplineYY_0_17_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_19_2_COURSE = 
  TESTING_DIR + "SineSplineYY_19_2_Course.mat";

const std::string NRWrapperTest::SINE_XX_FINE = TESTING_DIR + "SineXXFine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_0_FINE = 
  TESTING_DIR + "SineSplineYY_0_0_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_5_5_FINE = 
  TESTING_DIR + "SineSplineYY_5_5_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_17_FINE = 
  TESTING_DIR + "SineSplineYY_0_17_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_19_2_FINE = 
  TESTING_DIR + "SineSplineYY_19_2_Fine.mat";

void
NRWrapperTest::setUp() {
}

void
NRWrapperTest::tearDown() {

}

bool 
NRWrapperTest::IsBetween0And1( float x ) {
  bool isBetween = true;
  
  if( x>1.0 ) {
    isBetween = false;
  }
  
  if( x<0.0) {
    isBetween = false;
  }
  
  return isBetween;
}

bool 
NRWrapperTest::AreMatricesEqual( MATRIX *m1, MATRIX *m2, float tolerance, 
  int numberOfRows, int numberOfColumns) {
        
  bool areEqual = true;

  if( m1->rows >= numberOfRows &&
      m2->rows >= numberOfRows &&
      m1->cols >= numberOfColumns &&
      m2->cols >= numberOfColumns ) {

    for( int row=0; row<numberOfRows; row++ ) {
      for( int column=0; column<numberOfColumns; column++ ) {    

        int index1 = column + row * m1->cols;
        int index2 = column + row * m2->cols;
        
        float difference = fabs( m1->data[ index1 ] - m2->data[ index2 ] );
                
        if( difference > tolerance ) {
          areEqual = false;
        }

      }
    }
            
  } else {
    areEqual = false;    
  }
  return areEqual;
}

bool
NRWrapperTest::TestDFPMinHelper( int numberOfParameters, 
    float expectedFret, float * expectedP, 
    float ( *function )(float []),
    void ( *functionGradient )(float [], float []),
    void ( *stepFunction )( int itno, float sse, void *parms, float *p ), 
    void *params ) {
    
  int iter = 0;
  float fret = 0;
  float ftol = .1;
  
  float *p = vector( 1, numberOfParameters );
  
  dfpmin( p, numberOfParameters, ftol, &iter, &fret, function, functionGradient,
    stepFunction, params);

  bool gotExpectedValues = true;
  
  if( stepFunction != NULL ) {
    bool *isCalled = (bool *)params;
    
    if( !isCalled[0] ) {
      gotExpectedValues = false;
    }
  }

  const float equalsTolerance = 1e-4;
  gotExpectedValues = AreNByNArraysEqual( p, expectedP, numberOfParameters, 
    equalsTolerance );
    
  if (fret != expectedFret) gotExpectedValues = false;
    
  free_vector( p, 1, numberOfParameters );
  
  return gotExpectedValues;
}

void 
NRWrapperTest::SetUpMinimizationTestResults() {
  mFunction1Size = 1;
  mExpectedFret1 = -4.0;
  mExpectedP1 = new float[2];
  mExpectedP1[1] = 3.0;

  mFunction2Size = 2;
  mExpectedFret2 = -10.0;
  mExpectedP2 = new float[3];
  mExpectedP2[1] = 2.0;
  mExpectedP2[2] = -2.0;
}

void 
NRWrapperTest::TearDownMinimizationTestResults() {
  delete []mExpectedP1;
  mExpectedP1 = NULL;
  
  delete []mExpectedP2;
  mExpectedP2 = NULL;
}

void 
NRWrapperTest::TestDFPMin() {
  SetUpMinimizationTestResults();
  
  bool isSuccess = TestDFPMinHelper( mFunction1Size, mExpectedFret1, 
    mExpectedP1, TestFunction1, TestFunction1Gradient, NULL, NULL );
  CPPUNIT_ASSERT( isSuccess );

  isSuccess = TestDFPMinHelper( mFunction2Size, mExpectedFret2, mExpectedP2, 
    TestFunction2, TestFunction2Gradient, NULL, NULL );
  CPPUNIT_ASSERT( isSuccess );

  bool* isCalled = new bool[1];
  isSuccess = TestDFPMinHelper( mFunction2Size, mExpectedFret2, mExpectedP2, 
    TestFunction2, TestFunction2Gradient, StepFunction, isCalled );
  CPPUNIT_ASSERT( isSuccess );
  delete isCalled;
  
  TearDownMinimizationTestResults();
}

// static
void 
NRWrapperTest::StepFunction( int itno, float sse, void *parms, float *p ) {
  bool* isCalled = (bool *)parms;
  isCalled[0] = true;
}

// static
float 
NRWrapperTest::TestFunction1( float *p ) {
  float x = p[1];
  float y = ( x - 1 ) * ( x - 5 );
  
  return y;
}

void
NRWrapperTest::TestFunction1Gradient( float *p, float *g ) {
  float x = p[1];
  g[1] = 2*x - 6;
}

// static
float 
NRWrapperTest::TestFunction2( float *p ) {
  float Ax[2];
  Ax[0] = 3 * p[1] + 2 * p[2];
  Ax[1] = 2 * p[1] + 6 * p[2];
  
  float xTAx = p[1] * Ax[0] + p[2] * Ax[1];
  
  float bx = 2 * p[1] + (-8 * p[2]);

  return .5 * xTAx - bx;
}

void
NRWrapperTest::TestFunction2Gradient( float *p, float *g ) {
  float x1 = p[1];
  float x2 = p[2];
  
  g[1] = ( 3 * x1 ) + ( 2 * x2 ) - 2;
  g[2] = ( 6 * x2 ) + ( 2 * x1 ) + 8;
}

bool
NRWrapperTest::AreEqualWithinTolerance( const float expected, 
                                        const float actual,
                                        const float tolerance ) {
  float difference = std::fabs( expected - actual );
  return difference <= tolerance;                                            
}

void 
NRWrapperTest::TestPowell() {
  SetUpMinimizationTestResults();

  CPPUNIT_ASSERT( TestPowellHelper( mFunction1Size, mExpectedFret1, mExpectedP1, 
                                    TestFunction1 ) );

  CPPUNIT_ASSERT( TestPowellHelper( mFunction2Size, mExpectedFret2, mExpectedP2, 
                                    TestFunction2));

  TearDownMinimizationTestResults();
}

bool
NRWrapperTest::TestPowellHelper( const int numberOfParameters, 
                                 float expectedFret,
                                 float * expectedP, 
                                 float (*function)(float []) ) {                                  
  const float tolerance = 1e-5;
  float fret = 0;
  int iter = 0;
  
  float *p = vector( 1, numberOfParameters );
  float **xi = matrix( 1, numberOfParameters, 1, numberOfParameters );

  for (int i = 0; i < numberOfParameters; i++) {
    p[i + 1] = 0.0;
  }
  
  for (int i = 0; i < numberOfParameters; i++) {
    for (int j = 0; j < numberOfParameters; j++) {
      if (i == j) {
        xi[i + 1][j + 1] = 1.0;
      } else {
        xi[i + 1][j + 1] = 0.0;
      }
    }
  }
  
  powell( p, xi, numberOfParameters, tolerance, &iter, &fret, function );
    
  const float equalsTolerance = 1e-4;

  bool gotExpectedValues = AreNByNArraysEqual( p, expectedP, numberOfParameters,
    equalsTolerance );
    
  if (fret != expectedFret) gotExpectedValues = false;
  
  free_vector( p, 1, numberOfParameters );  
  free_matrix( xi, 1, numberOfParameters, 1, numberOfParameters );
  
  return gotExpectedValues;
}

bool
NRWrapperTest::AreNByNArraysEqual( float *expected, float *actual, int n,
  float tolerance ) {
  
  bool areEqual = true;
  for (int i = 0; i < n; i++) {
    if ( !AreEqualWithinTolerance( actual[i + 1], 
                                   expected[i + 1], 
                                   tolerance) ) {
      areEqual = false;
    }
  }
  
  return areEqual;
}

MATRIX*
NRWrapperTest::ReconstructFromSVD( MATRIX* u, VECTOR* wVector, MATRIX* v ) {
  int numberOfRows = u->rows;
  int numberOfColumns = v->cols;
  
  MATRIX *w = MatrixAlloc( numberOfRows, numberOfColumns, MATRIX_REAL );
  
  for( int row=0, column=0; 
       row < numberOfRows && column < numberOfColumns; 
       row++, column++ ) {
        
    w->rptr[ row+1 ][ column+1 ] = wVector->rptr[ 1 ][ column+1 ];
    
  }
    
  MATRIX* vTranspose = MatrixTranspose( v, NULL );      
  MATRIX* uw = MatrixMultiply( u, w, NULL );
  MATRIX* result = MatrixMultiply( uw, vTranspose, NULL );
  
  MatrixFree( &w );
  MatrixFree( &vTranspose );
  MatrixFree( &uw );

  return result;
}

int
NRWrapperTest::TestSVDcmpHelper( std::string matrixFile,
                                 int numberOfRows=-1, 
                                 int numberOfColumns=-1 ) {
    
  int status;
  
  MATRIX *x = MatrixRead( (char*) ( matrixFile.c_str() ) );
  
  if( numberOfRows == -1 ) {
    numberOfRows = x->rows;
  }
  
  if( numberOfColumns == -1 ) {
    numberOfColumns = x->cols;
  }
  
  MATRIX *u = MatrixCopy( x, NULL );
  VECTOR *w = RVectorAlloc( numberOfColumns, MATRIX_REAL );
  MATRIX *v = MatrixAlloc( numberOfColumns, numberOfColumns, MATRIX_REAL );
  
  int isError = NO_ERROR;
  
  // this will become the new svdcmp
  if( IS_USING_VNL ) {
    isError = OpenSvdcmp( u, numberOfRows, numberOfColumns, w, v );
  } else {
    isError = svdcmp( u->rptr, numberOfRows, numberOfColumns, 
                          w->rptr[1], v->rptr );
  }
    
  if( isError == NO_ERROR ) {
    
    // resize
    if( u->rows < u->cols ) {
      // TODO: maybe this should be ported over to the replacement algorithm too
      MATRIX *tmp = MatrixCopyRegion( u, NULL, 1, 1, u->rows, u->rows, 1, 1 );
      MatrixFree( &u );
      u = tmp;
    }

    MATRIX* result = ReconstructFromSVD( u, w, v );    
    bool areEqual = AreMatricesEqual( x, result, EPSILON, 
      numberOfRows, numberOfColumns );
                
    if( areEqual ) {
      status = MATRICES_EQUAL;
    } else {
      status = MATRICES_NOT_EQUAL;
    }    
    MatrixFree( &result );
  } else {
    status = MATRICES_ERROR;
  }
  
  MatrixFree( &x );
  MatrixFree( &u );
  VectorFree( &w );
  MatrixFree( &v );
  
  return status;
}

void 
NRWrapperTest::TestSVDcmp() {
  int status;

  status = TestSVDcmpHelper( PASCAL_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

  status = TestSVDcmpHelper( ZEROES_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

  status = TestSVDcmpHelper( IDENTITY_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

  status = TestSVDcmpHelper( SINGULAR_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

  status = TestSVDcmpHelper( ONE_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

  status = TestSVDcmpHelper( ONE_SMALL_MATRIX );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

// TODO: using vxl, this test case will fail--meaning the matrix is valid, which
// is a feature actually...
// TODO: with more columns than rows, this used to fail, but will now return
// a u of 5x11, rather than 5x5.  We should make a note of this.  I believe that
// the way that svdcmp has been used is to copy the 5x11 matrix into u initially
  status = TestSVDcmpHelper( RANDOM_5_11 );
  if( IS_USING_VNL ) {
    CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );
  } else {
    CPPUNIT_ASSERT_EQUAL( MATRICES_ERROR, status );
  }
  
  status = TestSVDcmpHelper( RANDOM_61_2, 61, 2 );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );
}

void 
NRWrapperTest::TestRan1() {
  const int numberOfBins = 10;
  int bins[numberOfBins];
  
  for(int i=0; i<numberOfBins; i++) {
    bins[i] = 0;
  }
  
  const int numberOfRuns = 1000;  
  long seed = -1L * (long)( abs( (int)time(NULL) ) );
  
  float randomNumber;
  
  for(int i=0; i<numberOfRuns; i++) {
    // TODO: done
    if( IS_USING_VNL ) {
      randomNumber = OpenRan1( &seed );
    } else {
      randomNumber = ran1( &seed );
    }
 
    CPPUNIT_ASSERT( IsBetween0And1( randomNumber ) );

    int bin = (int)floor(randomNumber * numberOfBins);
    bins[bin]++;
  }
    
  float binMean = (float)numberOfRuns/(float)numberOfBins;
    
  float x2 = 0;
  for( int i=0; i<numberOfBins; i++ ) {
    float n = (float)bins[i];
    float top = ( n - binMean );
    x2 += top * top / binMean;
  }  
  
  int x2Mean = numberOfBins - 1;
  float lowerBound = 2.0;
  float upperBound = x2Mean*3;
  
  // if x2 is close to 0, it means that the distribution it very uniform -- bad!  
  CPPUNIT_ASSERT( x2 > lowerBound );

  // if x2 is way above the mean, then it's too non-uniform
  CPPUNIT_ASSERT( x2 < upperBound );  
}

bool
NRWrapperTest::AreSplinesEqual(MATRIX* x, MATRIX* y, MATRIX* xx,
    std::string yyFile, float derivative1, float derivativeN) {    

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
NRWrapperTest::TestSplineAndSplInt() {
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
