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

extern "C" {
  #include "nr_wrapper.h"
  #include "matrix.h"
  #include "error.h"
}
  
class NRWrapperTest : public CppUnit::TestFixture {

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( NRWrapperTest );
  
//TODO: tested by test_matrix::inverse
//   CPPUNIT_TEST( TestLUDecomp );
//   CPPUNIT_TEST( TestLUMatrixInverse );

    CPPUNIT_TEST( TestDFPMin );
    CPPUNIT_TEST( TestPowell );
    CPPUNIT_TEST( TestSVDcmp );

// TODO: tested by test_matrix::eigensystem
//    CPPUNIT_TEST( TestTred2 );
//    CPPUNIT_TEST( TestTQLI );

    CPPUNIT_TEST( TestRan1 );
    CPPUNIT_TEST( TestSplineAndSplInt );

  CPPUNIT_TEST_SUITE_END();

private:
  MATRIX* mIdentityMatrix;
  MATRIX* mZeroMatrix;
  MATRIX* mWellFormedMatrix;
  MATRIX* mPascalMatrix;

public:
  static const int MATRICES_NOT_EQUAL;
  static const int MATRICES_EQUAL;
  static const int MATRICES_ERROR;
  
  static const float EPSILON = 5e-5;
  
  static const std::string TESTING_DIR;

//  static const std::string WELL_FORMED_MATRIX;
//  static const std::string WELL_FORMED_LUDCMP;

  static const std::string PASCAL_MATRIX;
//  static const std::string PASCAL_LUDCMP;
//  static const std::string PASCAL_INVERSE;
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

  static float testFunction1(float *p);
  static float testFunction2(float *p);

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

  MATRIX* ReconstructFromSVD( MATRIX* u, VECTOR* wVector, MATRIX* v );

  int TestSVDcmpHelper( std::string matrixFile, 
      int numberOfRows, int nubmerOfColumns );
      
  bool TestPowellHelper( const int numberOfParameters, float expectedFret,
                         float * expectedP, float (*function)(float []) );
      
  bool isBetween0And1( float x );
  
  void TestLUDecomp();
  void TestDFPMin();
  void TestPowell();
  void TestSVDcmp();
  void TestTred2();
  void TestTQLI();
  void TestRan1();
  void TestSplineAndSplInt();
  void TestLUMatrixInverse();

};

const int NRWrapperTest::MATRICES_NOT_EQUAL = 0;
const int NRWrapperTest::MATRICES_EQUAL = 1;
const int NRWrapperTest::MATRICES_ERROR = 2;

const std::string NRWrapperTest::TESTING_DIR = "test_nr_wrapper_data/";

//const std::string NRWrapperTest::WELL_FORMED_MATRIX = TESTING_DIR + "WellFormed.mat";
//const std::string NRWrapperTest::WELL_FORMED_LUDCMP = TESTING_DIR + "WellFormedLUDCMP.mat";

const std::string NRWrapperTest::PASCAL_MATRIX = TESTING_DIR + "Pascal.mat";
//const std::string NRWrapperTest::PASCAL_LUDCMP = TESTING_DIR + "PascalLUDCMP.mat";
//const std::string NRWrapperTest::PASCAL_INVERSE = TESTING_DIR + "PascalInverse.mat";

const std::string NRWrapperTest::ZEROES_MATRIX = TESTING_DIR + "Zeroes.mat";
const std::string NRWrapperTest::IDENTITY_MATRIX = TESTING_DIR + "Identity.mat";
const std::string NRWrapperTest::SINGULAR_MATRIX = TESTING_DIR + "Singular.mat";
const std::string NRWrapperTest::ONE_MATRIX = TESTING_DIR + "One.mat";
const std::string NRWrapperTest::ONE_SMALL_MATRIX = TESTING_DIR + "OneSmall.mat";

/**
 * This is actually a 61 x 61 matrix that has zeros outside of the 61 x 2 
 * bounds. 
 **/
const std::string NRWrapperTest::RANDOM_61_2 = TESTING_DIR + "Random_61_2.mat";
const std::string NRWrapperTest::RANDOM_5_11 = TESTING_DIR + "Random_5_11.mat";

const std::string NRWrapperTest::SINE_X = TESTING_DIR + "SineX.mat";
const std::string NRWrapperTest::SINE_Y = TESTING_DIR + "SineY.mat";

const std::string NRWrapperTest::SINE_XX_COURSE = TESTING_DIR + "SineXXCourse.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_0_COURSE = TESTING_DIR + "SineSplineYY_0_0_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_5_5_COURSE = TESTING_DIR + "SineSplineYY_5_5_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_17_COURSE = TESTING_DIR + "SineSplineYY_0_17_Course.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_19_2_COURSE = TESTING_DIR + "SineSplineYY_19_2_Course.mat";

const std::string NRWrapperTest::SINE_XX_FINE = TESTING_DIR + "SineXXFine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_0_FINE = TESTING_DIR + "SineSplineYY_0_0_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_5_5_FINE = TESTING_DIR + "SineSplineYY_5_5_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_0_17_FINE = TESTING_DIR + "SineSplineYY_0_17_Fine.mat";
const std::string NRWrapperTest::SINE_SPLINE_YY_19_2_FINE = TESTING_DIR + "SineSplineYY_19_2_Fine.mat";

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
NRWrapperTest::isBetween0And1( float x ) {
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

// static
float 
NRWrapperTest::testFunction1(float *p) {
  float x = p[1];
  float y = ( x - 1 ) * ( x - 5 );
  
  return y;
}

// static
float 
NRWrapperTest::testFunction2(float *p) {
  float Ax[2];
  Ax[0] = 3 * p[1] + 2 * p[2];
  Ax[1] = 2 * p[1] + 6 * p[2];
  
  float xTAx = p[1] * Ax[0] + p[2] * Ax[1];
  
  float bx = 2 * p[1] + (-8 * p[2]);

  return .5 * xTAx - bx;
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
  float expectedFret = -4.0;
  float *expectedP = new float[2];
  expectedP[1] = 3.0;
  CPPUNIT_ASSERT(TestPowellHelper(1, expectedFret, expectedP, testFunction1));
  

  delete [] expectedP;
  expectedP = new float[3];
  expectedFret = -10.0;
  expectedP[1] = 2.0;
  expectedP[2] = -2.0;
  CPPUNIT_ASSERT(TestPowellHelper(2, expectedFret, expectedP, testFunction2));
  
  delete [] expectedP;
}

bool
NRWrapperTest::TestPowellHelper( const int numberOfParameters, float expectedFret,
                                float * expectedP, float (*function)(float []) ) {
                                  
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
  bool gotExpectedValues = true;
  for (int i = 0; i < numberOfParameters; i++) {
    if ( !AreEqualWithinTolerance( p[i + 1], expectedP[i + 1], equalsTolerance) ) {
      gotExpectedValues = false;
    }
  }
  
  if (fret != expectedFret) gotExpectedValues = false;
  
  free_vector( p, 1, numberOfParameters );  
  free_matrix( xi, 1, numberOfParameters, 1, numberOfParameters );
  
  return gotExpectedValues;
}

MATRIX*
NRWrapperTest::ReconstructFromSVD( MATRIX* u, VECTOR* wVector, MATRIX* v ) {
  int numberOfRows = u->rows;
  int numberOfColumns = v->cols;
  
  MATRIX *w = MatrixAlloc( numberOfRows, numberOfColumns, MATRIX_REAL );
  
  for (int i=0; i < numberOfColumns; i++) {
    int index = i + numberOfColumns * i;
    w->data[index] = wVector->data[i];
  }
    
  MATRIX *vTranspose = MatrixTranspose( v, NULL );  
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
  
  MATRIX* x = MatrixRead( (char*) ( matrixFile.c_str() ) );

  if( numberOfRows == -1 ) {
    numberOfRows = x->rows;
  }
  
  if( numberOfColumns == -1 ) {
    numberOfColumns = x->cols;
  }
  
  MATRIX *u = MatrixCopy( x, NULL );
  VECTOR *w = RVectorAlloc( numberOfColumns, MATRIX_REAL );
  MATRIX *v = MatrixAlloc( numberOfColumns, numberOfColumns, MATRIX_REAL );
  
  int isError = svdcmp( u->rptr, numberOfRows, numberOfColumns, 
                        w->rptr[1], v->rptr );
      
  if( isError == NO_ERROR ) {
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
  int status = TestSVDcmpHelper( PASCAL_MATRIX );
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

  status = TestSVDcmpHelper( RANDOM_5_11 );
  CPPUNIT_ASSERT_EQUAL( MATRICES_ERROR, status );

  status = TestSVDcmpHelper( RANDOM_61_2, 61, 2 );
  CPPUNIT_ASSERT_EQUAL( MATRICES_EQUAL, status );

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
  const int numberOfBins = 10;
  int bins[numberOfBins];
  
  for(int i=0; i<numberOfBins; i++) {
    bins[i] = 0;
  }
  
  const int numberOfRuns = 1000;  
  long seed = -1L * (long)( abs( (int)time(NULL) ) );
  for(int i=0; i<numberOfRuns; i++) {
    float randomNumber = ran1( &seed );
 
    CPPUNIT_ASSERT( isBetween0And1( randomNumber ) );

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
  float lowerBound = 4.0;
  float upperBound = x2Mean*3;
  
  // if x2 is close to 0, it means that the distribution it very uniform -- bad!  
  CPPUNIT_ASSERT( x2 > lowerBound );

  // if x2 is way above the mean, then it's too non-uniform
  CPPUNIT_ASSERT( x2 < upperBound );
  
  // TODO: I think we need a test case to make sure that the numbers don't
  //       follow some period (try Unbinned KS testing).
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

