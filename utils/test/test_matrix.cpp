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
  #include "matrix.h"
}

using namespace std;
  
class MatrixTest : public CppUnit::TestFixture {

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( MatrixTest );  
    CPPUNIT_TEST( TestMatrixInverse );
    CPPUNIT_TEST( TestMatrixDeterminant );
    CPPUNIT_TEST( TestMatrixEigenSystem );
  CPPUNIT_TEST_SUITE_END();

private:
  /** a well-formed matrix */
  MATRIX* mPascalMatrix;
  
  /** matrix full of zeroes*/
  MATRIX* mZeroesMatrix;
  
  /** identity matrix */
  MATRIX* mIdentityMatrix;
  
  /** sparse adjacency matrix via the buckminster fuller's geodesic dome */
  MATRIX* mBuckyMatrix;
  
  /** singular matrix */
  MATRIX* mSingularMatrix;

  /** non-square matrix */  
  MATRIX* mNonSquareMatrix;

  /** 1x1 matrix of a 5 */  
  MATRIX* mOneMatrix;
  
  /** 1x1 matrix of a small number */
  MATRIX* mOneSmallMatrix;
  
  bool AreMatricesEqual(MATRIX *m1, MATRIX *m2, float tolerance);
  bool AreInversesEqual(MATRIX *matrix, const string inverseFile);
  void DeleteMatrix(MATRIX *m);
  bool DoesCreateValidEigenSystem( MATRIX* matrix );
  
public:  
  static const string TESTING_DIR;

  static const string PASCAL_MATRIX;
  static const string PASCAL_INVERSE;

  static const string BUCKY_MATRIX;
  static const string BUCKY_INVERSE;

  static const string ZEROES_MATRIX;

  static const string IDENTITY_MATRIX;
  
  static const string SINGULAR_MATRIX;
  
  static const string NON_SQUARE_MATRIX;

  static const string ONE_MATRIX;
  static const string ONE_INVERSE;

  static const string ONE_SMALL_MATRIX;
  static const string ONE_SMALL_INVERSE;
  
  // setUp is called automatically before each test
  void setUp();
  
  // tearDown is called automatically after each test
  void tearDown();
  
  void TestMatrixInverse();
  void TestMatrixDeterminant();
  void TestMatrixEigenSystem();

};

const string MatrixTest::TESTING_DIR = "test_matrix_data/";

const string MatrixTest::PASCAL_MATRIX = TESTING_DIR + "Pascal.mat";
const string MatrixTest::PASCAL_INVERSE = TESTING_DIR + "PascalInverse.mat";

const string MatrixTest::BUCKY_MATRIX = TESTING_DIR + "Bucky.mat";
const string MatrixTest::BUCKY_INVERSE = TESTING_DIR + "BuckyInverse.mat";

const string MatrixTest::ZEROES_MATRIX = TESTING_DIR + "Zeroes.mat";

const string MatrixTest::IDENTITY_MATRIX = TESTING_DIR + "Identity.mat";

const string MatrixTest::SINGULAR_MATRIX = TESTING_DIR + "Singular.mat";

const string MatrixTest::NON_SQUARE_MATRIX = TESTING_DIR + "NonSquare.mat";

const string MatrixTest::ONE_MATRIX = TESTING_DIR + "One.mat";
const string MatrixTest::ONE_INVERSE = TESTING_DIR + "OneInverse.mat";

const string MatrixTest::ONE_SMALL_MATRIX = TESTING_DIR + "OneSmall.mat";
const string MatrixTest::ONE_SMALL_INVERSE = TESTING_DIR + "OneSmallInverse.mat";

void
MatrixTest::setUp() {
  mPascalMatrix =    MatrixRead( (char*) ( PASCAL_MATRIX.c_str() ) );
  mBuckyMatrix =     MatrixRead( (char*) ( BUCKY_MATRIX.c_str() ) );
  mZeroesMatrix =    MatrixRead( (char*) ( ZEROES_MATRIX.c_str() ) );
  mIdentityMatrix =  MatrixRead( (char*) ( IDENTITY_MATRIX.c_str() ) );
  mSingularMatrix =  MatrixRead( (char*) ( SINGULAR_MATRIX.c_str() ) );
  mNonSquareMatrix = MatrixRead( (char*) ( NON_SQUARE_MATRIX.c_str() ) );
  mOneMatrix =       MatrixRead( (char*) ( ONE_MATRIX.c_str() ) );
  mOneSmallMatrix =  MatrixRead( (char*) ( ONE_SMALL_MATRIX.c_str() ) );
}

void
MatrixTest::tearDown() {
  MatrixFree( &mPascalMatrix );
  MatrixFree( &mBuckyMatrix );  
  MatrixFree( &mZeroesMatrix );
  MatrixFree( &mIdentityMatrix );
  MatrixFree( &mSingularMatrix );
  MatrixFree( &mNonSquareMatrix );
  MatrixFree( &mOneMatrix );
  MatrixFree( &mOneSmallMatrix );
}

bool 
MatrixTest::AreMatricesEqual( MATRIX *m1, MATRIX *m2, float tolerance=0.0 ) {
  bool areEqual = true;

  if( m1->rows == m2->rows &&
      m2->cols == m2->cols ) {  

    int numberOfRows =  m1->rows;
    int numberOfColumns = m1->cols;

    for( int row=0; row<numberOfRows; row++ ) {
      for( int column=0; column<numberOfColumns; column++ ) {    

        int index = column + row * numberOfColumns;        
        float difference = fabs( m1->data[ index ] - m2->data[ index ] );
                
        if( difference > tolerance ) {
          cerr << "diff: " << difference << endl;
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
MatrixTest::AreInversesEqual(MATRIX *matrix, const string inverseFile) {
  MATRIX *expectedInverse = MatrixRead((char*)( inverseFile.c_str() ));
  MATRIX *actualInverse = MatrixInverse(matrix, NULL);

  CPPUNIT_ASSERT( expectedInverse != NULL );
  CPPUNIT_ASSERT( actualInverse != NULL );
  
  bool areEqual = AreMatricesEqual( expectedInverse, actualInverse );
  
  DeleteMatrix(expectedInverse);  
  DeleteMatrix(actualInverse);  
  
  return areEqual;
}

void
MatrixTest::DeleteMatrix(MATRIX* m) {
  MatrixFree( &m );
  m = NULL;
}

void
MatrixTest::TestMatrixInverse() {
  // test a well-formed pascal matrix
  CPPUNIT_ASSERT( AreInversesEqual(mPascalMatrix, PASCAL_INVERSE) );

  // test a zero matrix
  CPPUNIT_ASSERT( MatrixInverse(mZeroesMatrix, NULL) == NULL );
  
  // test an identity matrix
  CPPUNIT_ASSERT( AreInversesEqual(mIdentityMatrix, IDENTITY_MATRIX) );
    
  // test a 1 x 1 matrix
  CPPUNIT_ASSERT( AreInversesEqual(mOneMatrix, ONE_INVERSE) );

  // test a singular matrix
  CPPUNIT_ASSERT( MatrixInverse(mSingularMatrix, NULL) == NULL );

  // test a sparse matrix
  CPPUNIT_ASSERT( AreInversesEqual(mBuckyMatrix, BUCKY_INVERSE) );  

  // test non-square matrix
  CPPUNIT_ASSERT( MatrixInverse(mNonSquareMatrix, NULL) == NULL );

  // test a 1x1 matrix with a small number
  CPPUNIT_ASSERT( AreInversesEqual(mOneSmallMatrix, ONE_SMALL_INVERSE) );  
}

void
MatrixTest::TestMatrixDeterminant() {
  double tolerance = 1e-6;
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mPascalMatrix),
                                1.0,
                                tolerance);
                        
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mZeroesMatrix),
                                0.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mIdentityMatrix), 
                                1.0,
                                tolerance );
  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mOneMatrix), 
                                5.0, 
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mSingularMatrix),
                                0.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mBuckyMatrix), 
                                2985984.0,
                                tolerance );

  // non-square matrices will have a determinant of 0 for us
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mNonSquareMatrix),
                                0.0,
                                tolerance );
                                
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mOneSmallMatrix), 
                                2.e-20, 
                                tolerance );
                                
}

bool
MatrixTest::DoesCreateValidEigenSystem( MATRIX* matrix ) {
  const float TOLERANCE = 1e-5;
  
  bool isValid = false;
  
  int size = matrix->rows;
  
  float *eigenValues = new float[ size ];
  MATRIX *eigenVectors = MatrixAlloc( size, size, MATRIX_REAL );

  MatrixEigenSystem( matrix, eigenValues, eigenVectors );

  MATRIX* eigenValuesMatrix = MatrixAlloc( size, size, MATRIX_REAL );
  
  for( int i=0; i<size; i++ ) {
    // index of the diagonal
    int index = i + i * size;
    eigenValuesMatrix->data[index] = eigenValues[i];
  }
  
  MATRIX* xv = MatrixMultiply( matrix, eigenVectors, NULL );
  MATRIX* vd = MatrixMultiply( eigenVectors, eigenValuesMatrix, NULL );

  isValid = AreMatricesEqual( xv, vd, TOLERANCE );

  delete []eigenValues;
  DeleteMatrix( eigenVectors );
  DeleteMatrix( eigenValuesMatrix );
  
  return isValid;
}

void
MatrixTest::TestMatrixEigenSystem() {
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mPascalMatrix ) );
  
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mBuckyMatrix ) );
  
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mZeroesMatrix ) );
  
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mIdentityMatrix ) );
  
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mSingularMatrix ) );
  
  // TODO: should this gracefully exit with a non-square matrix rather than seg fault?
//  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mNonSquareMatrix ) );

  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mOneMatrix ) );
  
  CPPUNIT_ASSERT( DoesCreateValidEigenSystem( mOneSmallMatrix ) );
}

int main ( int argc, char** argv ) {
  
  int SUCCESS = 0;
  int FAIL = 1;
  
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( MatrixTest::suite() );
  
  if( runner.run() ) {
    exit ( SUCCESS );
  }
  
  exit( FAIL );
}

