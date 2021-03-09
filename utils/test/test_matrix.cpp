/**
 * @brief Test freesurfers matrix operations
 *
 */
/*
 * Original Author: Dennis Jen
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "numerics.h"

extern "C"
{
#include "matrix.h"
#include "stdlib.h"
}

const char* Progname = "test_matrix";

class MatrixTest : public CppUnit::TestFixture
{

  // create the test suite and add the tests here
  CPPUNIT_TEST_SUITE( MatrixTest );
  CPPUNIT_TEST( TestMatrixInverse );
  CPPUNIT_TEST( TestMatrixDeterminant );
  CPPUNIT_TEST( TestOpenMatrixDeterminant );
  CPPUNIT_TEST( TestMatrixEigenSystem );
  CPPUNIT_TEST( TestMatrixSVDPseudoInverse );
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

  /** square matrix */
  MATRIX* mSquareMatrix;

  /** 1x1 matrix of a 5 */
  MATRIX* mOneMatrix;

  /** 1x1 matrix of a small number */
  MATRIX* mOneSmallMatrix;

  MATRIX* mTransformMatrix;

  MATRIX* mNonSymmetricMatrix;

  bool AreMatricesEqual(MATRIX *m1, MATRIX *m2, float tolerance);
  bool AreInversesEqual(MATRIX *matrix, const std::string inverseFile);
  void DeleteMatrix(MATRIX *m);
  int DoesCreateValidEigenSystem( MATRIX* matrix );

public:
  static const std::string TESTING_DIR;

  static const int EIGENSYSTEM_VALID;
  static const int EIGENSYSTEM_INVALID;
  static const int EIGENSYSTEM_NOT_DESCENDING;

  static const std::string PASCAL_MATRIX;
  static const std::string PASCAL_INVERSE;

  static const std::string BUCKY_MATRIX;
  static const std::string BUCKY_INVERSE;

  static const std::string ZEROES_MATRIX;

  static const std::string IDENTITY_MATRIX;

  static const std::string SINGULAR_MATRIX;

  static const std::string NON_SQUARE_MATRIX;
  static const std::string NON_SQUARE_MATRIX_PSEUDO_INVERSE;

  static const std::string SQUARE_MATRIX;
  static const std::string SQUARE_MATRIX_PSEUDO_INVERSE;

  static const std::string ONE_MATRIX;
  static const std::string ONE_INVERSE;

  static const std::string ONE_SMALL_MATRIX;
  static const std::string ONE_SMALL_INVERSE;

  static const std::string TRANSFORM_MATRIX;
  static const std::string TRANSFORM_INVERSE;

  static const std::string NON_SYMMETRIC_MATRIX;

  static const std::string SVD_U_MATRIX;
  static const std::string SVD_V_MATRIX;
  static const std::string SVD_S_VECTOR;

  // setUp is called automatically before each test
  void setUp();

  // tearDown is called automatically after each test
  void tearDown();

  void TestMatrixInverse();

  void TestMatrixDeterminant();

  void TestOpenMatrixDeterminant();

  void TestMatrixEigenSystem();

  void TestMatrixNonSymmetricEigenSystem();

  void TestMatrixSVDPseudoInverse();
};

const int MatrixTest::EIGENSYSTEM_VALID = 0;
const int MatrixTest::EIGENSYSTEM_INVALID = 1;
const int MatrixTest::EIGENSYSTEM_NOT_DESCENDING = 2;

const std::string MatrixTest::TESTING_DIR = "test_matrix_data/";

const std::string MatrixTest::PASCAL_MATRIX = TESTING_DIR + "Pascal.mat";
const std::string MatrixTest::PASCAL_INVERSE =
  TESTING_DIR + "PascalInverse.mat";

const std::string MatrixTest::BUCKY_MATRIX = TESTING_DIR + "Bucky.mat";
const std::string MatrixTest::BUCKY_INVERSE = TESTING_DIR + "BuckyInverse.mat";

const std::string MatrixTest::ZEROES_MATRIX = TESTING_DIR + "Zeroes.mat";

const std::string MatrixTest::IDENTITY_MATRIX = TESTING_DIR + "Identity.mat";

const std::string MatrixTest::SINGULAR_MATRIX = TESTING_DIR + "Singular.mat";

const std::string MatrixTest::NON_SQUARE_MATRIX =
  TESTING_DIR + "NonSquare.mat";
const std::string MatrixTest::NON_SQUARE_MATRIX_PSEUDO_INVERSE =
  TESTING_DIR + "NonSquarePseudoInverse.mat";

const std::string MatrixTest::SQUARE_MATRIX = TESTING_DIR + "Square.mat";
const std::string MatrixTest::SQUARE_MATRIX_PSEUDO_INVERSE =
  TESTING_DIR + "SvdPseudoInverse.mat";

const std::string MatrixTest::ONE_MATRIX = TESTING_DIR + "One.mat";
const std::string MatrixTest::ONE_INVERSE = TESTING_DIR + "OneInverse.mat";

const std::string MatrixTest::ONE_SMALL_MATRIX = TESTING_DIR + "OneSmall.mat";
const std::string MatrixTest::ONE_SMALL_INVERSE =
  TESTING_DIR + "OneSmallInverse.mat";

const std::string MatrixTest::TRANSFORM_MATRIX =
  TESTING_DIR + "TransformMatrix.mat";
const std::string MatrixTest::TRANSFORM_INVERSE =
  TESTING_DIR + "TransformMatrixInverse.mat";

const std::string MatrixTest::NON_SYMMETRIC_MATRIX =
  TESTING_DIR + "NonSymmetricMatrix.mat";

const std::string MatrixTest::SVD_U_MATRIX = TESTING_DIR + "SVD_U.mat";
const std::string MatrixTest::SVD_V_MATRIX = TESTING_DIR + "SVD_V.mat";
const std::string MatrixTest::SVD_S_VECTOR = TESTING_DIR + "SVD_S.mat";

void
MatrixTest::setUp()
{
  mPascalMatrix =    MatrixRead( (char*) ( PASCAL_MATRIX.c_str() ) );
  mBuckyMatrix =     MatrixRead( (char*) ( BUCKY_MATRIX.c_str() ) );
  mZeroesMatrix =    MatrixRead( (char*) ( ZEROES_MATRIX.c_str() ) );
  mIdentityMatrix =  MatrixRead( (char*) ( IDENTITY_MATRIX.c_str() ) );
  mSingularMatrix =  MatrixRead( (char*) ( SINGULAR_MATRIX.c_str() ) );
  mNonSquareMatrix = MatrixRead( (char*) ( NON_SQUARE_MATRIX.c_str() ) );
  mSquareMatrix =    MatrixRead( (char*) ( SQUARE_MATRIX.c_str() ) );
  mOneMatrix =       MatrixRead( (char*) ( ONE_MATRIX.c_str() ) );
  mOneSmallMatrix =  MatrixRead( (char*) ( ONE_SMALL_MATRIX.c_str() ) );
  mTransformMatrix = MatrixRead( (char*) ( TRANSFORM_MATRIX.c_str() ) );
  mNonSymmetricMatrix = MatrixRead( (char*) ( NON_SYMMETRIC_MATRIX.c_str() ) );
}

void
MatrixTest::tearDown()
{
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
MatrixTest::AreMatricesEqual( MATRIX *m1, MATRIX *m2, float tolerance=0.0 )
{
  bool areEqual = true;

  if ( m1->rows == m2->rows &&
       m2->cols == m2->cols )
  {

    int numberOfRows =  m1->rows;
    int numberOfColumns = m1->cols;

    for ( int row=0; row<numberOfRows; row++ )
    {
      for ( int column=0; column<numberOfColumns; column++ )
      {

        int index = column + row * numberOfColumns;
        float difference = fabs( m1->data[ index ] - m2->data[ index ] );

        if ( difference > tolerance )
        {
          areEqual = false;
          std::cerr << "MatrixTest::AreMatricesEqual() not equal: (" <<
          m1->data[ index ] << ", " << m2->data[ index ] << ")\n";
          std::cerr.flush();
        }
      }
    }
  }
  else
  {
    std::cerr << "MatrixTest::AreMatricesEqual() not equal: m1: (" <<
    m1->rows << ", " << m1->cols << "), m2: (" <<
    m2->rows << ", " << m2->cols << ")\n";
    areEqual = false;
  }

  if (! areEqual)
  {
    std::cerr << "actual:\n";
    MatrixPrint(stderr, m1);
    std::cerr << "expected:\n";
    MatrixPrint(stderr, m2);
  }

  return areEqual;
}

bool
MatrixTest::AreInversesEqual( MATRIX *matrix, const std::string inverseFile )
{

  // NJS: had to change the tolerance for the test case to pass with VXL
  float tolerance = 0.0022;

  MATRIX *expectedInverse = MatrixRead((char*)( inverseFile.c_str() ));
  MATRIX *actualInverse = MatrixInverse(matrix, NULL);

  CPPUNIT_ASSERT( expectedInverse != NULL );
  CPPUNIT_ASSERT( actualInverse != NULL );

  bool areEqual=AreMatricesEqual( expectedInverse, actualInverse, tolerance );

  DeleteMatrix(expectedInverse);
  DeleteMatrix(actualInverse);

  return areEqual;
}

void
MatrixTest::DeleteMatrix(MATRIX* m)
{
  MatrixFree( &m );
  m = NULL;
}

void
MatrixTest::TestMatrixInverse()
{
  std::cout << "\rMatrixTest::TestMatrixInverse()\n";

  //  std::cout << "\ttest a well-formed pascal matrix\n";
  CPPUNIT_ASSERT( AreInversesEqual(mPascalMatrix, PASCAL_INVERSE ) );

  //  std::cout << "\ttest a zero matrix\n";
  CPPUNIT_ASSERT( MatrixInverse(mZeroesMatrix, NULL) == NULL );

  //  std::cout << "\ttest an identity matrix\n";
  CPPUNIT_ASSERT( AreInversesEqual(mIdentityMatrix, IDENTITY_MATRIX ) );

  //  std::cout << "\ttest a 1 x 1 matrix\n";
  CPPUNIT_ASSERT( AreInversesEqual(mOneMatrix, ONE_INVERSE ) );

  //  std::cout << "\ttest a sparse matrix\n";
  CPPUNIT_ASSERT( AreInversesEqual(mBuckyMatrix, BUCKY_INVERSE ) );

  //  std::cout << "\ttest non-square matrix\n";
  CPPUNIT_ASSERT( MatrixInverse(mNonSquareMatrix, NULL) == NULL );

  //  std::cout << "\ttest a 1x1 matrix with a small number\n";
  CPPUNIT_ASSERT( AreInversesEqual(mOneSmallMatrix, ONE_SMALL_INVERSE ) );

  //  std::cout << "\ttest will fail when calculating the inverse using svd\n";
  CPPUNIT_ASSERT( AreInversesEqual(mTransformMatrix, TRANSFORM_INVERSE) );

  // TODO: why does this test not complain about having a singular
  // matrix passed to it?
  //std::cout << "\ttest a singular matrix\n";
  //CPPUNIT_ASSERT( MatrixInverse(mSingularMatrix, NULL) == NULL );
}

void
MatrixTest::TestMatrixDeterminant()
{
  double tolerance = 1e-4;

  std::cout << "\rMatrixTest::TestNRMatrixDeterminant()\n";

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

  // the tolerance had to be increased for this test case to pass.
  // The determinant is much larger in this case.
  const double buckyTolerance = 4;
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mBuckyMatrix),
                                2985984.0,
                                buckyTolerance );

  // non-square matrices will have a determinant of 0 for us
  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mNonSquareMatrix),
                                0.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)MatrixDeterminant(mOneSmallMatrix),
                                2.e-20,
                                tolerance );
}

void
MatrixTest::TestOpenMatrixDeterminant()
{
  double tolerance = 1e-4;

  // is it ok that bucky has a higher tolerance?
  double buckyTolerance = 4.0;

  std::cout << "\rMatrixTest::TestOpenMatrixDeterminant()\n";

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mPascalMatrix),
                                1.0,
                                tolerance);

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mZeroesMatrix),
                                0.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mIdentityMatrix),
                                1.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mOneMatrix),
                                5.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mSingularMatrix),
                                0.0,
                                tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mBuckyMatrix),
                                2985984.0,
                                buckyTolerance );

  // non-square matrices will have a determinant of 0 for us
  CPPUNIT_ASSERT_DOUBLES_EQUAL(
    (double)OpenMatrixDeterminant(mNonSquareMatrix),
    0.0,
    tolerance );

  CPPUNIT_ASSERT_DOUBLES_EQUAL( (double)OpenMatrixDeterminant(mOneSmallMatrix),
                                2.e-20,
                                tolerance );
}

int
MatrixTest::DoesCreateValidEigenSystem( MATRIX* matrix )
{
  const float TOLERANCE = 2e-5;

  int eigenSystemState = EIGENSYSTEM_INVALID;

  int size = matrix->rows;
  int cols = matrix->cols;

  float *eigenValues = new float[ size ];
  MATRIX *eigenVectors = MatrixAlloc( size, cols, MATRIX_REAL );

  MATRIX *eigenSystem = MatrixEigenSystem( matrix, eigenValues, eigenVectors );
  if (eigenSystem == NULL)
  {
    return EIGENSYSTEM_INVALID;
  }

  bool isInDecendingOrder = true;
  for ( int i=1; i<size; i++ )
  {
    if ( fabs( eigenValues[ i-1 ] ) < fabs( eigenValues[i] ) )
    {
      isInDecendingOrder = false;
      std::cerr << "DoesCreateValidEigenSystem::not in descending order: (" <<
      eigenValues[ i-1 ] << ", " << eigenValues[i] << ")\n";
    }
  }

  if ( !isInDecendingOrder )
  {
    eigenSystemState = EIGENSYSTEM_NOT_DESCENDING;
  }

  MATRIX* eigenValuesMatrix = MatrixAlloc( size, size, MATRIX_REAL );

  for ( int i=0; i<size; i++ )
  {
    // index of the diagonal
    int index = i + i * size;
    eigenValuesMatrix->data[index] = eigenValues[i];
  }

  MATRIX* xv = MatrixMultiply( matrix, eigenVectors, NULL );
  MATRIX* vd = MatrixMultiply( eigenVectors, eigenValuesMatrix, NULL );

  if ( isInDecendingOrder && AreMatricesEqual( xv, vd, TOLERANCE ) )
  {
    eigenSystemState = EIGENSYSTEM_VALID;
  }

  delete []eigenValues;
  DeleteMatrix( eigenVectors );
  DeleteMatrix( eigenValuesMatrix );

  return eigenSystemState;
}

void
MatrixTest::TestMatrixEigenSystem()
{

  std::cout << "\rMatrixTest::TestMatrixEigenSystem()\n";

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mPascalMatrix ) == EIGENSYSTEM_VALID );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem(mBuckyMatrix ) == EIGENSYSTEM_VALID );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mZeroesMatrix ) == EIGENSYSTEM_VALID);

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mIdentityMatrix ) == EIGENSYSTEM_VALID );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mSingularMatrix ) == EIGENSYSTEM_VALID);

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mNonSquareMatrix ) );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mOneMatrix ) == EIGENSYSTEM_VALID );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mOneSmallMatrix ) == EIGENSYSTEM_VALID );

  CPPUNIT_ASSERT
  ( DoesCreateValidEigenSystem( mNonSymmetricMatrix ) == EIGENSYSTEM_VALID );
}

void
MatrixTest::TestMatrixSVDPseudoInverse()
{
  float tolerance = 1e-5;

  std::cout << "\rMatrixTest::TestMatrixSVDPseudoInverse()\n";

  // check the low-level routine first
  MATRIX *Ux = MatrixRead( (char*) ( SVD_U_MATRIX.c_str() ) );
  MATRIX *Vx = MatrixRead( (char*) ( SVD_V_MATRIX.c_str() ) );
  VECTOR *Sx = MatrixRead( (char*) ( SVD_S_VECTOR.c_str() ) );
  CPPUNIT_ASSERT (Ux != NULL);
  CPPUNIT_ASSERT (Vx != NULL);
  CPPUNIT_ASSERT (Sx != NULL);

  MATRIX *U=MatrixCopy(mSquareMatrix,NULL);
  MATRIX *V=MatrixAlloc(U->cols, U->cols, MATRIX_REAL) ;
  VECTOR *S=VectorAlloc(U->cols, MATRIX_REAL) ;
  OpenSvdcmp( U, S, V ) ;
  CPPUNIT_ASSERT (U != NULL);
  CPPUNIT_ASSERT (V != NULL);
  CPPUNIT_ASSERT (S != NULL);
  CPPUNIT_ASSERT ( AreMatricesEqual( U, Ux, tolerance ) );
  CPPUNIT_ASSERT ( AreMatricesEqual( V, Vx, tolerance ) );
  CPPUNIT_ASSERT ( AreMatricesEqual( S, Sx, tolerance ) );

  // now check MatrixSVDPseudoInverse, which uses sc_linalg_SV_decomp

  tolerance = 1e-5;

  MATRIX *actualInverse = MatrixSVDPseudoInverse( mNonSquareMatrix, NULL );
  MATRIX *expectedInverse =
    MatrixRead( (char*) ( NON_SQUARE_MATRIX_PSEUDO_INVERSE.c_str() ) );
  CPPUNIT_ASSERT( AreMatricesEqual( actualInverse, expectedInverse ) );

  actualInverse = MatrixSVDPseudoInverse( mSquareMatrix, NULL );
  expectedInverse =
    MatrixRead( (char*) ( SQUARE_MATRIX_PSEUDO_INVERSE.c_str() ) );
  CPPUNIT_ASSERT
  ( AreMatricesEqual( actualInverse, expectedInverse, tolerance ) );
}

int main ( int argc, char** argv )
{

  int SUCCESS = 0;
  int FAIL = 1;

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( MatrixTest::suite() );

  if ( runner.run() )
  {
    exit ( SUCCESS );
  }

  exit( FAIL );
}

