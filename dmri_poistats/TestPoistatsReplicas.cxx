#include <iostream>

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "datamodel/data/PoistatsReplicas.h" // dmri_poistats

/** This is needed by the freesurfer utils library */
char *Progname;

class TestPoistatsReplicas : public CppUnit::TestFixture
{      
  // the test suite
  CPPUNIT_TEST_SUITE( TestPoistatsReplicas );

    CPPUNIT_TEST( TestRethreadPath );
  
    CPPUNIT_TEST( TestCopyCurrentToPreviousEnergy );

    CPPUNIT_TEST( TestSortArray );

    CPPUNIT_TEST( TestGetMinimumCurrentEnergy );
    
    CPPUNIT_TEST( TestGetCurrentMeanOfEnergies );

    CPPUNIT_TEST( TestGetNormalizedMeanCurrentPreviousEnergiesDifference );

    CPPUNIT_TEST( TestGetRandomSortedFirstSecondReplicas );

    CPPUNIT_TEST( TestCoolTemperatures );

    CPPUNIT_TEST( TestSetInitialPoints );

    CPPUNIT_TEST( TestGetBestTrialPaths );

  CPPUNIT_TEST_SUITE_END();
  
private:
  PoistatsModel *m_PoistatsModel;
  
  PoistatsReplicas *m_Replicas;

public:
  // these are required. setUp is called before every test and tearDown after  
  void setUp();
  void tearDown();

  void TestRethreadPath();  
    
  void TestCopyCurrentToPreviousEnergy();

  void TestSortArray();
  
  void TestGetMinimumCurrentEnergy();
  
  void TestGetCurrentMeanOfEnergies();
  
  void TestGetNormalizedMeanCurrentPreviousEnergiesDifference();
  
  void TestGetRandomSortedFirstSecondReplicas();
  
  void TestCoolTemperatures();
  
  void TestSetInitialPoints();

  void TestGetBestTrialPaths();

};

void TestPoistatsReplicas::setUp() {

  m_PoistatsModel = new PoistatsModel();
  
  m_Replicas = new PoistatsReplicas( m_PoistatsModel, 5 );
  
  m_PoistatsModel->SetNumberOfControlPoints( 2 );
  m_Replicas->SetNumberOfSteps( 50 );
  
  for( int cReplica=0; cReplica<m_Replicas->GetNumberOfReplicas(); cReplica++ ){
    m_Replicas->SetCurrentMeanEnergy( cReplica, cReplica );
    m_Replicas->SetTemperature( cReplica, cReplica * 2 );
  }
  m_Replicas->FillPreviousMeanEnergies( 11 );  

}

void TestPoistatsReplicas::tearDown() {

  delete m_PoistatsModel;
  m_PoistatsModel = NULL;
  
  delete m_Replicas;
  m_Replicas = NULL;

}

void TestPoistatsReplicas::TestRethreadPath() {
  
  std::cerr << "TestRethreadPath" << std::endl;

  m_PoistatsModel->SetNumberOfControlPoints( 1 );
  
  typedef itk::Array2D< double > MatrixType;

  const double actualRethreadShortPath[][3] = {
    {65, 46, 37},
    {68.9159, 60.3864, 39.0484},
    {64, 73, 36}
  };
  
  MatrixType originalPath( 2, 3 );
  
  originalPath[ 0 ][ 0 ] = 65;
  originalPath[ 0 ][ 1 ] = 46;
  originalPath[ 0 ][ 2 ] = 37;

  originalPath[ 1 ][ 0 ] = 64;
  originalPath[ 1 ][ 1 ] = 73;
  originalPath[ 1 ][ 2 ] = 36;

  int nSamples = 3;

  MatrixType *resultRethreadedPath = 
    m_PoistatsModel->RethreadPath( &originalPath, nSamples );
    
  // compare each element with what it should be now
  const double tolerance = 0.005;
  for( int cSample = 0; cSample < nSamples; cSample++ ) {
    for( int cColumn = 0; cColumn < 3; cColumn++ ) {
    
      const double actual = ( *resultRethreadedPath )[ cSample ][ cColumn ];
      const double result = actualRethreadShortPath[ cSample ][ cColumn ];

      CPPUNIT_ASSERT_DOUBLES_EQUAL( actual, result, tolerance );    
    }
  }

  delete resultRethreadedPath;  
  resultRethreadedPath = NULL;

  m_PoistatsModel->SetNumberOfControlPoints( 2 );

  const double actualRethreadPath[][3] = {
    {0,         0,         0},
    {0.0335,    0.0335,    0.0335},
    {0.0670,    0.0670,    0.0670},
    {0.1004,    0.1004,    0.1006},
    {0.1338,    0.1338,    0.1342},
    {0.1672,    0.1672,    0.1680},
    {0.2005,    0.2005,    0.2019},
    {0.2337,    0.2337,    0.2359},
    {0.2668,    0.2668,    0.2701},
    {0.2999,    0.2999,    0.3045},
    {0.3327,    0.3327,    0.3392},
    {0.3655,    0.3655,    0.3742},
    {0.3981,    0.3981,    0.4094},
    {0.4305,    0.4305,    0.4450},
    {0.4626,    0.4626,    0.4810},
    {0.4946,    0.4946,    0.5174},
    {0.5263,    0.5263,    0.5542},
    {0.5577,    0.5577,    0.5916},
    {0.5887,    0.5887,    0.6295},
    {0.6194,    0.6194,    0.6680},
    {0.6497,    0.6497,    0.7071},
    {0.6794,    0.6794,    0.7470},
    {0.7087,    0.7087,    0.7877},
    {0.7374,    0.7374,    0.8292},
    {0.7654,    0.7654,    0.8716},
    {0.7926,    0.7926,    0.9150},
    {0.8189,    0.8189,    0.9595},
    {0.8442,    0.8442,    1.0051},
    {0.8684,    0.8684,    1.0520},
    {0.8913,    0.8913,    1.1001},
    {0.9126,    0.9126,    1.1497},
    {0.9323,    0.9323,    1.2006},
    {0.9499,    0.9499,    1.2530},
    {0.9653,    0.9653,    1.3067},
    {0.9782,    0.9782,    1.3618},
    {0.9883,    0.9883,    1.4180},
    {0.9954,    0.9954,    1.4751},
    {0.9993,    0.9993,    1.5329},
    {0.9998,    0.9998,    1.5909},
    {0.9970,    0.9970,    1.6487},
    {0.9909,    0.9909,    1.7061},
    {0.9817,    0.9817,    1.7626},
    {0.9696,    0.9696,    1.8180},
    {0.9549,    0.9549,    1.8722},
    {0.9379,    0.9379,    1.9250},
    {0.9189,    0.9189,    1.9764},
    {0.8980,    0.8980,    2.0263},
    {0.8756,    0.8756,    2.0749},
    {0.8518,    0.8518,    2.1222},
    {0.8268,    0.8268,    2.1682},
    {0.8008,    0.8008,    2.2130},
    {0.7738,    0.7738,    2.2567},
    {0.7461,    0.7461,    2.2994},
    {0.7176,    0.7176,    2.3412},
    {0.6885,    0.6885,    2.3821},
    {0.6589,    0.6589,    2.4222},
    {0.6288,    0.6288,    2.4616},
    {0.5982,    0.5982,    2.5003},
    {0.5673,    0.5673,    2.5384},
    {0.5360,    0.5360,    2.5759},
    {0.5044,    0.5044,    2.6129},
    {0.4725,    0.4725,    2.6494},
    {0.4404,    0.4404,    2.6855},
    {0.4081,    0.4081,    2.7212},
    {0.3756,    0.3756,    2.7566},
    {0.3429,    0.3429,    2.7916},
    {0.3100,    0.3100,    2.8264},
    {0.2770,    0.2770,    2.8609},
    {0.2439,    0.2439,    2.8952},
    {0.2108,    0.2108,    2.9293},
    {0.1775,    0.1775,    2.9632},
    {0.1441,    0.1441,    2.9970},
    {0.1107,    0.1107,    3.0306},
    {0.0773,    0.0773,    3.0642},
    {0.0438,    0.0438,    3.0978},
    {0.0103,    0.0103,    3.1313},
    {-0.0232,   -0.0232,    3.1648},
    {-0.0566,   -0.0566,    3.1983},
    {-0.0901,   -0.0901,    3.2318},
    {-0.1235,   -0.1235,    3.2654},
    {-0.1569,   -0.1569,    3.2992},
    {-0.1902,   -0.1902,    3.3330},
    {-0.2235,   -0.2235,    3.3670},
    {-0.2566,   -0.2566,    3.4011},
    {-0.2897,   -0.2897,    3.4355},
    {-0.3226,   -0.3226,    3.4701},
    {-0.3554,   -0.3554,    3.5049},
    {-0.3880,   -0.3880,    3.5401},
    {-0.4205,   -0.4205,    3.5756},
    {-0.4527,   -0.4527,    3.6114},
    {-0.4848,   -0.4848,    3.6477},
    {-0.5165,   -0.5165,    3.6844},
    {-0.5480,   -0.5480,    3.7216},
    {-0.5792,   -0.5792,    3.7593},
    {-0.6100,   -0.6100,    3.7976},
    {-0.6404,   -0.6404,    3.8366},
    {-0.6703,   -0.6703,    3.8762},
    {-0.6997,   -0.6997,    3.9166},
    {-0.7286,   -0.7286,    3.9579},
    {-0.7568,   -0.7568,    4.0000}
  };
  
  const int nInputMatrixRows = 41;
  const int nInputMatrixColumns = 3;
          
  originalPath = MatrixType( nInputMatrixRows, nInputMatrixColumns);
  
  // set up our 3d path
  const double step = 0.1;
  double t = 0.0;
  for( int cRow = 0; cRow < nInputMatrixRows; cRow++ ) {
    
    double y = sin( t );
    int cColumn = 0;

    originalPath[ cRow ][ cColumn++ ] = y;      
    originalPath[ cRow ][ cColumn++ ] = y;      
    originalPath[ cRow ][ cColumn ] = t;
    
    t = t + step;
  }
  
  nSamples = 100;

  resultRethreadedPath = 
    m_PoistatsModel->RethreadPath( &originalPath, nSamples );  

  CPPUNIT_ASSERT_EQUAL( 
    nSamples, static_cast< int >( resultRethreadedPath->rows() ) );
  CPPUNIT_ASSERT_EQUAL( 
    nInputMatrixColumns, static_cast< int >( resultRethreadedPath->cols() ) );
   
  // compare each element with what it should be now
  for( int cSample = 0; cSample < nSamples; cSample++ ) {
    for( int cColumn = 0; cColumn < nInputMatrixColumns; cColumn++ ) {
    
      const double actual = ( *resultRethreadedPath )[ cSample ][ cColumn ];
      const double result = actualRethreadPath[ cSample ][ cColumn ];

      CPPUNIT_ASSERT_DOUBLES_EQUAL( actual, result, tolerance );    
    }
  }
     
  delete resultRethreadedPath;  
  resultRethreadedPath = NULL;
    
  originalPath = MatrixType( 4, 3 );
  originalPath[ 0 ][ 0 ] = 14.985982468709299;
  originalPath[ 0 ][ 1 ] = 63.3912989661489;
  originalPath[ 0 ][ 2 ] = 55.26209683302739;
  
  originalPath[ 1 ][ 0 ] = 13.679996356326065;
  originalPath[ 1 ][ 1 ] = 67.29527869522356;
  originalPath[ 1 ][ 2 ] = 46.53398359385062;

  originalPath[ 2 ][ 0 ] = 18.14736944328925;
  originalPath[ 2 ][ 1 ] = 71.9866493308684;
  originalPath[ 2 ][ 2 ] = 56.99132269186519;
  
  originalPath[ 3 ][ 0 ] = 23.518503643248224;
  originalPath[ 3 ][ 1 ] = 77.32140716429525;
  originalPath[ 3 ][ 2 ] = 55.944703364353195;  

  delete resultRethreadedPath;  
  resultRethreadedPath = NULL;
  
}

void TestPoistatsReplicas::TestCopyCurrentToPreviousEnergy() {

  std::cerr << "TestCopyCurrentToPreviousEnergy" << std::endl;  
  
  for( int cReplica=0; cReplica<m_Replicas->GetNumberOfReplicas(); cReplica++ ){
    const double current = m_Replicas->GetCurrentMeanEnergy( cReplica );
    const double previous = m_Replicas->GetPreviousMeanEnergy( cReplica );
    CPPUNIT_ASSERT( current != previous );
  }

  m_Replicas->CopyCurrentToPreviousEnergies();
  
  for( int cReplica=0; cReplica<m_Replicas->GetNumberOfReplicas(); cReplica++ ){
    const double current = m_Replicas->GetCurrentMeanEnergy( cReplica );
    const double previous = m_Replicas->GetPreviousMeanEnergy( cReplica );
    CPPUNIT_ASSERT_EQUAL( current, previous );
  }
    
}

void TestPoistatsReplicas::TestSortArray() {

  std::cerr << "TestSortArray" << std::endl;  

  const int arraySize = 5;
  
  int expectedIndices[] = { 0, 2, 4, 1, 3 };

  itk::Array< double > unsortedArray( arraySize );
  unsortedArray[ 0 ] = 1;
  unsortedArray[ 1 ] = 8;
  unsortedArray[ 2 ] = 2;
  unsortedArray[ 3 ] = 10;
  unsortedArray[ 4 ] = 5;
  
  itk::Array< int > sortedIndices( arraySize );
  
  PoistatsReplicas::SortArray( &unsortedArray, &sortedIndices );
  
  for( int cIndex=0; cIndex<arraySize; cIndex++ ) {
      int expected = expectedIndices[ cIndex ];
      int result = sortedIndices[ cIndex ];
      
      CPPUNIT_ASSERT_EQUAL( expected, result );    
  }
  
}

void TestPoistatsReplicas::TestGetMinimumCurrentEnergy() {
 
  std::cerr << "TestGetMinimumCurrentEnergy" << std::endl;  

  const double actual = m_Replicas->GetMinimumCurrentEnergy();
  const double expected = 0;
  
  CPPUNIT_ASSERT_EQUAL( expected, actual );
}

void TestPoistatsReplicas::TestGetCurrentMeanOfEnergies() {
  
  std::cerr << "TestGetCurrentMeanOfEnergies" << std::endl;  

  const double actual = m_Replicas->GetCurrentMeanOfEnergies();  
  const double expected = 2.0;
  
  CPPUNIT_ASSERT_EQUAL( expected, actual );  
  
}

void 
TestPoistatsReplicas::TestGetNormalizedMeanCurrentPreviousEnergiesDifference() {

  std::cerr << "TestGetNormalizedMeanCurrentPreviousEnergiesDifference" << std::endl;  
  
  const double actual = 
    m_Replicas->GetNormalizedMeanCurrentPreviousEnergiesDifference();  
  const double expected = -4.5;
  
  CPPUNIT_ASSERT_EQUAL( expected, actual );  
  
}

void TestPoistatsReplicas::TestGetRandomSortedFirstSecondReplicas() {

  std::cerr << "TestGetRandomSortedFirstSecondReplicas" << std::endl;  
  
  int firstActual;
  const int firstExpected = 2;

  int secondActual;
  const int secondExpected = 3;

  m_PoistatsModel->SetRandomSeed( 0 );

  m_Replicas->GetRandomSortedFirstSecondReplicas( 
    firstActual, secondActual );
  
  CPPUNIT_ASSERT_EQUAL( firstExpected, firstActual );  
  CPPUNIT_ASSERT_EQUAL( secondExpected, secondActual );    
  
}

void TestPoistatsReplicas::TestCoolTemperatures() {
  std::cerr << "TestCoolTemperatures" << std::endl;  

  m_Replicas->CoolTemperatures( 0.5 );
  
  for( int cReplica=0; cReplica<m_Replicas->GetNumberOfReplicas(); cReplica++ ){
    const double actualTemperature = m_Replicas->GetTemperature( cReplica );
    const double expectedTemperature = cReplica;
    CPPUNIT_ASSERT_EQUAL( expectedTemperature, actualTemperature );
  }
  
}

void TestPoistatsReplicas::TestSetInitialPoints() {
  
  std::cerr << "TestSetInitialPoints" << std::endl;    
  
  PoistatsModel::MatrixType initialPoints( 4, 3 );
  
  for( unsigned int cRow=0; cRow<initialPoints.rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<initialPoints.cols(); cCol++ ) {
      initialPoints[ cRow ][ cCol ] = cRow + cCol;
    }
  }
  
  m_Replicas->SetInitialPoints( &initialPoints );
  
  PoistatsModel::MatrixPointer expectedBasePath = 
    m_PoistatsModel->RethreadPath( &initialPoints, 4 );
    
  PoistatsModel::MatrixPointer expectedPreviousPath =
    m_PoistatsModel->RethreadPath( expectedBasePath, 50 );
      
  for( int cReplica=0; cReplica<m_Replicas->GetNumberOfReplicas(); cReplica++ ){

    PoistatsModel::MatrixPointer actualBasePath = 
      m_Replicas->GetBasePath( cReplica );

    for( unsigned int cRow=0; cRow<expectedBasePath->rows(); cRow++ ) {
      for( unsigned int cCol=0; cCol<expectedBasePath->cols(); cCol++ ) {
        
        const double expectedBase = ( *expectedBasePath )[ cRow ][ cCol ];
        const double actualBase = ( *actualBasePath )[ cRow ][ cCol ];
        
        CPPUNIT_ASSERT_EQUAL( expectedBase, actualBase );
        
      }
    }

    PoistatsReplica::MatrixPointer actualPreviousPath = 
      m_Replicas->GetPreviousTrialPath( cReplica );

    PoistatsReplica::MatrixPointer actualBestPath = 
      m_Replicas->GetBestTrialPath( cReplica );

    PoistatsReplica::MatrixPointer actualCurrentPath = 
      m_Replicas->GetCurrentTrialPath( cReplica );
    
    for( unsigned int cRow=0; cRow<expectedPreviousPath->rows(); cRow++ ) {
      for( unsigned int cCol=0; cCol<expectedPreviousPath->cols(); cCol++ ) {
        
        const double expectedPrevious = ( *expectedPreviousPath )[ cRow ][ cCol ];
        const double actualPrevious = ( *actualPreviousPath )[ cRow ][ cCol ];
        
        CPPUNIT_ASSERT_EQUAL( expectedPrevious, actualPrevious );

        const double actualBest = ( *actualBestPath )[ cRow ][ cCol ];
        CPPUNIT_ASSERT_EQUAL( 0.0, actualBest );

        const double actualCurrent = ( *actualCurrentPath )[ cRow ][ cCol ];
        CPPUNIT_ASSERT_EQUAL( 0.0, actualCurrent );

      }
    }
    
  }
  
  delete expectedPreviousPath;
  delete expectedBasePath;
}

void TestPoistatsReplicas::TestGetBestTrialPaths() {
  
  std::cerr << "TestGetBestTrialPaths" << std::endl;    
  
  PoistatsReplicas::MatrixListType bestPaths= m_Replicas->GetBestTrialPaths();
  
  const int expected = m_Replicas->GetNumberOfReplicas();
  const int actual = bestPaths.size();
  CPPUNIT_ASSERT_EQUAL( expected, actual );
  
}

int main ( int argc, char** argv ) {

  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  const int SUCCESS = 0;
  const int FAIL = 1;

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( TestPoistatsReplicas::suite() );

  if ( runner.run() ) {
    exit ( SUCCESS );
  }

  exit( FAIL );
}

