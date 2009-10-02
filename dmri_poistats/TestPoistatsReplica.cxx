#include <iostream>

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "datamodel/data/PoistatsReplica.h" // dmri_poistats

/** This is needed by the freesurfer utils library */
char *Progname;

class TestPoistatsReplica : public CppUnit::TestFixture
{      
  // the test suite
  CPPUNIT_TEST_SUITE( TestPoistatsReplica );
  
    CPPUNIT_TEST( TestGenerateUnitSphereRandom );

    CPPUNIT_TEST( TestGenerateConstrainedRandomPoint3D );

    CPPUNIT_TEST( TestCopyPath );
    
    CPPUNIT_TEST( TestFoundBestPath );    
    
  CPPUNIT_TEST_SUITE_END();
  
private:
  PoistatsModel *m_PoistatsModel;
  
  PoistatsReplica *m_Replica;

public:
  // these are required. setUp is called before every test and tearDown after  
  void setUp();
  void tearDown();
  
  void TestGenerateUnitSphereRandom();

  void TestGenerateConstrainedRandomPoint3D();

  void TestCopyPath();
  
  void TestFoundBestPath();
  
};


void TestPoistatsReplica::setUp() {
  m_PoistatsModel = new PoistatsModel();
  m_Replica = new PoistatsReplica( m_PoistatsModel );
}

void TestPoistatsReplica::tearDown() {
  delete m_Replica;
  m_Replica = NULL;
  
  delete m_PoistatsModel;
  m_PoistatsModel = NULL;
}

void TestPoistatsReplica::TestGenerateUnitSphereRandom() {

  std::cerr << "TestGenerateUnitSphereRandom" << std::endl;

  const double expectedRandomUnitSphere[][3] = {
    {-0.325104, 0.64072, 0.695547},
    {0.565207, -0.44426, -0.695107},
    {0.279902, 0.792168, 0.542333},
    {0.336483, 0.886751, -0.316942},
    {0.9177, 0.353203, 0.181866},
    {-0.274809, -0.958605, 0.0745433},
    {-0.106031, 0.841314, -0.530046},
    {-0.823139, -0.383543, 0.418732},
    {0.550587, -0.83409, 0.0338662},
    {0.324134, -0.698436, 0.638064},
    {-0.303495, -0.186503, 0.934402},
    {-0.545803, -0.833996, -0.080929},
    {0.839996, 0.422594, 0.340325},
    {-0.603853, 0.489283, 0.629256},
    {-0.886406, -0.0640528, -0.458455},
    {0.122969, 0.682238, -0.720715},
    {0.595194, 0.790339, 0.145291},
    {-0.465842, -0.872374, -0.148172},
    {0.442118, 0.21351, 0.871175},
    {-0.070468, 0.0236782, -0.997233},
    {0.0827501, 0.414049, -0.906485},
    {0.19256, -0.506928, -0.840205},
    {-0.736059, 0.409707, 0.538848},
    {-0.991326, 0.131315, 0.00532057},
    {0.466488, 0.692202, 0.550677},
    {0.583513, -0.778768, -0.230289},
    {-0.188957, -0.176297, -0.96603},
    {-0.507248, 0.242276, -0.827044},
    {0.604357, 0.794717, -0.056362},
    {0.457852, 0.372752, -0.807111},
    {-0.535797, -0.83549, -0.121978},
    {-0.701336, 0.653151, 0.285521},
    {0.940734, -0.104309, 0.322705},
    {-0.769364, 0.373685, -0.518111},
    {0.238121, -0.234705, -0.94245},
    {0.0241008, -0.895855, -0.443693},
    {0.324319, 0.339063, 0.883093},
    {-0.490603, 0.252394, 0.83403},
    {0.954552, -0.0671617, -0.290379},
    {0.275043, 0.939448, 0.204424},
    {0.0846503, -0.56852, -0.818303},
    {-0.31682, -0.797988, -0.512679},
    {0.681602, -0.622167, -0.385132},
    {0.763594, 0.261624, -0.590319},
    {-0.420984, -0.882777, -0.20851},
    {-0.142314, -0.796207, -0.588049},
    {-0.239678, 0.509192, -0.826606},
    {0.533088, 0.768713, -0.35341},
    {0.0195064, 0.889132, 0.457235},
    {0.282546, -0.643905, 0.711024},
    {-0.816064, 0.527527, 0.236123},
    {-0.493313, -0.866754, -0.0733469},
    {0.199724, 0.268474, -0.942354},
    {0.0890264, 0.352722, -0.931483},
    {-0.68264, 0.331689, 0.651141},
    {0.146557, 0.13132, -0.980447},
    {0.411658, -0.156015, -0.897884},
    {0.818925, -0.326835, -0.471742},
    {0.135153, -0.549321, 0.824609},
    {0.265159, -0.808856, 0.524826},
    {-0.461971, -0.401313, -0.790905},
    {-0.0397155, -0.374678, -0.926304},
    {0.506233, 0.226471, 0.832129},
    {0.115004, -0.334185, 0.935465},
    {0.268745, 0.961123, 0.063399},
    {0.611666, 0.122459, -0.78158},
    {0.987461, 0.00192603, -0.157852},
    {-0.511267, 0.406938, 0.756973},
    {0.575657, 0.637076, -0.512594},
    {0.346071, -0.938054, 0.0169896},
    {-0.738647, 0.511228, 0.439371},
    {-0.10269, -0.571761, -0.813968},
    {0.527004, -0.633733, -0.56626},
    {-0.824698, 0.0229299, -0.565109},
    {-0.511788, 0.217579, 0.831103},
    {0.255733, 0.917494, -0.30464},
    {0.681492, -0.542539, -0.491142},
    {-0.665437, -0.727732, -0.166132},
    {-0.505926, -0.344531, -0.790783},
    {-0.317099, -0.336594, -0.886653},
    {-0.230165, 0.863757, 0.448272},
    {-0.114041, 0.816784, -0.56556},
    {-0.430175, -0.0299895, -0.902247},
    {0.952986, -0.299839, 0.0437597},
    {0.731178, -0.298405, -0.61346},
    {0.702433, -0.103449, 0.704192},
    {0.249527, 0.943325, 0.2188},
    {0.723046, -0.218025, 0.655491},
    {-0.837047, 0.48149, -0.259844},
    {-0.332394, -0.743191, -0.580673},
    {-0.948613, 0.295655, 0.112787},
    {-0.600378, 0.558226, -0.572652},
    {-0.0097279, 0.893123, -0.449708},
    {0.898567, -0.404973, -0.16904},
    {0.0469486, -0.336879, 0.940377},
    {0.3201, 0.848696, 0.421012},
    {0.222737, 0.966395, 0.128328},
    {-0.73637, 0.238522, -0.633141},
    {0.266053, -0.85653, 0.442235},
    {-0.734796, -0.607653, -0.301385}
  };

  const int numberOfPoints = 100;
  const int spatialDimensions = 3;
  
  const long seed = 10;
  m_PoistatsModel->SetRandomSeed( seed );

  PoistatsReplica::MatrixType resultRandomUnitSphere(
    numberOfPoints, spatialDimensions );
  
  m_Replica->GenerateUnitSphereRandom( &resultRandomUnitSphere );
      
  const double tolerance = 0.00005;
  for( unsigned int cRow=0; cRow<resultRandomUnitSphere.rows(); cRow++ ) {    
    for( unsigned int cColumn=0; cColumn<resultRandomUnitSphere.cols(); cColumn++ ) {
      double expected = expectedRandomUnitSphere[ cRow ][ cColumn ];
      double result = resultRandomUnitSphere[ cRow ][ cColumn ];
      
      CPPUNIT_ASSERT_DOUBLES_EQUAL( expected, result, tolerance );    
    }
  }
      
}

void TestPoistatsReplica::TestGenerateConstrainedRandomPoint3D() {

  std::cerr << "TestGenerateConstrainedRandomPoint3D" << std::endl;

  const int spatialDimensions = 3;    
  const int numberOfPoints = 3;    

  vnl_vector< double > currentPoint( spatialDimensions );
  currentPoint[ 0 ] = 1;
  currentPoint[ 1 ] = 2;
  currentPoint[ 2 ] = 3;
  
  PoistatsReplica::MatrixType possibleNewPoints( numberOfPoints, spatialDimensions );
  possibleNewPoints[ 0 ][ 0 ] = 1;
  possibleNewPoints[ 0 ][ 1 ] = 2;
  possibleNewPoints[ 0 ][ 2 ] = 3;

  possibleNewPoints[ 1 ][ 0 ] = 1;
  possibleNewPoints[ 1 ][ 1 ] = 1;
  possibleNewPoints[ 1 ][ 2 ] = 1;

  possibleNewPoints[ 2 ][ 0 ] = 3.5;
  possibleNewPoints[ 2 ][ 1 ] = 3.1;
  possibleNewPoints[ 2 ][ 2 ] = 2.0;
    
  const double sigma = 2.5;
  long seed = 10;
  m_PoistatsModel->SetRandomSeed( seed );

  PoistatsReplica::ArrayType newRandomPoint( spatialDimensions );

  m_Replica->GenerateConstrainedRandomPoint3D( currentPoint, 
    &possibleNewPoints, sigma, &newRandomPoint );

  const double tolerance = 0.5;

  // with a seed of 10, we expect the random point to be the one at 0
  int expectedNewPointIndex = 0;
  for( int cColumn=0; cColumn<spatialDimensions; cColumn++ ) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( newRandomPoint[cColumn],
      possibleNewPoints[ expectedNewPointIndex ][ cColumn ], tolerance );
  }
  
  seed = 9;
  m_PoistatsModel->SetRandomSeed( seed );
  m_Replica->GenerateConstrainedRandomPoint3D( currentPoint, 
    &possibleNewPoints, sigma, &newRandomPoint );
  expectedNewPointIndex = 1;
  for( int cColumn=0; cColumn<spatialDimensions; cColumn++ ) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( newRandomPoint[cColumn],
      possibleNewPoints[ expectedNewPointIndex ][ cColumn ], tolerance );
  }
  
}
  
void TestPoistatsReplica::TestCopyPath() {

  std::cerr << "TestCopyPath" << std::endl;
  
  const double points1[][ 3 ] =
    {
      {1, 2, 3},
      {4, 5, 6},
      {1, 5, -2},
      {5, 1, 0}
    };
  
  vnl_matrix< double > vnlPath1( *points1, 4, 3 );  
  PoistatsReplica::MatrixType path1( vnlPath1 );
  
  PoistatsReplica::MatrixType destinationPath( 4, 3 );
  PoistatsReplica::CopyPath( &path1, &destinationPath );
  
  for( unsigned int cRow=0; cRow<path1.rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<path1.cols(); cCol++ ) {
      
      const double source = path1[ cRow ][ cCol ];
      const double destination = destinationPath[ cRow ][ cCol ];
      
      CPPUNIT_ASSERT_EQUAL( source, destination );
    }
  }
  
  vnl_matrix< double > vnlPath2( *points1, 4, 3 );  
  PoistatsReplica::MatrixType path2( vnlPath2 );
  PoistatsReplica::CopyPath( &path2, &destinationPath );
  
  for( unsigned int cRow=0; cRow<path2.rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<path2.cols(); cCol++ ) {
      
      const double source = path2[ cRow ][ cCol ];
      const double destination = destinationPath[ cRow ][ cCol ];
      
      CPPUNIT_ASSERT_EQUAL( source, destination );
    }
  }
  
}

void TestPoistatsReplica::TestFoundBestPath() {

  std::cerr << "TestFoundBestPath" << std::endl;
  
  const double trialPathPoints[][ 3 ] =
    {
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6},
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2},
      {4, 5, 6}
    };
  vnl_matrix< double > vnlTrialPath( *trialPathPoints, 50, 3 );  
  PoistatsReplica::MatrixType trialPath( vnlTrialPath );
  
  m_Replica->SetCurrentTrialPath( &trialPath );
    
  PoistatsReplica::MatrixType zeroes( 4, 3 );
  zeroes.Fill( 0.0 );
  m_Replica->SetBasePath( &zeroes );
  
  zeroes = PoistatsReplica::MatrixType( 50, 3 );
  zeroes.Fill( 0.0 );
  m_Replica->SetBestTrialPath( &zeroes );

  const double basePoints[][ 3 ] =
    {
      {4, 5, 6},
      {1, 2, 3},
      {5, 1, 0},
      {1, 5, -2}
    };
  vnl_matrix< double > vnlBasePath( *basePoints, 4, 3 );  
  PoistatsReplica::MatrixType basePath( vnlBasePath );
  m_Replica->FoundBestPath( &basePath );
  
  for( unsigned int cRow=0; cRow<basePath.rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<basePath.cols(); cCol++ ) {
      
      PoistatsReplica::MatrixPointer actualBasePath = m_Replica->GetBasePath();
      const double expected = basePath[ cRow ][ cCol ];
      const double actual = ( *actualBasePath )[ cRow ][ cCol ];

      CPPUNIT_ASSERT_EQUAL( expected, actual );
      
    }
  }
  
  for( unsigned int cRow=0; cRow<trialPath.rows(); cRow++ ) {
    for( unsigned int cCol=0; cCol<trialPath.cols(); cCol++ ) {
      
      PoistatsReplica::MatrixPointer actualBasePath = 
        m_Replica->GetBestTrialPath();
      const double expected = trialPath[ cRow ][ cCol ];
      const double actual = ( *actualBasePath )[ cRow ][ cCol ];

      CPPUNIT_ASSERT_EQUAL( expected, actual );
      
    }
  }
  

}

int main ( int argc, char** argv ) {

  // this is needed by the freesurfer utils library
  Progname = (char*)"TestPoistatsReplica";

  const int SUCCESS = 0;
  const int FAIL = 1;

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( TestPoistatsReplica::suite() );

  if ( runner.run() ) {
    exit ( SUCCESS );
  }

  exit( FAIL );
}
