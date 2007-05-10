#include "PoistatsModel.h"

#include <iostream>
#include <cmath>

PoistatsModel::PoistatsModel() {

  std::srand( ( unsigned ) time( 0 ) );
  const long seed = std::rand();

  this->m_RandomNumberGenerator = vnl_random( seed );
  
  this->Init();  
}

PoistatsModel::PoistatsModel( const long seed ){
  this->m_RandomNumberGenerator = vnl_random( seed );
  this->Init();
}

PoistatsModel::~PoistatsModel() {  
//  this->FreeVector( this->m_SeedValues );
}

void PoistatsModel::SetRandomSeed( const long seed ) {
  this->m_RandomNumberGenerator.reseed( seed );
}

double PoistatsModel::GetRandomNumber() {

  const double rangeMin = 0.0;
  const double rangeMax = 1.0;
  
  const double randomNumber = 
    this->m_RandomNumberGenerator.drand64( rangeMin, rangeMax );

  return randomNumber;
}

double PoistatsModel::GetNormallyDistributedRandomNumber() {

  const double randomNumber = 
    this->m_RandomNumberGenerator.normal64();
  return randomNumber;
}

int PoistatsModel::GetRandomInt( const int floor, const int ceiling ) {
  int randomInt  = this->m_RandomNumberGenerator.lrand32( floor, ceiling );
  
  return randomInt;
}

// TODO: this is wrong...  The range will be between 0 and 1
double PoistatsModel::GetRandomNumberWithoutRange() {
  return m_RandomNumberGenerator.drand64();
}

bool PoistatsModel::IsDebug() {
  return m_IsDebug;
}

void PoistatsModel::SetDebug( const bool isDebug ) {
  m_IsDebug = isDebug;
}

bool PoistatsModel::IsUsingPathInitialization() {
  return m_IsUsingPathInitialization;
}

void PoistatsModel::SetUsingPathInitialization( const bool isUsing ) {
  m_IsUsingPathInitialization = isUsing;
}

MRI* 
PoistatsModel::GetEigenVectors() {
  return m_EigenVectors;
}

void 
PoistatsModel::SetEigenVectors( MRI* vectors ) {
  m_EigenVectors = vectors;
}
  
MRI*
PoistatsModel::GetSeedVolume() {
  return m_SeedVolume;
}

void 
PoistatsModel::SetSeedVolume( MRI* volume ) {
  m_SeedVolume = volume;
}
  
std::vector< int >* 
PoistatsModel::GetSeedValues() {
  return m_SeedValues;
}

void 
PoistatsModel::SetSeedValues( std::vector< int >* values ) {
  m_SeedValues = values;
}

void 
PoistatsModel::Init() {

  this->m_IsDebug = false;
  
  this->m_IsUsingPathInitialization = false;

  this->m_EigenVectors = NULL;
  this->m_SeedVolume = NULL;
  this->m_SeedValues = NULL;

  // set up cubic spline filter
  m_CubicSplineFilter = CubicSplineFilterType::New();
  OutputImageType::SizeType size;  
//  size.Fill( 11 );
  size.Fill( 128 );
  m_CubicSplineFilter->SetSize( size );
  
  OutputImageType::PointType origin;
  origin.Fill( 0.0 );
  m_CubicSplineFilter->SetOrigin( origin );
  
  OutputImageType::SpacingType spacing;
  spacing.Fill( 0.1 );

  m_CubicSplineFilter->SetSpacing( spacing );

  m_CubicSplineFilter->SetSplineOrder( 3 );    
  
  // TODO: try to adjust the levels to speed up the algorithm
//  m_CubicSplineFilter->SetNumberOfLevels( 20 );
//  m_CubicSplineFilter->SetNumberOfLevels( 15 );
  m_CubicSplineFilter->SetNumberOfLevels( 8 );
  
  m_CubicSplineFilter->SetGenerateOutputImage( false );

}

void 
PoistatsModel::FreeVector( std::vector< int > *v ) {

  if( v != NULL ) {
    v->clear();
    delete v;
  }
  
}

/**
 * Returns the path with evenly spaced sample points.
 */
PoistatsModel::MatrixPointer
PoistatsModel::RethreadPath(
  MatrixPointer originalPath, const int nNewSamples ) {
    
  // create evenly spaced parametric points for the original path
  const double gridFloor = 0.0;
  const double gridCeiling = 1.0;
  ArrayType originalPathGrid( originalPath->rows() );
  PoistatsModel::SpaceEvenly( &originalPathGrid, gridFloor, gridCeiling );
  
  // interploate the spline
  MatrixPointer rethreadedPath = this->CubicSplineInterpolation( 
    originalPath, &originalPathGrid, nNewSamples );

  // calculate the path length of the new path at each point
  ArrayType magnitude( nNewSamples-1 );
  PoistatsModel::CalculatePathMagnitude( rethreadedPath, &magnitude );

  // create an array of the total path lengths as you progress along the path
  ArrayType cumulativeSum( nNewSamples );
  PoistatsModel::CalculateCumulativeSum( &magnitude, &cumulativeSum );

  // we want to calculate a new spline that is evenly spaces in image
  // coordinates.  We do this by calculating a new spline based on the old one
  // and setting the parametric coordinates based on the path length
  ArrayType normalizedCumulativeSum( nNewSamples );
  double pathLength = cumulativeSum[ cumulativeSum.size() - 1 ];
    
  for( unsigned int cRow=0; cRow<cumulativeSum.size(); cRow++ ) {
    normalizedCumulativeSum[ cRow ] = cumulativeSum[ cRow ] / pathLength;
  }

  MatrixPointer reRethreadedPath = this->CubicSplineInterpolation( 
    rethreadedPath, &normalizedCumulativeSum, nNewSamples );

  delete rethreadedPath;
  rethreadedPath = NULL;

  return reRethreadedPath;
}

/**
 * Returns the cubic spline interpolation given a 3D path and parametric values
 * for that path.
 */
PoistatsModel::MatrixPointer
PoistatsModel::CubicSplineInterpolation( 
  MatrixPointer originalPath, ArrayPointer originalPathGrid, 
  const int nNewSamples ) {
    
  const int spatialDimension = 3;

  PointSetType::Pointer pointSet = PointSetType::New();
  for( unsigned int cRow=0; cRow<originalPath->rows(); cRow++ ) {
    unsigned long i = pointSet->GetNumberOfPoints();

    PointSetType::PointType point;
    point[ 0 ]  = ( *originalPathGrid )[ cRow ];
    
    pointSet->SetPoint(i, point);        
    
    VectorType V;
    
    for( int cColumn = 0; cColumn<spatialDimension; cColumn++ ) {
      V[ cColumn ] = ( *originalPath )[ cRow ][ cColumn ];
    }

    pointSet->SetPointData( i, V );
  }

  CubicSplineFilterType::ArrayType nControlPoints;

  // this is the number of control points that the spline will use to
  // interpolate -- not completely related to the number of control points that
  // poistats uses for the monte carlo part of the algorithm.  Sometimes with 
  // too few points, it still can't find a good fit and fails, so we're 
  // increasing it a bit.
  const int nTotalControlPoints = this->GetNumberOfControlPoints() + 3;
  nControlPoints.Fill( nTotalControlPoints );
  m_CubicSplineFilter->SetNumberOfControlPoints( nControlPoints );

  m_CubicSplineFilter->SetInput( pointSet );

  m_CubicSplineFilter->Update();

  MatrixPointer rethreadedPath = 
    new MatrixType( nNewSamples, spatialDimension );

  const double resampledStepSize = 1.0 / 
    ( static_cast< double >( nNewSamples ) - 1.0 );
        
  // evaluate at regular steps
  for( int cRow=0; cRow<nNewSamples; cRow++ ) {
    
    double t = static_cast< double >( cRow ) * resampledStepSize;
    const double maxT = 1.0;
    if( t > maxT ) {
      t = maxT;
    }
    
    PointSetType::PointType point;
    point[ 0 ] = t;
    VectorType V; 
    m_CubicSplineFilter->EvaluateAtPoint( point, V );
    
    for( int cColumn = 0; cColumn<spatialDimension; cColumn++ ) {
      ( *rethreadedPath )[ cRow ][ cColumn ] = V[ cColumn ];
    }
    
  }
    
  return rethreadedPath;  
}

void PoistatsModel::SetNumberOfControlPoints( const int nPoints ) {
  m_NumberOfControlPoints = nPoints;
}

int PoistatsModel::GetNumberOfControlPoints() {
  return m_NumberOfControlPoints;
}

void
PoistatsModel::SpaceEvenly( ArrayPointer outputArray, const double floor, 
  const double ceiling ) {
  
  const double spacing = ( ceiling - floor ) / ( outputArray->size()-1 );
  for( unsigned int cOutputArray=0; cOutputArray<outputArray->size(); cOutputArray++ ) {
    ( *outputArray )[ cOutputArray ] = floor + cOutputArray * spacing;
  }  
}

void
PoistatsModel::CalculatePathMagnitude( 
  MatrixPointer path, ArrayPointer magnitude ) {

  const int spatialDimension = path->columns();
  
  itk::Array2D< double > pathDifference( path->rows()-1, spatialDimension );

  CalculatePathVectors( path, &pathDifference );
  CalculateMagnitude( &pathDifference, magnitude );  
}

/**
 * Given a 1D array, calculate at an index and all the values preceeding it.
 */
void
PoistatsModel::CalculateCumulativeSum(
  ArrayPointer inputArray, ArrayPointer cumulativeSum ) {
                            
  ( *cumulativeSum )[ 0 ] = 0.0;
  for( unsigned int cRow=1; cRow<cumulativeSum->size(); cRow++ ) {
    ( *cumulativeSum )[ cRow ] = 0.0;
    ( *cumulativeSum )[ cRow ] = ( *cumulativeSum )[ cRow-1 ]
                                  + ( *inputArray )[ cRow-1 ] ;
  }
  
}

/**
 * Given path, gets the vector pointing from one point to the following point.
 */
void
PoistatsModel::CalculatePathVectors( 
  MatrixPointer path, MatrixPointer vectors ) {

  const int spatialDimension = path->columns();

  for( int cRow=0; cRow < static_cast< int >( path->rows() )-1; cRow++ ) {
    
    for( int cColumn=0; cColumn<spatialDimension; cColumn++ ) {

      double difference = ( *path )[ cRow+1 ][ cColumn ]- 
                          ( *path )[ cRow ][ cColumn ];

      ( *vectors )[ cRow ][ cColumn ] = difference;
    }
  }
                            
}

void
PoistatsModel::CalculateMagnitude( 
  MatrixPointer path, ArrayPointer magnitude ) {

  const int spatialDimension = path->columns();
  for( int cRow=0; cRow < static_cast< int >( path->rows() ); cRow++ ) {
    
    double rowSumOfDifferenceSquared = 0.0;
    for( int cColumn=0; cColumn<spatialDimension; cColumn++ ) {
      
      rowSumOfDifferenceSquared += ( *path )[ cRow ][ cColumn ] * 
                                   ( *path )[ cRow ][ cColumn ];
      
    }
    ( *magnitude )[ cRow ] = sqrt( rowSumOfDifferenceSquared );
  }
  
}
