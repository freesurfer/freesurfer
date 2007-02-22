#include "PoistatsReplica.h"

#include <iostream>
#include <cmath>

PoistatsReplica::PoistatsReplica() {
  Init();
}
  
PoistatsReplica::PoistatsReplica( PoistatsModel* model ) {
  
  this->SetModel( model );
  
  Init();  
}

PoistatsReplica::~PoistatsReplica() {  
  DeletePathIfNotNull( m_BasePath );
  DeletePathIfNotNull( m_PreviousTrialPath );
  DeletePathIfNotNull( m_CurrentTrialPath );  
  DeletePathIfNotNull( m_BestTrialPath );
  DeleteArrayIfNotNull( m_BestTrialPathProbabilities );
}

void 
PoistatsReplica::Init() {
  m_BasePath = NULL;
  m_PreviousTrialPath = NULL;
  m_CurrentTrialPath = NULL;
  m_BestTrialPath = NULL;
  m_BestTrialPathProbabilities = NULL;
}

void
PoistatsReplica::SetModel( PoistatsModel* model ) {
  m_PoistatsModel = model;
}

double 
PoistatsReplica::GetCurrentMeanEnergy() const {
  return m_CurrentMeanEnergy;
}

void 
PoistatsReplica::SetCurrentMeanEnergy( const double energy ) {
  m_CurrentMeanEnergy = energy;
}
  
double 
PoistatsReplica::GetPreviousMeanEnergy() const {
  return m_PreviousMeanEnergy;
}

void 
PoistatsReplica::SetPreviousMeanEnergy( const double energy ) {
  m_PreviousMeanEnergy = energy;
}

double 
PoistatsReplica::GetTemperature() {
  return m_Temperature;
}

void
PoistatsReplica::SetTemperature( const double temperature ) {
  m_Temperature = temperature;  
}

void 
PoistatsReplica::CopyCurrentToPreviousEnergy() {
  m_PreviousMeanEnergy = m_CurrentMeanEnergy;
}

void PoistatsReplica::CoolTemperature( const double coolingFactor ) {
  m_Temperature *= coolingFactor;
}

bool
PoistatsReplica::ShouldUpdateEnergy() {

  const double delta  = ( m_CurrentMeanEnergy - m_PreviousMeanEnergy ) / 
    m_Temperature;

  double updateProbability = exp( -delta );
  
  static const double maxProbablity = 1.0;
  
  if( updateProbability > maxProbablity ) {
    updateProbability = maxProbablity;
  }
  
  const double randomNumber = m_PoistatsModel->GetRandomNumber();
  const bool shouldUpdate = ( randomNumber <= updateProbability );

  return shouldUpdate;  
}

void PoistatsReplica::ResetCurrentToPreviousEnergy() {
  m_CurrentMeanEnergy = m_PreviousMeanEnergy;
}

PoistatsReplica::MatrixPointer PoistatsReplica::GetBasePath() {
  return m_BasePath;
}

void PoistatsReplica::SetBasePath( const MatrixPointer basePath ) {
  DeletePathIfNotNull( m_BasePath );
  m_BasePath = new MatrixType( *basePath );
}

PoistatsReplica::MatrixPointer PoistatsReplica::GetPreviousTrialPath() {
  return m_PreviousTrialPath;
}

void PoistatsReplica::SetPreviousTrialPath( const MatrixPointer path ) {
  DeletePathIfNotNull( m_PreviousTrialPath );
  m_PreviousTrialPath = new MatrixType( *path );
}

void
PoistatsReplica::SpaceEvenly( ArrayPointer outputArray, const double floor, 
  const double ceiling ) {
  
  const double spacing = ( ceiling - floor ) / ( outputArray->size()-1 );
  for( unsigned int cOutputArray=0; cOutputArray<outputArray->size(); cOutputArray++ ) {
    ( *outputArray )[ cOutputArray ] = floor + cOutputArray * spacing;
  }  
}

void
PoistatsReplica::CalculatePathMagnitude( 
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
PoistatsReplica::CalculateCumulativeSum(
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
PoistatsReplica::CalculatePathVectors( 
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
PoistatsReplica::CalculateMagnitude( 
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

PoistatsReplica::MatrixPointer
PoistatsReplica::GetCurrentTrialPath() {
  return m_CurrentTrialPath;
}

void
PoistatsReplica::SetCurrentTrialPath( MatrixPointer path) {
  DeletePathIfNotNull( m_CurrentTrialPath );
  m_CurrentTrialPath = new MatrixType( *path );
}

PoistatsReplica::MatrixPointer
PoistatsReplica::GetBestTrialPath() {
  return m_BestTrialPath;
}

void
PoistatsReplica::SetBestTrialPath( MatrixPointer path) {
  DeletePathIfNotNull( m_BestTrialPath );
  m_BestTrialPath = new MatrixType( *path );
}

void PoistatsReplica::DeletePathIfNotNull( MatrixPointer path ) {
  if( path != NULL ) {
    delete path;
  }
  path = NULL;
}

void PoistatsReplica::DeleteArrayIfNotNull( ArrayPointer array ) {
  if( array != NULL ) {
    delete array;
  }
  array = NULL;
}

void PoistatsReplica::GetPerturbedBasePath( MatrixPointer perturbedPath, 
  const double sigma, const MatrixPointer possibleStartPoints, 
  const MatrixPointer possibleEndPoints) {
    
  // perturb = sigma*diag(rand(ncontrolpoints,1))*sphererand(ncontrolpoints);
  MatrixType randomUnitSphere( m_BasePath->rows(), m_BasePath->cols() );

  this->GenerateUnitSphereRandom( &randomUnitSphere );
  
  const int nBasePoints = m_BasePath->rows();
  vnl_matrix< double > randomPoints( nBasePoints, nBasePoints );
  randomPoints.fill( 0.0 );
  
  // fill the diagonals with random numbers
  for( int cControlPoint=0; cControlPoint<nBasePoints; cControlPoint++ ) {
    randomPoints[ cControlPoint ][ cControlPoint ] = m_PoistatsModel->GetRandomNumber();
  }
        
  vnl_matrix< double > perturb = ( randomPoints * randomUnitSphere ) * sigma;
    
  MatrixPointer currentBasePath = m_BasePath;
  
  // lowtrialpath(2:end-1,:) = basepath{i}(2:end-1,:) + perturb;
  for( unsigned int cRow=1; cRow<currentBasePath->rows()-1; cRow++ ) {
    for( unsigned int cColumn=0; cColumn<currentBasePath->cols(); cColumn++ ) {
      const double currentPerturb = perturb[ cRow-1 ][ cColumn ];
      ( *perturbedPath )[ cRow ][ cColumn ] = 
        ( *currentBasePath )[ cRow ][ cColumn ] + currentPerturb;
    }
  }
  
  // TODO: I can probably get the row using vnl rather than copying
  const int nSpatialDimensions = 3;

  // lowtrialpath(1,:) = ...
  //   constrainedrnd(basepath{i}(1,:), seedstart, sigma);
  ArrayType newRandomPoint( nSpatialDimensions );
  
  this->GenerateConstrainedRandomPoint3D( currentBasePath->get_row( 0 ), 
    possibleStartPoints, sigma, &newRandomPoint );
  
  for( int cDimension=0; cDimension<nSpatialDimensions; cDimension++ ) {
    ( *perturbedPath )[ 0 ][ cDimension ] = newRandomPoint[ cDimension ];
  }

  // lowtrialpath(end,:) = ...
  //   constrainedrnd(basepath{i}(end,:), seedend, sigma);
  this->GenerateConstrainedRandomPoint3D( 
    currentBasePath->get_row( currentBasePath->rows()-1 ), 
    possibleEndPoints, sigma, &newRandomPoint );
  for( int cDimension=0; cDimension<nSpatialDimensions; cDimension++ ) {
    ( *perturbedPath )[ perturbedPath->rows()-1 ][ cDimension ] = 
      newRandomPoint[ cDimension ];
  }
  
}

/**
 * Generate a random point on sphere
 * Muller, M. E. "A Note on a Method for Generating Points Uniformly on 
 * N-Dimensional Spheres" Comm. Assoc. Comput. Mach. 2, 19-20, 1959. 
 * Marsaglia, G. "Choosing a Point from the Surface of a Sphere." Ann. Math. 
 * Stat. 43, 645-646, 1972.
 */
void
PoistatsReplica::GenerateUnitSphereRandom( MatrixPointer randomUnitSphere ) {
    
  const int numberOfPoints = randomUnitSphere->rows();

  for( unsigned int cRow=0; cRow<randomUnitSphere->rows(); cRow++ ) {
    for( unsigned int cColumn=0; cColumn<randomUnitSphere->cols(); cColumn++ ) {

      ( *randomUnitSphere )[ cRow ][ cColumn ] = 
        m_PoistatsModel->GetNormallyDistributedRandomNumber();
    }
  }
  
  ArrayType magnitude( numberOfPoints );  
  CalculateMagnitude( randomUnitSphere, &magnitude );
  
  for( unsigned int cRow=0; cRow<randomUnitSphere->rows(); cRow++ ) {
    
    const double inverseMagnitude = ( 1.0 / magnitude[ cRow ] );
    for( unsigned int cColumn=0; cColumn<randomUnitSphere->cols(); cColumn++ ) {

      ( *randomUnitSphere )[ cRow ][ cColumn ] = 
        ( *randomUnitSphere )[ cRow ][ cColumn ] * inverseMagnitude;

    }
    
  }
  
}

void
PoistatsReplica::GenerateConstrainedRandomPoint3D(

  vnl_vector< double > currentPoint,
  MatrixPointer possibleNewPoints,
  const double sigma,
  ArrayPointer newRandomPoint ) {
    
  std::vector< int > indicesOfPointsWithinRadius;

  bool isPointFound = FindPointsWithinRadius( &currentPoint, possibleNewPoints, 
    sigma, &indicesOfPointsWithinRadius );
  
  if( isPointFound ) {
        
    int min = 0;
    int max = indicesOfPointsWithinRadius.size()-1;    
    int randomIndex = m_PoistatsModel->GetRandomInt( min, max );
    
    for( unsigned int cColumn=0; cColumn<newRandomPoint->size(); cColumn++ ) {
      
      double perturb = m_PoistatsModel->GetRandomNumberWithoutRange() - 0.5; 

      int randomPointIndex = indicesOfPointsWithinRadius[randomIndex];

      ( *newRandomPoint )[ cColumn ] = 
        ( *possibleNewPoints )[ randomPointIndex ][ cColumn ] + perturb;
        
    }
  } else {

    // if the point isn't found, then we perturb the current point a little
    itk::Array2D< double > randomUnitSphere( 1, newRandomPoint->size() );
    GenerateUnitSphereRandom( &randomUnitSphere );

    for( unsigned int cColumn=0; cColumn<newRandomPoint->size(); cColumn++ ) {
      
      double perturb = sigma * m_PoistatsModel->GetRandomNumberWithoutRange()
        * randomUnitSphere[ 1 ][ cColumn ];
      
      ( *newRandomPoint )[ cColumn ] = currentPoint[ cColumn ] + perturb;
      
    }
    
  }
  
}

/**
 * Returns the indices of the points that are within the radius of the center
 * point.
 */
bool PoistatsReplica::FindPointsWithinRadius(
  vnl_vector< double >* center,
  MatrixPointer points,
  const double sigma,
  std::vector< int > *indices ) {
    
  bool isPointFound = false;

  // create a difference matrix between the center and the other points
  MatrixType differenceMatrix( points->rows(), points->cols() );
  for( unsigned int cRow=0; cRow<points->rows(); cRow++ ) {
    for( unsigned int cColumn=0; cColumn<points->cols(); cColumn++ ) {
      differenceMatrix[ cRow ][ cColumn ] = ( *points )[ cRow ][ cColumn ]
                                            - ( *center )[ cColumn ];
    }  
  }
  
  ArrayType distances( points->rows() );
  PoistatsReplica::CalculateMagnitude( &differenceMatrix, &distances );
      
  for( unsigned int cRow=0; cRow<distances.size(); cRow++ ) {
    
    if( distances[ cRow ] < sigma ) {
        
      indices->insert( indices->begin(), cRow );

      // this causes problems...
//      indices->push_back( cRow );
      isPointFound = true;
    }
    
  }
  
  return isPointFound;
}

void
PoistatsReplica::CopyPath( 
  const MatrixPointer source,
  MatrixPointer destination ) {

  for( unsigned int cRow = 0; cRow < source->rows(); cRow++ ) {
    for( unsigned int cColumn = 0; cColumn < source->cols(); cColumn++ ) {
    
      ( *destination )[ cRow ][ cColumn] = ( *source )[ cRow ][ cColumn];
      
    }
  }

}

void PoistatsReplica::FoundBestPath( const MatrixPointer basePath ) {  
  CopyPath( basePath, m_BasePath );
  CopyPath( m_CurrentTrialPath, m_BestTrialPath );
}

void PoistatsReplica::CopyCurrentToPreviousTrialPath() {
  CopyPath( m_CurrentTrialPath, m_PreviousTrialPath );
}

void 
PoistatsReplica::SetBestTrialPathProbabilities( ArrayPointer probabilities ) {
  DeleteArrayIfNotNull( m_BestTrialPathProbabilities );
  m_BestTrialPathProbabilities = new ArrayType( *probabilities );
}
  
PoistatsReplica::ArrayPointer 
PoistatsReplica::GetBestTrialPathProbabilities() {
  return this->m_BestTrialPathProbabilities;
}
