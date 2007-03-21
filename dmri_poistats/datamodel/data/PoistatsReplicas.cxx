#include <iostream>

#include "PoistatsReplicas.h"

PoistatsReplicas::PoistatsReplicas() 
{
}

PoistatsReplicas::PoistatsReplicas( PoistatsModel *model, const int nReplicas ) 
{

  m_PoistatsModel = model;

  this->SetNumberOfReplicas( nReplicas );
  
  m_InitialPoints = NULL;
  
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
  spacing.Fill( 0.01 );

  m_CubicSplineFilter->SetSpacing( spacing );

  m_CubicSplineFilter->SetSplineOrder( 3 );    
  
  // TODO: try to adjust the levels to speed up the algorithm
//  m_CubicSplineFilter->SetNumberOfLevels( 20 );
//  m_CubicSplineFilter->SetNumberOfLevels( 15 );
  m_CubicSplineFilter->SetNumberOfLevels( 8 );
  
  m_CubicSplineFilter->SetGenerateOutputImage( false );
  
}

PoistatsReplicas::~PoistatsReplicas() {

  if( m_Replicas != NULL ) {
    delete []m_Replicas;
  }
  m_Replicas = NULL;
  
  if( m_InitialPoints != NULL ) {
    delete m_InitialPoints;
  }
  m_InitialPoints = NULL;
}

int
PoistatsReplicas::GetNumberOfReplicas() {
  return m_NumberOfReplicas;
}

void
PoistatsReplicas::SetNumberOfReplicas( const int nReplicas ) {

  m_NumberOfReplicas = nReplicas;

//  if( m_Replicas != NULL ) {
//    delete[] m_Replicas;
//  }
  
  m_Replicas = new PoistatsReplica[ m_NumberOfReplicas ];

  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetModel( m_PoistatsModel );
  }
  
}


double
PoistatsReplicas::GetMinimumCurrentEnergy() {
  
  double min = m_Replicas[ 0 ].GetCurrentMeanEnergy();
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    const double energy = m_Replicas[ cReplica ].GetCurrentMeanEnergy();
    if( min > energy ) {
      min = energy;
    }
  }
  
  return min;
}

void 
PoistatsReplicas::FillCurrentMeanEnergies( const double energy ) {
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetCurrentMeanEnergy( energy );
  }
}

void 
PoistatsReplicas::FillPreviousMeanEnergies( const double energy ) {
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetPreviousMeanEnergy( energy );
  }
}

double
PoistatsReplicas::GetCurrentMeanOfEnergies() const {
  double sum = 0.0;
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    sum += m_Replicas[ cReplica ].GetCurrentMeanEnergy();
  }
  
  const double mean = sum / static_cast< double >( m_NumberOfReplicas );
  
  return mean;
}

double
PoistatsReplicas::GetPreviousMeanOfEnergies() const {
  double sum = 0.0;
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    sum += m_Replicas[ cReplica ].GetPreviousMeanEnergy();
  }
  
  const double mean = sum / static_cast< double >( m_NumberOfReplicas );
  
  return mean;
}

double 
PoistatsReplicas::GetNormalizedMeanCurrentPreviousEnergiesDifference() const {
  const double currentMean = GetCurrentMeanOfEnergies();
  const double previousMean = GetPreviousMeanOfEnergies();
  const double difference = currentMean - previousMean;
    
  const double normalizedDifference = difference / currentMean;

  return normalizedDifference;
}

double
PoistatsReplicas::GetCurrentMeanEnergy( const int replica ) {
  return m_Replicas[ replica ].GetCurrentMeanEnergy();
}

void PoistatsReplicas::SetCurrentMeanEnergy( 
  const int replica, const double energy ) {

  m_Replicas[ replica ].SetCurrentMeanEnergy( energy );

}

double
PoistatsReplicas::GetPreviousMeanEnergy( const int replica ) {
  return m_Replicas[ replica ].GetPreviousMeanEnergy();
}

void PoistatsReplicas::GetRandomSortedFirstSecondReplicas( 
  int &first, int &second ) {
    
  itk::Array< double >energies( m_NumberOfReplicas );
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    energies[ cReplica ]= m_Replicas[ cReplica ].GetCurrentMeanEnergy();
  }

  itk::Array< int > sortedReplicaIndices( m_NumberOfReplicas );
  SortArray( &energies, &sortedReplicaIndices );
        
  const int randomIndex = 
    static_cast< int >( floor( 
      static_cast< double >( m_NumberOfReplicas - 1 )
        * this->m_PoistatsModel->GetRandomNumber() ) );

  first = sortedReplicaIndices[ randomIndex ];
  second = sortedReplicaIndices[ randomIndex + 1 ];
    
}

void
PoistatsReplicas::SortArray(
  itk::Array< double >* unsorted,
  itk::Array< int >* sortedIndices ) {

  itk::Array< double >sorted( *unsorted );

  std::sort( sorted.begin(), sorted.end() );
  
  const int invalidIndex = -1;
  sortedIndices->Fill( invalidIndex );

  // create the list of indices mapping the unsorted values into ascending order
  for( unsigned int cSorted=0; cSorted<sorted.size(); cSorted++ ) {
    double nextLowest = sorted[ cSorted ];
    
    for( unsigned int cUnsorted=0; cUnsorted<unsorted->size(); cUnsorted++ ) {
      if( nextLowest == ( *unsorted )[ cUnsorted ] ) {
        ( *sortedIndices )[ cSorted ] = cUnsorted;
      }
    }
    
  }

}

void 
PoistatsReplicas::CopyCurrentToPreviousEnergies() {
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].CopyCurrentToPreviousEnergy();
  }
}

void 
PoistatsReplicas::CopyCurrentToPreviousEnergy( int replica ) {
  m_Replicas[ replica ].CopyCurrentToPreviousEnergy();
}

double PoistatsReplicas::GetTemperature( const int replica ) {
  return m_Replicas[ replica ].GetTemperature();
}

void PoistatsReplicas::SetTemperature( 
  const int replica, const double temperature ) {

  m_Replicas[ replica ].SetTemperature( temperature );

}

void PoistatsReplicas::CoolTemperatures( const double coolingFactor ) {
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].CoolTemperature( coolingFactor );
  }
} 

bool PoistatsReplicas::ShouldUpdateEnergy( const int replica ) {

  return m_Replicas[ replica ].ShouldUpdateEnergy();
} 

void PoistatsReplicas::ResetCurrentToPreviousEnergy( const int replica ) {
  m_Replicas[ replica ].ResetCurrentToPreviousEnergy();
}

double PoistatsReplicas::CalculateProbablityExchange( const int replica1,
  const int replica2 ) {

  // Dbeta = 1/temp(idx) - 1/temp(jdx);        
  // Denergy = replicaenergies(idx) - replicaenergies(jdx);
  // Delta = -Dbeta * Denergy;
  // exchprob = min(1,exp(-Delta));
  
  const double temperature1 = m_Replicas[ replica1 ].GetTemperature();
  const double energy1 = m_Replicas[ replica1 ].GetCurrentMeanEnergy();
  
  const double temperature2 = m_Replicas[ replica2 ].GetTemperature();
  const double energy2 = m_Replicas[ replica2 ].GetCurrentMeanEnergy();
                               
  const double betaDifference = 1.0/temperature1 - 1.0/temperature2;

  const double energyDifference = energy1 - energy2;
    
  const double delta = -betaDifference * energyDifference;

  const double minProbablity = 1.0;
  double exchangeProbablity = exp( -delta );
  if( minProbablity < exchangeProbablity ) {
    exchangeProbablity = minProbablity;
  }
  
  return exchangeProbablity;
}

void PoistatsReplicas::ExchangeTemperatures( const int replica1,
  const int replica2 ) {

  const double temperature1 = m_Replicas[ replica1 ].GetTemperature();
  const double temperature2 = m_Replicas[ replica2 ].GetTemperature();
  
  m_Replicas[ replica1 ].SetTemperature( temperature2 );    
  m_Replicas[ replica2 ].SetTemperature( temperature1 );
  
}

void PoistatsReplicas::SetBasePaths( const MatrixPointer basePath ) {

  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetBasePath( basePath );
  }

}
  
PoistatsReplicas::MatrixPointer 
PoistatsReplicas::GetBasePath( const int replica ) {
  return m_Replicas[ replica ].GetBasePath();
}

PoistatsReplicas::MatrixPointer 
PoistatsReplicas::GetPreviousTrialPath( const int replica ) {
  return m_Replicas[ replica ].GetPreviousTrialPath();
}

void PoistatsReplicas::SetInitialPoints( const MatrixPointer points ) {

  if( m_InitialPoints != NULL ) {
    delete m_InitialPoints;
  }
  
  m_InitialPoints = new MatrixType( *points );
  
  const int nEndPoints = 2;
  const int nTotalControlPoints = GetNumberOfControlPoints() + nEndPoints;

  MatrixPointer originalPath = this->RethreadPath( 
    m_InitialPoints, nTotalControlPoints );
  
  SetBasePaths( originalPath );
  
  MatrixPointer previousPath = 
    this->RethreadPath( originalPath, GetNumberOfSteps() );    

  SetPreviousTrialPaths( previousPath );
  
  const int nSpatialDimensions = 3;
  MatrixType zeroes( GetNumberOfSteps(), nSpatialDimensions );
  zeroes.Fill( 0.0 );
  
  SetCurrentTrialPaths( &zeroes );
  SetBestTrialPaths( &zeroes );  
  
  delete previousPath;
  delete originalPath;  
}

void PoistatsReplicas::SetNumberOfControlPoints( const int nPoints ) {
  m_NumberOfControlPoints = nPoints;
}

int PoistatsReplicas::GetNumberOfControlPoints() {
  return m_NumberOfControlPoints;
}

void PoistatsReplicas::SetPreviousTrialPaths( const MatrixPointer path ) {

  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetPreviousTrialPath( path );
  }

}

void PoistatsReplicas::SetNumberOfSteps( const int nSteps ) {
  m_NumberOfSteps = nSteps;
}

int PoistatsReplicas::GetNumberOfSteps() {
  return m_NumberOfSteps;
}

void PoistatsReplicas::SetCurrentTrialPaths( const MatrixPointer path ) {

  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetCurrentTrialPath( path );
  }

}

PoistatsReplicas::MatrixPointer 
PoistatsReplicas::GetCurrentTrialPath( const int replica ) {
  return m_Replicas[ replica ].GetCurrentTrialPath();
}

void PoistatsReplicas::SetBestTrialPaths( const MatrixPointer path ) {

  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    m_Replicas[ cReplica ].SetBestTrialPath( path );
  }

}

PoistatsReplicas::MatrixPointer 
PoistatsReplicas::GetBestTrialPath( const int replica ) {
  return m_Replicas[ replica ].GetBestTrialPath();
}

PoistatsReplicas::MatrixListType 
PoistatsReplicas::GetBestTrialPaths() {
  MatrixListType bestPaths;
  for( int cReplica=0; cReplica<this->GetNumberOfReplicas(); cReplica++ ) {
    bestPaths.push_back( GetBestTrialPath( cReplica ) );
  }
  
  return bestPaths;
}

void PoistatsReplicas::SetBestTrialPathProbabilities( const int replica,
  ArrayPointer probabilities ){
  
  m_Replicas[ replica ].SetBestTrialPathProbabilities( probabilities );
  
}

void
PoistatsReplicas::GetBestTrialPathsProbabilities( MatrixPointer probabilities ){

  for( int cReplica=0; cReplica<this->GetNumberOfReplicas(); cReplica++ ) {

    ArrayPointer replicaProbabilities = 
      m_Replicas[ cReplica ].GetBestTrialPathProbabilities();
      
    for( unsigned int cStep=0; cStep<replicaProbabilities->size(); cStep++ ) {
      ( *probabilities )[ cReplica ][ cStep ] = 
        ( *replicaProbabilities )[ cStep ];
    }    
    
  }
  
}

PoistatsReplicas::ArrayPointer 
PoistatsReplicas::GetBestTrialPathProbabilities( const int replica) {
  return m_Replicas[ replica ].GetBestTrialPathProbabilities();
}


void PoistatsReplicas::SpaceTemperaturesEvenly( 
  const double floorTemp, const double ceilingTemp ) {

  if( this->GetNumberOfReplicas() > 1 ) {

    ArrayType temperatures( this->GetNumberOfReplicas() );  
    PoistatsReplica::SpaceEvenly( &temperatures, floorTemp, ceilingTemp );
    for( int cReplica=0; cReplica<this->GetNumberOfReplicas(); cReplica++ ) {
      const double temperature = temperatures[ cReplica ];
      this->SetTemperature( cReplica, temperature );
    }
   
  } else {
    // if we just have one replica, then we'll just take the center temp
    this->SetTemperature( 0, ( ceilingTemp-floorTemp ) / 2.0 );
  }
    
}

void PoistatsReplicas::GetPerturbedBasePath( const int replica, 
  MatrixPointer perturbedPath, 
  const double sigma, const MatrixPointer possibleStartPoints, 
  const MatrixPointer possibleEndPoints) {
    
  m_Replicas[ replica ].GetPerturbedBasePath( perturbedPath, sigma, 
    possibleStartPoints, possibleEndPoints);
    
}

void PoistatsReplicas::PerturbCurrentTrialPath( const int replica, 
  MatrixPointer lowTrialPath, const int nSteps ) {
    
  MatrixPointer perturbedTrialPath = this->RethreadPath( lowTrialPath, nSteps );
  
  m_Replicas[ replica ].SetCurrentTrialPath( perturbedTrialPath );
  
  delete perturbedTrialPath;
    
}


void PoistatsReplicas::FoundBestPath( 
  const int replica, const MatrixPointer basePath ) {
  
  m_Replicas[ replica ].FoundBestPath( basePath );

}

void PoistatsReplicas::CopyCurrentToPreviousTrialPath( const int replica ) {
  m_Replicas[ replica ].CopyCurrentToPreviousTrialPath();
}

/**
 * Returns the path with evenly spaced sample points.
 */
PoistatsReplicas::MatrixPointer
PoistatsReplicas::RethreadPath(
  MatrixPointer originalPath, const int nNewSamples ) {
    
  // create evenly spaced parametric points for the original path
  const double gridFloor = 0.0;
  const double gridCeiling = 1.0;
  ArrayType originalPathGrid( originalPath->rows() );
  PoistatsReplica::SpaceEvenly( &originalPathGrid, gridFloor, gridCeiling );
  
  // interploate the spline
  MatrixPointer rethreadedPath = this->CubicSplineInterpolation( 
    originalPath, &originalPathGrid, nNewSamples );

  // calculate the path length of the new path at each point
  ArrayType magnitude( nNewSamples-1 );
  PoistatsReplica::CalculatePathMagnitude( rethreadedPath, &magnitude );

  // create an array of the total path lengths as you progress along the path
  ArrayType cumulativeSum( nNewSamples );
  PoistatsReplica::CalculateCumulativeSum( &magnitude, &cumulativeSum );

  // we want to calculate a new spline that is evenly spaces in imag
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
PoistatsReplicas::MatrixPointer
PoistatsReplicas::CubicSplineInterpolation( 
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
