
#include "PoistatsReplicas.h"

#include <iostream>

#include "dmri_poistats/datamodel/utils/EigenVectorInitPathStrategy.h"

PoistatsReplicas::PoistatsReplicas() 
{
}

PoistatsReplicas::PoistatsReplicas( PoistatsModel *model, const int nReplicas ) 
{

  m_PoistatsModel = model;

  this->SetNumberOfReplicas( nReplicas );
  
  m_InitialPoints = NULL;
    
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
  const double temperature2 = m_Replicas[ replica2 ].GetTemperature();
  // TODO: should this be the absolute value of the difference?
  const double betaDifference = 1.0/temperature1 - 1.0/temperature2;

  const double energy1 = m_Replicas[ replica1 ].GetCurrentMeanEnergy();
  const double energy2 = m_Replicas[ replica2 ].GetCurrentMeanEnergy();
  // TODO: should this be the absolute value of the energy differene?
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

// TODO: I guess you could move this to a strategy too
void PoistatsReplicas::SetInitialPoints( const MatrixPointer points ) {

  if( m_InitialPoints != NULL ) {
    delete m_InitialPoints;
  }
  
  m_InitialPoints = new MatrixType( *points );
  
  const int nEndPoints = 2;
  const int nTotalControlPoints = m_PoistatsModel->GetNumberOfControlPoints() + nEndPoints;

  MatrixPointer originalPath = m_PoistatsModel->RethreadPath( 
    m_InitialPoints, nTotalControlPoints );
  
  SetBasePaths( originalPath );
  
  MatrixPointer previousPath = 
    m_PoistatsModel->RethreadPath( originalPath, GetNumberOfSteps() );    

  SetPreviousTrialPaths( previousPath );
  
  const int nSpatialDimensions = 3;
  MatrixType zeroes( GetNumberOfSteps(), nSpatialDimensions );
  zeroes.Fill( 0.0 );
  
  SetCurrentTrialPaths( &zeroes );
  SetBestTrialPaths( &zeroes );  
  
  delete previousPath;
  delete originalPath;  
}

void PoistatsReplicas::InitializePaths() {
  
  const int nEndPoints = 2;
  const int nTotalControlPoints = m_PoistatsModel->GetNumberOfControlPoints() + nEndPoints;
  
  InitializePath *initializePath = m_PoistatsModel->GetPathInitializer();

  // for each replica create a different initial path  
  for( int cReplica=0; cReplica<m_NumberOfReplicas; cReplica++ ) {
    
    // calculate a new path
    initializePath->CalculateInitialPath();

    // get the new path
    MatrixPointer initialPath = initializePath->GetInitialPath();

    // rethread the path to the base path
    MatrixPointer base = m_PoistatsModel->RethreadPath( initialPath, nTotalControlPoints );    
    m_Replicas[ cReplica ].SetBasePath( base );
    delete base;
    
    // rethread to the previous path
    MatrixPointer previous = m_PoistatsModel->RethreadPath( initialPath, 
      GetNumberOfSteps() );
    m_Replicas[ cReplica ].SetPreviousTrialPath( previous );    
    delete previous;
    
    // there's no need to delete initialPath.  InitializePath will do it 
      
  }
  
  // set all the trial paths to zero
  const int nSpatialDimensions = 3;
  MatrixType zeroes( GetNumberOfSteps(), nSpatialDimensions );
  zeroes.Fill( 0.0 );
  
  SetCurrentTrialPaths( &zeroes );
  SetBestTrialPaths( &zeroes );  
  
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
    PoistatsModel::SpaceEvenly( &temperatures, floorTemp, ceilingTemp );
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
    
  MatrixPointer perturbedTrialPath = m_PoistatsModel->RethreadPath( lowTrialPath, nSteps );
  
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
