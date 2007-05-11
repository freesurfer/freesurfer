#include "EigenVectorInitPathStrategy.h"

// for MRIextractDistanceMap
#include "fastmarching.h"

// for random number generator
#include "numerics.h"

#include <iostream>

EigenVectorInitPathStrategy::EigenVectorInitPathStrategy() {
  m_EigenVectors = NULL;
  m_GradientVolumes.clear();  

  m_IsGradientCached = false;
}

EigenVectorInitPathStrategy::EigenVectorInitPathStrategy( PoistatsModel *model ) {

  m_PoistatsModel = model;

  this->SetEigenVectors( m_PoistatsModel->GetEigenVectors() );
  m_GradientVolumes.clear();
  
  this->SetSeedVolume( m_PoistatsModel->GetSeedVolume() );
  this->SetSeedValues( m_PoistatsModel->GetSeedValues() );
  
  m_IsGradientCached = false;

}

EigenVectorInitPathStrategy::~EigenVectorInitPathStrategy() {
  this->DeleteGradientVolumes();
}

void 
EigenVectorInitPathStrategy::SetEigenVectors( MRI *eigenVectors ) {
  m_EigenVectors = eigenVectors;
}

void
EigenVectorInitPathStrategy::CacheGradientVolumes() {
  
  std::vector< int > *seedValues = m_PoistatsModel->GetSeedValues();
  
  MRI *seedVolume = m_PoistatsModel->GetSeedVolume();
    
  // make sure that the seed volume and seed values have been set
  if( seedVolume != NULL && seedValues != NULL ) {
    
    if( !m_GradientVolumes.empty() ) {
      // deallocate the gradient volumes if they've been allocated already
      this->DeleteGradientVolumes();
    }
    
    // there will be one less gradient volumes as seed values
    const int nGradientVolumes = seedValues->size() - 1;
      
    for( int cSeeds = 0; cSeeds < nGradientVolumes; cSeeds++ ) {
    
      // using the next seed as the destination
      const int destinationSeedLabel = ( *seedValues )[ cSeeds+1 ];
      
      // create the distance transform, all the values should be set before running
      MRI *distanceTransform = this->GetDistanceVolume( destinationSeedLabel );
      
      // now that we have the transform to get to the second seed, we want to take
      // the derivatives of the distance transform to see what direction we need
      // to take to get to the second seed
      MRI *gradientsVolume = MRIsobel( distanceTransform, NULL, NULL );
      
      m_GradientVolumes.push_back( gradientsVolume );
      
      if( distanceTransform != NULL ) {
        MRIfree( &distanceTransform );
      }
      
    }
    
    m_IsGradientCached = true;
  
  } else {
    std::cerr << "EigenVectorInitPathStrategy::CacheGradientVolumes() -- seed volume is null" << std::endl;
  }
    
}

void 
EigenVectorInitPathStrategy::CalculateInitialPath() {
  
  if( !m_IsGradientCached ) {
    this->CacheGradientVolumes();
  }
  
  std::vector< int > *seedValues = m_PoistatsModel->GetSeedValues();
  
  // this will be a vector with 3 columns
  std::vector< double* > path;
  
  for( unsigned int cSeeds = 0; cSeeds < ( seedValues->size() - 1 ); cSeeds++ ) {
    
    // using the next seed as the destination
    const int destinationSeedLabel = ( *seedValues )[ cSeeds+1 ];
    
    // now that we have the transform to get to the second seed, we want to take
    // the derivatives of the distance transform to see what direction we need
    // to take to get to the second seed
    MRI *gradientsVolume = m_GradientVolumes[ cSeeds ];
    
    // pick a random starting point within the starting seed region
    const int startSeedLabel = ( *seedValues )[ cSeeds ];
    int currentPointInt[3];
    this->GetRandomSeedPoint( currentPointInt, startSeedLabel );

    // pick a random end point.  This will be used for determining the initial
    // previous point so that the correct orientation of the eigenvector can be
    // obtained.
    int endPoint[3];
    this->GetRandomSeedPoint( endPoint, destinationSeedLabel );  
    
    // there values are only use when iterating and saving double values, rather
    // than rounding to the int
    double currentPointDouble[3];
    double previousPointDouble[3];
    for( int i=0; i<3; i++ ) {
      currentPointDouble[i] = currentPointInt[i];
      
      // the previous point is set to this initially so that when the orientation
      // of the eigenvector is determined, we can find the one that points most
      // toward the end point initially
      previousPointDouble[i] = currentPointInt[i] - endPoint[i];
    }
    
    // TODO: we can speed things up by caching the seed points
    // starting at our random seed point, we want to move it toward the second
    // seed region and also following the eigenvector  
    std::vector< int* > destinationSeeds = this->GetSeedPoints( destinationSeedLabel );
        
    // keep iterating until you've reached your destination region
    while( !this->IsInRegion( currentPointInt, destinationSeeds ) ) {

      // get the eigenvector
      double eigenVector[3];
      this->GetEigenVector( currentPointInt, eigenVector );
      
      // get the eigenvector or flip it based on the previous point
      if( this->ShouldFlipEigenVector( previousPointDouble, currentPointDouble, eigenVector ) ) {
        for( int cDim=0; cDim<3; cDim++ ) {
          eigenVector[cDim] *= -1;
        }
      }
      
      // get the gradients
      double gradients[3];
      this->GetGradient( gradientsVolume, currentPointInt, gradients );
  
      // this random number between 0 and 1 will be the percentage that the point
      // is moved along the eigenvector or toward the destination point
      const float pixelJump = 1;
      const float randomNumber = OpenRan1( &m_RandomTimeSeed ) * pixelJump;
      
      // get the next point in the path
      for( int cDim=0; cDim<3; cDim++ ) {
  
        currentPointDouble[cDim] = currentPointDouble[cDim] - randomNumber * gradients[cDim]
          + ( pixelJump - randomNumber ) * eigenVector[cDim];
        
      }
      
      // make sure that the point is within the bounds of the image
      this->EnsureWithinBounds( currentPointDouble );

      // add point to be saved
      double* point = new double[ 3 ];
      path.push_back( point );
      
      // copy the current point in
      for( int cDim=0; cDim<3; cDim++ ) {
        point[ cDim ] = currentPointDouble[ cDim ];
        
        previousPointDouble[cDim] = currentPointDouble[cDim];
        currentPointInt[ cDim ] = static_cast< int >( round( currentPointDouble[cDim] ) );
      }
      
    }
    
    // clean up the data
    this->FreeVector( destinationSeeds );
  
  }

  // copy the path to the matrix output
  this->CopyPathToOutput( path );
  this->FreeVector( path );
  
}

MRI* 
EigenVectorInitPathStrategy::GetDistanceVolume( const int label ) {
  
  // TODO: this might need to be changed, but I don't think letting it default
  // will guarentee something larger than the maximum distance all the time
  const int maxDistance = 999;

  // from mri_distance_transform  mode 1 is outside
  const int mode = 1;
  
  MRI* seedVolume = m_PoistatsModel->GetSeedVolume();
  
  // from fastmarching
  MRI* distanceTransform = MRIextractDistanceMap( seedVolume, NULL, label, maxDistance, mode);
  MRIcopyHeader( seedVolume, distanceTransform);
                          
  return distanceTransform;
}

void 
EigenVectorInitPathStrategy::GetGradient( MRI* gradientsVol, int *index, double *oGradient ) {

  for( int cDim=0; cDim<3; cDim++ ) {
    // get the gradient at the current point
    oGradient[ cDim ] = MRIgetVoxVal( gradientsVol, index[ 0 ], index[ 1 ], 
      index[ 2 ], cDim );
      
  }
  
}

void 
EigenVectorInitPathStrategy::GetEigenVector( int *index, double *oVector ) {

  for( int cDim=0; cDim<3; cDim++ ) {
    // get the eigenvector at the current point
    oVector[ cDim ] = MRIgetVoxVal( m_EigenVectors, index[ 0 ], index[ 1 ], 
      index[ 2 ], cDim );
  }
  
}

bool 
EigenVectorInitPathStrategy::ShouldFlipEigenVector( const double* const previousPoint, 
  const double* const currentPoint, const double* const eigenVector ) {
  
  bool shouldFlip = false;
  
  // the vector pointing from the previous to the current point
  double differenceVector[3];
  for( int cDim=0; cDim<3; cDim++ ) {
    differenceVector[cDim] = currentPoint[cDim] - previousPoint[cDim];
  }
  
  double dotProduct = 0.0;
  for( int cDim=0; cDim<3; cDim++ ) {
    dotProduct += differenceVector[ cDim ] * eigenVector[ cDim ];
  }
  
  // if the dot project is less than 0, then the flipped eigenVector is a better
  // fit -- it doesn't move move the pathway in a largely different angle.
  if( dotProduct < 0 ) {
    shouldFlip = true;
  }
  
  return shouldFlip;
  
}

void 
EigenVectorInitPathStrategy::DeleteGradientVolumes() {
  
  // go through all the gradient volumes and deallocate them
  for( std::vector< MRI* >::iterator it = m_GradientVolumes.begin(); it != m_GradientVolumes.end(); it++ ) {
    MRI* volume = ( *it );
    if( volume != NULL ) {
      MRIfree( &volume );
    }
  }
  
  m_GradientVolumes.clear();
  
}

