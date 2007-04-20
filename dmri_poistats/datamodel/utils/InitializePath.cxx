#include "InitializePath.h"

// for MRIextractDistanceMap
#include "fastmarching.h"

// for random number generator
#include "numerics.h"

#include <iostream>

InitializePath::InitializePath() {
  m_EigenVectors = NULL;
  m_SeedVolume = NULL;
  m_SeedValues = NULL;
  m_InitialPath = NULL;
  
  // this will be the seed value we use for getting random numbers
  std::srand( ( unsigned ) time( 0 ) );
  m_RandomTimeSeed = std::rand();
  
}

InitializePath::~InitializePath() {

  if( m_InitialPath != NULL ) {
      MatrixFree( &m_InitialPath );
      m_InitialPath = NULL;
  }
  
}

void 
InitializePath::SetEigenVectors( MRI *eigenVectors ) {
  m_EigenVectors = eigenVectors;
}

void 
InitializePath::SetSeedVolume( MRI *volume ) {
  m_SeedVolume = volume;
}

void 
InitializePath::SetSeedValues( std::vector< int > *values ) {
  m_SeedValues = values;
}

void 
InitializePath::CalculateInitialPath() {
  
  std::cerr << "CalculateInitialPath" << std::endl;

  // TODO: using the second seed as the destination
  const int destinationSeedLabel = ( *m_SeedValues )[ 1 ];
  
  // create the distance transform, all the values should be set before running
  MRI *distanceTransform = this->GetDistanceVolume( destinationSeedLabel );
  
  MRIwrite( distanceTransform, "/autofs/space/heraclitus_001/users/dsjen/data/InitPath/DistanceTransform.mgz" );
  
  // now that we have the transform to get to the second seed, we want to take
  // the derivatives of the distance transform to see what direction we need
  // to take to get to the second seed
  MRI *gradients = MRIsobel( distanceTransform, NULL, NULL );
  
  // pick a random starting point within the starting seed region
  const int startSeedLabel = ( *m_SeedValues )[ 0 ];
  int currentPoint[3];
  this->GetRandomSeedPoint( currentPoint, startSeedLabel );
  
  // starting at our random seed point, we want to move it toward the second
  // seed region and also following the eigenvector
  
  std::vector< int* > destinationSeeds = this->GetSeedPoints( destinationSeedLabel );
  
  // this will be a vector with 3 columns
  std::vector< int* > path;
  
  // keep iterating until you've reached your destination region
  while( !this->IsInRegion( currentPoint, destinationSeeds ) ) {

  // TODO: remove these lines
//  for( int i=0; i<100; i++ ) {
//    std::cerr << "  i: " << i << std::endl;
        
    // this random number between 0 and 1 will be the percentage that the point
    // is moved along the eigenvector or toward the destination point
    const float pixelJump = 1;
    const float randomNumber = OpenRan1( &m_RandomTimeSeed ) * pixelJump;
    for( int cDim=0; cDim<3; cDim++ ) {
      
      // TODO: you have to make sure that the dimensions are correct here
      
      // get the eigenvector at the current point
      const float eigenVector = MRIFseq_vox( m_EigenVectors, currentPoint[ 0 ], 
        currentPoint[ 1 ], currentPoint[ 2 ], cDim );

      // get the gradient at the current point
      const float gradient = MRIFseq_vox( gradients, currentPoint[ 0 ], 
        currentPoint[ 1 ], currentPoint[ 2 ], cDim );
        
      const float newPoint = currentPoint[ cDim ] + randomNumber * -gradient
        + ( pixelJump - randomNumber ) * eigenVector;
        
      // TODO: maybe I should be storing doubles instead of ints...
      currentPoint[ cDim ] = static_cast< int >( round( newPoint ) );
      
    }
                
    // add point to be saved
    int* point = new int[ 3 ];
    path.push_back( point );
    
    // copy the current point in
    for( int cDim=0; cDim<3; cDim++ ) {
      point[ cDim ] = currentPoint[ cDim ];
      
      std::cerr << point[ cDim ] << "  ";
    }
    std::cerr << std::endl;
    
  }
  
  std::cerr << "  clean up path data" << std::endl;
    
  // clean up the data
  this->FreeVector( destinationSeeds );
  this->FreeVector( path );
  
  std::cerr << "  free gradient" << std::endl;
  if( gradients != NULL ) {
    MRIfree( &gradients );
  }

  std::cerr << "  free distance transform" << std::endl;
  if( distanceTransform != NULL ) {
    MRIfree( &distanceTransform );
  }
  
  std::cerr << "end CalculateInitialPath" << std::endl;
  
}

MATRIX* 
InitializePath::GetInitialPath() {
  return m_InitialPath;
}

MRI* 
InitializePath::GetDistanceVolume( const int label ) {
  
  // TODO: this might need to be changed, but this is the default in mri_distance_transform
  const int maxDistance = 10;

  // from mri_distance_transform  mode : 1 = outside , mode : 2 = inside , mode : 3 = both
  const int mode = 1;
  
  // from fastmarching -- doesn't work properly
  // MRI* distanceTransform = MRIextractDistanceMap( m_SeedVolume, NULL, label, maxDistance, mode);

  // this is from mri.h and seems to work
  MRI* distanceTransform = MRIdistanceTransform( m_SeedVolume, NULL, label, maxDistance, mode );
                          
  return distanceTransform;
}

void 
InitializePath::GetRandomSeedPoint( int *startSeed, const int label ) {
  
  std::vector< int* > seeds = this->GetSeedPoints( label );
    
  // pick a random index
  // this will be a number between 0 and 1
  const float randomNumber = OpenRan1( &m_RandomTimeSeed );  
  const int randomIndex = static_cast< int >( floor( randomNumber * seeds.size() ) );

  // return the random index  
  for( int i=0; i<3; i++ ) {
    startSeed[ i ] = seeds[ randomIndex ][ i ];
  }
  
  // clean up the data
  this->FreeVector( seeds );
  
}

std::vector< int* > 
InitializePath::GetSeedPoints( const int label ) {

  // the seed points of a label
  std::vector< int* > seeds;
  
  // get the indices with this label
  for( int x = 0; x < m_SeedVolume->width; x++ ) {
    for( int y = 0; y < m_SeedVolume->height; y++ ) {
      for( int z = 0; z < m_SeedVolume->depth; z++ ) {
        
        // get the current value in the volume
        const float currentVoxel = MRIFseq_vox( m_SeedVolume, x, y, z, 0 );
                
        // if the value in the volume is the same as the label we're looking
        // for, then save the location
        if( label == currentVoxel ) {
          int *seed = new int[ 3 ];
          seed[ 0 ] = x;
          seed[ 1 ] = y;
          seed[ 2 ] = z;          
          seeds.push_back( seed );
        }
        
      }
    }
  }
  
  return seeds;
  
}

bool 
InitializePath::IsInRegion( int *point, std::vector< int* > region ) {

  bool isInRegion = false;

  // go through all the seeds and determine if point is in the region
  for( std::vector< int* >::iterator itRegion = region.begin(); !isInRegion && itRegion != region.end(); itRegion++ ) {
    
    // is all dimensions are the same, then we've got a match and can quit
    bool isSame = true;
    for( int i=0; i<3; i++ ) {
      // if the point differs, then we need to keep searching
      if( point[ i ] != ( *itRegion )[ i ] ) {
        isSame = false;
      }
    }
    
    // if we went through all the dimensions and they're the same, then it's
    // a match
    if( isSame ) {
     isInRegion = true; 
    }
    
  }
  
  return isInRegion;
  
}

void 
InitializePath::FreeVector( std::vector< int* > v ) {

  // clean up the data
  for( std::vector< int* >::iterator it = v.begin(); it != v.end(); it++ ) {
    delete[] ( *it );
  }
  v.clear();
  
}
