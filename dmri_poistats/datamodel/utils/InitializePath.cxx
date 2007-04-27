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
  
  // TODO: using the second seed as the destination
  const int destinationSeedLabel = ( *m_SeedValues )[ 1 ];
  
  // create the distance transform, all the values should be set before running
  MRI *distanceTransform = this->GetDistanceVolume( destinationSeedLabel );
  
  // now that we have the transform to get to the second seed, we want to take
  // the derivatives of the distance transform to see what direction we need
  // to take to get to the second seed
  MRI *gradientsVolume = MRIsobel( distanceTransform, NULL, NULL );
  
  // pick a random starting point within the starting seed region
  const int startSeedLabel = ( *m_SeedValues )[ 0 ];
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
    
  // starting at our random seed point, we want to move it toward the second
  // seed region and also following the eigenvector  
  std::vector< int* > destinationSeeds = this->GetSeedPoints( destinationSeedLabel );
  
  // this will be a vector with 3 columns
  std::vector< double* > path;

  std::cerr << "starting..." << std::endl;
  
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
        + (pixelJump - randomNumber ) * eigenVector[cDim];

      // TODO: this is only to follow the princliple eigendirection
//      currentPointDouble[cDim] = previousPointDouble[cDim] + eigenVector;
      
    }
                
    // add point to be saved
    double* point = new double[ 3 ];
    path.push_back( point );
    
    // copy the current point in
    for( int cDim=0; cDim<3; cDim++ ) {
      point[ cDim ] = currentPointDouble[ cDim ];
      
      previousPointDouble[cDim] = currentPointDouble[cDim];
      currentPointInt[ cDim ] = static_cast< int >( round( currentPointDouble[cDim] ) );
      std::cerr << point[ cDim ] << "  ";
    }
    std::cerr << std::endl;
        
  }
  
  // copy the path to the matrix output
  this->CopyPathToOutput( path );
  
  // clean up the data
  this->FreeVector( destinationSeeds );
  this->FreeVector( path );
  
  if( gradientsVolume != NULL ) {
    MRIfree( &gradientsVolume );
  }

  if( distanceTransform != NULL ) {
    MRIfree( &distanceTransform );
  }
  
}

MATRIX* 
InitializePath::GetInitialPath() {
  return m_InitialPath;
}

MRI* 
InitializePath::GetDistanceVolume( const int label ) {
  
  // TODO: this might need to be changed, but this is the default in mri_distance_transform
  const int maxDistance = 999;

  // from mri_distance_transform  mode : 1 = outside , mode : 2 = inside , mode : 3 = both
  const int mode = 1;
  
  // from fastmarching
  MRI* distanceTransform = MRIextractDistanceMap( m_SeedVolume, NULL, label, maxDistance, mode);
  MRIcopyHeader( m_SeedVolume, distanceTransform);
                          
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

void 
InitializePath::FreeVector( std::vector< double* > v ) {

  // clean up the data
  for( std::vector< double* >::iterator it = v.begin(); it != v.end(); it++ ) {
    delete[] ( *it );
  }
  v.clear();
  
}

void 
InitializePath::CopyPathToOutput( std::vector< double* > path ) {

  // just in case something was in it
  if( m_InitialPath != NULL ) {
      MatrixFree( &m_InitialPath );
      m_InitialPath = NULL;
  }
  
  // allocate the matrix
  m_InitialPath = MatrixAlloc( path.size(), 3, MATRIX_REAL );
  
  int row = 0;
  // copy the elements
  for( std::vector< double* >::iterator it = path.begin(); it != path.end(); it++, row++ ) {
    for( int col=0; col<3; col++ ) {      
      m_InitialPath->rptr[ row+1 ][ col+1 ] = ( *it )[ col ];
    }
  }
  
}

void 
InitializePath::GetGradient( MRI* gradientsVol, int *index, double *oGradient ) {

  for( int cDim=0; cDim<3; cDim++ ) {
    // get the gradient at the current point
    oGradient[ cDim ] = MRIgetVoxVal( gradientsVol, index[ 0 ], index[ 1 ], 
      index[ 2 ], cDim );
      
  }
  
}

void 
InitializePath::GetEigenVector( int *index, double *oVector ) {

  for( int cDim=0; cDim<3; cDim++ ) {
    // get the eigenvector at the current point
    oVector[ cDim ] = MRIgetVoxVal( m_EigenVectors, index[ 0 ], index[ 1 ], 
      index[ 2 ], cDim );
  }
  
}

bool 
InitializePath::ShouldFlipEigenVector( const double* const previousPoint, 
  const double* const currentPoint, const double* const eigenVector ) {
  
  bool shouldFlip = false;
  
  // the vector pointing from the previous to the current point
  double differenceVector[3];
  for( int cDim=0; cDim<3; cDim++ ) {
    differenceVector[cDim] = currentPoint[cDim] - previousPoint[cDim];
  }
  
  double dotProject = 0.0;
  for( int cDim=0; cDim<3; cDim++ ) {
    dotProject += differenceVector[ cDim ] * eigenVector[ cDim ];
  }
  
  // if the dot project is less than 0, then the flipped eigenVector is a better
  // fit -- it doesn't move move the pathway in a largely different angle.
  if( dotProject < 0 ) {
    shouldFlip = true;
  }
  
  return shouldFlip;
  
}
