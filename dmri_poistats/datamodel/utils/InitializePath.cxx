#include "InitializePath.h"

// for MRIextractDistanceMap
#include "fastmarching.h"

// for random number generator
#include "numerics.h"

#include <iostream>

InitializePath::InitializePath() {
  m_PoistatsModel = NULL;
  m_SeedVolume = NULL;
  m_SeedValues = NULL;
  m_InitialPath = NULL;
  
  // this will be the seed value we use for getting random numbers
  std::srand( ( unsigned ) time( 0 ) );
  m_RandomTimeSeed = std::rand();

}

InitializePath::InitializePath( PoistatsModel *model ) {

  m_PoistatsModel = model;

  this->SetSeedVolume( m_PoistatsModel->GetSeedVolume() );
  this->SetSeedValues( m_PoistatsModel->GetSeedValues() );
  
  m_InitialPath = NULL;

  // this will be the seed value we use for getting random numbers
  std::srand( ( unsigned ) time( 0 ) );
  m_RandomTimeSeed = std::rand();

}

InitializePath::~InitializePath() {
  
  if( m_InitialPath != NULL ) {
    delete m_InitialPath;
  }
  
}

void 
InitializePath::SetSeedVolume( MRI *volume ) {
  m_SeedVolume = volume;
}

void 
InitializePath::SetSeedValues( std::vector< int > *values ) {
  m_SeedValues = values;
}

itk::Array2D< double >* 
InitializePath::GetInitialPath() {
  return m_InitialPath;
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

  // copy the output to the itk matrix
  if( m_InitialPath == NULL ) {
    m_InitialPath = new itk::Array2D< double >;
  } else {
    // if there's already something there, erase it
    m_InitialPath->SetSize( 0, 0 );
  }
  
  m_InitialPath->SetSize( path.size(), 3 );    
    
  int row = 0;
  // copy the elements
  for( std::vector< double* >::iterator it = path.begin(); it != path.end(); it++, row++ ) {
    for( int col=0; col<3; col++ ) {      
      ( *m_InitialPath )( row, col ) = ( *it )[ col ];
    }
  }
  
}

void
InitializePath::EnsureWithinBounds( double* point ) {

  // bounds of the volume
  const double bounds[] = {
    m_SeedVolume->width - 1,
    m_SeedVolume->height - 1,
    m_SeedVolume->depth - 1
  };
  
  const double minimum = 0.0;
  
  // go through all the dimensions and make sure the point is in the bounds
  for( int cDim=0; cDim<3; cDim++ ) {
    
    // if it is less than zero, then it's out of bounds
    if( point[ cDim ] < minimum ) {
      point[ cDim ] = minimum;
    } else if( point[ cDim ] > bounds[ cDim ] ) {
      // if greater than the bounds of this dimension, then out of bounds
      point[ cDim ] = bounds[ cDim ];
    }
    
  }
  
}
