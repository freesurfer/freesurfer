#ifndef __InitializePath_h
#define __InitializePath_h

extern "C" {
#include "mri.h"
#include "matrix.h"
}

#include <vector>

class InitializePath
{
public:

  InitializePath();
  ~InitializePath();
  
  /**
   * The eigenvectors that will be used in determining an initialization.
   */
  void SetEigenVectors( MRI *eigenVectors );

  /**
   * The seed volume that has the seed regions to be connected.
   */
  void SetSeedVolume( MRI *volume );
  
  /**
   * These are the initial seed values in the seed volume that are to be 
   * connected.
   */
  void SetSeedValues( std::vector< int > *values );

  /**
   * This does the work of finding the initial path.
   */
  void CalculateInitialPath();
  
  /**
   * Returns the intial path that was found.
   */
  MATRIX* GetInitialPath();  

private:

  MRI *m_EigenVectors;

  MRI *m_SeedVolume;
  
  std::vector< int > *m_SeedValues; 
  
  MATRIX *m_InitialPath;
  
  long m_RandomTimeSeed;
  
  /**
   * Gets the distance transform for label.
   */
  MRI* GetDistanceVolume( const int label ); 
  
  /**
   * Get a random seed point at index location at a label.
   */
  void GetRandomSeedPoint( int *startSeed, const int label );
  
  /**
   * Returns the image coordinates of the seed points for a label.
   */
  std::vector< int* > GetSeedPoints( const int label );
  
  /**
   * Returns true if the point is within the region, consisting of points
   */
  bool IsInRegion( int *point, std::vector< int* > region );

  /**
   * Frees the memory stored in the vector and clears the vector.
   */ 
  void FreeVector( std::vector< int* > v );

  void FreeVector( std::vector< double* > v );
  
  void CopyPathToOutput( std::vector< double* > path );
  
};

#endif
