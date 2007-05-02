#ifndef __InitializePath_h
#define __InitializePath_h

extern "C" {
#include "mri.h"
}

#include <vector>

#include <itkArray2D.h>

#include "dmri_poistats/datamodel/PoistatsModel.h"

class InitializePath
{
public:

  InitializePath();
  InitializePath( PoistatsModel* model );
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
   * Calculate and cache the gradient volumes.  This will save time when
   * creating multiple initial paths for replicas.
   */
  void CacheGradientVolumes();

  /**
   * This does the work of finding the initial path.
   */
  void CalculateInitialPath();
  
  /**
   * Returns the intial path that was found.
   */
  itk::Array2D< double >* GetInitialPath();  

private:

  PoistatsModel *m_PoistatsModel;

  MRI *m_EigenVectors;

  MRI *m_SeedVolume;
  
  std::vector< int > *m_SeedValues; 
  
  /** Initial path in itk 2d array form */
  itk::Array2D< double > *m_InitialPath;
  
  long m_RandomTimeSeed;
  
  /**
   * Vector of gradient volumes.  The gradients will be one less than the
   * number of seed values.  The will be calculated from the distances
   * transform of the second seed value, third, etc.
   **/
  std::vector< MRI * > m_GradientVolumes;
  
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
  
  /**
   * Gets the gradient and returns it in oGradient.
   */
  void GetGradient( MRI* gradientsVol, int *index, double *oGradient );

  /**
   * Gets the eigenvector at index and returns it in oVector.
   */
  void GetEigenVector( int *index, double *oVector );
  
  /**
   * Determines if the eigenvector should be flipped (because they're symmetric)
   * when using it to find the next position along the initial path.
   */
  bool ShouldFlipEigenVector( const double* const previousPoint, const double* const currentPoint, const double* const eigenVector );
  
  /**
   * Deallocates memory of the gradient volumes.
   */
  void DeleteGradientVolumes();
  
  /**
   * Changes the point if it is out of bounds to be at the extent of the volume.
   */
  void EnsureWithinBounds( double* point );
};

#endif
