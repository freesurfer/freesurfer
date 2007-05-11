#ifndef __EigenVectorInitPathStrategy_h
#define __EigenVectorInitPathStrategy_h

#include "InitializePath.h"

class EigenVectorInitPathStrategy : public InitializePath
{
public:

  EigenVectorInitPathStrategy();
  EigenVectorInitPathStrategy( PoistatsModel* model );
  ~EigenVectorInitPathStrategy();  

  /**
   * This does the work of finding the initial path.
   */
  void CalculateInitialPath();

  /**
   * The eigenvectors that will be used in determining an initialization.
   */
  void SetEigenVectors( MRI *eigenVectors );

protected:

  /**
   * Calculate and cache the gradient volumes.  This will save time when
   * creating multiple initial paths for replicas.
   */
  void CacheGradientVolumes();

  /**
   * Gets the distance transform for label.
   */
  MRI* GetDistanceVolume( const int label ); 

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
  
private:

  MRI *m_EigenVectors;

  bool m_IsGradientCached;

  /**
   * Vector of gradient volumes.  The gradients will be one less than the
   * number of seed values.  The will be calculated from the distances
   * transform of the second seed value, third, etc.
   **/
  std::vector< MRI * > m_GradientVolumes;

};

#endif 
