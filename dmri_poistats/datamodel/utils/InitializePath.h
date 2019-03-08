#ifndef __InitializePath_h
#define __InitializePath_h

 
#include "mri.h"


#include <vector>

#include <itkArray2D.h>

#include "datamodel/PoistatsModel.h" // dmri_poistats

class InitializePath
{
public:

  InitializePath();
  InitializePath( PoistatsModel* model );
  virtual ~InitializePath();
  
  /**
   * The seed volume that has the seed regions to be connected.
   */
  virtual void SetSeedVolume( MRI *volume );
  
  /**
   * These are the initial seed values in the seed volume that are to be 
   * connected.
   */
  virtual void SetSeedValues( std::vector< int > *values );

  /**
   * This does the work of finding the initial path.
   */
  virtual void CalculateInitialPath() = 0;
  
  /**
   * Returns the intial path that was found.
   */
  itk::Array2D< double >* GetInitialPath();  
  
protected:

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
   * Changes the point if it is out of bounds to be at the extent of the volume.
   */
  void EnsureWithinBounds( double* point );

  PoistatsModel *m_PoistatsModel;

  MRI *m_SeedVolume;
  
  std::vector< int > *m_SeedValues; 
  
  /** Initial path in itk 2d array form */
  itk::Array2D< double > *m_InitialPath;
  
  long m_RandomTimeSeed;

private:

    
};

#endif
