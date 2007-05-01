#ifndef __PoistatsModel_h
#define __PoistatsModel_h

#include <vnl/vnl_random.h>
#include <vector>

extern "C" {
#include "mri.h"
}

class PoistatsModel
{
public:

  PoistatsModel();
  PoistatsModel( const long seed );
  ~PoistatsModel();

  void SetRandomSeed( const long seed );
      
  double GetRandomNumber();
  
  double GetNormallyDistributedRandomNumber();

  int GetRandomInt( const int floor, const int ceiling );
  
  double GetRandomNumberWithoutRange();
  
  bool IsDebug();

  void SetDebug( const bool isDebug );
  
  bool IsUsingPathInitialization();

  void SetUsingPathInitialization( const bool isUsing );
  
  /**
   * Returns the eigenvectors.
   */
  MRI* GetEigenVectors();
  void SetEigenVectors( MRI* vectors );   
  
  /**
   * Returns the volume with the seed labels.
   */
  MRI* GetSeedVolume();
  void SetSeedVolume( MRI* volume );   
  
  /**
   * Returns the seed values to be used.
   */
  std::vector< int >* GetSeedValues();
  void SetSeedValues( std::vector< int >* values);
    
private:

  vnl_random m_RandomNumberGenerator;
  
  bool m_IsDebug;
  
  bool m_IsUsingPathInitialization;
  
  MRI* m_EigenVectors;
  MRI* m_SeedVolume;
  std::vector< int >* m_SeedValues;
  
  /**
   * Initializes the class.
   */
  void Init();

};

#endif

