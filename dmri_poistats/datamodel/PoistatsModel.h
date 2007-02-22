#ifndef __PoistatsModel_h
#define __PoistatsModel_h

#include <vnl/vnl_random.h>

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
  
private:

  vnl_random m_RandomNumberGenerator;
  
  bool m_IsDebug;

};

#endif

