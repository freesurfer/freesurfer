#include "PoistatsModel.h"

#include <iostream>
#include <cmath>

PoistatsModel::PoistatsModel() {
  std::srand( ( unsigned ) time( 0 ) );
  const long seed = std::rand();

  this->m_RandomNumberGenerator = vnl_random( seed );
}

PoistatsModel::PoistatsModel( const long seed ){
  this->m_RandomNumberGenerator = vnl_random( seed );
}

PoistatsModel::~PoistatsModel() {
  
}

void PoistatsModel::SetRandomSeed( const long seed ) {
  this->m_RandomNumberGenerator.reseed( seed );
}

double PoistatsModel::GetRandomNumber() {

  const double rangeMin = 0.0;
  const double rangeMax = 1.0;
  
  const double randomNumber = 
    this->m_RandomNumberGenerator.drand64( rangeMin, rangeMax );

  return randomNumber;
}

double PoistatsModel::GetNormallyDistributedRandomNumber() {

  const double randomNumber = 
    this->m_RandomNumberGenerator.normal64();
  return randomNumber;
}

int PoistatsModel::GetRandomInt( const int floor, const int ceiling ) {
  int randomInt  = this->m_RandomNumberGenerator.lrand32( floor, ceiling );
  
  return randomInt;
}

double PoistatsModel::GetRandomNumberWithoutRange() {
  return m_RandomNumberGenerator.drand64();
}

bool PoistatsModel::IsDebug() {
  return m_IsDebug;
}

void PoistatsModel::SetDebug( const bool isDebug ) {
  m_IsDebug = isDebug;
}
