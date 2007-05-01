#include "PoistatsModel.h"

#include <iostream>
#include <cmath>

PoistatsModel::PoistatsModel() {

  std::srand( ( unsigned ) time( 0 ) );
  const long seed = std::rand();

  this->m_RandomNumberGenerator = vnl_random( seed );
  
  this->Init();  
}

PoistatsModel::PoistatsModel( const long seed ){
  this->m_RandomNumberGenerator = vnl_random( seed );
  this->Init();
}

PoistatsModel::~PoistatsModel() {  
  this->FreeVector( this->m_SeedValues );
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

// TODO: this is wrong...  The range will be between 0 and 1
double PoistatsModel::GetRandomNumberWithoutRange() {
  return m_RandomNumberGenerator.drand64();
}

bool PoistatsModel::IsDebug() {
  return m_IsDebug;
}

void PoistatsModel::SetDebug( const bool isDebug ) {
  m_IsDebug = isDebug;
}

bool PoistatsModel::IsUsingPathInitialization() {
  return m_IsUsingPathInitialization;
}

void PoistatsModel::SetUsingPathInitialization( const bool isUsing ) {
  m_IsUsingPathInitialization = isUsing;
}

MRI* 
PoistatsModel::GetEigenVectors() {
  return m_EigenVectors;
}

void 
PoistatsModel::SetEigenVectors( MRI* vectors ) {
  m_EigenVectors = vectors;
}
  
MRI*
PoistatsModel::GetSeedVolume() {
  return m_SeedVolume;
}

void 
PoistatsModel::SetSeedVolume( MRI* volume ) {
  m_SeedVolume = volume;
}
  
std::vector< int >* 
PoistatsModel::GetSeedValues() {
  return m_SeedValues;
}

void 
PoistatsModel::SetSeedValues( std::vector< int >* values ) {
  m_SeedValues = values;
}

void 
PoistatsModel::Init() {

  this->m_IsDebug = false;
  
  this->m_IsUsingPathInitialization = false;

  this->m_EigenVectors = NULL;
  this->m_SeedVolume = NULL;
  this->m_SeedValues = NULL;
}

void 
PoistatsModel::FreeVector( std::vector< int > *v ) {

  v->clear();
  delete v;
  
}
