#include <stdexcept>
#include "ProgressDisplayManager.h"

ProgressDisplayManager* ProgressDisplayManager::sManager = NULL;

void
ProgressDisplayManager::SetManager ( ProgressDisplayManager* iManager ) {
  
  sManager = iManager;
}


ProgressDisplayManager& 
ProgressDisplayManager::GetManager() {

  if( NULL == sManager ) {
    throw std::runtime_error( "No manager set yet." );
  }

  return *sManager;
}

