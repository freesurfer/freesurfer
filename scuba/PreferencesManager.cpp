#include <iostream>
#include <string>
#include <list>
#include "DataManager.h"

using namespace std;

PreferencesManager::PreferencesManager() {

}

PreferencesManager& 
PreferencesManager::GetManager() {

  static PreferencesManager sManager;

  return sManager;
}


void 
PreferencesManager::UseFile( std::string ifnPrefs ) {
  
}

void 
PreferencesManager::SetHeader( std::string isHeader ) {

}

template <class T> void 
PreferencesManager::RegisterValue( std::string isKeyName,
				   std::string isDescription,
				   T iDefaultValue ) {

}

template <class T> T 
PreferencesManager::GetValue( std::string isKeyName ) {

}

template <class T> void 
PreferencesManager::SetValue( std::string isKeyName, T iValue ) {

}
    
