//
// PreferencesManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2003/10/09 23:07:15 $
// Revision       : $Revision: 1.1 $

#ifndef PreferencesManager_h
#define PreferencesManager_h

#include <stdlib.h>
#include <map>
#include <string>

class PreferencesManager {
 public:

  static PreferencesManager& GetManager();

  void UseFile( std::string ifnPrefs );

  void SetHeader( std::string isHeader );

  template <class T> void RegisterValue( std::string isKeyName,
					 std::string isDescription,
					 T iDefaultValue );

  template <class T> T GetValue( std::string isKeyName );

  template <class T> void SetValue( std::string isKeyName, T iValue );
    
 protected:
  PreferencesManager();

  class PrefsValue {
    friend class PreferencesManager;
    std::string msKeyName;
    std::string msDescription;
    std::string msDefaultValue;
    std::string msValue;
  };

  std::string mfnPrefs;
  std::string msHeader;
  std::map<std::string,PrefsValue> mPrefValues;
};

#endif
