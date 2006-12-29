/**
 * @file  PreferencesManager.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


//
// PreferencesManager.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/12/29 02:09:14 $
// Revision       : $Revision: 1.6 $

#ifndef PreferencesManager_h
#define PreferencesManager_h

#include <stdlib.h>
#include <map>
#include "string_fixed.h"
#include "DebugReporter.h"

class PreferencesManager : public DebugReporter {
public:

  // An interface for a preference value. The PreferencesManager needs
  // string accessors to get and set the value. You can implement a
  // subclass for this however you'd like. Three simple versions are
  // included in PreferencesManager.
  class SimplePreferenceValue {
  public:
    SimplePreferenceValue() {}
    virtual ~SimplePreferenceValue() {}
    virtual void SetFromString( std::string const ) = 0;
    virtual std::string ValueToString() = 0;
  };


  // Static manager accessor.
  static PreferencesManager& GetManager();

  // Specify a file name. If a relative file name, will look in the
  // normal search path to find a file. If found, this will read and
  // parse the file.
  void UseFile( std::string const ifnPrefs );

  std::string GetFileName() {
    return mfnPrefs;
  }

  // Sets the header description that will go at the top of the prefs
  // file.
  void SetHeader( std::string const isHeader );


  // Register a preferences value. If not already in the file, the
  // value will be set as the default value. If already found in the
  // file, nothing will be changed.
  void RegisterValue( std::string const isKeyName,
                      std::string const isDescription,
                      SimplePreferenceValue& iValue );

  // Sets a preference value.
  void SetValue( std::string const isKeyName,
                 SimplePreferenceValue& iValue );

  // Gets a preference value.
  std::string GetValue( std::string const isKeyName );


  void Clear();


  void SaveFile() {
    WriteFile();
  }


  // These are three simple subclasses for SimplePreferenceValue that
  // implement int, float, and string values.
class IntPrefValue : public SimplePreferenceValue {
  public:
    IntPrefValue( int const i );
    IntPrefValue( std::string const i );
    virtual ~IntPrefValue() {}
    void SetFromString( std::string const isValue );
    std::string ValueToString();
    int GetValue() const;
  protected:
    int mValue;
  };

class FloatPrefValue : public SimplePreferenceValue {
  public:
    FloatPrefValue( float const i );
    FloatPrefValue( std::string const i );
    virtual ~FloatPrefValue() {}
    void SetFromString( std::string const isValue );
    std::string ValueToString();
    float GetValue() const;
  protected:
    float mValue;
  };

class StringPrefValue : public SimplePreferenceValue {
  public:
    StringPrefValue( std::string const i );
    virtual ~StringPrefValue() {}
    void SetFromString( std::string const isValue );
    std::string ValueToString();
    std::string GetValue() const;
  protected:
    std::string msValue;
  };


protected:
  PreferencesManager();

  void ReadFile();
  void WriteFile();

  class PreferenceValue {
    friend class PreferencesManager;
    std::string msKeyName;
    std::string msDescription;
    std::string msDefaultValue;
    std::string msValue;
  };

  int mVersion;
  std::string mfnPrefs;
  std::string msHeader;

  bool mbPrefsFileDirty;

  typedef std::map<std::string,PreferenceValue*> PreferenceValueMap;
  PreferenceValueMap mPrefValues;
};

#endif
