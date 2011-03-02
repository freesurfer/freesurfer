/**
 * @file  PreferencesManager.cpp
 * @brief Assigns, reads, and writes preferences values
 *
 * This object manages the setting and accessing of preferences values
 * of different data types and writing and reading them to a file. Use
 * RegisterValue() to create a value and its default value, and
 * Set/GetValue() to access those values. UseFile() specifies a file
 * from which to read values and write them when done.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.18 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <iostream>
#include <fstream>
#include "string_fixed.h"
#include <list>
#include <stdio.h>
#include <stdexcept>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "PreferencesManager.h"

using namespace std;

PreferencesManager::PreferencesManager()
    : DebugReporter() {

  mVersion = 1;
}

PreferencesManager&
PreferencesManager::GetManager() {

  static PreferencesManager sManager;

  return sManager;
}


void
PreferencesManager::UseFile( string const ifnPrefs ) {

  string fnPrefs = ifnPrefs;
  struct stat info;
  int rStat;

  // If we're already using this file and it's not dirty, don't reread it.
  if ( mfnPrefs == fnPrefs &&
       mbPrefsFileDirty ) {
    return;
  }

  // If this is an absolute file name, try it.
  if ( fnPrefs[0] == '/' ) {
    rStat = stat( fnPrefs.c_str(), &info );

    // If that didn't work, at least make sure we can write it. If
    // not, return an error.
    if ( !S_ISREG(info.st_mode) ) {
      fstream testFile( fnPrefs.c_str(), ios::out );
      if ( !testFile.good() ) {
        stringstream ssError;
        ssError << "Cannot find or create specified preferences file: "
        << fnPrefs;
        throw runtime_error( ssError.str() );
      }
    }
  }


  // If this is not an absolute file name, try a few search paths.
  if ( fnPrefs[0] != '/' ) {

    // Local file first, then in home dir, then in usr/share.
    fnPrefs = "./" + ifnPrefs;

    info.st_mode = (mode_t) 0; // Keep valgrind from complaining.
    rStat = stat( fnPrefs.c_str(), &info );
    if ( !S_ISREG(info.st_mode) ) {
      char* homeDir = getenv("HOME");
      if ( NULL != homeDir ) {
        fnPrefs = string(homeDir) + "/" + ifnPrefs;
        rStat = stat( fnPrefs.c_str(), &info );
      }
    }

    if ( !S_ISREG(info.st_mode) ) {
      fnPrefs = "/usr/share" + ifnPrefs;
      rStat = stat( fnPrefs.c_str(), &info );
    }

    // If still nothing, we'll write it in their home dir.
    char* homeDir = getenv("HOME");
    if ( NULL != homeDir ) {
      fnPrefs = string(homeDir) + "/" + ifnPrefs;
    } else {
      fnPrefs = "./" +  ifnPrefs;
    }
  }


  mfnPrefs = fnPrefs;
  mbPrefsFileDirty = false;

  DebugOutput( << "Using prefs file " << mfnPrefs );

  ReadFile();
}

void
PreferencesManager::SetHeader( string const isHeader ) {

  msHeader = isHeader;
}

void
PreferencesManager::RegisterValue( string const isKeyName,
                                   string const isDescription,
                                   SimplePreferenceValue& iValue ) {

  PreferenceValueMap::iterator tPref = mPrefValues.find(isKeyName);
  if ( tPref != mPrefValues.end() ) {

    PreferenceValue pref = *mPrefValues[isKeyName];

  } else {

    PreferenceValue* pref = new PreferenceValue();
    pref->msKeyName = isKeyName;
    pref->msDescription = isDescription;
    pref->msDefaultValue = iValue.ValueToString();
    pref->msValue = iValue.ValueToString();

    mPrefValues[isKeyName] = pref;
  }
}

void
PreferencesManager::SetValue( std::string const isKeyName,
                              SimplePreferenceValue& iValue ) {

  PreferenceValueMap::iterator tPref = mPrefValues.find(isKeyName);
  if ( tPref != mPrefValues.end() ) {

    PreferenceValue* pref = mPrefValues[isKeyName];
    pref->msValue = iValue.ValueToString();

    mbPrefsFileDirty = true;

  } else {
    throw runtime_error( isKeyName + ": Value not found." );
  }
}

string
PreferencesManager::GetValue( std::string const isKeyName )  {

  PreferenceValueMap::iterator tPref = mPrefValues.find(isKeyName);
  if ( tPref != mPrefValues.end() ) {

    PreferenceValue* pref = mPrefValues[isKeyName];
    return pref->msValue;

  } else {
    throw runtime_error( isKeyName + ": Value not found." );
  }
}

void
PreferencesManager::Clear() {
  mPrefValues.clear();
}


void
PreferencesManager::ReadFile() {

  ifstream fPrefs( mfnPrefs.c_str(), ios::in );
  if ( !fPrefs || fPrefs.bad() ) {
    //    throw runtime_error("Can't open prefs file");
    return;
  }

  string sKeyword;
  while ( !fPrefs.eof() ) {
    getline( fPrefs, sKeyword );

    if ( sKeyword == "begin-version" ) {
      fPrefs >> mVersion;
      fPrefs >> sKeyword;
      if ( sKeyword != "end-version" ) {
        stringstream sError;
        sError << "Bad prefs file: expected end-version, got " << sKeyword;
        throw runtime_error( sError.str() );
      }
      DebugOutput( << "Reading prefs file version " << mVersion );

    } else if ( sKeyword == "begin-header" ) {
      stringstream sHeader;
      while ( !fPrefs.eof() ) {
        getline( fPrefs, sKeyword );
        if ( sKeyword == "end-header" ) {
          msHeader = sHeader.str();
          break;
        } else {
          sHeader << sKeyword;
        }
      }

      /* Just read in but ignore the timestamp. */
    } else if ( sKeyword == "begin-timestamp" ) {
      while ( !fPrefs.eof() ) {
        getline( fPrefs, sKeyword );
        if ( sKeyword == "end-timestamp" ) {
          break;
        }
      }

    } else if ( sKeyword == "begin-pref" ) {

      PreferenceValue* pref = new PreferenceValue();

      while ( !fPrefs.eof() ) {
        getline( fPrefs, sKeyword );

        if ( sKeyword == "begin-description" ) {
          stringstream sDescription;
          while ( !fPrefs.eof() ) {
            getline( fPrefs, sKeyword );
            if ( sKeyword == "end-description" ) {
              pref->msDescription = sDescription.str();
              break;
            } else {
              sDescription << sKeyword;
            }
          }

        } else if ( sKeyword == "begin-name" ) {
          fPrefs >> pref->msKeyName;
          fPrefs >> sKeyword;
          if ( sKeyword != "end-name" ) {
            stringstream sError;
            sError << "Bad prefs file: expected end-name, got " << sKeyword;
            throw runtime_error( sError.str() );
          }

        } else if ( sKeyword == "begin-value" ) {
          stringstream sValue;
          while ( !fPrefs.eof() ) {
            getline( fPrefs, sKeyword );
            if ( sKeyword == "end-value" ) {
              pref->msValue = sValue.str();
              break;
            } else {
              sValue << sKeyword;
            }
          }

        } else if ( sKeyword == "end-pref" ) {
          break;
        }
      }

      mPrefValues[pref->msKeyName] = pref;
    }
  }

  fPrefs.close();
}

void
PreferencesManager::WriteFile() {

  ofstream fPrefs( mfnPrefs.c_str(), ios::out );
  if ( fPrefs.bad() ) {
    throw runtime_error( "Can't open prefs file" );
  }

  time_t curTime;
  time( &curTime );
  fPrefs << "begin-version" << endl;
  fPrefs << mVersion << endl;
  fPrefs << "end-version" << endl << endl;

  fPrefs << "begin-timestamp" << endl;
  fPrefs << "# Scuba preferences file written " << ctime(&curTime);
  fPrefs << "# $Id: PreferencesManager.cpp,v 1.18 2011/03/02 00:04:36 nicks Exp $" << endl;
  fPrefs << "end-timestamp" << endl << endl;

  fPrefs << "begin-header" << endl;
  fPrefs << msHeader << endl << endl;
  fPrefs << "end-header" << endl << endl;

  PreferenceValueMap::iterator tPref;
  for ( tPref = mPrefValues.begin(); tPref != mPrefValues.end(); ++tPref ) {

    PreferenceValue* pref = (*tPref).second;

    fPrefs << "begin-pref" << endl;
    fPrefs << "begin-description" << endl
    << pref->msDescription << endl << "end-description" << endl;
    fPrefs << "begin-name" << endl
    << pref->msKeyName << endl << "end-name" << endl;
    fPrefs << "begin-value" << endl
    << pref->msValue << endl << "end-value" << endl;
    fPrefs << "end-pref" << endl << endl;
  }

  fPrefs.close();

  mbPrefsFileDirty = false;
}





PreferencesManager::IntPrefValue::IntPrefValue( int const i ) {
  mValue = i;
}

PreferencesManager::IntPrefValue::IntPrefValue( string const i ) {
  mValue = atoi(i.c_str());
  mValue = strtol(i.c_str(), (char**)NULL, 10);
}

void
PreferencesManager::IntPrefValue::SetFromString( string const isValue ) {
  mValue = atoi(isValue.c_str());
}

string
PreferencesManager::IntPrefValue::ValueToString() {
  char sValue[256] = "";
  sprintf( sValue, "%d", mValue );
  return string(sValue);
}

int
PreferencesManager::IntPrefValue::GetValue() const {
  return mValue;
}



PreferencesManager::FloatPrefValue::FloatPrefValue( float const i ) {
  mValue = i;
}

PreferencesManager::FloatPrefValue::FloatPrefValue( string const i ) {
  mValue = atof(i.c_str());
}

void
PreferencesManager::FloatPrefValue::SetFromString( string const isValue ) {
  mValue = atof(isValue.c_str());
}

string
PreferencesManager::FloatPrefValue::ValueToString() {
  char sValue[256] = "";
  sprintf( sValue, "%f", mValue );
  return string(sValue);
}

float
PreferencesManager::FloatPrefValue::GetValue() const {
  return mValue;
}




PreferencesManager::StringPrefValue::StringPrefValue( string const i ) {
  msValue = i;
}

void
PreferencesManager::StringPrefValue::SetFromString( string const isValue ) {
  msValue = isValue;
}

string
PreferencesManager::StringPrefValue::ValueToString() {
  return msValue;
}

string
PreferencesManager::StringPrefValue::GetValue() const {
  return msValue;
}

