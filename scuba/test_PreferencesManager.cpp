/**
 * @file  test_PreferencesManager.cpp
 * @brief test PreferencesManager class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.11 $
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

#include <stdlib.h>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "PreferencesManager.h"
#include "Scuba-impl.h"

#define AssertS(x,s)   if(!(x)) { throw (char const*) s; }
#define Assert(x,e)   if(!(x)) { throw e; }

using namespace std;

const char* Progname = "test_PreferencesManager";


// Try to register a value.
void TestRegister( PreferencesManager::SimplePreferenceValue& iValue,
                   string isKeyName,
                   string isError ) {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.RegisterValue( isKeyName, "A sample value.", iValue );
}

// Try to write a value.
void TestW( PreferencesManager::SimplePreferenceValue& iValue,
            string                isKeyName,
            string                isError ) {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  prefsMgr.SetValue( isKeyName, iValue );
}

// Try to read a value and compare it to the argument.
void TestR( PreferencesManager::SimplePreferenceValue& iValue,
            string                isKeyName,
            string                isError ) {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();
  string value( prefsMgr.GetValue( isKeyName ) );
  Assert( (iValue.ValueToString() == value), logic_error(isError.c_str()) );
}

// Write and read a few times with different values.
void TestWR( PreferencesManager::SimplePreferenceValue& iValue,
             PreferencesManager::SimplePreferenceValue& iValue2,
             string                isKeyName,
             string                isError ) {

  PreferencesManager::GetManager();

  TestW( iValue, isKeyName, isError );
  TestR( iValue, isKeyName, isError );

  TestW( iValue2, isKeyName, isError );
  TestR( iValue2, isKeyName, isError );
}


// Run lots of test cycles on a certain pref value type.
template <class TPrefValue, class TValue>
void TestType ( int icValues,
                string(*iGetName)(int n),
                TValue(*iGetNewValue)() ) {

  PreferencesManager& prefsMgr = PreferencesManager::GetManager();

  // Clear the prefs manager and reload our file.
  prefsMgr.Clear();
  prefsMgr.UseFile( "/tmp/testPrefs" );

  // Create our values using the iGetNewValue function we were passed.
  map<string,TPrefValue*> aValues;
  for ( int nValue = 0; nValue < icValues; nValue++ ) {
    TPrefValue* value = new TPrefValue( iGetNewValue() );
    aValues[iGetName(nValue)] = value;
  }

  // TestRegister and TestWR all our values.
  typename map<string,TPrefValue*>::iterator tValue;
  for ( tValue = aValues.begin(); tValue != aValues.end(); ++tValue ) {

    string sName = (*tValue).first;
    TPrefValue* value = (*tValue).second;

    TestRegister( *value, sName, "value didn't match" );
    TestWR( *value, *value, sName, "value didn't match" );
  }

  // Write the file, close it, and reload it.
  prefsMgr.SaveFile();
  prefsMgr.Clear();
  prefsMgr.UseFile( "/tmp/testPrefs" );

  // Read the values again and make sure they were what we wrote
  // before saving and reloading the file.
  for ( tValue = aValues.begin(); tValue != aValues.end(); ++tValue ) {

    string sName = (*tValue).first;
    TPrefValue* value = (*tValue).second;

    TestR( *value, sName, "didn't read correct value" );
  }
}

// Generate names and random values for our test values.
string GetIntName( int inValue ) {
  stringstream sName;
  sName << "intValue-" << inValue;
  return sName.str();
}

int GetIntValue() {
  return random();
}

string GetFloatName( int inValue ) {
  stringstream sName;
  sName << "floatValue-" << inValue;
  return sName.str();
}

float GetFloatValue() {
  float divisor = random();
  while ( divisor == 0 ) divisor = random();
  return (float)random() / divisor;
}

string GetStringName( int inValue ) {
  stringstream sName;
  sName << "stringValue-" << inValue;
  return sName.str();
}

string GetStringValue() {
  int cChars = random() % 500;
  stringstream sValue;
  for ( int nChar = 0; nChar < cChars; nChar++ ) {
    if ( nChar % 20 == 0 )
      sValue << " ";
    else
      sValue << (char)( (int)'a' + random() % 26 );
  }
  return sValue.str();
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  string sHeader = "This is the header.";

  try {

    PreferencesManager& prefsMgr = PreferencesManager::GetManager();
    prefsMgr.DisableOutput();

    // Try filenames. Punching in a relative filename should end up in
    // the home dir. A full filename should be the same.
    string fnPrefs = "TestPreferenceFile";
    prefsMgr.UseFile( fnPrefs );
    string fnPrefsTest = prefsMgr.GetFileName();
    char* sHomeDir = getenv("HOME");
    if ( sHomeDir ) {
      string fnHomeVersion = string(sHomeDir) + "/" + fnPrefs;
      stringstream sError;
      sError << "Filenames don't match: should be" << fnHomeVersion
      << ", got " << fnPrefs;
      Assert( (fnHomeVersion == fnPrefsTest), logic_error(sError.str()) );
    }

    // Full filename should be the same.
    fnPrefs = "/tmp/TestPreferenceFile";
    prefsMgr.UseFile( fnPrefs );
    fnPrefsTest = prefsMgr.GetFileName();
    stringstream sError;
    sError << "Filenames don't match: should be" << fnPrefs
    << ", got " << fnPrefsTest;
    Assert( (fnPrefs == fnPrefsTest), logic_error(sError.str()) );

    // Set the header.
    prefsMgr.SetHeader( sHeader );

    // Run a bunch of tests on our prefs types.
    int const kzValuesToTest = 100;
    TestType<PreferencesManager::IntPrefValue,int>
    ( kzValuesToTest, GetIntName, GetIntValue );
    TestType<PreferencesManager::FloatPrefValue,float>
    ( kzValuesToTest, GetFloatName, GetFloatValue );
    TestType<PreferencesManager::StringPrefValue,string>
    ( kzValuesToTest, GetStringName, GetStringValue );
  } catch ( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }
  catch ( logic_error& e ) {
    cerr << "failed: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed." << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

