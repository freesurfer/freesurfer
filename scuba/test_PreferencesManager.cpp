#include <stdlib.h>
#include <string>
#include <iostream>
#include <tcl.h>
#include "PreferencesManager.h"

#define Assert(x,s)   if(!(x)) { throw (char const*) s; }

using namespace std;

char* Progname = "test_PreferencesManager";

template <class T>
void TestType( T iValue, string isKeyName, string isError ) {

  T origValue = iValue;
  T value = iValue;
  prefsMgr.RegisterValue<T>( isKeyName, "A sample value.",
			     origValue );
  
  value = prefsMgr.GetValue<T>( isKeyName );
  
  Assert( (origValue == value), isError );
}


int main ( int argc, char** argv ) {

  string fnPrefs = "TestPreferenceFile";
  string sHeader = "This is the header.";

  try { 
 
    PreferencesManager& prefsMgr = PreferencesManager::GetManager();

    prefsMgr.UseFile( fnPrefs );

    prefsMgr.SetHeader( sHeader );

    TestType<int>( 5, "intValue", "int value didn't match" );
    TestType<float>( 5.5, "floatValue", "float value didn't match" );
    TestType<string>( "a string", "stringValue", "string value didn't match" );

  }
  catch( char const* iMsg ) {
    cerr << "failed: " << iMsg << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}


