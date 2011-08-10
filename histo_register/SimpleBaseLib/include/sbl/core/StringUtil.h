#ifndef _SBL_STRING_UTIL_H_
#define _SBL_STRING_UTIL_H_
#include <sbl/core/String.h>
#include <sbl/core/Array.h>
namespace sbl {


/*! \file StringUtil.h
    \brief The StringUtil module includes functions for working with strings
    (represented by String objects).
*/


// register commands, etc. defined in this module
void initStringUtil();


//-------------------------------------------
// STRING UTILITIES
//-------------------------------------------


/// create a formatting string in the manner of sprintf()
String sprintF( const char *format, ... );


/// python-compatible string hash
int strHash( const String &s );
int strHash( const char *cstr );
int strHash( const char *bytes, int length );


//-------------------------------------------
// STRING ARRAY UTILITIES
//-------------------------------------------


/// concatenate the given strings using the given joiner between each one
String join( const Array<String> &strArr, const String &joiner );


/// load a file into an array of strings (one per line)
Array<String> loadStrings( const String &fileName, bool includeComments );


/// sort the array of strings (case sensitive, ascii ordered, using strcmp)
Array<String> sort( const Array<String> &strArr );


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns formatted memory quantity;
/// note: returns pointer to internal memory
const char *memString( int byteCount );


/// returns bool formatted a string "yes" or "no"
String toStrYesNo( bool val );


} // end namespace sbl
#endif // _SBL_STRING_UTIL_H_

