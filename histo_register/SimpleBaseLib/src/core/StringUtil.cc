// Licensed under MIT license; see license.txt.

#include <sbl/core/StringUtil.h>
#include <sbl/core/File.h>
#include <sbl/core/UnitTest.h>
#include <stdarg.h>
namespace sbl {


//-------------------------------------------
// STRING UTILITIES
//-------------------------------------------


/// allocate a string and fill in using sprintf()
String sprintF( const char *format, ... ) {
    char buf[ 100000 ]; // fix!!!

    // do printf
    va_list argList;
    va_start( argList, format );
    vsprintf( buf, format, argList ); 

    // allocate new string
    return String( buf );
}


/// python string hash
// (from string_hash in stringobject.c)
int strHash( const String &s ) {
    const char *cstr = s.c_str();
    int count = s.length();
    int x = *cstr << 7;
    while (--count >= 0)
        x = (1000003*x) ^ *cstr++;
    x ^= s.length();
    if (x == -1)
        x = -2;
    return x;
}


/// python string hash
// (from string_hash in stringobject.c)
int strHash( const char *cstr ) {
    int len = (int) strlen( cstr ); // assume 32-bit length
    int count = len;
    int x = *cstr << 7;
    while (--count >= 0)
        x = (1000003*x) ^ *cstr++;
    x ^= len;
    if (x == -1)
        x = -2;
    return x;
}


/// python string hash
// (from string_hash in stringobject.c)
int strHash( const char *bytes, int length ) {
    int count = length;
    int x = *bytes << 7;
    while (--count >= 0)
        x = (1000003*x) ^ *bytes++;
    x ^= length;
    if (x == -1)
        x = -2;
    return x;
}


//-------------------------------------------
// STRING ARRAY UTILITIES
//-------------------------------------------


/// concatenate the given strings using the given joiner between each one
String join( const Array<String> &strArr, const String &joiner ) {

    // compute length of new string
    int joinerLength = joiner.length();
    int newLength = 0;
    for (int i = 0; i < strArr.count(); i++) {
        newLength += strArr[ i ].length();
        if (i + 1 < strArr.count())
            newLength += joinerLength;
    }

    // allocate string to hold result
    unsigned short *newStr = new unsigned short[ newLength + 10 ];
    assertDebug( newStr );

    // build joined string
    int pos = 0;
    for (int i = 0; i < strArr.count(); i++) {

        // add array item
        const String &item = strArr[ i ];
        for (int j = 0; j < item.length(); j++)
            newStr[ pos++ ] = item.str()[ j ];

        // add joiner
        if (i + 1 < strArr.count())
            for (int j = 0; j < joiner.length(); j++) 
                newStr[ pos++ ] = joiner.str()[ j ];
    }
    newStr[ pos ] = 0;
    String s( newStr );
    delete [] newStr;
    return s;
}


/// load a file into an array of strings (one per line)
Array<String> loadStrings( const String &fileName, bool includeComments ) {
    File file( fileName, FILE_READ, FILE_TEXT );
    Array<String> strArray;
    while (file.endOfFile() == false) {
        String line = file.readLine();
        if (line.length() && (line.startsWith( "#" ) == false || includeComments)) {
            strArray.append( new String( line ));
        }
    }
    return strArray;
}


// comparison function used by string array sorting function
// fix(later): should use unicode comparison
int stringCompare( const void *v1, const void *v2 ) {
    const String *s1 = *((const String **) v1);
    const String *s2 = *((const String **) v2);
    return strcmp( s1->c_str(), s2->c_str() );
}


/// sort the array of strings (case sensitive, ascii ordered, using strcmp)
Array<String> sort( const Array<String> &strArr ) {
    Array<String> sorted;
    int count = strArr.count();
    String const ** data = new String const*[ count ]; // non-const array of pointers to const data
    for (int i = 0; i < count; i++) 
        data[ i ] = &(strArr[ i ]);
    qsort( data, count, sizeof(String*), stringCompare );
    for (int i = 0; i < count; i++) 
        sorted.append( new String( *(data[ i ]) ) );    
    delete [] data;
    return sorted;
}


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns formatted memory quantity 
/// (note: returns pointer to internal memory)
const char *memString( int byteCount ) {
    static char memString[1000];
    if (byteCount < 1024)
        sprintf( memString, "%d", byteCount );
    else if (byteCount < 1024 * 1024)
        sprintf( memString, "%dkb", (byteCount - 1) / 1024 + 1 );
    else 
        sprintf( memString, "%dmb", (byteCount - 1) / (1024 * 1024) + 1 );
    return memString;
}


/// returns bool formatted a string "yes" or "no"
String yesNo( bool val ) {
    if (val)
        return "yes";
    else
        return "no";
}


//-------------------------------------------
// TESTING
//-------------------------------------------


// test string array sorting
bool testSortStringArrayCase( const String &list, const String &listGoal ) {
    Array<String> arr = list.split( ";" );
    Array<String> arrGoal = listGoal.split( ";" );
    unitAssert( arr.count() == arrGoal.count() );
    Array<String> sorted = sort( arr );
    unitAssert( sorted.count() == arrGoal.count() )
    for (int i = 0; i < arr.count(); i++) 
        unitAssert( arrGoal[ i ] == sorted[ i ] );
    return true;
}


// test string array sorting
bool testSortStringArray() {
    unitAssert( testSortStringArrayCase( "z;a;d", "a;d;z" ) );
    unitAssert( testSortStringArrayCase( "z;d;a", "a;d;z" ) );
    unitAssert( testSortStringArrayCase( "a;d;z", "a;d;z" ) );
    unitAssert( testSortStringArrayCase( "5;1;4;3;2", "1;2;3;4;5" ) );
    unitAssert( testSortStringArrayCase( "05;01;04;03;02", "01;02;03;04;05" ) );
    unitAssert( testSortStringArrayCase( "aab;aaa", "aaa;aab" ) );
    unitAssert( testSortStringArrayCase( "x", "x" ) );
    return true;
}


// test string split
bool testSplitCase( const char *str, const char *splitChars, int count, const char *firstItem ) {
    disp( 3, "%s:", str );
    Array<String> strArray = String( str ).split( splitChars );
    for (int i = 0; i < strArray.count(); i++) 
        disp( 4, "[%s]", strArray[ i ].c_str());
    unitAssert( strArray.count() == count );
    unitAssert( strArray[ 0 ] == firstItem );
    return true;
}


// test string split
bool testSplit() {
    unitAssert( testSplitCase( "abc-def-ghi", "-", 3, "abc" ) );
    unitAssert( testSplitCase( "a-b-c", "-", 3, "a" ) );
    unitAssert( testSplitCase( "a-b", "-", 2, "a" ) );
    unitAssert( testSplitCase( "abc", "-", 1, "abc" ) );
    return true;
}


// test string join
bool testJoin() {
    Array<String> strArr;
    String join0 = join( strArr, "-" );
    unitAssert( join0 == "" );
    strArr.append( new String( "abc" ));
    String join1 = join( strArr, "-" );
    unitAssert( join1 == "abc" );
    strArr.append( new String( "def" ));
    strArr.append( new String( "123" ));
    String join2 = join( strArr, "-" );
    unitAssert( join2 == "abc-def-123" );
    String join3 = join( strArr, "" );
    unitAssert( join3 == "abcdef123" );
    return true;
}


// register commands, etc. defined in this module
void initStringUtil() {
    registerUnitTest( testSplit );
    registerUnitTest( testJoin );
    registerUnitTest( testSortStringArray );
}


} // end namespace sbl

