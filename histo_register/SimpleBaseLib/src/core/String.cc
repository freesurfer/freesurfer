// Licensed under MIT license; see license.txt.

#include <sbl/core/String.h>
#include <sbl/core/Display.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> // for atoi/atof
namespace sbl {


/// unicode strlen
int strLen( const unsigned short *str ) {
    assertDebug( str );
    int i = 0;
    while (str[ i ])
        i++;
    return i;
}


/// sub-string constructor
String::String( const String &s, int startPosition, int length ) {

    // check args
    assertDebug( length >= 0 );
    assertDebug( startPosition >= 0 );
    assertDebug( startPosition + length <= s.length() );

    // allocate new string
    m_length = length;
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );

    // copy input string 
    for (int i = 0; i < m_length; i++) {
        m_string[ i ] = s.m_string[ startPosition + i ];
        m_cstring[ i ] = s.m_cstring[ startPosition + i ];
    }
    
    // add terminating zero
    m_string[ m_length ] = 0;
    m_cstring[ m_length ] = 0;
}


/// copy constructor
String::String( const String &s ) {
    m_length = s.m_length;
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );

    // copy input string (including terminating 0)
    for (int i = 0; i < m_length + 1; i++) {
        m_string[ i ] = s.m_string[ i ];
        m_cstring[ i ] = s.m_cstring[ i ];
    }
}


/// construct from c-style string
String::String( const char *cstr ) {
    assertDebug( cstr );
    m_length = (int) strlen( cstr ); // assume 32-bit length
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );

    // copy input string (including terminating 0)
    for (int i = 0; i < m_length + 1; i++) {
        m_string[ i ] = cstr[ i ];
        m_cstring[ i ] = cstr[ i ];
    }
}


/// construct from unicode string
String::String( const unsigned short *str ) {
    assertDebug( str );
    m_length = strLen( str );
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );

    // copy input string (including terminating 0)
    for (int i = 0; i < m_length + 1; i++) {
        m_string[ i ] = str[ i ];
        m_cstring[ i ] = (char) str[ i ];
    }
}


/// construct empty string
String::String() {
    m_length = 0;
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );
    m_string[ 0 ] = 0;
    m_cstring[ 0 ] = 0;
}


/// construct repeating string
String::String( unsigned short val, int repeatCount ) {
    m_length = repeatCount;
    m_string = new unsigned short[ m_length + 1 ];
    m_cstring = new char[ m_length + 1 ];
    assertDebug( m_string );
    assertDebug( m_cstring );
    for (int i = 0; i < m_length; i++) {
        m_string[ i ] = val;
        m_cstring[ i ] = (char) val;
    }
    m_string[ m_length ] = 0;
    m_cstring[ m_length ] = 0;
}


// deallocate string
String::~String() {
    delete [] m_string;
    delete [] m_cstring;
}


/// bytes used by this object
int String::memUsed() const {
    return sizeof( String ) + (m_length + 1) * (sizeof( short ) + sizeof( char ));
}


/// compare with another string
bool String::operator==( const String &s ) const {
    if (m_length != s.m_length)
        return false;
    for (int i = 0; i < m_length; i++)
        if (m_string[ i ] != s.m_string[ i ])
            return false;
    return true;
}


/// lower case version of this string
String String::lower() const {
    String lower( *this );
    for (int i = 0; i < m_length; i++) {
        int c = m_string[ i ];
        if (c >= 'A' && c <= 'Z') {
            lower.m_string[ i ] = c - 'A' + 'a';
            lower.m_cstring[ i ] = c - 'A' + 'a';
        }
    }
    return lower;
}


/// returns true if all characters lower case a-z (assuming basic ascii)
bool String::isLower() const {
    if (m_length == 0)
        return false;
    for (int i = 0; i < m_length; i++) {
        int c = m_string[ i ];
        if (c < 'a' || c > 'z') {
            return false;
        }    
    }
    return true;
}


/// returns true if all characters upper case A-Z (assuming basic ascii)
bool String::isUpper() const {
    if (m_length == 0)
        return false;
    for (int i = 0; i < m_length; i++) {
        int c = m_string[ i ];
        if (c < 'A' || c > 'Z') {
            return false;
        }    
    }
    return true;
}


/// returns string with whitespace removed from beginning and end
String String::strip() const {
    // fix(clean): clean up
    if (m_length == 1 && m_cstring[ 0 ] <= 32)
        return "";
    int startIndex = 0;
    while (startIndex < m_length - 1 && m_cstring[ startIndex ] <= 32)
        startIndex++;
    int endIndex = m_length - 1;
    while (endIndex > 0 && m_cstring[ endIndex ] <= 32)
        endIndex--;
    if (endIndex < startIndex)
        return "";
    return subString( startIndex, endIndex - startIndex + 1 );
}


/// returns string with stripChar instances removed from beginning and end
String String::strip( unsigned short stripChar ) const {
    // fix(clean): clean up
    if (m_length == 1 && m_string[ 0 ] == stripChar)
        return "";
    int startIndex = 0;
    while (m_string[ startIndex ] == stripChar && startIndex < m_length - 1)
        startIndex++;
    int endIndex = m_length - 1;
    while (m_string[ endIndex ] == stripChar && endIndex > 0)
        endIndex--;
    if (endIndex < startIndex)
        return "";
    return subString( startIndex, endIndex - startIndex + 1 );
}


/// split a string at any of the given split characters; see python split() examples 
// fix(later): error checking
Array<String> String::split( const char *splitChars, bool quoteAware ) const {
    Array<String> strArray;
    int len = length();
    if (len == 0)
        return strArray;
    int splitCharCount = (int) strlen( splitChars ); // assume 32-bit length
    int lastSplit = -1;
    bool insideQuotes = false;
    for (int i = 0; i < len; i++) {

        // check for one of the split characters
        bool foundSplit = false;
        for (int j = 0; j < splitCharCount; j++) {
            if (m_string[ i ] == splitChars[ j ]) {
                foundSplit = true;
                break;
            }
        }

        // if quote aware, ignore splits inside quotes
        if (quoteAware) {
            if (m_string[ i ] == '"')
                insideQuotes = !insideQuotes;
        }

        // if split here
        if (foundSplit && insideQuotes == false) {
            strArray.appendCopy( subString( lastSplit + 1, i - lastSplit - 1 ));
            lastSplit = i;
        }
    }
    strArray.appendCopy( subString( lastSplit + 1, len - lastSplit - 1 ));
    return strArray;
}


/// convert string to numeric value
int String::toInt() const { 
    return atoi( m_cstring ); 
}


/// convert string to numeric value
float String::toFloat() const { 
    return (float) atof( m_cstring ); 
}


/// convert string to numeric value
double String::toDouble() const { 
    return atof( m_cstring ); 
}


/// append a string to self
void String::append( const String &s ) {
    int lengthNew = m_length + s.length();
    unsigned short *stringNew = new unsigned short[ lengthNew + 1 ];
    char *cstringNew = new char[ lengthNew + 1 ];
    assertDebug( stringNew );
    assertDebug( cstringNew );

    // copy original 
    for (int i = 0; i < m_length; i++) {
        stringNew[ i ] = m_string[ i ];
        cstringNew[ i ] = m_cstring[ i ];
    }

    // copy appended string
    for (int i = 0; i < s.length() + 1; i++) {
        stringNew[ i + m_length ] = s.m_string[ i ];
        cstringNew[ i + m_length ] = s.m_cstring[ i ];
    }

    // replace original
    delete [] m_string;
    delete [] m_cstring;
    m_length = lengthNew;
    m_string = stringNew;
    m_cstring = cstringNew;
}


/// return position of first instance of character; -1 if not found
int String::firstCharPos( unsigned short c ) const {
    for (int i = 0; i < m_length; i++) 
        if (m_string[ i ] == c)
            return i;
    return -1;
}


/// return position of last instance of character; -1 if not found
int String::lastCharPos( unsigned short c ) const {
    for (int i = m_length - 1; i >= 0; i--) 
        if (m_string[ i ] == c)
            return i;
    return -1;
}


/// replace each instance of the specified character
void String::replaceInPlace( unsigned short find, unsigned short replace ) {
    for (int i = 0; i < m_length; i++) {
        if (m_string[ i ] == find) {
            m_string[ i ] = replace;
            m_cstring[ i ] = (char) replace;
        }
    }
}


/// replace each instance of the specified character
String String::replace( unsigned short find, unsigned short replace ) const {
    String s( *this );
    s.replaceInPlace( find, replace );
    return s;
}


/// set an individual character
void String::set( int index, unsigned short c ) { 
    assertDebug( index >= 0 && index < m_length ); 
    m_string[ index ] = c; 
    m_cstring[ index ] = (char) c; 
}


/// get an individual character
unsigned short String::get( int index ) const { 
    assertDebug( index >= 0 && index < m_length ); 
    return m_string[ index ];
}


/// true if starts with given string
bool String::startsWith( const String &s ) const {
    if (s.length() == 0)
        return true;
    if (s.length() > length())
        return false;
    String subStr = subString( 0, s.length() );
    return s == subStr;
}


/// true if ends with given string
bool String::endsWith( const String &s ) const {
    if (s.length() == 0)
        return true;
    if (s.length() > length())
        return false;
    String subStr = subString( length() - s.length(), s.length() );
    return s == subStr;
}


/// returns true if the string contains the given string as a sub-string
// fix(later): should handle unicode substrings
bool String::contains( const String &s ) const {
    return strstr( c_str(), s.c_str() ) ? true : false;
}


/// basic assignment operator
String &String::operator=( const String &s ) {
    if (&s != this) {
        delete [] m_string;
        delete [] m_cstring;
        m_length = s.m_length;
        m_string = new unsigned short[ m_length + 1 ];
        m_cstring = new char[ m_length + 1 ];
        assertDebug( m_string );
        assertDebug( m_cstring );

        // copy input string (including terminating 0)
        for (int i = 0; i < m_length + 1; i++) {
            m_string[ i ] = s.m_string[ i ];
            m_cstring[ i ] = s.m_cstring[ i ];
        }
    }
    return *this;
}


} // end namespace sbl

