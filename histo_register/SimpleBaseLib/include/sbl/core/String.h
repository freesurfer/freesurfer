#ifndef _SBL_STRING_H_
#define _SBL_STRING_H_
#include <sbl/core/Array.h>
namespace sbl {


// conventions:
// a generic String variable is called "s"
// a generic 8-bit C string variable is called "cstr"
// a generic 16-bit unicode string variable is called "str"


/// The String class represents unicode strings and has python-inspired methods.
class String {
public:

    // basic constructors
    String( const String &s );
    String( const char *cstr );
    explicit String( const unsigned short *str );
    String();
    ~String();

    /// construct repeating string
    String( unsigned short val, int repeatCount );

    /// sub-string constructor
    String( const String &s, int startPosition, int length );

    /// basic assignment operator
    String &operator=( const String &s );

    /// bytes used by this object
    int memUsed() const;

    //-------------------------------------------
    // ACCESS STRING CONTENTS
    //-------------------------------------------

    /// get length of string
    inline int length() const { return m_length; }

    /// 8-bit version of string
    inline const char *c_str() const { return m_cstring; }

    /// unicode version of string
    inline const unsigned short *str() const { return m_string; }

    /// set an individual character
    void set( int index, unsigned short c );

    /// get an individual character
    unsigned short get( int index ) const;

    //-------------------------------------------
    // STRING SEARCH / COMPARISON
    //-------------------------------------------

    /// return position of first/last instance of character; -1 if not found
    int firstCharPos( unsigned short c ) const;
    int lastCharPos( unsigned short c ) const;
    
    /// true if starts with given string
    bool startsWith( const String &s ) const;

    /// true if ends with given string
    bool endsWith( const String &s ) const;

    /// returns true if the string contains the given string as a sub-string
    bool contains( const String &s ) const;
    bool contains( unsigned short c ) const { return firstCharPos( c ) >= 0; }

    /// compare with another string
    bool operator==( const String &s ) const;
    inline bool operator!=( const String &s ) const { return !operator==( s ); }

    /// returns true if all characters lower case a-z (assuming basic ascii)
    bool isLower() const;

    /// returns true if all characters upper case A-Z(assuming basic ascii)
    bool isUpper() const;

    //-------------------------------------------
    // SUB-STRINGS
    //-------------------------------------------

    /// create sub-string (wrapper for sub-string constructor)
    inline String subString( int startPosition, int length ) const { return String( *this, startPosition, length ); }

    /// get string left/right of given position
    inline String leftOf( int position ) const { return subString( 0, position ); }
    inline String rightOf( int position ) const { return subString( position + 1, m_length - position - 1 ); }

    /// get string left/right of first/last instance of the specified character
    inline String rightOfFirst( unsigned short c ) const { int pos = firstCharPos( c ); if (pos >= 0) return rightOf( pos ); return *this; }
    inline String rightOfLast( unsigned short c ) const { int pos = lastCharPos( c ); if (pos >= 0) return rightOf( pos ); return *this; }
    inline String leftOfFirst( unsigned short c ) const { int pos = firstCharPos( c ); if (pos >= 0) return leftOf( pos ); return *this; }
    inline String leftOfLast( unsigned short c ) const { int pos = lastCharPos( c ); if (pos >= 0) return leftOf( pos ); return *this; }

    /// split a string at any of the given split characters; see python split() examples 
    Array<String> split( const char *splitChars, bool quoteAware = false ) const;

    //-------------------------------------------
    // STRING MANIPULATION / OPERATIONS
    //-------------------------------------------

    /// append a string to self
    void append( const String &s );

    /// return a string with new string appended
    inline void operator+=( const String &s ) { append( s ); } // fix
    inline const String operator+( const String &s ) const { String sNew( *this ); sNew.append( s ); return sNew; }

    /// replace each instance of the specified character
    void replaceInPlace( unsigned short find, unsigned short replace );
    String replace( unsigned short find, unsigned short replace ) const;

    /// returns string with whitespace removed from beginning and end
    String strip() const;

    /// returns string with stripChar instances removed from beginning and end
    String strip( unsigned short stripChar ) const;

    //-------------------------------------------
    // STRING CONVERSION
    //-------------------------------------------

    /// lower case version of this string
    String lower() const;

    /// convert string to numeric value
    int toInt() const;
    float toFloat() const;
    double toDouble() const;
    inline bool toBool() const { return toInt() ? true : false; }

private:

    // 8-bit (latin) string value
    char *m_cstring; // hack: put m_cstring first so that printf of string object (but not pointer to string object) works (at least when compiled by Visual C++)

    // unicode string value
    unsigned short *m_string;

    // length of string (not including null terminator)
    int m_length;
};


} // end namespace sbl
#endif // _SBL_STRING_H_

