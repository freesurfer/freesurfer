#ifndef _SBL_TAGGED_FILE_H_
#define _SBL_TAGGED_FILE_H_
#include <sbl/core/File.h>
namespace sbl {


// all array tags should be above this value
#define TAG_OFFSET_ARRAY 10000


// a marker for the end of an item
#define TAG_END_SECTION 20000


/*
    A tagged file consists of a sequence of items.  Each item has one of two forms:
    1. a single value (int, float, string, int array, float array, etc.)
        - a tag
        - the size of the data (in bytes)
        - followed by the data
    2. an array of sub-items (variable size objects, as opposed to an array of floats or an array of characters)
        - a tag
        - the number of items
        - the sequence of sub-items 
            - each sub-item consists of one or more items
            - each sub-item ends with a TAG_END_SECTION 
*/


/// The TaggedFile class provides facilities for reading/write tag-based file formats.
class TaggedFile : public File {
public:

    /// basic constructor; see File constructor for more info
    TaggedFile( const String &fileName, FileOpenMode mode, FileOpenType type )
                   : File( fileName, mode, type ) {
        assertAlways( type == FILE_BINARY );
    }
    virtual ~TaggedFile() {}

    /// write basic data types with size info
    inline void writeBool( bool val ) { writeTagInfo( sizeof( int ) ); writeInt( (int) val ); }
    inline void writeUChar( unsigned char val ) { writeTagInfo( sizeof( unsigned char ) ); fwrite( &val, sizeof(unsigned char), 1, m_file ); }
    inline void writeUShort( unsigned short val ) { writeTagInfo( sizeof( unsigned short ) ); fwrite( &val, sizeof(unsigned short), 1, m_file ); }
    inline void writeInt( int val ) { writeTagInfo( sizeof( int ) ); fwrite( &val, sizeof(int), 1, m_file ); }
    inline void writeFloat( float val ) { writeTagInfo( sizeof( float ) ); fwrite( &val, sizeof(float), 1, m_file ); }
    inline void writeDouble( double val ) { writeTagInfo( sizeof( double ) ); fwrite( &val, sizeof(double), 1, m_file ); }
    void writeString( const String &s );

    /// read basic data types with size info
    inline bool readBool() { assertAlways( readTagInfo() == sizeof( int ) ); return readInt() ? true : false; }
    inline unsigned char readUChar() { assertAlways( readTagInfo() == sizeof( unsigned char ) ); unsigned char val = 0; CHECK_READ( fread( &val, sizeof(unsigned char), 1, m_file ) ); return val; }
    inline unsigned short readUShort() { assertAlways( readTagInfo() == sizeof( unsigned short ) ); unsigned short val = 0; CHECK_READ( fread( &val, sizeof(unsigned short), 1, m_file ) ); return val; }
    inline int readInt() { assertAlways( readTagInfo() == sizeof( int ) ); int val = 0; CHECK_READ( fread( &val, sizeof(int), 1, m_file ) ); return val; }
    inline float readFloat() { assertAlways( readTagInfo() == sizeof( float ) ); float val = 0; CHECK_READ( fread( &val, sizeof(float), 1, m_file ) ); return val; }
    inline double readDouble() { assertAlways( readTagInfo() == sizeof( double ) ); double val = 0; CHECK_READ( fread( &val, sizeof(double), 1, m_file ) ); return val; }
    String readString();

    /// write a vector
    // fix(later): use writeRawData
    template<typename T> void writeVector( const Vector<T> &v ) {
        writeTagInfo( v.length() * sizeof( T ) );
        writeBlock( v.dataPtr(), v.length() * sizeof( T ) );
    }

    /// read a vector
    // fix(later): use readRawData
    template<typename T> Vector<T> readVector() {
        Vector<T> v;
        int size = readTagInfo();
        int len = size / sizeof( T );
        if (len < 0) {
            warning( "invalid vector in file" );
        } else { 
            v.setLength( len );
            readBlock( v.dataPtr(), v.length() * sizeof( T ) );
        }
        return v;
    }

    /// write array; assumes type T has a .save() method
    template<typename T> void writeArray( const Array<T> &a ) {
        writeTagInfo( a.count() );
        for (int i = 0; i < a.count(); i++) {
            a[ i ].save( *this );
            writeTagInfo( TAG_END_SECTION );
        }
    }

    /// read array; assumes type T has a file load constructor
    template<typename T> void readArray( Array<T> &a ) {
        int count = readTagInfo();
        if (count < 0) {
            warning( "invalid object count in file" );
        } else {
            for (int i = 0; i < count; i++) 
                a.append( new T( *this ) );
        }
    }

    /// write string array (same file structure as other arrays)
    // note: could have used above readArray/writeArray but don't want to add load/save methods to String class
    template<typename T> void writeStrings( const Array<T> &a ) {
        writeTagInfo( a.count() );
        for (int i = 0; i < a.count(); i++) {
            writeString( a[ i ] ); // note: this is tagged
            writeTagInfo( TAG_END_SECTION );
        }
    }

    /// read string array (same file structure as other arrays)
    template<typename T> void readStrings( Array<T> &a ) {
        int count = readTagInfo();
        if (count < 0) {
            warning( "invalid object count in file" );
        } else {
            for (int i = 0; i < count; i++) {
                a.appendCopy( readString() ); // note: this is tagged
                int tag = readTagInfo();
                if (tag != TAG_END_SECTION) {
                    warning( "invalid string array in file" );
                    break;
                }
            }
        }
    }

    /// write raw data of a particular type (length is number of elements not number of bytes)
    template<typename T> void writeRawData( const T *data, int length ) {
        writeTagInfo( length * sizeof( T ) );
        writeBlock( data, length * sizeof( T ) );
    }

    /// write raw data of a particular type (length is number of elements not number of bytes)
    template<typename T> T *readRawData( int &length ) {
        T *data = NULL;
        int size = readTagInfo();
        int len = size / sizeof( T );
        if (len < 0) {
            warning( "invalid raw data in file" );
            length = 0;
        } else { 
            data = new T[ len ];
            readBlock( data, len * sizeof( T ) );
            length = len;
        }
        return data;
    }

    /// write a tag and return file reference, allowing this syntax: file.tag( TAG ).writeString( "foo" )
    inline TaggedFile &tag( int tagVal ) { writeTagInfo( tagVal ); return *this; }

    /// read and discard a value from a tagged file
    void skipTag( int tag );

    /// read a tag
    int readTag() { if (endOfFile()) return TAG_END_SECTION; else return readTagInfo(); }

    /// read/write data size for tagged file
    inline int readTagInfo() { int size = 0; CHECK_READ( fread( &size, sizeof(int), 1, m_file ) ); return size; }
    inline void writeTagInfo( int size ) { fwrite( &size, sizeof(int), 1, m_file ); }

private:

    // disable copy constructor and assignment operator
    TaggedFile( const TaggedFile &x );
    TaggedFile &operator=( const TaggedFile &x );
};


} // end namespace sbl
#endif // _SBL_TAGGED_FILE_H_

