// Licensed under MIT license; see license.txt.

#include <sbl/core/File.h>
#include <sbl/core/Display.h>
#include <sbl/core/UnitTest.h>
#include <string.h>
#include <stdarg.h>
#ifdef USE_ZLIB
    #include <zlib.h>
#endif
namespace sbl {


// fix(clean): pay better attention to fscanf/fread return values
void readError() {
}


//-------------------------------------------
// FILE CLASS
//-------------------------------------------


// open a file
File::File( const String &fileName, FileOpenMode mode, FileOpenType type ) {
    m_binary = false;
    m_file = NULL;
    m_gzFile = NULL;
    if (fileName.length())
        open( fileName, mode, type );
}


// close the file
File::~File() {
    if (m_file)
        close();
}


/// open a file
void File::open( const String &fileName, FileOpenMode mode, FileOpenType type ) {
    assertAlways( m_file == NULL );
    m_binary = (type == FILE_BINARY || type == FILE_GZIP_BINARY);
    const char *modeStr = NULL;
    if (mode == FILE_APPEND)
        modeStr = type == FILE_TEXT ? "a" : "ab";
    else if (mode == FILE_WRITE)
        modeStr = type == FILE_TEXT ? "w" : "wb";
    else
        modeStr = type == FILE_TEXT ? "r" : "rb";
    m_file = fopen( fileName.c_str(), modeStr );
#ifdef USE_ZLIB
    if (type == FILE_GZIP_BINARY && m_file) {
        m_gzFile = gzopen( fileName.c_str(), modeStr ); // fix(clean): don't also open as m_file
    }
#endif
}


/// close the file
void File::close() { 
#ifdef USE_ZLIB
    if (m_gzFile) {
        gzclose( m_gzFile );
        m_file = NULL;
    }
#endif
    if (m_file) {
        fclose( m_file ); 
        m_file = NULL;
    }
}


/// flush writes to disk
void File::flush() { 
    if (m_gzFile) {
#ifdef USE_ZLIB
        gzflush( m_gzFile, Z_FINISH ); // note: this flush should be used sparingly
#endif
    } else {
        fflush( m_file ); 
    }
} 


/// returns true if end of file reached on last read; returns true if file not successfully opened
bool File::endOfFile() { 
    bool eof = true;
    if (m_file != NULL) {
        if (m_gzFile) {
#ifdef USE_ZLIB
            eof = gzeof( m_gzFile ) ? true : false;
#endif
        } else {
            eof = feof( m_file ) ? true : false;
        }
    }
    return eof; 
} 


/// seek to the given position
void File::seek( int offset, bool relative ) {
    if (m_gzFile) {
#ifdef USE_ZLIB
        gzseek( m_gzFile, offset, relative ? SEEK_CUR : SEEK_SET ); 
#endif
    } else {
        fseek( m_file, offset, relative ? SEEK_CUR : SEEK_SET ); 
    }
}


/// returns the current position (assumes file position fits in int)
int File::tell() {
    int pos = 0;
    if (m_gzFile) {
#ifdef USE_ZLIB
        pos = gztell( m_gzFile ); 
#endif
    } else {
        pos = ftell( m_file ); 
    }
    return pos;
}


/// write a formatted string to a file
void File::writeF( const char *str, ... ) {
    assertDebug( m_file );
    char buf[ 10000 ];
    va_list argList;
    va_start( argList, str );
    vsprintf( buf, str, argList ); 
    int len = strlen( buf );
    assertAlways( m_binary == false );
    fwrite( buf, len, 1, m_file ); 
}


/// write length and raw bytes, regardless of whether text or binary format
void File::writeString( const String &s ) {
    assertDebug( m_file );
    writeInt( s.length() );
    writeBlock( s.c_str(), s.length() );
}


/// write raw bytes without writing length
void File::writeRawString( const String &s ) {
    assertDebug( m_file );
    writeBlock( s.c_str(), s.length() );
}


/// read up to end of line; strip trailing end-of-line characters
String File::readLine() {
    assertDebug( m_file );

    // read line into buffer
    char buf[ 10000 ];
    buf[ 0 ] = 0;
    CHECK_READ( fgets( buf, 10000, m_file ) );

    // strip trailing CRLF
    int len = (int) strlen( buf ); // assumes 32-bit string length
    while (len && (buf[ len - 1 ] == 10 || buf[ len - 1 ] == 13)) {
        buf[ len - 1 ] = 0;
        len--;
    }

    // return string object copy of buffer
    return String( buf );
}


/// reads length and raw bytes, regardless of whether text or binary format;
/// assumes string contains no 0 values
String File::readString() {
    assertDebug( m_file );
    int len = readInt();
    if (len < 0 || len > 100000000) {
        warning( "invalid string" );
        return String();
    }

    // read string data
    char *data = new char[ len + 1 ];
    if (m_binary == false)
        readBlock( data, 1 );
    readBlock( data, len );
    data[ len ] = 0;
    String s( data );
    delete [] data;
    return s;
}


/// write block of data (assumes binary)
void File::writeBlock( const void *data, int byteCount ) { 
    if (m_gzFile) {
#ifdef USE_ZLIB
        gzwrite( m_gzFile, data, byteCount );
#endif
    } else {
        fwrite( data, byteCount, 1, m_file ); 
    }
}


/// read block of data (assumes binary)
void File::readBlock( void *data, int byteCount ) { 
#ifdef USE_ZLIB
    if (m_gzFile)
        gzread( m_gzFile, data, byteCount );
    else
#endif
        CHECK_READ( fread( data, byteCount, 1, m_file ) );
}


//-------------------------------------------
// TESTING
//-------------------------------------------


// test the file class
bool testFile( FileOpenType type ) {
    String fileName = "test.dat";

    // write data
    File outFile( fileName, FILE_WRITE, type );
    unitAssert( outFile.openSuccess() );
    outFile.writeBool( false );
    outFile.writeUChar( 'A' );
    outFile.writeUShort( 7777 );
    outFile.writeInt( 123456789 );
    outFile.writeFloat( 3.14f );
    outFile.writeDouble( 3.1415 );
    outFile.close();

    // read data
    File inFile( fileName, FILE_READ, type );
    unitAssert( inFile.openSuccess() );
    unitAssert( inFile.readBool() == false );
    unitAssert( inFile.readUChar() == 'A' );
    unitAssert( inFile.readUShort() == 7777 );
    unitAssert( inFile.readInt() == 123456789 );
    unitAssert( inFile.readFloat() == 3.14f );
    unitAssert( inFile.readDouble() == 3.1415 );
    return true;
}


// test the file class
bool testFile() {
    unitAssert( testFile( FILE_BINARY ) );
    unitAssert( testFile( FILE_TEXT ) );
    unitAssert( testFile( FILE_GZIP_BINARY ) );
    return true;
}


// register commands, etc. defined in this module
void initFile() {
    registerUnitTest( testFile );
}


} // end namespace sbl

