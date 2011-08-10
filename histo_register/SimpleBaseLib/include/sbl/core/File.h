#ifndef _SBL_FILE_H_
#define _SBL_FILE_H_
#include <sbl/core/String.h>
#include <sbl/math/Vector.h>
#include <stdio.h>
namespace sbl {


// register commands, etc. defined in this module
void initFile();


/// file open modes
enum FileOpenMode {
    FILE_READ,
    FILE_WRITE,
    FILE_APPEND
};


/// file open types
enum FileOpenType {
    FILE_TEXT,
    FILE_BINARY,
    FILE_GZIP_BINARY,
    FILE_TEXT_FROM_BINARY
};


// fix(clean): pay better attention to fscanf/fread return values
void readError();
#define CHECK_READ( x ) if ((x) == 0) readError()


/// The File class is a fast and simple wrapper for the stdio FILE objects.
class File {
public:

    /// open a file
    explicit File( const String &fileName = "", FileOpenMode mode = FILE_READ, FileOpenType type = FILE_TEXT );

    /// close the file
    virtual ~File();

    /// open a file
    void open( const String &fileName, FileOpenMode mode, FileOpenType type );

    /// close the file
    void close();

    /// returns true if file was opened for binary i/o
    inline bool binary() const { return m_binary; }

    /// returns true if file was opened successfully
    inline bool openSuccess() const { return m_file != NULL || m_gzFile != NULL; }

    /// flush writes to disk
    void flush();

    /// returns true if end of file reached on last read; returns true if file not successfully opened
    bool endOfFile();

    /// seek to the given position
    void seek( int offset, bool relative );

    /// returns the current position (assumes file position fits in int)
    int tell();

    /// write various basic types
    virtual inline void writeBool( bool val ) { writeInt( (int) val ); }
    virtual inline void writeUChar( unsigned char val ) { if (m_binary) writeBlock( &val, sizeof(unsigned char) ); else fprintf( m_file, "%d ", val ); }
    virtual inline void writeUShort( unsigned short val ) { if (m_binary) writeBlock( &val, sizeof(unsigned short) ); else fprintf( m_file, "%d ", val ); }
    virtual inline void writeInt( int val ) { if (m_binary) writeBlock( &val, sizeof(int) ); else fprintf( m_file, "%d ", val ); }
    virtual inline void writeFloat( float val ) { if (m_binary) writeBlock( &val, sizeof(float) ); else fprintf( m_file, "%f ", val ); }
    virtual inline void writeDouble( double val ) { if (m_binary) writeBlock( &val, sizeof(double) ); else fprintf( m_file, "%lf ", val ); }
    virtual void writeString( const String &s );
    void writeRawString( const String &s );

    /// write formatted string; equivalent of fprintf (printf for files)
    void writeF( const char *str, ... );

    /// read various basic types
    virtual inline bool readBool() { return readInt() ? true : false; }
    virtual inline unsigned char readUChar() { unsigned char val = 0; if (m_binary) readBlock( &val, sizeof(unsigned char) ); else { int valInt = 0; CHECK_READ( fscanf( m_file, "%d", &valInt ) ); val = (unsigned char) valInt; } return val; } // fix(clean): what is correct format specifier? not %c, not %d, not %u?
    virtual inline unsigned short readUShort() { unsigned short val = 0; if (m_binary) readBlock( &val, sizeof(unsigned short) ); else CHECK_READ( fscanf( m_file, "%hu", &val ) ); return val; }
    virtual inline int readInt() { int val = 0; if (m_binary) readBlock( &val, sizeof(int) ); else CHECK_READ( fscanf( m_file, "%d", &val ) ); return val; }
    virtual inline float readFloat() { float val = 0; if (m_binary) readBlock( &val, sizeof(float) ); else CHECK_READ( fscanf( m_file, "%f", &val ) ); return val; }
    virtual inline double readDouble() { double val = 0; if (m_binary) readBlock( &val, sizeof(double) ); else CHECK_READ( fscanf( m_file, "%lf", &val ) ); return val; }
    virtual String readString();

    /// read up to end of line; strip trailing end-of-line characters
    String readLine();

    /// read/write block of data (assumes binary)
    void writeBlock( const void *data, int byteCount );
    void readBlock( void *data, int byteCount );

    /// write array (assumes type T has a .save() method)
    template<typename T> void writeArray( const Array<T> &a ) {
        writeInt( a.count() );
        for (int i = 0; i < a.count(); i++) 
            a[ i ].save( *this );
    }

    /// read array (assumes type T has a file load constructor)
    template<typename T> void readArray( Array<T> &a ) {
        int count = readInt();
        if (count < 0) {
            warning( "invalid object count in file" );
        } else {
            for (int i = 0; i < count; i++) 
                a.append( new T( *this ) );
        }
    }

    /// write a vector, with length
    template<typename T> void writeVector( const Vector<T> &v ) {
        assertAlways( binary() );
        writeInt( v.length() );
        writeBlock( v.dataPtr(), v.length() * sizeof(T) );
    }

    /// read a vector, with length (returns copy)
    template<typename T> Vector<T> readVector() {
        assertAlways( binary() );
        Vector<T> v;
        int len = readInt();
        if (len < 0) {
            warning( "invalid vector in file" );
        } else { 
            v.setLength( len );
            readBlock( v.dataPtr(), v.length() * sizeof(T) );
        }
        return v;
    }

    /// read a vector, with length (avoids copy)
    template<typename T> void readVector( Vector<T> &v ) {
        assertAlways( binary() );
        int len = readInt();
        if (len < 0) {
            warning( "invalid vector in file" );
        } else { 
            v.setLength( len );
            readBlock( v.dataPtr(), v.length() * sizeof(T) );
        }
    }

protected:

    // the stdio FILE handle
    FILE *m_file;

    // zlib file handle (if any)
    void *m_gzFile;

private:

    // true if file opened in binary mode
    bool m_binary;

    // disable copy constructor and assignment operator
    File( const File &x );
    File &operator=( const File &x );
};


} // end namespace sbl
#endif // _SBL_FILE_H_

