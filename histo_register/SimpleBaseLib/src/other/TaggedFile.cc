// Licensed under MIT license; see license.txt.

#include <sbl/other/TaggedFile.h>
namespace sbl {


/// write length and raw bytes
void TaggedFile::writeString( const String &s ) {
    assertDebug( m_file );
    writeTagInfo( s.length() );
    writeBlock( s.c_str(), s.length() );
}


/// reads length and raw bytes
String TaggedFile::readString() {
    assertDebug( m_file );
    int len = readTagInfo();
    if (len < 0 || len > 100000000) {
        warning( "invalid string" );
        return String();
    }

    // read string data
    char *data = new char[ len + 1 ];
    readBlock( data, len );
    data[ len ] = 0;
    String s( data );
    delete [] data;
    return s;
}


/// read and discard a value from a tagged file
void TaggedFile::skipTag( int tag ) { 

    // if end of section, don't need to skip anything
    if (tag == TAG_END_SECTION)
        return;

    // get the length of data to skip over
    int len = readTagInfo(); 
    if (tag > TAG_OFFSET_ARRAY)
        fatalError( "need to implement skip array" );

    // read the skipped data
    char *data = new char[ len ]; 
    readBlock( data, len ); 
    delete [] data; 
}


} // end namespace sbl

