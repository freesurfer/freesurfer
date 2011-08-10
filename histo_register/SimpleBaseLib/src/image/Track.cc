// Licensed under MIT license; see license.txt.

#include <sbl/image/Track.h>
namespace sbl {


//-------------------------------------------
// TRACK CLASS
//-------------------------------------------


/// load a track from a file
Track::Track( File &file ) {
    m_startFrameIndex = file.readInt();
    m_x = file.readVector<float>();
    m_y = file.readVector<float>();
}


/// set the track position, assuming frameIndex is at or after startFrameIndex
void Track::setPosition( int frameIndex, float x, float y ) {
    assertAlways( frameIndex >= m_startFrameIndex );
    int index = frameIndex - m_startFrameIndex;
    while (index >= m_x.length()) {
        m_x.append( 0 );
        m_y.append( 0 );
    }
    m_x[ index ] = x;
    m_y[ index ] = y;    
}


/// save track to file
void Track::save( File &file ) const {
    file.writeInt( m_startFrameIndex );
    file.writeVector( m_x );
    file.writeVector( m_y );
}


//-------------------------------------------
// TRACK SET CLASS
//-------------------------------------------


/// load tracks from binary file
void TrackSet::load( const String &fileName ) {
    File file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess()) {
        m_tracks.reset();
        file.readArray( m_tracks );
    }
}


/// save tracks to binary file
void TrackSet::save( const String &fileName ) const {
    File file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess()) {
        file.writeArray( m_tracks );
    }
}


} // end namespace sbl

