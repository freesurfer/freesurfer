#ifndef _SBL_TRACK_H_
#define _SBL_TRACK_H_
#include <sbl/core/File.h>
namespace sbl {


//-------------------------------------------
// TRACK CLASS
//-------------------------------------------


/// The Track class represents the trajectory of a feature point in an image sequence.
class Track {
public:

    /// create a track that will start at the given frame
    explicit Track( int startFrameIndex ) { m_startFrameIndex = startFrameIndex; }

    /// load a track from a file
    explicit Track( File &file );

    /// get position
    inline float x( int frameIndex ) const { return m_x[ frameIndex + m_startFrameIndex ]; }
    inline float y( int frameIndex ) const { return m_y[ frameIndex + m_startFrameIndex ]; }

    /// get current start and end frame indices
    inline int startFrameIndex() const { return m_startFrameIndex; }
    inline int endFrameIndex() const { return m_startFrameIndex + m_x.length() - 1; }

    /// get length of track (in frames)
    inline int length() const { return m_x.length(); }

    /// returns true if track is defined for this frame
    inline int inBounds( int frameIndex ) const { return frameIndex >= m_startFrameIndex && frameIndex <= endFrameIndex(); }

    /// set the track position, assuming frameIndex is at or after startFrameIndex
    void setPosition( int frameIndex, float x, float y );

    /// save track to file
    void save( File &file ) const;

private:

    // the track data
    int m_startFrameIndex;
    VectorF m_x;
    VectorF m_y;

    // disable copy constructor and assignment operator
    Track( const Track &x );
    Track &operator=( const Track &x );
};


//-------------------------------------------
// TRACK SET CLASS
//-------------------------------------------


/// The TrackSet represents a set of tracks (assumed to be for a single image sequence).
class TrackSet {
public:

    // create empty track set 
    TrackSet() {}

    /// the number of tracks
    inline int count() const { return m_tracks.count(); }

    /// access a track
    inline Track &track( int index ) { return m_tracks[ index ]; }
    inline const Track &track( int index ) const { return m_tracks[ index ]; }

    /// add a track to the set; takes ownership of pointer (don't deallocate outside this class)
    inline void add( Track *track ) { m_tracks.append( track ); }

    /// load tracks from binary file
    void load( const String &fileName );

    /// save tracks to binary file
    void save( const String &fileName ) const;

private:

    // the tracks
    Array<Track> m_tracks;

    // disable copy constructor and assignment operator
    TrackSet( const TrackSet &x );
    TrackSet &operator=( const TrackSet &x );
};


} // end namespace sbl
#endif // _SBL_TRACK_H_

