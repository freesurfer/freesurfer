#ifndef _SIMPLE_PARTICLE_DATA_H_
#define _SIMPLE_PARTICLE_DATA_H_
#include <sbl/core/Array.h>
#include <sbl/core/File.h>
#include <sbl/core/Config.h>
#include <sbl/core/Pointer.h>
#include <sbl/math/Vector.h>
#include <sbl/image/Image.h> 
#include <sbl/system/Timer.h>
using namespace sbl;
namespace pvl {


/// The SimpleParticle class represents a single point trajectory.
class SimpleParticle {
public:

    // constructor / destructor
    SimpleParticle( int frameIndex, int x, int y, int channelCount );

    // access position
    inline float x( int frameIndex ) const { return m_x[ frameIndex - m_startFrame ]; };
    inline float y( int frameIndex ) const { return m_y[ frameIndex - m_startFrame ]; };
    void setPosition( int frameIndex, float x, float y );

    // extend by one frame
    void extend();

    // prune at the given frame (so particle terminates at frameIndex - 1)
    void prune( int frameIndex );

    // access frame start/end
    inline int startFrame() const { return m_startFrame; }
    inline int endFrame() const { return m_endFrame; }
    inline int length() const { return m_endFrame - m_startFrame + 1; }
    inline bool active( int frameIndex ) const { return frameIndex >= m_startFrame && frameIndex <= m_endFrame; }
    inline bool active( int startFrame, int endFrame ) const { return startFrame >= m_startFrame && endFrame <= m_endFrame; }

    // get/set particle channels
    inline float channel( int frameIndex, int channelIndex ) const { return m_channel[ channelIndex ][ frameIndex - m_startFrame ]; }
    inline void setChannel( int frameIndex, int channelIndex, float val ) { m_channel[ channelIndex ][ frameIndex - m_startFrame ] = val; }
    void setChannels( int frameIndex, const Array<ImageGrayF> &chan, const VectorF *channelBlurKernel );
    inline int channelCount() const { return m_channel.count(); }

    /// position of particle in (an abstract) cluster space
    inline const VectorF &clusterSpacePosition() const { return m_clusterSpacePosition; }

    /// set the position of the particle in (an abstract) cluster space
    inline void setClusterSpacePosition( const VectorF &clusterSpacePosition ) { m_clusterSpacePosition = clusterSpacePosition; }

    /// access temporally-blurred channel values
    float channelBlur( int frameIndex, int channelIndex ) const;

    /// load particle from file
    bool load( File &file );

    /// save particle to file
    void save( File &file ) const;

    /// save particle info in simplified text format
    void saveSimple( File &file ) const;

private:

    // particle position
    VectorF m_x;
    VectorF m_y;

    // start/end frame index
    int m_startFrame;
    int m_endFrame;

    // position of particle in cluster space
    VectorF m_clusterSpacePosition;

    // particle appearance in each frame (observed channel values for the particle)
    Array<VectorF> m_channel; 

    // blurred channel values
    VectorF m_channelBlur;

    // disable copy constructor and assignment operator
    SimpleParticle( const SimpleParticle &x );
    SimpleParticle &operator=( const SimpleParticle &x );
};


/// The SimpleParticleCluster class represents cluster membership within a SimpleParticleSet.
class SimpleParticleCluster {
public:

    // basic constructors
    SimpleParticleCluster( const VectorF &mean, const VectorI &particleIndex );
    SimpleParticleCluster() {}

    // access members
    inline int memberCount() const { return m_particleIndex.length(); };
    inline int memberIndex( int index ) const { return m_particleIndex[ index ]; };

    // cluster center in cluster space
    inline const VectorF &mean() const { return m_mean; };

    // load/save particle cluster
    bool load( File &file );
    void save( File &file ) const;

private:

    // indices of particles in this cluster
    VectorI m_particleIndex;

    // cluster center in cluster space
    VectorF m_mean;

    // disable copy constructor and assignment operator
    SimpleParticleCluster( const SimpleParticleCluster &x );
    SimpleParticleCluster &operator=( const SimpleParticleCluster &x );
};


// kinds of timers stored in particle set class
enum SimpleParticleTimerType {
    TIMER_OPT,
    TIMER_SOLVER,
    TIMER_COUNT
};


/// The SimpleParticleSet class represents a collection of particle trajectories.
class SimpleParticleSet {
public:

    // basic constructor 
    SimpleParticleSet() {}

    /// access particles
    inline SimpleParticle &ref( int index ) { return m_simpleParticleSet[ index ]; };
    inline const SimpleParticle &ref( int index ) const { return m_simpleParticleSet[ index ]; };
    inline int count() const { return m_simpleParticleSet.count(); };
    int activeCount( int frameIndex ) const;
    int activeCount( int startFrame, int endFrame ) const;

    /// min start frame across set of particles
    int minStartFrame() const;
    
    /// max end frame across set of particles
    int maxEndFrame() const;

    /// access clusters
    inline SimpleParticleCluster &clusterRef( int index ) { return m_simpleParticleClusterSet[ index ]; };
    inline const SimpleParticleCluster &clusterRef( int index ) const { return m_simpleParticleClusterSet[ index ]; };
    inline int clusterCount() const { return m_simpleParticleClusterSet.count(); };
    void clearClusters();

    /// add a particle (takes ownership of pointer; do not deallocate outside)
    inline void append( SimpleParticle *p ) { m_simpleParticleSet.append( p ); }

    /// add a particle cluster (takes ownership of pointer; do not deallocate outside)
    inline void appendCluster( SimpleParticleCluster *pc ) { m_simpleParticleClusterSet.append( pc ); }

    /// access other particle set properties    
    inline Timer &timeRef( SimpleParticleTimerType spt ) { return m_timeSum[ spt ]; };

    /// file name used for last load/save (empty if none)
    inline const String &fileName() const { return m_fileName; };

    /// load particle set from file; returns true on success
    bool load( const String &fileName );

    /// save particle set to file
    void save( const String &fileName );
    void save() const;

    /// saves the particle set in simple text format
    void saveSimple( const String &fileName ) const;

    /// draw the particles on (a copy of) the given image
    aptr<ImageColorU> draw( int frameIndex, const ImageColorU &frame );

    /// display summary statistics about particles
    void dispStats( int indent );

private:

    // the particles
    Array<SimpleParticle> m_simpleParticleSet;

    // the clusters (if any)
    Array<SimpleParticleCluster> m_simpleParticleClusterSet;

    // current file name (empty if none)
    String m_fileName;

    // algorithm timing info
    Timer m_timeSum[ TIMER_COUNT ];

    // disable copy constructor and assignment operator
    SimpleParticleSet( const SimpleParticleSet &x );
    SimpleParticleSet &operator=( const SimpleParticleSet &x );
};


} // end namespace pvl
#endif // _SIMPLE_PARTICLE_DATA_H_

