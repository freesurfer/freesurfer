#include <pvl/SimpleParticleData.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MatrixUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageDraw.h> // for SimpleParticleSet::draw()
namespace pvl {


//-------------------------------------------
// SIMPLE PARTICLE CLASS
//-------------------------------------------


// basic constructor
SimpleParticle::SimpleParticle( int frameIndex, int x, int y, int channelCount ) {
    m_startFrame = frameIndex;
    m_endFrame = frameIndex;
    m_x.append( (float) x );
    m_y.append( (float) y );

    // allocate channel data
    for (int i = 0; i < channelCount; i++) {
        VectorF *channel = new VectorF( 1 );
        channel->data( 0 ) = 0;
        m_channel.append( channel );
    }
    m_channelBlur.setLength( channelCount );
}


/// set position of particle in given frame
void SimpleParticle::setPosition( int frameIndex, float x, float y ) {
    if (frameIndex < m_startFrame || frameIndex > m_endFrame)
        fatalError( "SimpleParticle::SetPosition: invalid frame index" );
    m_x[ frameIndex - m_startFrame ] = x;
    m_y[ frameIndex - m_startFrame ] = y;
}


/// prune at the given frame (so particle terminates at frameIndex - 1)
void SimpleParticle::prune( int frameIndex ) {
    if (frameIndex > m_endFrame || frameIndex <= m_startFrame)
        fatalError( "SimpleParticle::Prune: invalid frame index" );
    m_endFrame = frameIndex - 1;
}


/// get time-blurred channel value 
float SimpleParticle::channelBlur( int frameIndex, int channelIndex ) const {
    if (frameIndex <= m_startFrame || frameIndex > m_endFrame)
        fatalError( "Particle::ChannelBlur: frame index out of bounds" );
    return m_channelBlur[ channelIndex ];
}


/// store channel values from image into particle
void SimpleParticle::setChannels( int frameIndex, const Array<ImageGrayF> &chan, const VectorF *channelBlurKernel ) {
    float xp = x( frameIndex ), yp = y( frameIndex );
    if (chan[ 0 ].inBounds( xp, yp ))
        for (int i = 0; i < m_channel.count(); i++) {
            float val = chan[ i ].interp( xp, yp );
            setChannel( frameIndex, i, val );
            m_channelBlur[ i ] = val;
        }

    // compute blurred channel values (including values just set)
    if (channelBlurKernel) {
        int tableRadius = (channelBlurKernel->length() - 1) / 2;
        int length = m_channel[ 0 ].length();
        if (length == 0)
            fatalError( "SimpleParticle::SetChannels: length is 0" );
        if (frameIndex < m_startFrame || frameIndex > m_endFrame)
            fatalError( "SimpleParticle::SetChannels: frame index out of bounds" );

        // loop over channels, applying blur kernel
        for (int i = 0; i < m_channel.count(); i++) {
            float sum = 0, sumWt = 0;
            for (int j = -tableRadius; j <= tableRadius; j++) {
                int srcInd = frameIndex + j - m_startFrame; // index into channel vectors
                if (srcInd >= 0 && srcInd < length) { 
                    float wt = channelBlurKernel->data( j + tableRadius );
                    sum += wt * m_channel[ i ][ srcInd ];
                    sumWt += wt;
                }
            }
            m_channelBlur[ i ] = sum / sumWt;
        }
    }
}


/// extend by one frame
void SimpleParticle::extend() {
    m_x.append( x( m_endFrame ) );
    m_y.append( y( m_endFrame ) );
    for (int i = 0; i < m_channel.count(); i++) 
        m_channel[ i ].append( 0 );
    m_endFrame++;
}


/// load particle from file
bool SimpleParticle::load( File &file ) {
    // fix(later): check for errors
    m_startFrame = file.readInt();
    m_endFrame = file.readInt();
    m_x = file.readVector<float>();
    m_y = file.readVector<float>();
    m_clusterSpacePosition = file.readVector<float>();
    m_channel.reset();
    int channelCount = file.readInt();
    for (int i = 0; i < channelCount; i++) 
        m_channel.appendCopy( file.readVector<float>() );
    m_channelBlur = file.readVector<float>();
    return true;
}


/// save particle to file
void SimpleParticle::save( File &file ) const {
    file.writeInt( m_startFrame );
    file.writeInt( m_endFrame );
    file.writeVector( m_x );
    file.writeVector( m_y );
    file.writeVector( m_clusterSpacePosition );
    file.writeInt( m_channel.count() );
    for (int i = 0; i < m_channel.count(); i++) 
        file.writeVector( m_channel[ i ] );
    file.writeVector( m_channelBlur );
}


/// save particle info in simplified text format
void SimpleParticle::saveSimple( File &file ) const {
    file.writeInt( m_startFrame );
    file.writeInt( m_endFrame );
    int length = m_endFrame - m_startFrame + 1;
    file.writeInt( length );
    for (int i = 0; i < length; i++) {
        file.writeFloat( m_x[ i ] );
        file.writeFloat( m_y[ i ] );
    }
    file.writeRawString( "\n" );
}


//-------------------------------------------
// SIMPLE PARTICLE CLUSTER CLASS
//-------------------------------------------


/// basic constructor
SimpleParticleCluster::SimpleParticleCluster( const VectorF &mean, const VectorI &particleIndex ) {
    m_mean = mean;
    m_particleIndex = particleIndex;
}


/// load cluster from file
bool SimpleParticleCluster::load( File &file ) {
    // fix(later): check for errors
    m_particleIndex = file.readVector<int>();
    m_mean = file.readVector<float>();
    return true;
}


/// save cluster to file
void SimpleParticleCluster::save( File &file ) const {
    file.writeVector( m_particleIndex );
    file.writeVector( m_mean );
}


//-------------------------------------------
// SIMPLE PARTICLE SET CLASS
//-------------------------------------------


/// number of particles active in given frame
int SimpleParticleSet::activeCount( int frameIndex ) const {
    int activeCount = 0;
    for (int i = 0; i < m_simpleParticleSet.count(); i++) 
        if (m_simpleParticleSet[ i ].active( frameIndex ))
            activeCount++;
    return activeCount;
}


/// number of particles active over entire range of frames
int SimpleParticleSet::activeCount( int startFrame, int endFrame ) const {
    int activeCount = 0;
    for (int i = 0; i < m_simpleParticleSet.count(); i++) 
        if (m_simpleParticleSet[ i ].active( startFrame, endFrame ))
            activeCount++;
    return activeCount;
}


/// min start frame across set of particles
int SimpleParticleSet::minStartFrame() const {
    int minStartFrame = -1;
    for (int i = 0; i < m_simpleParticleSet.count(); i++) {
        const SimpleParticle &p = m_simpleParticleSet[ i ];
        if (p.startFrame() < minStartFrame || minStartFrame == -1)
            minStartFrame = p.startFrame();
    }
    return minStartFrame;
}


/// max end frame across set of particles
int SimpleParticleSet::maxEndFrame() const {
    int maxEndFrame = -1;
    for (int i = 0; i < m_simpleParticleSet.count(); i++) {
        const SimpleParticle &p = m_simpleParticleSet[ i ];
        if (p.endFrame() > maxEndFrame)
            maxEndFrame = p.endFrame();
    }
    return maxEndFrame;
}


/// deallocate any existing clusters
void SimpleParticleSet::clearClusters() {
    m_simpleParticleClusterSet.reset();
}


/// load particle set from file; returns true on success
bool SimpleParticleSet::load( const String &fullFileName ) {
    File file( fullFileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess()) {

        // read header
        char header[10];
        file.readBlock( header, 4 );
        if (strcmp( header, "sps" )) {
            warning( "invalid particle set header" );
            return false;
        }

        // clear existing data
        m_simpleParticleSet.reset();
        m_simpleParticleClusterSet.reset();
        m_fileName = fullFileName;

        // read particles
        int particleCount = file.readInt();
        if (particleCount < 0 || particleCount > 1000000) {
            warning( "invalid particle count" );
            return false;
        }
        for (int i = 0; i < particleCount; i++) {
            SimpleParticle *sp = new SimpleParticle( 0, 0, 0, 0 ); // fix(clean): create simple constructor
            if (sp->load( file )) {
                m_simpleParticleSet.append( sp );
            } else {
                warning( "invalid particle" );
                delete sp;
                return false;
            }
        }

        // read clusters
        int clusterCount = file.readInt();
        if (clusterCount < 0 || clusterCount > 1000000) {
            warning( "invalid cluster count" );
            return false;
        }
        for (int j = 0; j < clusterCount; j++) {
            SimpleParticleCluster *spc = new SimpleParticleCluster();
            if (spc->load( file )) {
                m_simpleParticleClusterSet.append( spc );
            } else {
                warning( "invalid particle cluster" );
                delete spc;
                return false;
            }
        }
    }
    return true;
}


/// saves the particle set and stores the specified file name
void SimpleParticleSet::save( const String &fullFileName ) {
    m_fileName = fullFileName;
    save();
}


/// saves the particle set to previously used file name
void SimpleParticleSet::save() const {
    if (m_fileName.length()) {
        File file( m_fileName, FILE_WRITE, FILE_BINARY );
        if (file.openSuccess()) {

            // write header
            file.writeBlock( "sps", 4 );

            // write particles
            file.writeInt( m_simpleParticleSet.count() );
            for (int i = 0; i < m_simpleParticleSet.count(); i++)
                m_simpleParticleSet[ i ].save( file );

            // write clusters
            file.writeInt( m_simpleParticleClusterSet.count() );
            for (int i = 0; i < m_simpleParticleClusterSet.count(); i++)
                m_simpleParticleClusterSet[ i ].save( file );
        }
    }
}


/// saves the particle set in simple text format
void SimpleParticleSet::saveSimple( const String &fileName ) const {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {

        // write particles
        file.writeInt( m_simpleParticleSet.count() );
        file.writeRawString( "\n" );
        for (int i = 0; i < m_simpleParticleSet.count(); i++)
            m_simpleParticleSet[ i ].saveSimple( file );
    }
}


/// draw the particles on (a copy of) the given image
aptr<ImageColorU> SimpleParticleSet::draw( int frameIndex, const ImageColorU &frame ) {
    aptr<ImageColorU> vis( new ImageColorU( frame ));
    for (int i = 0; i < m_simpleParticleSet.count(); i++) {
        const SimpleParticle &particle = m_simpleParticleSet[ i ];
        if (particle.active( frameIndex )) {
            drawCircleFilled( *vis, round( particle.x( frameIndex ) ), round( particle.y( frameIndex ) ), 2, 255, 0, 0 );
        }
    }
    return vis;
}


/// display summary statistics about particles
void SimpleParticleSet::dispStats( int indent ) {
    disp( indent, "particle count: %d", count() ); 
    disp( indent, "cluster count: %d", clusterCount() );
    disp( indent, "min start frame: %d", minStartFrame() );
    disp( indent, "max end frame: %d", maxEndFrame() );
    VectorD length;
    for (int i = 0; i < m_simpleParticleSet.count(); i++) {
        const SimpleParticle &sp = m_simpleParticleSet[ i ];
        length.append( sp.length() );
    }
    disp( indent, "particle lifetime min: %3.1f, mean: %5.3f, max: %3.1f", 
        length.min(), length.mean(), length.max() );
    if (clusterCount()) {
        VectorD clusterSize;
        for (int i = 0; i < clusterCount(); i++) {
            const SimpleParticleCluster &pc = clusterRef( i );
            clusterSize.append( pc.memberCount() );
        }
        disp( indent, "particle cluster size min: %3.1f, mean: %5.3f, max: %3.1f", 
            clusterSize.min(), clusterSize.mean(), clusterSize.max() );
    }
}


} // end namespace pvl

