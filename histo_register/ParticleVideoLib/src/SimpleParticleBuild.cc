#include <pvl/SimpleParticleBuild.h>
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/math/KMeans.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/image/MotionFieldUtil.h>
#include <pvl/VarMotionUtil.h>
#include <pvl/SimpleParticleOptimize.h>
#include <math.h>
namespace pvl {


//-------------------------------------------
// ESTIMATE PARTICLE POSITIONS
//-------------------------------------------


/// compute the seperation between particles according to the flow gradient
int particleSeparation( float grad, float maxFlowGradient, int minSeparation, int maxSeparation ) {
    if (grad > maxFlowGradient)
        grad = maxFlowGradient;
    return minSeparation + round( (float) (maxSeparation - minSeparation) * (1.0f - grad / maxFlowGradient) );
}


/// add new particles in gaps between existing particles
void simpleParticleAdd( SimpleParticleSet &particleSet, int frameIndex, const Array<ImageGrayF> &chan, const MotionField &mf, Config &spConf ) {
    int width = mf.width(), height = mf.height();

    // get config parameters
    float maxFlowGradient = spConf.readFloat( "maxFlowGradient" );
    int flowGradientBlur = spConf.readInt( "flowGradientBlur" );
    int minSeparation = spConf.readInt( "minSeparation" );
    int maxSeparation = spConf.readInt( "maxSeparation" );
    int frameBorder = spConf.readInt( "frameBorder" );

    // compute motion gradient
    aptr<ImageGrayF> gradMagSqd = motionGradMagSqd( mf );

    // blur motion gradient
    aptr<ImageGrayF> gradBlur = blurBox( *gradMagSqd, flowGradientBlur );

    // mark existing particles
    ImageGrayU partMask( width, height );
    partMask.clear( 0 );
    for (int i = 0; i < particleSet.count(); i++) {
        const SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex )) {
            int x = round( p.x( frameIndex ));
            int y = round( p.y( frameIndex ));
            if (gradBlur->inBounds( x, y )) {
                int sep = particleSeparation( gradBlur->data( x, y ), maxFlowGradient, minSeparation, maxSeparation );
                drawCircleFilled( partMask, x, y, sep, sep );
            }
        }
    }

    // add new particles
    for (int y = frameBorder; y < height - frameBorder; y++) {
        for (int x = frameBorder; x < width - frameBorder; x++) {
            if (partMask.data( x, y ) == 0) {
                SimpleParticle *p = new SimpleParticle( frameIndex, x, y, chan.count() );
                p->setChannels( frameIndex, chan, NULL );
                particleSet.append( p );
                int sep = particleSeparation( gradBlur->data( x, y ), maxFlowGradient, minSeparation, maxSeparation );
                drawCircleFilled( partMask, x, y, sep, sep );
            }
        }
    }
}


/// prune particles for which current channel values do not match blurred channel values
void simpleParticlePrune( SimpleParticleSet &particleSet, int frameIndex, Config &spConf ) {
    float pruneChannelThresh = spConf.readFloat( "pruneChannelThresh" );
    int pruneCount = 0;
    
    // prune if image appearance doesn't match history
    for (int i = 0; i < particleSet.count(); i++) {
        SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex - 1, frameIndex )) {

            // compute sum of absolute difference over channels
            int count = p.channelCount();
            float sumDiff = 0;
            for (int channelIndex = 0; channelIndex < count; channelIndex++) {
                float diff = p.channel( frameIndex, channelIndex ) - p.channelBlur( frameIndex, channelIndex );
                if (diff < 0) diff = -diff;
                sumDiff += diff;
            }

            // if mean channel diff greater than threshold, prune
            if (sumDiff > pruneChannelThresh * count) {
                p.prune( frameIndex );
                pruneCount++;
            }
        }
    }
    disp( 2, "pruned %d particles", pruneCount ); 
}


/// update particle properties
void simpleParticleUpdate( SimpleParticleSet &particleSet, int frameIndex, const Array<ImageGrayF> &chan, Config &spConf ) {
    float channelTimeSigma = spConf.readFloat( "channelTimeSigma" );

    // compute kernel used for channel smoothing in data term
    VectorF channelBlurKernel = gaussKernel( round( channelTimeSigma * 3 ), channelTimeSigma );

    // update each particle
    for (int i = 0; i < particleSet.count(); i++) {
        SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex - 1, frameIndex ))
            p.setChannels( frameIndex, chan, &channelBlurKernel );
    }
}


/// propagate particles into current frame using motion field;
/// does not propagate particles in occluded regions or regions with large negative divergences
void simpleParticlePropagate( SimpleParticleSet &particleSet, int frameIndex, const MotionField &mfPrev, Config &spConf ) {
    int occCount = 0, boundCount = 0, propCount = 0;
    int width = mfPrev.width(), height = mfPrev.height();

    // get config parameters
    int frameBorder = spConf.readInt( "frameBorder" );
    int divBlurSize = spConf.readInt( "divBlurSize" );
    float divThresh = spConf.readFloat( "divThresh" );

    // compute flow divergence
    aptr<ImageGrayF> div = motionDivergence( mfPrev );

    // blur flow divergence
    aptr<ImageGrayF> divBlur = blurBox( *div, divBlurSize );

    // propagate each particle into current frame
    for (int i = 0; i < particleSet.count(); i++) {
        SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex - 1 )) {
            float xPrev = p.x( frameIndex - 1 );
            float yPrev = p.y( frameIndex - 1 );
            int xPrevRound = round( xPrev );
            int yPrevRound = round( yPrev );
            
            // if in bounds (in previous frame) and is not occluded
            if (mfPrev.inBounds( xPrevRound, yPrevRound )) {
                if (mfPrev.occ( xPrevRound, yPrevRound ) > 128 || divBlur->data( xPrevRound, yPrevRound ) < divThresh) {
                    occCount++;
                } else {
                    float xNew = xPrev + mfPrev.u( xPrevRound, yPrevRound );
                    float yNew = yPrev + mfPrev.v( xPrevRound, yPrevRound );

                    // if in bounds (in dest), propagate
                    if (xNew > frameBorder && xNew < width - frameBorder && yNew > frameBorder && yNew < height - frameBorder) {
                        p.extend();
                        p.setPosition( frameIndex, xNew, yNew );
                        propCount++;
                    } else {
                        boundCount++;
                    }
                }
            }
        }
    }
    disp( 2, "propagated particles: %d (%d occluded, %d out of bounds)", propCount, occCount, boundCount );
}


/// extend current particle set by one frame;
/// to build a particle set, call this function on consecutive frames (starting with frameIndex 0);
/// mf maps frameIndex to frameIndex + 1;
/// mfPrev maps frameIndex - 1 to frameIndex
void simpleParticleFrame( SimpleParticleSet &particleSet, int frameIndex, const ImageColorU &frame, 
                          const MotionField *mfPrev, const MotionField &mf, Config &spConf ) {

    // check frame index
    if (frameIndex != particleSet.maxEndFrame() + 1) {
        warning( "SimpleParticleFrame: invalid frame index: %d (expected next consecutive frame: %d)", frameIndex, particleSet.maxEndFrame() + 1 );
        return;
    }

    // get image channels
    Array<ImageGrayF> channels = extractChannels( frame, spConf );

    // if not first frame
    if (frameIndex) {

        // propagate particles into current frame using motion field
        simpleParticlePropagate( particleSet, frameIndex, *mfPrev, spConf );

        // optimizes particle positions using blurred channel values
        simpleParticleOptimize( particleSet, frameIndex, channels, spConf );
            
        // prune particles for which current channel values do not match blurred channel values
        simpleParticlePrune( particleSet, frameIndex, spConf );

        // update particle properties
        simpleParticleUpdate( particleSet, frameIndex, channels, spConf );
    }

    // add new particles in gaps between existing particles
    simpleParticleAdd( particleSet, frameIndex, channels, mf, spConf );
}


//-------------------------------------------
// CLUSTER PARTICLES
//-------------------------------------------


/// compute distance between particles according to space, appearance, and motion
float particleDistanceSqd( const SimpleParticle &p1, const SimpleParticle &p2, float maxDistSqd,
                           float spaceSigmaSqd, float chanSigmaSqd, float motionSigmaSqd ) {
    int commonStartFrame = max( p1.startFrame(), p2.startFrame() );
    int commonEndFrame = min( p1.endFrame(), p2.endFrame() );
    int commonLength = commonEndFrame - commonStartFrame + 1;
    if (commonLength <= 1) 
        return maxDistSqd;
    float spaceDistSqd = 0, chanDistSqd = 0, motionDistSqd = 0;
    int channelCount = min( p1.channelCount(), p2.channelCount() );

    // loop over frames in common
    for (int frameIndex = commonStartFrame; frameIndex <= commonEndFrame; frameIndex++) {

        // accumulate space distance
        float xDiff = p1.x( frameIndex ) - p2.x( frameIndex );
        float yDiff = p1.y( frameIndex ) - p2.y( frameIndex );
        spaceDistSqd += xDiff * xDiff + yDiff * yDiff;

        // accumulate channel distance
        for (int channelIndex = 0; channelIndex < channelCount; channelIndex++) {
            float chanDiff = p1.channel( frameIndex, channelIndex ) - p2.channel( frameIndex, channelIndex );
            chanDistSqd += chanDiff * chanDiff;
        }

        // accumulate motion distance
        if (frameIndex > commonStartFrame) {
            float uDiff = (p1.x( frameIndex ) - p1.x( frameIndex - 1 ))
                        - (p2.x( frameIndex ) - p2.x( frameIndex - 1 ));
            float vDiff = (p1.y( frameIndex ) - p1.y( frameIndex - 1 ))
                        - (p2.y( frameIndex ) - p2.y( frameIndex - 1 ));
            motionDistSqd += uDiff * uDiff + vDiff * vDiff;
        }
    }

    // scale each distance
    spaceDistSqd /= spaceSigmaSqd * (float) commonLength;
    chanDistSqd /= chanSigmaSqd * (float) commonLength;
    motionDistSqd /= motionSigmaSqd * (float) (commonLength - 1);

    // compute overlap fraction
    int maxLength = max( p1.endFrame() - p1.startFrame() + 1, p2.endFrame() - p2.startFrame() + 1 );
    float overlapFrac = ((float) commonLength - 1.0f) / (float) maxLength; 

    // combine distances (weighted average of distance components, weighting motionDist by overlap)
    return (spaceDistSqd + chanDistSqd + motionDistSqd * overlapFrac) / (2 + overlapFrac);
}


/// display affinity matrix values along particle links
void visualizeParticleAffinity( const SimpleParticleSet &particleSet, int frameIndex, int width, int height, 
                                const MatrixF &affinity ) {
    ImageColorU vis( width, height );
    vis.clear( 255, 255, 255 );

    // draw each link
    aptr<MatrixI> links = createLinks( particleSet, frameIndex, false );
    for (int i = 0; i < links->rows(); i++) {
        int ind1 = links->data( i, 0 ), ind2 = links->data( i, 1 );
        const SimpleParticle &p1 = particleSet.ref( ind1 );
        const SimpleParticle &p2 = particleSet.ref( ind2 );
        if (p1.active( frameIndex ) && p2.active( frameIndex )) {
            int r = 0, g = 0, b = 0;
            colorize( 1.0f - affinity.data( ind1, ind2 ), 0, 1, r, g, b );
            drawLine( vis, round( p1.x( frameIndex )), round( p1.y( frameIndex )), round( p2.x( frameIndex )), round( p2.y( frameIndex )), r, g, b, true );
        }
    }
    dispImage( vis );
}


/// display particle clusters; each cluster has a pre-defined color
void visualizeParticleClusters( const SimpleParticleSet &particleSet, int frameIndex, int width, int height ) {
    ImageColorU vis( width, height );
    vis.clear( 255, 255, 255 );

    // loop over clusters
    for (int k = 0; k < particleSet.clusterCount(); k++) {
        int r = 0, g = 0, b = 0;
        colorizeDiscrete( k, r, g, b );
        const SimpleParticleCluster &pc = particleSet.clusterRef( k );

        // loop over particles in cluster
        for (int i = 0; i < pc.memberCount(); i++) {
            const SimpleParticle &p = particleSet.ref( pc.memberIndex( i ));
            if (p.active( frameIndex ))
                drawCircleFilled( vis, round( p.x( frameIndex )), round( p.y( frameIndex )), 4, r, g, b );
        }
    }
    dispImage( vis );
}


/// create particle clusters using a particle-pair distance measure
void buildParticleClusters( SimpleParticleSet &particleSet, int width, int height, Config &spConf, bool visualize ) {

    // get clustering parameters
    float spaceSigma = spConf.readFloat( "clustSpaceSigma" );
    float chanSigma = spConf.readFloat( "clustChanSigma" );
    float motionSigma = spConf.readFloat( "clustMotionSigma" );
    float combinedSigma = spConf.readFloat( "clustCombinedSigma" );
    float maxParticleDist = spConf.readFloat( "maxParticleDist" );
    int clusterCount = spConf.readInt( "clusterCount" );
    float maxParticleDistSqd = maxParticleDist * maxParticleDist;
    float spaceSigmaSqd = spaceSigma * spaceSigma;
    float chanSigmaSqd = chanSigma * chanSigma;
    float motionSigmaSqd = motionSigma * motionSigma;

    // compute distance matrix
    int particleCount = particleSet.count(), fillCount = 0;
    if (particleCount == 0) 
        return;
    disp( 1, "computing %d by %d distance matrix...", particleCount, particleCount );
    checkCommandEvents();
    float factor = gaussFactor( combinedSigma );
    MatrixF affinity( particleCount, particleCount );
    for (int i = 0; i < particleCount; i++) {
        affinity.data( i, i ) = 0; // set diagonal to 0
        for (int j = i + 1; j < particleCount; j++) {
            const SimpleParticle &p1 = particleSet.ref( i );
            const SimpleParticle &p2 = particleSet.ref( j );
            float distSqd = particleDistanceSqd( p1, p2, maxParticleDistSqd, spaceSigmaSqd, chanSigmaSqd, motionSigmaSqd );
            if (distSqd < maxParticleDistSqd) {
                affinity.data( j, i ) = affinity.data( i, j ) = gauss( distSqd, factor );
                fillCount++;
            } else {
                affinity.data( j, i ) = affinity.data( i, j ) = 0;
            }
        }
    }
    if (visualize) 
        visualizeParticleAffinity( particleSet, 0, width, height, affinity );
    disp( 2, "filled %f entries per particle", (float) fillCount / (float) particleCount );

    // compute space for clustering 
    disp( 1, "running spectral transform to find %d-dimensional space...", clusterCount );
    checkCommandEvents();
    aptr<MatrixF> clusterSpace = spectralTransform( affinity, clusterCount, true );
    for (int i = 0; i < particleCount; i++) {
        SimpleParticle &p = particleSet.ref( i );
        p.setClusterSpacePosition( clusterSpace->row( i ) );
    }

    // remove previous clusters
    particleSet.clearClusters();

    // run k-means to create clusters
    disp( 1, "running k-means to find %d clusters...", clusterCount );
    checkCommandEvents();
    aptr<MatrixF> means;
    aptr<VectorI> assign;
    kMeans( *clusterSpace, clusterCount, means, assign );
    for (int j = 0; j < clusterCount; j++) {
        VectorF clusterMean = means->row( j );
        VectorI members;
        for (int i = 0; i < particleCount; i++)
            if (assign->data( i ) == j)
                members.append( i );
        disp( 2, "cluster %d: %d members", j, members.length() );
        SimpleParticleCluster *pc = new SimpleParticleCluster( clusterMean, members ); 
        particleSet.appendCluster( pc );
    }
    if (visualize)
        visualizeParticleClusters( particleSet, 0, width, height );
    disp( 1, "done" );
}


} // end namespace pvl

