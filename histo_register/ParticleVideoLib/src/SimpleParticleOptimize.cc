#include <pvl/SimpleParticleOptimize.h>
#include <sbl/core/Command.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/Triangulation.h>
#include <pvl/SparseSystem.h>
#include <pvl/VarMotionUtil.h>
#include <math.h>
namespace pvl {


// offsets to index X and Y coordinates (for variables and equations)
#define _X_ 0
#define _Y_ 1


/// returns a matrix in which each row is a (particleIndex1, particleIndex2) indicating a link between a pair of particle;
/// if activePrev set, requires that particle be active in previous frame
aptr<MatrixI> createLinks( const SimpleParticleSet &particleSet, int frameIndex, bool activePrev ) {
    int minFrameIndex = activePrev ? frameIndex - 1 : frameIndex;

    // get particles that are active in this frame
    int activeCount = particleSet.activeCount( minFrameIndex, frameIndex );
    if (activeCount < 2)
        fatalError( "CreateLinks: not enough active particles" );
    int activeIndex = 0;
    MatrixF points( activeCount, 2 );
    VectorI index( activeCount ); // maps active index to particle index
    for (int i = 0; i < particleSet.count(); i++) {
        const SimpleParticle &p = particleSet.ref( i );
        if (p.active( minFrameIndex, frameIndex )) {
            points.data( activeIndex, 0 ) = p.x( frameIndex );
            points.data( activeIndex, 1 ) = p.y( frameIndex );
            index.data( activeIndex ) = i;
            activeIndex++;
        }
    }

    // perform triangulation (if enough particles)
    aptr<MatrixI> links = triangulateIndex( points );

    // transform active index to particle index
    for (int i = 0; i < links->rows(); i++) {
        links->data( i, 0 ) = index.data( links->data( i, 0 ) );
        links->data( i, 1 ) = index.data( links->data( i, 1 ) );
    }
    return links; 
}


/// draw links and particles 
void plotLinks( const SimpleParticleSet &particleSet, const MatrixI &links, int frameIndex, int width, int height ) {

    // draw each link
    ImageColorU vis( width, height );
    vis.clear( 255, 255, 255 );
    for (int i = 0; i < links.rows(); i++) {
        int index1 = links.data( i, 0 );
        int index2 = links.data( i, 1 );
        if (index1 == index2)
            disp( 1, "self link" );
        if (index1 < 0 || index1 >= particleSet.count() || index2 < 0 || index2 >= particleSet.count())
            disp( 1, "invalid link index" );
        const SimpleParticle &p1 = particleSet.ref( index1 );
        const SimpleParticle &p2 = particleSet.ref( index2 );
        drawLine( vis, round( p1.x( frameIndex )), round( p1.y( frameIndex )), round( p2.x( frameIndex )), round( p2.y( frameIndex )), 0, 0, 255, false );
    }

    // draw each particle
    for (int i = 0; i < particleSet.count(); i++) {
        const SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex ))
            drawCross( vis, round( p.x( frameIndex )), round( p.y( frameIndex )), 2, 255, 0, 0, true );
    }

    // display image
    dispImage( vis );
}


/// returns map from point index (active index) to particle index
aptr<VectorI> pointMap( const SimpleParticleSet &particleSet, int frameIndex ) {
    aptr<VectorI> ptMap( new VectorI( 0 ) );
    for (int i = 0; i < particleSet.count(); i++) {
        const SimpleParticle &p = particleSet.ref( i );
        if (p.active( frameIndex - 1, frameIndex ))
            ptMap->append( i );
    }
    return ptMap;
}


/// returns map from particle index to point index (active index)
aptr<VectorI> pointMapInverse( const SimpleParticleSet &particleSet, int frameIndex ) {
    int particleCount = particleSet.count();
    aptr<VectorI> ptMapInv( new VectorI( particleCount ) );
    ptMapInv->clear( -1 );
    int index = 0;
    for (int i = 0; i < particleCount; i++) {
        const SimpleParticle &p = particleSet.ref( i );        
        if (p.active( frameIndex - 1, frameIndex ))
            ptMapInv->data( i ) = index++;
    }
    return ptMapInv;
}


/// check point map and inverse point map used to construct sparse system 
void checkPointMap( const SimpleParticleSet &particleSet, int frameIndex, const VectorI &ptMap, const VectorI &ptMapInv ) {

    // check forward map
    for (int i = 0; i < ptMap.length(); i++) {
        int particleIndex = ptMap.data( i );
        if (particleIndex < 0 || particleIndex >= particleSet.count())
            fatalError( "CheckPointMap: invalid forward particle index" );
        if (particleSet.ref( particleIndex ).active( frameIndex - 1, frameIndex ) == false)
            fatalError( "CheckPointMap: inactive particle" );
        if (ptMapInv[ particleIndex ] != i)
            fatalError( "CheckPointMap: forward doesn't match inverse" );
    }

    // check inverse map
    if (ptMapInv.length() != particleSet.count())
        fatalError( "CheckPointMap: incorrect inverse length" );
    for (int i = 0; i < ptMapInv.length(); i++) {
        int pointIndex = ptMapInv[ i ];
        const SimpleParticle &p = particleSet.ref( i );
        if (pointIndex == -1) {
            if (p.active( frameIndex - 1, frameIndex ))
                fatalError( "CheckPointMap: supposed to be inactive" );
        } else {
            if (p.active( frameIndex - 1, frameIndex ) == false)
                fatalError( "CheckPointMap: supposed to be active" );
            if (ptMap[ pointIndex ] != i)
                fatalError( "CheckPointMap: inverse doesn't match forward" );
        }
    }
}


/// applies a position offset (dx, dy) to each particle specified in the ptMap
void updateParticlePositions( SimpleParticleSet &particleSet, int frameIndex, const VectorI &ptMap, const VectorF &dx, const VectorF &dy, float maxMove ) {
    int pointCount = ptMap.length();
    if (pointCount != dx.length() || pointCount != dy.length())
        fatalError( "UpdateParticlePositions" );
    for (int i = 0; i < pointCount; i++) {
        SimpleParticle &p = particleSet.ref( ptMap[ i ] );
        float xNew = p.x( frameIndex ) + bound( dx[ i ], -maxMove, maxMove );
        float yNew = p.y( frameIndex ) + bound( dy[ i ], -maxMove, maxMove );
        p.setPosition( frameIndex, xNew, yNew );
    }
}


/// adds an equation of the form: factor * (var[varIndex1] - var[varIndex2]) == 0
void addPairEquation( SparseSystem &sparseSystem, int eqnIndex, int varIndex1, int varIndex2, float factor ) {
    if (eqnIndex < 0 || varIndex1 < 0 || varIndex2 < 0)
        fatalError( "AddPairEquation: invalid index" );

    // add x terms
    if (sparseSystem.moreTermsAllowed( eqnIndex + _X_ ) >= 2) {
        sparseSystem.addTerm( eqnIndex + _X_, varIndex1 + _X_, factor );
        sparseSystem.addTerm( eqnIndex + _X_, varIndex2 + _X_, -factor );
    }

    // add y terms
    if (sparseSystem.moreTermsAllowed( eqnIndex + _Y_ ) >= 2) {
        sparseSystem.addTerm( eqnIndex + _Y_, varIndex1 + _Y_, factor );
        sparseSystem.addTerm( eqnIndex + _Y_, varIndex2 + _Y_, -factor );
    }
}


/// adds data equation to sparse system
void addDataEquation( SparseSystem &sparseSystem, int eqnIndex, 
                      float ixix, float ixiy, float iyiy, float ixiz, float iyiz, 
                      float factor, float condFactor ) {

    sparseSystem.addTerm( eqnIndex + _X_, eqnIndex + _X_, factor * ixix + condFactor );
    sparseSystem.addTerm( eqnIndex + _X_, eqnIndex + _Y_, factor * ixiy );
    sparseSystem.addToB( eqnIndex + _X_, -factor * ixiz ); // negate because subtract from both sides

    sparseSystem.addTerm( eqnIndex + _Y_, eqnIndex + _Y_, factor * iyiy + condFactor );
    sparseSystem.addTerm( eqnIndex + _Y_, eqnIndex + _X_, factor * ixiy );
    sparseSystem.addToB( eqnIndex + _Y_, -factor * iyiz ); // negate because subtract from both sides
}


/// robust distance function from Brox et al.
float spPsi( float diffSqd ) {
    return sqrtf( diffSqd + 1e-6f );
}


/// derivative of robust distance function from Brox et al.
float spPsiDeriv( float diffSqd ) {
    return 0.5f / sqrtf( diffSqd + 1e-6f );
}


/// visualize offsets from inner loop
void visualizeOffsets( const SimpleParticleSet &particleSet, const VectorI &ptMap, int frameIndex, const Array<ImageGrayF> &chanScaled, const VectorF &dx, const VectorF &dy ) {
    float arrowScale = 10;

    // store frame
    aptr<ImageGrayU> frameGray = toUChar( chanScaled.ref( 0 ) );
    aptr<ImageColorU> frame = toColor( *frameGray );

    // store stats
    int count = 0;
    float dxSum = 0, dySum = 0;
    float dxMin = 1e8, dxMax = -1e8;
    float dyMin = 1e8, dyMax = -1e8;

    // loop over points
    for (int i = 0; i < ptMap.length(); i++) {
        int particleIndex = ptMap[ i ];
        const SimpleParticle &p = particleSet.ref( particleIndex );

        // gather stats
        float dxVal = dx[ i ];
        float dyVal = dy[ i ];
        dxSum += dxVal;
        dySum += dyVal;
        if (dxVal > dxMax) dxMax = dxVal;
        if (dxVal < dxMin) dxMin = dxVal;
        if (dyVal > dyMax) dyMax = dyVal;
        if (dyVal < dyMin) dyMin = dyVal;
        count++;

        // draw scaled arrow
        float x = p.x( frameIndex );
        float y = p.y( frameIndex );
        drawArrow( *frame, round( x ), round( y ), round( x + dxVal * arrowScale ), round( y + dyVal * arrowScale ), 255, 0, 0, true );
    }

    // display stats
    if (count) {
        float dxMean = dxSum / (float) count;
        float dyMean = dySum / (float) count;
        disp( 3, "%d points: dx: (%f, %f, %f), dy: (%f, %f, %f)", count, 
            dxMin, dxMean, dxMax, dyMin, dyMean, dyMax );
    }

    // display visualization
    dispImage( *frame );
}


/// build and solve sparse linear system to update particle positions
void solveSystem( SimpleParticleSet &particleSet, int frameIndex, MatrixI &links, 
                  VectorI &ptMap, VectorI &ptMapInv, Config &spConf, 
                  VectorF &ixix, VectorF &ixiy, VectorF &iyiy, VectorF &ixiz, VectorF &iyiz, 
                  VectorF &dxVect, VectorF &dyVect ) {
    int pointCount = ptMap.length();
    int linkCount = links.rows();
    int X = 0;
    int Y = 1;

    // get config parameters
    int solverIters = spConf.readInt( "solverIters" );
    float dataFactor = spConf.readFloat( "dataFactor" );
    float smoothFactor = spConf.readFloat( "smoothFactor" );
    float sorFactor = spConf.readFloat( "sorFactor" );
    float condFactor = spConf.readFloat( "condFactor" );

    // initialize sparse linear system
    SparseSystem sparseSystem( pointCount * 2 );
    sparseSystem.setMaxIter( solverIters );
    sparseSystem.setSORFactor( sorFactor );
    sparseSystem.setDisplayIndent( 3 );
    VectorF init( pointCount * 2 );

    // add data term for each point
    for (int pointIndex = 0; pointIndex < pointCount; pointIndex++) {
        int varIndex = pointIndex * 2;

        // store initial value
        init[ varIndex + _X_ ] = dxVect[ pointIndex ];
        init[ varIndex + _Y_ ] = dyVect[ pointIndex ];

        // set data terms
        sparseSystem.setB( varIndex + _X_, 0 );
        sparseSystem.setB( varIndex + _Y_, 0 );
        addDataEquation( sparseSystem, varIndex, ixix[ pointIndex ], ixiy[ pointIndex ], 
                         iyiy[ pointIndex ], ixiz[ pointIndex ], iyiz[ pointIndex ], 
                         dataFactor, condFactor );
    }

    // add smoothness terms
    if (smoothFactor) {
        for (int linkIndex = 0; linkIndex < linkCount; linkIndex++) {
            int partIndex = links.data( linkIndex, 0 );
            int otherPartIndex = links.data( linkIndex, 1 );
            if (partIndex < 0 || partIndex >= ptMapInv.length() || otherPartIndex < 0 || otherPartIndex >= ptMapInv.length() )
                fatalError( "SolveSystem: invalid particle index for smoothness" );
            int varIndex = ptMapInv[ partIndex ] * 2;
            int otherVarIndex = ptMapInv[ otherPartIndex ] * 2;
            if (varIndex < 0 || otherVarIndex < 0)
                fatalError( "SolveSystem: invalid var index for smoothness" );

            // add equation to minimize (dx_i - dx_j)
            addPairEquation( sparseSystem, varIndex, varIndex, otherVarIndex, smoothFactor );
        }
    }

    // set initial value for iterative solver
    sparseSystem.setInit( init );

    // solve the system
    particleSet.timeRef( TIMER_SOLVER ).start();
    VectorF result = sparseSystem.solve();
    particleSet.timeRef( TIMER_SOLVER ).stop();

    // copy solution from result vector into dx and dy (returned results)
    for (int pointIndex = 0; pointIndex < pointCount; pointIndex++) {
        dxVect[ pointIndex ] = result[ pointIndex * 2 + X ];
        dyVect[ pointIndex ] = result[ pointIndex * 2 + Y ];
    }
}


/// optimizes particle positions using blurred channel values
void simpleParticleOptimize( SimpleParticleSet &particleSet, int frameIndex, 
                             const Array<ImageGrayF> &chan, Config &spConf ) {
    int width = chan[ 0 ].width(), height = chan[ 0 ].height();
    aptr<MatrixI> links = createLinks( particleSet, frameIndex, true );
    plotLinks( particleSet, *links, frameIndex, width, height );

    // start timer
    particleSet.timeRef( TIMER_OPT ).start();

    // get config parameters
    int outerIters = spConf.readInt( "outerIters" );
    float maxMove = spConf.readFloat( "maxMove" );

    // get variable maps
    aptr<VectorI> ptMap = pointMap( particleSet, frameIndex );
    aptr<VectorI> ptMapInv = pointMapInverse( particleSet, frameIndex );
    checkPointMap( particleSet, frameIndex, *ptMap, *ptMapInv );
    int pointCount = ptMap->length();

    // compute gradients of image channels
    Array<ImageGrayF> chanDx = dx( chan );
    Array<ImageGrayF> chanDy = dy( chan );

    // perform outer loop
    int outerIter = 0;
    for (; outerIter < outerIters; outerIter++) {

        // computation for data term
        int channelCount = chan.count();
        VectorF ixix( pointCount ), ixiy( pointCount ), iyiy( pointCount ), ixiz( pointCount ), iyiz( pointCount );
        ixix.clear( 0 );
        ixiy.clear( 0 );
        iyiy.clear( 0 );
        ixiz.clear( 0 );
        iyiz.clear( 0 );
        for (int i = 0; i < pointCount; i++) {
            int particleIndex = ptMap->data( i );
            if (particleIndex < 0 || particleIndex >= particleSet.count())
                fatalError( "particle index out of bounds" );
            const SimpleParticle &p = particleSet.ref( particleIndex );
            if (p.active( frameIndex - 1, frameIndex ) == false)
                fatalError( "particle not active" );
            float x = p.x( frameIndex );
            float y = p.y( frameIndex );
            float xPrev = p.x( frameIndex - 1 );
            float yPrev = p.y( frameIndex - 1 );

            // loop over channels; set data terms
            int boundary = 3;
            float ixixSum = 0, ixiySum = 0, iyiySum = 0, ixizSum = 0, iyizSum = 0;
            if (x - boundary > 0 && x + boundary < width && y - boundary > 0 && y + boundary < height
                    && xPrev - boundary > 0 && xPrev + boundary < width && yPrev - boundary > 0 && yPrev + boundary < height) {
                for (int j = 0; j < channelCount; j++) {
                    float ix = chanDx[ j ].interp( x, y );
                    float iy = chanDy[ j ].interp( x, y );
                    float iz = chan[ j ].interp( x, y ) 
                             - p.channelBlur( frameIndex, j );
                    float dPsi = 1; //_PsiDeriv( iz * iz ); // fix(later): what is correct term here?
                    ixixSum += ix * ix * dPsi;
                    ixiySum += ix * iy * dPsi;
                    iyiySum += iy * iy * dPsi;
                    ixizSum += ix * iz * dPsi;
                    iyizSum += iy * iz * dPsi;
                }
            }
            ixix[ i ] = ixixSum;
            ixiy[ i ] = ixiySum;
            iyiy[ i ] = iyiySum;
            ixiz[ i ] = ixizSum;
            iyiz[ i ] = iyizSum;
        }
        
        // start with zero point update
        VectorF dx( pointCount ), dy( pointCount );
        dx.clear( 0 );
        dy.clear( 0 );

        // solve for new dx, dy
        solveSystem( particleSet, frameIndex, *links, *ptMap, *ptMapInv, spConf, 
                     ixix, ixiy, iyiy, ixiz, iyiz, dx, dy );

        // visualize offsets from inner loop
        visualizeOffsets( particleSet, *ptMap, frameIndex, chan, dx, dy );

        // apply point update
        updateParticlePositions( particleSet, frameIndex, *ptMap, dx, dy, maxMove );

        // check for user cancel
        if (checkCommandEvents())
            break;
    }

    // display status
    disp( 3, "opt complete after %d outer iters", outerIter );
    
    // stop timer
    particleSet.timeRef( TIMER_OPT ).stop();
}


} // end namespace pvl

