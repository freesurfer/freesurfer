#include <pvl/VarMotionMultiRes.h>
#include <sbl/core/Command.h> // for checkCommandEvents()
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/MotionFieldUtil.h>
#include <pvl/SparseSystem.h>
#include <pvl/VarMotionUtil.h>
#include <pvl/Occlusion.h>
namespace pvl {


// This file provides an algorithm that uses a lower resolution for the motion field than for the image data, allowing faster optimization.
// This file duplicates much of the basic VarMotion algorithm, but is kept separate in order to avoid complicating the original implementation.


// dealing with edge cases:
// - the mask specifies which pixels are variables in the sparse system (to have flow estimated)
// - a flow vector that projects outside the image has no data penalty
// - the lower (y=0) and left (x=0) pixels have no data term and only a basic smoothness term


//-------------------------------------------
// VARIATIONAL OPTIMIZATION
//-------------------------------------------


/// solves for du, dv (for all pixels in the mask)
void varMultiResSolveSystem( const ImageGrayU &mask, int pixelCount, Config &varConf, const ImageGrayF &localSmoothness, 
                             const ImageGrayI &uIndex, const ImageGrayI &vIndex,
                             const ImageGrayF &uDx, const ImageGrayF &uDy, 
                             const ImageGrayF &vDx, const ImageGrayF &vDy, 
                             const Array<ImageGrayF> &Ix, const Array<ImageGrayF> &Iy, const Array<ImageGrayF> &Iz,
                             const Array<ImageGrayF> &dataFactorImg, const ImageGrayF &dPsiSmoothness,
                             ImageGrayF &du, ImageGrayF &dv ) {

    // get parameters
    float smoothness = varConf.readFloat( "smoothness" ); // aka alpha (in Brox papers)
    int solverItersPerMegapixel = varConf.readInt( "solverItersPerMegapixel" );
    float sorFactor = varConf.readFloat( "sorFactor" );
    float maxMove = varConf.readFloat( "maxMove" );
    int adaptIterMask = varConf.readInt( "adaptIterMask" );
    float adaptUpdateThresh = varConf.readFloat( "adaptUpdateThresh" );
    
    // add more iters if small 
    int solverIters = sbl::round( (float) solverItersPerMegapixel * (float) pixelCount / 1000000.0f);

    // allocate sparse linear system
    SparseSystem sparseSystem( pixelCount * 2 );
    VectorF init( pixelCount * 2 );

    // loop over pixels, adding equation for each variable (two variables per pixel)
    int width = mask.width(), height = mask.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y )) {

                // get indices of variables used in this pair of equation
                int uInd = uIndex.data( x, y );
                int vInd = vIndex.data( x, y );
                int uMxInd = x > 0 ? uIndex.data( x - 1, y ) : -1;
                int uPxInd = x < width - 1 ? uIndex.data( x + 1, y ) : -1;
                int uMyInd = y > 0 ? uIndex.data( x, y - 1 ) : -1;
                int uPyInd = y < height - 1 ? uIndex.data( x, y + 1 ) : -1;
                int vMxInd = x > 0 ? vIndex.data( x - 1, y ) : -1;
                int vPxInd = x < width - 1 ? vIndex.data( x + 1, y ) : -1;
                int vMyInd = y > 0 ? vIndex.data( x, y - 1 ) : -1;
                int vPyInd = y < height - 1 ? vIndex.data( x, y + 1 ) : -1;

                // check indices
                if (uInd < 0 || uInd >= init.length() || vInd < 0 || vInd >= init.length()) 
                    fatalError( "SolveInnerVarSystem: sanity check failed" );
    
                // store current value of du, dv as initial values for iterative solver
                init[ uInd ] = du.data( x, y );
                init[ vInd ] = dv.data( x, y );

                // these coefs are for the first equation (derivative w.r.t. du)
                float duCoef1 = 0, dvCoef1 = 0, b1 = 0;

                // these coefs are for the second equation (derivative w.r.t. dv)
                float duCoef2 = 0, dvCoef2 = 0, b2 = 0;

                // loop over channels
                int xLarger = x * 2, yLarger = y * 2;
                for (int yo = 0; yo < 2; yo++) {
                    for (int xo = 0; xo < 2; xo++) {
                        for (int k = 0; k < Ix.count(); k++) {

                            // prep for data coefs
                            float ix = Ix[ k ].data( xLarger + xo, yLarger + yo );
                            float iy = Iy[ k ].data( xLarger + xo, yLarger + yo );
                            float iz = Iz[ k ].data( xLarger + xo, yLarger + yo );
                            float dataFactor = dataFactorImg[ k ].data( xLarger + xo, yLarger + yo );

                            // compute data coefs
                            duCoef1 += dataFactor * ix * ix;
                            dvCoef1 += dataFactor * ix * iy;
                            b1 += -dataFactor * ix * iz; // negate because subtract from both sides
                            duCoef2 += dataFactor * iy * ix;
                            dvCoef2 += dataFactor * iy * iy;
                            b2 += -dataFactor * iy * iz; // negate because subtract from both sides
                        }
                    }
                }

                // prep for smoothness term
                float alpha = smoothness + localSmoothness.data( x, y );
                float smoothFactor = alpha * dPsiSmoothness.data( x, y );
                float smoothFactorPx = x < width - 1 ? alpha * dPsiSmoothness.data( x + 1, y ) : 0;
                float smoothFactorPy = y < height - 1 ? alpha * dPsiSmoothness.data( x, y + 1 ) : 0;

                // smoothness: u eqn coefs
                duCoef1 += 2.0f * smoothFactor;
                duCoef1 += smoothFactorPx;
                duCoef1 += smoothFactorPy;
                float duMxCoef = -smoothFactor;
                float duMyCoef = -smoothFactor;
                float duPxCoef = -smoothFactorPx;
                float duPyCoef = -smoothFactorPy;
                b1 -= (uDx.data( x, y ) + uDy.data( x, y )) * smoothFactor;
                if (x < width - 1) b1 += uDx.data( x + 1, y ) * smoothFactorPx;
                if (y < height - 1) b1 += uDy.data( x, y + 1 ) * smoothFactorPy;

                // smoothness: v eqn coefs
                dvCoef2 += 2.0f * smoothFactor;
                dvCoef2 += smoothFactorPx;
                dvCoef2 += smoothFactorPy;
                float dvMxCoef = -smoothFactor;
                float dvMyCoef = -smoothFactor;
                float dvPxCoef = -smoothFactorPx;
                float dvPyCoef = -smoothFactorPy;
                b2 -= (vDx.data( x, y ) + vDy.data( x, y )) * smoothFactor;
                if (x < width - 1) b2 += vDx.data( x + 1, y ) * smoothFactorPx;
                if (y < height - 1) b2 += vDy.data( x, y + 1 ) * smoothFactorPy;
                                       
                // smoothness: handle edge cases 
                if (uMxInd == -1) {
                    uMxInd = 0;
                    duMxCoef = 0;
                    vMxInd = 0;
                    dvMxCoef = 0;
                }
                if (uPxInd == -1) {
                    uPxInd = 0;
                    duPxCoef = 0;
                    vPxInd = 0;
                    dvPxCoef = 0;
                }
                if (uMyInd == -1) {
                    uMyInd = 0;
                    duMyCoef = 0;
                    vMyInd = 0;
                    dvMyCoef = 0;
                }
                if (uPyInd == -1) {
                    uPyInd = 0;
                    duPyCoef = 0;
                    vPyInd = 0;
                    dvPyCoef = 0;
                }

                // add equations
                if (x == 0) {
                    sparseSystem.addEquation( uInd, vInd, uMxInd, uPxInd, uMyInd, uPyInd, 
                                              1, 0, 0, -1, 0, 0, 0 );
                    sparseSystem.addEquation( vInd, uInd, vMxInd, vPxInd, vMyInd, vPyInd, 
                                              1, 0, 0, -1, 0, 0, 0 );
                } else if (y == 0) {
                    sparseSystem.addEquation( uInd, vInd, uMxInd, uPxInd, uMyInd, uPyInd, 
                                              1, 0, 0, 0, 0, -1, 0 );
                    sparseSystem.addEquation( vInd, uInd, vMxInd, vPxInd, vMyInd, vPyInd, 
                                              1, 0, 0, 0, 0, -1, 0 );
                } else {
                    sparseSystem.addEquation( uInd, vInd, uMxInd, uPxInd, uMyInd, uPyInd, 
                                              duCoef1, dvCoef1, duMxCoef, duPxCoef, duMyCoef, duPyCoef, b1 );
                    sparseSystem.addEquation( vInd, uInd, vMxInd, vPxInd, vMyInd, vPyInd, 
                                              dvCoef2, duCoef2, dvMxCoef, dvPxCoef, dvMyCoef, dvPyCoef, b2 );
                }
            }
        }
    }

    // set initial value for iterative solver
    sparseSystem.setInit( init );

    // solve the system
    sparseSystem.setMaxIter( solverIters );
    sparseSystem.setSORFactor( sorFactor );
    sparseSystem.enableAdaptive( adaptIterMask, adaptUpdateThresh );
    VectorF result = sparseSystem.solve();

    // copy solution from result vector into du and dv (returned results)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y )) {
                float duBounded = bound( result.data( uIndex.data( x, y ) ), -maxMove, maxMove );
                float dvBounded = bound( result.data( vIndex.data( x, y ) ), -maxMove, maxMove );
                du.data( x, y ) = duBounded;
                dv.data( x, y ) = dvBounded;
            }
        }
    }
}


/// incrementally computes du, dv
void varMultiResInnerLoop( const ImageGrayU &maskMotion, int pixelCount, Config &varConf, const ImageGrayF &localSmoothness,  
                           const ImageGrayI &uIndex, const ImageGrayI &vIndex,
                           const ImageGrayF &uDx, const ImageGrayF &uDy, 
                           const ImageGrayF &vDx, const ImageGrayF &vDy, 
                           const Array<ImageGrayF> &Ix, const Array<ImageGrayF> &Iy, const Array<ImageGrayF> &Iz,
                           const ImageGrayF &u, const ImageGrayF &v, 
                           ImageGrayF &du, ImageGrayF &dv ) {

    // get parameters
    int innerIters = varConf.readInt( "innerIters" );

    // allocate workspace (don't need to init)
    int mfWidth = maskMotion.width(), mfHeight = maskMotion.height();    
    int chanWidth = Ix[ 0 ].width(), chanHeight = Ix[ 0 ].height();
    int chanCount = Ix.count();
    Array<ImageGrayF> dPsiData = allocChannels( chanWidth, chanHeight, chanCount, 0 );
    ImageGrayF dPsiSmoothness( mfWidth, mfHeight );
    dPsiSmoothness.clear( 0 ); // necessary?

    // run inner loop
    for (int innerIter = 0; innerIter < innerIters; innerIter++) {

        // check for cancel
        if (checkCommandEvents())
            break;
    
        // compute dPsiData using current values of du, dv
        for (int y = 0; y < chanHeight; y++) {
            for (int x = 0; x < chanWidth; x++) {
                int xSmaller = sbl::round( x * 0.5f );
                int ySmaller = sbl::round( y * 0.5f );
                if (xSmaller >= mfWidth) xSmaller = mfWidth - 1;
                if (ySmaller >= mfHeight) ySmaller = mfHeight - 1;

                // compute dPsi terms for each pixel
                if (maskMotion.data( xSmaller, ySmaller )) {

                    // data term: compute dPsi for each channel
                    for (int k = 0; k < chanCount; k++) {
                        float diff = Ix[ k ].data( x, y ) * du.data( xSmaller, ySmaller ) 
                                   + Iy[ k ].data( x, y ) * dv.data( xSmaller, ySmaller ) 
                                   + Iz[ k ].data( x, y );
                        dPsiData[ k ].data( x, y ) = broxDistDeriv( diff * diff );
                    }
                }
            }
        }

        // compute dPsiSmoothness using current values of du, dv (skip x = 0 and y = 0 pixels)
        for (int y = 1; y < mfHeight; y++) {
            for (int x = 1; x < mfWidth; x++) {

                // compute dPsi terms for each pixel
                if (maskMotion.data( x, y )) {

                    // smoothness term (using gradient of (u + du, v + dv))
                    float uNewDx = uDx.data( x, y ) + dx( du, x, y );
                    float uNewDy = uDy.data( x, y ) + dy( du, x, y );
                    float vNewDx = vDx.data( x, y ) + dx( dv, x, y );
                    float vNewDy = vDy.data( x, y ) + dy( dv, x, y );
                    dPsiSmoothness.data( x, y ) = broxDistDeriv( uNewDx * uNewDx + uNewDy * uNewDy + vNewDx * vNewDx + vNewDy * vNewDy );
                }
            }
        }

        // solve the system
        varMultiResSolveSystem( maskMotion, pixelCount, varConf, localSmoothness, uIndex, vIndex, 
                                uDx, uDy, vDx, vDy, Ix, Iy, Iz, 
                                dPsiData, dPsiSmoothness, du, dv );
    }
}


/// computes an update (du, dv) to the flow field using a different resolution for the image vs. the motion field
void varMultiResIteration( MotionField &mf, const Array<ImageGrayF> &srcChanScaled, const Array<ImageGrayF> &destChanScaled, 
                           const ImageGrayU &maskMotion, Config &varConf ) {    

    // get config parameters
    float flatSigma = varConf.readFloat( "flatSigma" );
    float flatSmoothness = varConf.readFloat( "flatSmoothness" );

    // info for this iteration
    int chanCount = srcChanScaled.count(), chanWidth = srcChanScaled[ 0 ].width(), chanHeight = srcChanScaled[ 0 ].height();
    int mfWidth = mf.width(), mfHeight = mf.height();

    // compute index maps (indices of solver variables)
    aptr<ImageGrayI> uIndex, vIndex;
    int pixelCount = buildIndexMaps( maskMotion, uIndex, vIndex );

    // get derivatives of destination channels
    Array<ImageGrayF> destChanDx = dx( destChanScaled );
    Array<ImageGrayF> destChanDy = dy( destChanScaled );

    // get flow components for quick reference
    ImageGrayF &u = mf.uRef();
    ImageGrayF &v = mf.vRef();

    // allocate workspace 
    Array<ImageGrayF> Ix = allocChannels( chanWidth, chanHeight, chanCount, 0 );
    Array<ImageGrayF> Iy = allocChannels( chanWidth, chanHeight, chanCount, 0 );
    Array<ImageGrayF> Iz = allocChannels( chanWidth, chanHeight, chanCount, 0 );

    // use u and v to compute values needed by inner loop
    for (int y = 0; y < chanHeight; y++) {
        for (int x = 0; x < chanWidth; x++) {
            int xSmaller = sbl::round( x * 0.5f );
            int ySmaller = sbl::round( y * 0.5f );
            if (xSmaller >= mfWidth) xSmaller = mfWidth - 1;
            if (ySmaller >= mfHeight) ySmaller = mfHeight - 1;

            // if inside optimization region
            if (maskMotion.data( xSmaller, ySmaller ) > 150 && mf.occSafe( xSmaller, ySmaller ) < 128) {
                float xProj = x + u.data( xSmaller, ySmaller ) * 2.0f;
                float yProj = y + v.data( xSmaller, ySmaller ) * 2.0f;

                // if projects inside image (no penalty if project outside image)
                if (xProj > 0 && xProj < chanWidth - 1 && yProj > 0 && yProj < chanHeight - 1) {
                    for (int k = 0; k < chanCount; k++) {
                        Ix[ k ].data( x, y ) = destChanDx[ k ].interp( xProj, yProj );
                        Iy[ k ].data( x, y ) = destChanDy[ k ].interp( xProj, yProj );
                        Iz[ k ].data( x, y ) = destChanScaled[ k ].interp( xProj, yProj ) - srcChanScaled[ k ].data( x, y );
                    }
                }
            }
        }
    }

    // compute flow derivatives
    aptr<ImageGrayF> uDx = dx( u );
    aptr<ImageGrayF> uDy = dy( u );
    aptr<ImageGrayF> vDx = dx( v );
    aptr<ImageGrayF> vDy = dy( v );

    // compute local smoothness factor
    aptr<ImageGrayF> lSmoothness = localSmoothness( srcChanScaled[ 0 ], flatSigma, flatSmoothness );

    // initialize du, dv = 0
    ImageGrayF du( mfWidth, mfHeight ), dv( mfWidth, mfHeight );
    du.clear( 0 );
    dv.clear( 0 );

    // run inner loop to compute du, dv
    varMultiResInnerLoop( maskMotion, pixelCount, varConf, *lSmoothness, *uIndex, *vIndex, 
                          *uDx, *uDy, *vDx, *vDy, Ix, Iy, Iz, u, v, du, dv );

    // increment u and v according to du and dv
    for (int y = 0; y < mfHeight; y++) {
        for (int x = 0; x < mfWidth; x++) {
            if (maskMotion.data( x, y )) {
                u.data( x, y ) = u.data( x, y ) + du.data( x, y );
                v.data( x, y ) = v.data( x, y ) + dv.data( x, y );
            }
        }
    }
}


} // end namespace pvl

