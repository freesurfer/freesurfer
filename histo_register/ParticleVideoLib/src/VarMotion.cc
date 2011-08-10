#include <pvl/VarMotion.h>
#include <sbl/core/Command.h> // for checkCommandEvents
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/MotionFieldUtil.h>
#include <pvl/SparseSystem.h>
#include <pvl/VarMotionUtil.h>
#include <pvl/Occlusion.h>
#include <pvl/VarMotionMultiRes.h>
namespace pvl {


// dealing with edge cases:
// - the mask specifies which pixels are variables in the sparse system (to have flow estimated)
// - a flow vector that projects outside the image has no data penalty
// - the lower (y=0) and left (x=0) pixels have no data term and only a basic smoothness term


//-------------------------------------------
// DIAGNOSTIC UTILS
//-------------------------------------------


/// display statistics about the given image
void dispStats( const String &name, const ImageGrayF &img ) {
    int width = img.width(), height = img.height();
    float min = 1e8, max = -1e8, mean = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            float val = img.data( x, y );
            if (val < min)
                min = val;
            if (val > max)
                max = val;
            mean += val;
        }
    mean /= (float) (width * height);
    disp( 5, "[%s]: min: %f, mean: %f, max: %f", name.c_str(), min, mean, max );
}


//-------------------------------------------
// VARIATIONAL OPTIMIZATION
//-------------------------------------------


/// solves for du, dv (for all pixels in the mask)
void varSolveSystem( const ImageGrayU &mask, int pixelCount, Config &varConf, const ImageGrayF &localSmoothness, 
                     const ImageGrayI &uIndex, const ImageGrayI &vIndex,
                     const ImageGrayF &uDx, const ImageGrayF &uDy, 
                     const ImageGrayF &vDx, const ImageGrayF &vDy, 
                     const Array<ImageGrayF> &Ix, const Array<ImageGrayF> &Iy, const Array<ImageGrayF> &Iz,
                     const Array<ImageGrayF> &dataFactorImg, const ImageGrayF &dPsiSmoothness,
                     ImageGrayF &du, ImageGrayF &dv, bool fast ) {

    // get parameters
    float smoothness = varConf.readFloat( "smoothness" ); // aka alpha (in Brox papers)
    int solverItersPerMegapixel = varConf.readInt( "solverItersPerMegapixel" );
    float sorFactor = varConf.readFloat( "sorFactor" );
    float maxMove = varConf.readFloat( "maxMove" );
    int adaptIterMask = varConf.readInt( "adaptIterMask" );
    float adaptUpdateThresh = varConf.readFloat( "adaptUpdateThresh" );

    // add more iters if small
    int solverIters = sbl::round( (float) solverItersPerMegapixel * (float) pixelCount / 1000000.0f);
    if (fast) {
        adaptUpdateThresh *= 10.0f;
        solverIters /= 2;
    }

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
                for (int k = 0; k < Ix.count(); k++) {

                    // prep for data coefs
                    float ix = Ix[ k ].data( x, y );
                    float iy = Iy[ k ].data( x, y );
                    float iz = Iz[ k ].data( x, y );
                    float dataFactor = dataFactorImg[ k ].data( x, y );

                    // compute data coefs
                    duCoef1 += dataFactor * ix * ix;
                    dvCoef1 += dataFactor * ix * iy;
                    b1 += -dataFactor * ix * iz; // negate because subtract from both sides
                    duCoef2 += dataFactor * iy * ix;
                    dvCoef2 += dataFactor * iy * iy;
                    b2 += -dataFactor * iy * iz; // negate because subtract from both sides
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
void varInnerLoop( const ImageGrayU &mask, int pixelCount, Config &varConf, const ImageGrayF &localSmoothness,  
                   const ImageGrayI &uIndex, const ImageGrayI &vIndex,
                   const ImageGrayF &uDx, const ImageGrayF &uDy, 
                   const ImageGrayF &vDx, const ImageGrayF &vDy, 
                   const Array<ImageGrayF> &Ix, const Array<ImageGrayF> &Iy, const Array<ImageGrayF> &Iz,
                   const ImageGrayF &u, const ImageGrayF &v, 
                   ImageGrayF &du, ImageGrayF &dv, bool fast ) {

    // get parameters
    int innerIters = varConf.readInt( "innerIters" );
    int verbosity = varConf.readInt( "verbosity" );

    // allocate workspace (don't need to init)
    int width = mask.width(), height = mask.height();    
    int chanCount = Ix.count();
    Array<ImageGrayF> dPsiData = allocChannels( width, height, chanCount, 0 );
    ImageGrayF dPsiSmoothness( width, height );
    dPsiSmoothness.clear( 0 ); // necessary?

    // run inner loop
    for (int innerIter = 0; innerIter < innerIters; innerIter++) {

        // check for cancel
        if (checkCommandEvents())
            break;
    
        // compute dPsiData, dPsiSmoothness using current values of du, dv (skip x = 0 and y = 0 pixels)
        for (int y = 1; y < height; y++) {
            for (int x = 1; x < width; x++) {

                // compute dPsi terms for each pixel
                if (mask.data( x, y )) {

                    // data term: compute dPsi for each channel
                    for (int k = 0; k < chanCount; k++) {
                        float diff = Ix[ k ].data( x, y ) * du.data( x, y ) 
                                   + Iy[ k ].data( x, y ) * dv.data( x, y ) 
                                   + Iz[ k ].data( x, y );
                        dPsiData[ k ].data( x, y ) = broxDistDeriv( diff * diff );
                    }

                    // smoothness term (using gradient of (u + du, v + dv))
                    float uNewDx = uDx.data( x, y ) + dx( du, x, y );
                    float uNewDy = uDy.data( x, y ) + dy( du, x, y );
                    float vNewDx = vDx.data( x, y ) + dx( dv, x, y );
                    float vNewDy = vDy.data( x, y ) + dy( dv, x, y );
                    dPsiSmoothness.data( x, y ) = broxDistDeriv( uNewDx * uNewDx + uNewDy * uNewDy + vNewDx * vNewDx + vNewDy * vNewDy );
                }
            }
        }

        // diagnostics
        if (verbosity > 12) {
            for (int k = 0; k < chanCount; k++)
                dispStats( "dPsiData", dPsiData[ k ] );
            dispStats( "dPsiSmoothness", dPsiSmoothness );
        }

        // solve the system
        varSolveSystem( mask, pixelCount, varConf, localSmoothness, uIndex, vIndex, 
                        uDx, uDy, vDx, vDy, Ix, Iy, Iz, 
                        dPsiData, dPsiSmoothness, du, dv, fast );
    }
}


/// computes an update (du, dv) to the flow field
void varIteration( MotionField &mf, const Array<ImageGrayF> &srcChanScaled, const Array<ImageGrayF> &destChanScaled, 
                   const ImageGrayU &maskMotion, Config &varConf, bool fast ) {    

    // get config parameters
    float flatSigma = varConf.readFloat( "flatSigma" );
    float flatSmoothness = varConf.readFloat( "flatSmoothness" );
    int verbosity = varConf.readInt( "verbosity" );

    // info for this iteration    
    int chanCount = srcChanScaled.count(), width = mf.width(), height = mf.height();
    
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
    Array<ImageGrayF> Ix = allocChannels( width, height, chanCount, 0 );
    Array<ImageGrayF> Iy = allocChannels( width, height, chanCount, 0 );
    Array<ImageGrayF> Iz = allocChannels( width, height, chanCount, 0 );

    // compute local smoothness factor
    aptr<ImageGrayF> lSmoothness = localSmoothness( srcChanScaled[ 0 ], flatSigma, flatSmoothness );

    // use u and v to compute values needed by inner loop
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {

            // if inside optimization region (and unoccluded)
            if (maskMotion.data( x, y ) && mf.occSafe( x, y ) < 128) {
                float xProj = x + u.data( x, y );
                float yProj = y + v.data( x, y );

                // if projects inside image (no penalty if project outside image)
                if (xProj > 0 && xProj < width - 1 && yProj > 0 && yProj < height - 1) {
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

    // initialize du, dv = 0
    ImageGrayF du( width, height ), dv( width, height );
    du.clear( 0 );
    dv.clear( 0 );

    // run inner loop to compute du, dv
    varInnerLoop( maskMotion, pixelCount, varConf, *lSmoothness, *uIndex, *vIndex, 
                  *uDx, *uDy, *vDx, *vDy, Ix, Iy, Iz, u, v, du, dv, fast );

    // increment u and v according to du and dv
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (maskMotion.data( x, y )) {
                u.data( x, y ) = u.data( x, y ) + du.data( x, y );
                v.data( x, y ) = v.data( x, y ) + dv.data( x, y );
            }
        }
    }

    // display diagnostics
    if (verbosity >= 8) {
        dispStats( "Ix", Ix[ 0 ] );
        dispStats( "Iy", Iy[ 0 ] );
        dispStats( "Iz", Iz[ 0 ] );
        dispStats( "uDx", *uDx );
        dispStats( "uDy", *uDy );
        dispStats( "vDx", *vDx );
        dispStats( "vDy", *vDy );
        dispStats( "du", du );
        dispStats( "dv", dv );
        dispMotionStats( 3, mf );
    }
}


/// loops over scales; for each scale; executes outer loop
aptr<MotionField> varScaleLoop( const Array<ImageGrayF> &srcChan, const Array<ImageGrayF> &destChan, const MotionField *init, const ImageGrayU &mask, Config &varConf ) {

    // get config parameters
    float sigma = varConf.readFloat( "sigma" );
    float scaleFactor = varConf.readFloat( "scaleFactor" );
    float minScale = varConf.readFloat( "minScale" );
    int verbosity = varConf.readInt( "verbosity" );
    bool visualizeIntermediateFlow = varConf.readBool( "visualizeIntermediateFlow" );

    // configure alternate versions
    bool enableNormal = true;
    int multiScaleFactor = 1;
    bool enableMultiScale = varConf.readBool( "enableMultiScale" );
    bool finalPass = false;
    if (enableMultiScale) {
        multiScaleFactor = 2;
        enableNormal = false;
        finalPass = varConf.readBool( "enableFinalPass" );
    }

    // display channels
    if (verbosity >= 8)
        dispChannels( srcChan );

    // compute scales 
    VectorF scales = computeScales( scaleFactor, minScale * (float) multiScaleFactor );
    if (finalPass) {

        // prepend -1
        scales.append( 0 );
        for (int i = scales.length() - 1; i > 0; i--) 
            scales[ i ] = scales[ i - 1 ];
        scales[ 0 ] = -1;
    }
    int scaleCount = scales.length();
    if (verbosity >= 4)
        disp( 1, "iterating over %d scales (scale-factor = %f)", scaleCount, scaleFactor );

    // result motion field
    aptr<MotionField> mf;

    // loop over scales (from smallest to largest)
    for (int i = scaleCount - 1; i >= 0; i--) {
        float scale = scales[ i ];

        // if final pass, switch out of multi-res mode
        bool fast = false;
        if (scale < 0) {
            enableMultiScale = false;
            enableNormal = true;
            multiScaleFactor = 1;
            scale = 1;
            fast = true;
        }

        // compute motion and image dimensions
        int mfWidth = sbl::round( srcChan.ref( 0 ).width() * scale / (float) multiScaleFactor );
        int mfHeight = sbl::round( srcChan.ref( 0 ).height() * scale / (float) multiScaleFactor );
        int chanWidth = mfWidth * multiScaleFactor;
        int chanHeight = mfHeight * multiScaleFactor;
        if (verbosity >= 4)
            disp( 2, "scale %f: motion: %d by %d, image: %d by %d", scale, mfWidth, mfHeight, chanWidth, chanHeight );

        // rescale and pad the motion field
        if (mf.get() == NULL) {
            mf = createInitMotion( mfWidth, mfHeight, init );
        } else {
            mf->resize( mfWidth, mfHeight, true ); // uses linear interp
        }

        // shrink and blur the channels to current size
        Array<ImageGrayF> srcChanScaled = shrinkChannels( srcChan, chanWidth, chanHeight, sigma );
        Array<ImageGrayF> destChanScaled = shrinkChannels( destChan, chanWidth, chanHeight, sigma );
        aptr<ImageGrayU> maskScaled = resize( mask, chanWidth, chanHeight, true );
        Array<ImageGrayF> srcChanMotion = srcChanScaled;
        Array<ImageGrayF> destChanMotion = destChanScaled; 
        if (chanWidth != mfWidth) {
            srcChanMotion = shrinkChannels( srcChanScaled, mfWidth, mfHeight, sigma ); // fix(faster): don't need to resize all channels, just first
            destChanMotion = shrinkChannels( destChanScaled, mfWidth, mfHeight, sigma ); // fix(faster): don't need to resize all channels, just first
            maskScaled = resize( *maskScaled, mfWidth, mfHeight, true );
        }

        // run outer iteration (inner loop)
        if (enableMultiScale) 
            varMultiResIteration( *mf, srcChanScaled, destChanScaled, *maskScaled, varConf );
        if (enableNormal)
            varIteration( *mf, srcChanScaled, destChanScaled, *maskScaled, varConf, fast );

        // run occlusion processing after each inner loop
        if (mf->width() > 200 && checkCommandEvents() == false) {

            // run bilateral filter loop with occlusion detection
            filterMotionLoop( *mf, srcChanMotion[ 0 ], destChanMotion[ 0 ], varConf, fast );
        }

        // display intermediate flow field
        if (visualizeIntermediateFlow) {
            String caption = sprintF( "size: %d, %d", mf->width(), mf->height() );
            dispMotion( *mf, caption );
        }
    }
    
    // clean up and return computed motion field
    return mf;
}


//-------------------------------------------
// TOP-LEVEL VARIATIONAL ALGORITHM
//-------------------------------------------


/// top-level variational (modified Brox et al.) motion estimation algorithm;
/// mask (if any) specifies where the algorithm should estimate flow vectors;
/// init (if any) provides an initial motion estimation;
/// scaleFactor (if not 1) is used to estimate flow at a lower resolution (for example, at half resolution if scaleFactor == 2)
aptr<MotionField> varMotion( const ImageColorU &src, const ImageColorU &dest, const String &configFileName, const MotionField *init, const ImageGrayU *mask, int scaleFactor ) {
    Config varConf;
    varConf.load( configFileName );
    if (varConf.entryCount() == 0) 
        warning( "failed to load varMotion config: %s", configFileName.c_str() );

    // create small images (will be same size of scaleFactor == 1)
    int newWidth = src.width() / scaleFactor;
    int newHeight = src.height() / scaleFactor;
    aptr<ImageColorU> srcSmall = resize( src, newWidth, newHeight, true );
    aptr<ImageColorU> destSmall = resize( dest, newWidth, newHeight, true );
    aptr<ImageGrayU> maskSmall;
    if (mask) {
        maskSmall = resize( *mask, newWidth, newHeight, true );
    } else {
        maskSmall.reset( new ImageGrayU( newWidth, newHeight ) );
        maskSmall->clear( 255 ); // if mask not specified, estimate flow everywhere
    }

    // extact channels from source and dest images
    Array<ImageGrayF> srcChan = extractChannels( *srcSmall, varConf );
    Array<ImageGrayF> destChan = extractChannels( *destSmall, varConf );

    // loop over scales
    aptr<MotionField> mf = varScaleLoop( srcChan, destChan, init, *maskSmall, varConf );

    // resize flow to full image size
    mf->resize( src.width(), src.height(), true );
    return mf;
}


/// top-level variational (modified Brox et al.) motion estimation algorithm;
/// mask (if any) specifies where the algorithm should estimate flow vectors;
/// init (if any) provides an initial motion estimation;
/// scaleFactor (if not 1) is used to estimate flow at a lower resolution (for example, at half resolution if scaleFactor == 2)
aptr<MotionField> varMotion( const ImageGrayU &src, const ImageGrayU &dest, const String &configFileName, const MotionField *init, const ImageGrayU *mask, int scaleFactor ) {

    // convert input images to color (note: colorFactor should be zero when running on grayscale images)
    aptr<ImageColorU> srcColor = toColor( src );
    aptr<ImageColorU> destColor = toColor( dest );

    // run motion estimation
    return varMotion( *srcColor, *destColor, configFileName, init, mask, scaleFactor );
}


} // end namespace pvl

