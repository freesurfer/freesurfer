#include <pvl/Occlusion.h>
#include <sbl/core/Command.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/MotionFieldUtil.h>
namespace pvl {


//-------------------------------------------
// FLOW FILTERING
//-------------------------------------------


/// filter the motion field using a bilateral filter
// good param vals: distSigma = 4, flowSigma = 0.5, imageSigma = 7
void dualFilterMotion( MotionField &mf, const ImageGrayF &src, const ImageGrayU &filtMask, 
                       float distSigma, float flowSigma, float imageSigma, float minWeight ) {
    int width = src.width(), height = src.height();
    int filtMaskThresh = 10;

    // items computed from parameters
    int filtRadius = sbl::round( distSigma * 2.5f );
    int filtRadiusSqd = filtRadius * filtRadius;
    float distFactor = -0.5f / (distSigma * distSigma);
    float flowFactor = -0.5f / (flowSigma * flowSigma);
    float imageFactor = -0.5f / (imageSigma * imageSigma);

    // pointers to motion field data for quick reference
    ImageGrayF &u = mf.uRef();
    ImageGrayF &v = mf.vRef();
    ImageGrayU &occMask = mf.occRef();
    if (mf.occDefined() == false) 
        fatalError( "FilterMotion: occ mask not defined" );

    // loop over image
    for (int yc = 0; yc < height; yc++)
        for (int xc = 0; xc < width; xc++) 
            if (filtMask.data( xc, yc ) > filtMaskThresh) {

                // get info about center pixel
                float uc = u.data( xc, yc );
                float vc = v.data( xc, yc );
                float ic = src.data( xc, yc );

                // contains weighted average of filter inputs
                float uSum1 = 0, vSum1 = 0, weightSum1 = 0;
                float uSum2 = 0, vSum2 = 0, weightSum2 = 0;

                // compute filter bounds (outside loop)
                int xMin = -filtRadius;
                int xMax = filtRadius;
                int yMin = -filtRadius;
                int yMax = filtRadius;
                if (xc + xMin < 0)
                    xMin = xc;
                if (xc + xMax > width - 1)
                    xMax = width - 1 - xc;
                if (yc + yMin < 0)
                    yMin = yc;
                if (yc + yMax > height - 1)
                    yMax = height - 1 - yc;

                // loop over neighborhood
                for (int yo = yMin; yo <= yMax; yo++) {
                    for (int xo = xMin; xo <= xMax; xo++) {
                        if (yo * yo + xo * xo < filtRadiusSqd) {
                            int x = xc + xo;
                            int y = yc + yo;

                            // weight spatially
                            float weight1 = expf( (float) (xo * xo + yo * yo) * distFactor );

                            // weight by occlusion
                            weight1 *= 1.0f - (float) occMask.data( x, y ) / 255.0f;
                            float weight2 = weight1;

                            // weight by flow value
                            float uVal = u.data( x, y );
                            float vVal = v.data( x, y );
                            float uDiff = uc - uVal;
                            float vDiff = vc - vVal;
                            weight1 *= expf( (uDiff * uDiff + vDiff * vDiff) * flowFactor );

                            // weight by pixel value
                            float iDiff = ic - src.data( x, y );
                            weight2 *= expf( iDiff * iDiff * imageFactor );
                            
                            // accumulate filter values
                            uSum1 += uVal * weight1;
                            vSum1 += vVal * weight1;
                            weightSum1 += weight1;
                            uSum2 += uVal * weight2;
                            vSum2 += vVal * weight2;
                            weightSum2 += weight2;
                        }
                    }
                }

                // store filtered value
                if (weightSum1 + weightSum2 > minWeight) {
                    u.data( xc, yc ) = (uSum1 + uSum2) / (weightSum1 + weightSum2);
                    v.data( xc, yc ) = (vSum1 + vSum2) / (weightSum1 + weightSum2);
                }
            }
}


//-------------------------------------------
// OCCLUSION DETECTION
//-------------------------------------------


/// labels possibly occluded pixels using flow divergence and src/dest projection difference
void gaussianOcclusionMask( MotionField &mf, const ImageGrayF &src, const ImageGrayF &dest,
                            const ImageGrayU &occUpdateRegion, const ImageGrayF &div, 
                            Config &varConf, bool visualize ) {

    // get parameters
    float divSigma = varConf.readFloat( "occDivSigma" );
    float imgSigma = varConf.readFloat( "occImageSigma" );
    float divFactor = -0.5f / (divSigma * divSigma);
    float imgFactor = -0.5f / (imgSigma * imgSigma);
    
    // loop over image
    int width = src.width(), height = src.height();
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            if (occUpdateRegion.data( x, y )) {
                float u = mf.u( x, y );
                float v = mf.v( x, y );
                float xProj = x + u;
                float yProj = y + v;

                // if projects into dest
                if (dest.inBounds( xProj, yProj )) {
                    float imgDiff = src.data( x, y ) - dest.interp( xProj, yProj );
                    float divVal = div.data( x, y );
                    float negDiv = divVal < 0 ? divVal : 0;
                    float occVal = expf( negDiv * negDiv * divFactor )
                                 * expf( imgDiff * imgDiff * imgFactor );
                    int occInt = sbl::round( 255 * (1.0f - occVal) );
                    if (occInt > 255) 
                        occInt = 255;
                    mf.setOcc( x, y, occInt );

                // if projects out of image, set occluded
                } else {
                    mf.setOcc( x, y, 255 );
                }
            }
        }
    if (visualize)
        dispImage( mf.occRef(), "occ mask" );
}


/// labels possibly occluded pixels using flow divergence and src/dest projection difference;
/// also computes region in which to run bilateral filter
void identifyOcclusions( MotionField &mf, const ImageGrayF &src, 
                         const ImageGrayF &dest, Config &varConf, ImageGrayU &filtMask, bool visualize ) {

    // get parameters
    float divThresh = -0.2f;
    float gradMagThresh = 0.5f; // was 0.25f
    float gradMagThreshSqd = gradMagThresh * gradMagThresh;

    // enable occlusion mask
    if (mf.occDefined() == false)
        mf.enableOcc();

    // compute flow gradients
    aptr<ImageGrayF> uDx = xGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> uDy = xGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> vDx = yGrad( mf.vRef(), 3 );
    aptr<ImageGrayF> vDy = yGrad( mf.vRef(), 3 );

    // compute flow divergence and gradient magnitude 
    int width = mf.width(), height = mf.height();
    ImageGrayF gradMagSqd( width, height ), divergence( width, height );
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            float uDxVal = uDx->data( x, y );
            float uDyVal = uDy->data( x, y );
            float vDxVal = vDx->data( x, y );
            float vDyVal = vDy->data( x, y );
            gradMagSqd.data( x, y ) = uDxVal * uDxVal + uDyVal * uDyVal + vDxVal * vDxVal + vDyVal * vDyVal;
            divergence.data( x, y ) = uDxVal + vDyVal;
        }

    // bilateral filter region: union of prev occ and current motion boundaries
    ImageGrayU filtRegion( mf.occRef() );

    // occlusion update region: union of prev occ and current occlusion boundaries
    ImageGrayU occUpdateRegion( mf.occRef() );

    // compute filt region, occ update region, occ bound mask
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            if (gradMagSqd.data( x, y ) > gradMagThreshSqd) {
                filtRegion.data( x, y ) = 255;
                if (divergence.data( x, y ) < divThresh) {
                    occUpdateRegion.data( x, y ) = 255;
                }
            }
        }

    // expand occlusion update region
    aptr<ImageGrayU> occUpdateRegionBlur = blurBox( occUpdateRegion, 7 );

    // compute occlusion mask
    gaussianOcclusionMask( mf, src, dest, *occUpdateRegionBlur, divergence, varConf, visualize );

    // compute filter mask
    aptr<ImageGrayU> filtRegionBlur = blurBox( filtRegion, 5 );
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            filtMask.data( x, y ) = filtRegionBlur->data( x, y ); // fix(clean): use aptr to avoid copy
}


//-------------------------------------------
// TOP-LEVEL FILTER/OCCLUSION CODE
//-------------------------------------------


/// run bilateral filter with occlusion detection (single step)
void filterMotionStep( MotionField &mf, MotionField *mfRev, const ImageGrayF &src, const ImageGrayF &dest, Config &varConf ) {

    // get parameters
    float distSigma = varConf.readFloat( "distSigma" );
    float flowSigma = varConf.readFloat( "flowSigma" );
    float imageSigma = varConf.readFloat( "imageSigma" );
    float minWeight = 1.5f;
    bool visualize = false;

    // get filter mask and occlusion mask
    ImageGrayU filtMask( mf.width(), mf.height() );
    identifyOcclusions( mf, src, dest, varConf, filtMask, visualize );

    // run filtering 
    dualFilterMotion( mf, src, filtMask, distSigma, flowSigma, imageSigma, minWeight );

    // display motion
    if (visualize) {
        dispImage( filtMask, "filter mask" );
        dispMotion( mf, "filtered" );
    }
}


/// run bilateral filter loop with occlusion detection (multiple steps)
void filterMotionLoop( MotionField &mf, const ImageGrayF &src, const ImageGrayF &dest, Config &varConf, bool fast ) {

    // get parameters
    int filterIters = varConf.readInt( "filterIters" );
    if (fast)
        filterIters /= 2;

    // repeatedly apply filter
    for (int i = 0; i < filterIters; i++) {

        // filter forward motion
        filterMotionStep( mf, NULL, src, dest, varConf );

        // check for cancel
        if (checkCommandEvents())
            break;
    }
}


} // end namespace pvl

