#include <pvl/VarMotionUtil.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/MotionFieldUtil.h>
namespace pvl {


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns vector of scales for multi-res optimization
VectorF computeScales( float scaleFactor, float minScale ) {
    VectorF scales( 1 );
    scales[ 0 ] = 1;
    do {
        scales.append( scales.endValue() * scaleFactor );
    } while (scales.endValue() * scaleFactor > minScale);
    return scales;
}


/// allocate and build mapping of flow components to system variables;
/// each pixel in mask corresponds to 2 system variables (du and dv);
/// returns count of pixels in mask
int buildIndexMaps( const ImageGrayU &mask, aptr<ImageGrayI> &uIndex, aptr<ImageGrayI> &vIndex ) {

    // allocate and initialize index maps
    int width = mask.width(), height = mask.height();
    uIndex.reset( new ImageGrayI( width, height ) );
    vIndex.reset( new ImageGrayI( width, height ) );
    uIndex->clear( -1 );
    vIndex->clear( -1 );

    // fill in variable indices
    int index = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y )) {
                uIndex->data( x, y ) = index++;
                vIndex->data( x, y ) = index++;
            }
        }
    }

    // return count of pixels in index maps
    return index / 2;
}


/// create a motion field of the desired size;
/// if init is specified, use it (resize as appropriate)
aptr<MotionField> createInitMotion( int width, int height, const MotionField *init ) {
    aptr<MotionField> mf;
    if (init) {
        mf.reset( new MotionField( *init ) );
        mf->resize( width, height, true ); // uses linear interp
        disp( 3, "resize init field to %d by %d", width, height );
    } else {
        mf.reset( new MotionField( width, height ) );
        mf->clear();
    }
    return mf;
}


/// compute local smoothness factor (higher where image gradient is low)
aptr<ImageGrayF> localSmoothness( const ImageGrayF &src, float flatSigma, float flatSmoothness ) {
    int width = src.width(), height = src.height();
    aptr<ImageGrayF> flatness( new ImageGrayF( width, height ) );
    if (flatSmoothness) {
        float flatGaussFactor = gaussFactor( flatSigma * flatSigma );
        aptr<ImageGrayF> srcDx = xGrad( src, 3 );
        aptr<ImageGrayF> srcDy = yGrad( src, 3 );
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float srcDxVal = srcDx->data( x, y );
                float srcDyVal = srcDy->data( x, y );
                flatness->data( x, y ) = flatSmoothness * gauss( srcDxVal * srcDxVal + srcDyVal * srcDyVal, flatGaussFactor );
            }
        }
    } else {
        flatness->clear( 0 );
    }
    return flatness;
}


//-------------------------------------------
// DISCRETE DERIVATIVES
//-------------------------------------------


/// discrete derivative in x direction
float dx( const ImageGrayF &img, int x, int y ) {
    return img.data( x, y ) - img.data( x - 1, y );
}


/// discrete derivative in y direction
float dy( const ImageGrayF &img, int x, int y ) {
    return img.data( x, y ) - img.data( x, y - 1 );
}


/// apply discrete derivative in x direction to image
aptr<ImageGrayF> dx( const ImageGrayF &img ) {
    int width = img.width(), height = img.height();
    aptr<ImageGrayF> dx( new ImageGrayF( width, height ) );
    dx->clear( 0 );
    for (int y = 0; y < height; y++) {
        for (int x = 1; x < width; x++) {
            dx->data( x, y ) = img.data( x, y ) - img.data( x - 1, y );
        }
    }
    return dx;
}


/// apply discrete derivative in y direction to image
aptr<ImageGrayF> dy( const ImageGrayF &img ) {
    int width = img.width(), height = img.height();
    aptr<ImageGrayF> dy( new ImageGrayF( width, height ) );
    dy->clear( 0 );
    for (int y = 1; y < height; y++) {
        for (int x = 0; x < width; x++) {
            dy->data( x, y ) = img.data( x, y ) - img.data( x, y - 1 );
        }
    }
    return dy;
}


/// apply x discrete derivative to each channel
Array<ImageGrayF> dx( const Array<ImageGrayF> &chan ) {
    Array<ImageGrayF> chanDx;
    for (int i = 0; i < chan.count(); i++) 
        chanDx.append( dx( chan[ i ] ).release() );
    return chanDx;
}


/// apply y discrete derivative to each channel
Array<ImageGrayF> dy( const Array<ImageGrayF> &chan ) {
    Array<ImageGrayF> chanDy;
    for (int i = 0; i < chan.count(); i++) 
        chanDy.append( dy( chan[ i ] ).release() );
    return chanDy;
}


//-------------------------------------------
// HANDLE MULTI-CHANNEL IMAGES
//-------------------------------------------


/// given a video frame, extract a set of channels used for the data term
Array<ImageGrayF> extractChannels( const ImageColorU &img, Config &varConf ) {
    float gradWeight = varConf.readFloat( "gradWeight" );
    float colorWeight = varConf.readFloat( "colorWeight" );
    Array<ImageGrayF> channels;
    ImageGrayU *gray = toGray( img ).release();
    aptr<ImageColorU> imgMedian = median( img, 3 );

    // add gray scale channel
    ImageGrayF *grayFloat = toFloat( *gray, 1.0f / 255.0f ).release(); // scale to [0,1] range
    colorWeight = colorWeight / 255.0f;
    channels.append( grayFloat );

    // add gradient channels
    if (gradWeight) {
        ImageGrayF *grayDx = xGrad( *grayFloat, 3 ).release();
        ImageGrayF *grayDy = yGrad( *grayFloat, 3 ).release();
        multiply( *grayDx, gradWeight, *grayDx );
        multiply( *grayDy, gradWeight, *grayDy );
        channels.append( grayDx );
        channels.append( grayDy );
    }

    // add color channels
    if (colorWeight) {
        int width = grayFloat->width(), height = grayFloat->height();
        ImageGrayF *gr = new ImageGrayF( width, height );
        ImageGrayF *gb = new ImageGrayF( width, height );
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++) {
                float g = imgMedian->g( x, y );
                gr->data( x, y ) = (g - imgMedian->r( x, y )) * colorWeight; // green minus red
                gb->data( x, y ) = (g - imgMedian->b( x, y )) * colorWeight; // green minus blue
            }
        channels.append( gr );
        channels.append( gb );
    }

    // clean up
    delete gray;
    return channels;
}


/// display info about an image channel
void dispChannel( const ImageGrayF &img, const char *name ) {
    float min = 0, mean = 0, max = 0;
    imageStats( img, min, mean, max );
    String caption = sprintF( "%s [%f,%f]", name, min, max );
    aptr<ImageGrayU> vis = toUChar( img );
    dispImage( *vis, caption );
}


/// displays each defined channel
void dispChannels( const Array<ImageGrayF> &chan ) {
    if (chan.count() > 0) dispChannel( chan[ 0 ], "gray" );
    if (chan.count() > 1) dispChannel( chan[ 1 ], "grad-x" );
    if (chan.count() > 2) dispChannel( chan[ 2 ], "grad-y" );
    if (chan.count() > 3) dispChannel( chan[ 3 ], "color-gr" );
    if (chan.count() > 4) dispChannel( chan[ 4 ], "color-gb" );
}


/// scale and blur each channel
Array<ImageGrayF> shrinkChannels( const Array<ImageGrayF> &chan, int scaledWidth, int scaledHeight, float blurSigma ) {
    int chanCount = chan.count();
    Array<ImageGrayF> chanScaled;

    // apply to each channel
    for (int i = 0; i < chanCount; i++) {
        const ImageGrayF &img = chan[ i ];

        // shrink the image using linear interpolation and pre-blur
        aptr<ImageGrayF> imgScaled = resize( img, scaledWidth, scaledHeight, true );

        // blur the images to ease optimization
        aptr<ImageGrayF> imgBlur = blurGauss( *imgScaled, blurSigma );

        // add to output channel set
        chanScaled.append( imgBlur.release() );
    }
    return chanScaled;
}


/// apply a Gaussian blur to each channel
Array<ImageGrayF> blurChannels( const Array<ImageGrayF> &chan, float blurSigma ) {
    int chanCount = chan.count();
    Array<ImageGrayF> chanBlurred;
    for (int i = 0; i < chanCount; i++) 
        chanBlurred.append( blurGauss( chan[ i ], blurSigma ).release() );
    return chanBlurred;
}


/// allocate multi-channel image
Array<ImageGrayF> allocChannels( int width, int height, int count, float init ) {
    Array<ImageGrayF> chan;
    for (int i = 0; i < count; i++) {
        ImageGrayF *img = new ImageGrayF( width, height );
        img->clear( init );
        chan.append( img );
    }
    return chan;
}


//-------------------------------------------
// VARIATIONAL OBJECTIVE FUNCTION
//-------------------------------------------


/// this is the robust norm used in Brox et al. (denoted as Psi)
float broxDist( float diffSqd ) {
    return sqrtf( diffSqd + 1e-6f );
}


/// derivative of robust norm w.r.t. squared diff
float broxDistDeriv( float diffSqd ) {
    return 0.5f / sqrtf( diffSqd + 1e-6f );
}


/// the variational objective function
float varEnergy( const MotionField &mf, const ImageGrayU &mask, const Array<ImageGrayF> &srcChanBlurred, const Array<ImageGrayF> &destChanBlurred, 
                 Config &varConf, bool verbose, float &dataEnergyRet, float &smoothnessEnergyRet,
                 ImageGrayF *dataEnergyMap, ImageGrayF *smoothnessEnergyMap ) {

    // get config parameters
    float smoothness = varConf.readFloat( "smoothness" );
    float flatSigma = varConf.readFloat( "flatSigma" );
    float flatSmoothness = varConf.readFloat( "flatSmoothness" );

    // initialize energy maps
    if (dataEnergyMap) dataEnergyMap->clear( 0 );
    if (smoothnessEnergyMap) smoothnessEnergyMap->clear( 0 );

    // get flatness image
    aptr<ImageGrayF> lSmoothness = localSmoothness( srcChanBlurred[ 0 ], flatSigma, flatSmoothness );

    // loop over motion field
    int chanCount = srcChanBlurred.count();
    float sumDataEnergy = 0, sumSmoothnessEnergy = 0;
    int width = mf.width(), height = mf.height(), count = 0;
    for (int y = 1; y < height - 1; y++)
        for (int x = 1; x < width - 1; x++) {

            // if in mask
            if (mask.data( x, y )) {
                count++;

                // compute smoothness term using derivatives of flow field
                float uDx = dx( mf.uRef(), x, y );
                float uDy = dy( mf.uRef(), x, y );
                float vDx = dx( mf.vRef(), x, y );
                float vDy = dy( mf.vRef(), x, y );
                float smoothnessEnergy = (smoothness + lSmoothness->data( x, y )) * broxDist( uDx * uDx + uDy * uDy + vDx * vDx + vDy * vDy );
                sumSmoothnessEnergy += smoothnessEnergy;
                if (verbose)
                    disp( 1, "u: %f, v: %f, uDx: %f, uDy: %f, vDx: %f, vDy: %f", mf.u( x, y ), mf.v( x, y ), uDx, uDy, vDx, vDy );
                if (smoothnessEnergyMap) 
                    smoothnessEnergyMap->data( x, y ) = smoothnessEnergy;

                // compute data term (no penalty if project outside image)
                float dataEnergy = 0;
                float xProj = x + mf.u( x, y );
                float yProj = y + mf.v( x, y );
                if (destChanBlurred.ref( 0 ).inBounds( xProj, yProj )) {

                    // loop over channels, computing src/dest distance for each one
                    for (int i = 0; i < chanCount; i++) {
                        float diff = (destChanBlurred[ i ].interp( xProj, yProj ) - srcChanBlurred[ i ].data( x, y ));
                        dataEnergy += broxDist( diff * diff );
                    }
                }
                sumDataEnergy += dataEnergy;
                if (dataEnergyMap)
                    dataEnergyMap->data( x, y ) = dataEnergy;
            }
        }

    // compute mean of pixel energies
    float meanEnergy = 0;
    if (count) {
        float meanSmoothnessEnergy = sumSmoothnessEnergy / count;
        float meanDataEnergy = sumDataEnergy / count;
        meanEnergy = meanSmoothnessEnergy + meanDataEnergy;
        dataEnergyRet = meanDataEnergy;
        smoothnessEnergyRet = meanSmoothnessEnergy;
    }
    return meanEnergy;
}


/// higher-level version of above; computes variational objective value
float varEnergy( const MotionField &mf, const ImageGrayU &mask, const ImageColorU &src, const ImageColorU &dest, 
                 Config &varConf, float &dataEnergy, float &smoothnessEnergy,
                 ImageGrayF *dataEnergyMap, ImageGrayF *smoothnessEnergyMap ) {

    // get config parameters
    float blurSigma = varConf.readFloat( "sigma" );

    // extract channels
    Array<ImageGrayF> srcChan = extractChannels( src, varConf );
    Array<ImageGrayF> destChan = extractChannels( dest, varConf );
    Array<ImageGrayF> srcChanBlurred = blurChannels( srcChan, blurSigma );
    Array<ImageGrayF> destChanBlurred = blurChannels( destChan, blurSigma );

    // call inner energy function
    return varEnergy( mf, mask, srcChanBlurred, destChanBlurred, varConf, false, dataEnergy, smoothnessEnergy, dataEnergyMap, smoothnessEnergyMap );
}


/// higher-level version of above; computes variational objective value
float varEnergy( const MotionField &mf, const ImageGrayU &mask, 
                 const Array<ImageGrayF> &srcChanBlurred, const Array<ImageGrayF> &destChanBlurred,
                 Config &varConf, ImageGrayF &energyMap ) {

    // call inner energy function
    int width = mf.width(), height = mf.height();
    float dataEnergy = 0, smoothnessEnergy = 0;
    ImageGrayF dataEnergyMap( width, height ), smoothnessEnergyMap( width, height );
    float energy = varEnergy( mf, mask, srcChanBlurred, destChanBlurred, varConf, false, dataEnergy, smoothnessEnergy, &dataEnergyMap, &smoothnessEnergyMap );

    // combine the energy maps
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) 
            energyMap.data( x, y ) = dataEnergyMap.data( x, y ) + smoothnessEnergyMap.data( x, y );
    return energy;
}


/// computes a map of the variational objective function (including both data and smoothness terms)
aptr<ImageGrayF> varEnergyMap( MotionField &mf, ImageGrayU &mask, 
                               Array<ImageGrayF> &srcChan, Array<ImageGrayF> &destChan, 
                               Config &varConf ) {

    // get config parameters
    float blurSigma = varConf.readFloat( "sigma" );

    // blur the channels
    // fix(clean): move outside this function?
    Array<ImageGrayF> srcChanBlurred = blurChannels( srcChan, blurSigma );
    Array<ImageGrayF> destChanBlurred = blurChannels( destChan, blurSigma );

    // get data and smoothness maps 
    aptr<ImageGrayF> energyMap( new ImageGrayF( mf.width(), mf.height() ) );
    varEnergy( mf, mask, srcChanBlurred, destChanBlurred, varConf, *energyMap );
    return energyMap;
}


} // end namespace pvl

