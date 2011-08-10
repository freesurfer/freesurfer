// Licensed under MIT license; see license.txt.

#include <sbl/image/MotionFieldUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/Triangulation.h>
#include <sbl/image/ImageUtil.h>
namespace sbl {


//-------------------------------------------
// MOTION FIELD STATS
//-------------------------------------------


/// get min/mean/max of u and v components 
void motionFieldStats( const MotionField &mf, float &uMinRet, float &vMinRet, float &uMeanRet, float &vMeanRet, float &uMaxRet, float &vMaxRet ) {
    int width = mf.width(), height = mf.height();
    float uMin = mf.u( 0, 0 ), uMax = uMin;
    float vMin = mf.v( 0, 0 ), vMax = vMin;
    float vSum = 0, uSum = 0;
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) {
            float u = mf.u( x, y );
            float v = mf.v( x, y );
            if (u < uMin) uMin = u;
            if (u > uMax) uMax = u;
            if (v < vMin) vMin = v;
            if (v > vMax) vMax = v;
            uSum += u;
            vSum += v;
        }
    uMinRet = uMin;
    vMinRet = vMin;
    uMaxRet = uMax;
    vMaxRet = vMax;
    float size = (float) (width * height);
    uMeanRet = uSum / size;
    vMeanRet = vSum / size;
}


/// display statistics for motion field
void dispMotionStats( int indent, const MotionField &mf ) {
    float uMin = 0, vMin = 0, uMean = 0, vMean = 0, uMax = 0, vMax = 0;
    motionFieldStats( mf, uMin, vMin, uMean, vMean, uMax, vMax );
    disp( indent, "mf: min: %f, %f, mean: %f, %f, max: %f, %f", uMin, vMin, uMean, vMean, uMax, vMax );
    if (mf.occDefined()) {
        disp( indent, "occ count: %d", mf.occCount() );
    }
}


/// compute mean magnitude of motion vectors
float motionFieldMag( const MotionField &mf ) {
    int width = mf.width(), height = mf.height();
    float magSum = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float u = mf.u( x, y );
            float v = mf.v( x, y );
            magSum += sqrtf( u * u + v * v );
        }
    }
    return magSum / (float) (width * height);
}


//-------------------------------------------
// MOTION FIELD VISUALIZATION
//-------------------------------------------


/// create color visualization of motion field; different colors represent different motions;
/// gray corresponds to no motion; black indicates occluded regions
aptr<ImageColorU> colorizeMotion( const MotionField &mf ) {
    int width = mf.width(), height = mf.height();
    aptr<ImageColorU> display( new ImageColorU( width, height ) );

    // compute max motion
    float max = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            float val = mf.u(x, y);
            if (val < 0) val = -val;
            if (val > max) max = val;
            val = mf.v(x, y);
            if (val < 0) val = -val;
            if (val > max) max = val;
        }

    // create image
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) {
            int u = 128 + round( mf.u( x, y ) * 128.0f / max );
            if (u > 255) u = 255;
            if (u < 0) u = 0;
            int v = 128 + round( mf.v( x, y ) * 128.0f / max );
            if (v > 255) v = 255;
            if (v < 0) v = 0;
            display->setRGB(x, y, u, 128, v);
        }

    // black out occluded areas
    if (mf.occDefined()) {
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++) 
                if (mf.occ( x, y ) > 128)
                    display->setRGB( x, y, 0, 0, 0 );
    }

    // return motion plot
    return display; 
}


/// display motion field visualization
void dispMotion( const MotionField &mf, const String &caption ) {
    aptr<ImageColorU> vis = colorizeMotion( mf );
    dispImage( *vis, caption );
}


//-------------------------------------------
// MOTION FIELD FILE I/O
//-------------------------------------------


/// load motion field from binary data file
aptr<MotionField> loadMotionField( File &file ) {
    aptr<MotionField> mf;
    
    // read header
    char header[10];
    file.readBlock( header, 3 );
    if (strcmp( header, "mf" )) {
        warning( "invalid motion field header" );
        return mf;
    }
    int version = 0;
    int width = file.readInt();
    if (width == 0) {
        version = file.readInt();
        width = file.readInt();
    }
    int height = file.readInt();

    // check dimensions
    if (width <= 0 || height <= 0)
        return mf;

    // read dimensions
    mf.reset( new MotionField( width, height ) );

    // read meta data
    mf->setBuildTime( (float) file.readDouble() );
    file.readDouble(); // note: remove this for old files

    // read motion
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            mf->setU( x, y, file.readFloat() );
            mf->setV( x, y, file.readFloat() );
        }
    }

    // read occlusion mask
    int containsOccMask = file.readInt();
    if (containsOccMask) {
        mf->enableOcc();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                mf->setOcc( x, y, file.readUChar() );
            }
        }
    }
    
    // return loaded motion field
    return mf;
}


/// load motion field from binary data file
aptr<MotionField> loadMotionField( const String &fileName ) {
    aptr<MotionField> mf;
    File file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess()) 
        mf = loadMotionField( file );
    return mf;
}


/// load motion field from floating-point image file
aptr<MotionField> loadMotionFieldImage( const String &fileName ) {
    aptr<MotionField> mf;
    aptr<ImageColorF> image = load<ImageColorF>( fileName );
    if (image.get()) {
        int width = image->width(), height = image->height();
        mf.reset( new MotionField( width, height ) );
        mf->enableOcc(); // fix(later): should check whether required first
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float u = image->r( x, y ) - 100.0f;
                float v = image->g( x, y ) - 100.0f;
                int occ = round( image->b( x, y ) );
                mf->setU( x, y, u );
                mf->setV( x, y, v );
                mf->setOcc( x, y, occ );
            }
        }
    }
    return mf;
}


/// save motion field to binary data file
void saveMotionField( const MotionField &mf, File &file ) {

    // write header
    file.writeBlock( "mf", 3 );

    // write dimensions
    int width = mf.width(), height = mf.height();
    file.writeInt( 0 );
    file.writeInt( 1 );
    file.writeInt( width );
    file.writeInt( height );

    // write meta data
    file.writeDouble( mf.buildTime() );
    file.writeDouble( 0 ); // note: remove this for old files

    // write motion
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            file.writeFloat( mf.u( x, y ));
            file.writeFloat( mf.v( x, y ));
        }
    }

    // write occlusion mask (if any)
    file.writeInt( mf.occDefined() ? 1 : 0 );
    if (mf.occDefined()) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                file.writeUChar( mf.occ( x, y ));
            }
        }
    }
}


/// save motion field to binary data file
void saveMotionField( const MotionField &mf, const String &fileName ) {
    File file( fileName, FILE_WRITE, FILE_BINARY ); 
    if (file.openSuccess()) 
        saveMotionField( mf, file );
}


/// save motion field to floating-point image file
void saveMotionFieldImage( const MotionField &mf, const String &fileName ) {
    int width = mf.width(), height = mf.height();
    ImageColorF image( width, height );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float r = 100.0f + mf.u( x, y );
            float g = 100.0f + mf.v( x, y );
            float b = (float) mf.occSafe( x, y );
            image.setRGB( x, y, r, g, b );
        }
    }
    saveImage( image, fileName );
}


/// save motion field to text file
void saveMotionFieldText( const MotionField &mf, const String &fileName ) {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {

        // write dimensions
        int width = mf.width(), height = mf.height();
        file.writeInt( width );
        file.writeInt( height );
        file.writeRawString( "\n" );

        // write motion
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                file.writeFloat( mf.u( x, y ) );
                file.writeFloat( mf.v( x, y ) );
                file.writeInt( mf.occSafe( x, y ) );
            }
            file.writeRawString( "\n" );
        }
    }
}


//-------------------------------------------
// MISC. MOTION FIELD UTILS
//-------------------------------------------


/// compute motion field divergence (du/dx + dv/dy)
aptr<ImageGrayF> motionDivergence( const MotionField &mf ) {

    // compute flow gradients
    aptr<ImageGrayF> uDx = xGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> vDy = yGrad( mf.vRef(), 3 );

    // compute divergence
    int width = mf.width(), height = mf.height();
    aptr<ImageGrayF> div( new ImageGrayF( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
             div->data( x, y ) = uDx->data( x, y ) + vDy->data( x, y );
        }
    }
    return div;
}


/// compute motion field gradient magnitude squared [(du/dx)^2 + (du/dy)^2 + (dv/dx)^2 + (dv/dy)^2]
aptr<ImageGrayF> motionGradMagSqd( const MotionField &mf ) {

    // compute flow gradients
    aptr<ImageGrayF> uDx = xGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> uDy = yGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> vDx = xGrad( mf.vRef(), 3 );
    aptr<ImageGrayF> vDy = yGrad( mf.vRef(), 3 );

    // compute gradient magnitude squared
    int width = mf.width(), height = mf.height();
    aptr<ImageGrayF> gradMagSqd( new ImageGrayF( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float uDxVal = uDx->data( x, y );
            float uDyVal = uDy->data( x, y );
            float vDxVal = vDx->data( x, y );
            float vDyVal = vDy->data( x, y );
             gradMagSqd->data( x, y ) = uDxVal * uDxVal + uDyVal * uDyVal + vDxVal * vDxVal + vDyVal * vDyVal;
        }
    }
    return gradMagSqd;
}


/// compute mean motion field gradient magnitude
float meanGradMag( const MotionField &mf ) {

    // compute flow gradients
    aptr<ImageGrayF> uDx = xGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> uDy = yGrad( mf.uRef(), 3 );
    aptr<ImageGrayF> vDx = xGrad( mf.vRef(), 3 );
    aptr<ImageGrayF> vDy = yGrad( mf.vRef(), 3 );

    // compute gradient magnitude squared
    double sum = 0;
    int width = mf.width(), height = mf.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float uDxVal = uDx->data( x, y );
            float uDyVal = uDy->data( x, y );
            float vDxVal = vDx->data( x, y );
            float vDyVal = vDy->data( x, y );
             sum += (double) sqrtf( uDxVal * uDxVal + uDyVal * uDyVal + vDxVal * vDxVal + vDyVal * vDyVal );
        }
    }
    return (float) (sum / (double) (width * height));
}


/// interpolate (triangulate) offsets to create motion field
aptr<MotionField> computeDenseMotion( const VectorF &x, const VectorF &y, const VectorF &u, const VectorF &v, 
                                      int width, int height ) {
    aptr<MotionField> mf;
#ifndef USE_CDT
    fatalError( "requires CDT" );
#else

    // for efficiency, compute motion on a smaller resolution then rescale motion field
    float mfScaleFactor = 0.5f;

    // prepare matrix of match points
    int pointCount = x.size();
    MatrixF points( pointCount, 2 );
    for (int i = 0; i < pointCount; i++) {
        points.data( i, 0 ) = x[ i ] * mfScaleFactor;
        points.data( i, 1 ) = y[ i ] * mfScaleFactor;
    }

    // fix(clean): should not be needed, but current triangulationInterpolation requires non-const u and v
    VectorF uCopy( u ), vCopy( v );
    
    // triangulate the points and linearly intepolate on the triangles
    int mfWidth = round( width * mfScaleFactor );
    int mfHeight = round( height * mfScaleFactor );
    aptr<ImageGrayF> uMap = triangulationInterpolation( points, uCopy, mfWidth, mfHeight );
    aptr<ImageGrayF> vMap = triangulationInterpolation( points, vCopy, mfWidth, mfHeight );
    mf.reset( new MotionField( mfWidth, mfHeight ) );
    for (int y = 0; y < mfHeight; y++) {
        for (int x = 0; x < mfWidth; x++) {
            mf->setU( x, y, uMap->data( x, y ) );
            mf->setV( x, y, vMap->data( x, y ) );
        }
    }
    mf->resize( width, height, false );
#endif
    return mf;
}


/// compute inverse motion field using point sample triangulation;
/// sampleStep determines distance between samples
aptr<MotionField> invertMotionField( const MotionField &mf, const ImageGrayU &mask, int sampleStep ) {
    int width = mf.width(), height = mf.height();
    VectorF xVect, yVect, uVect, vVect;
    for (int y = 0; y < height; y += sampleStep) {
        for (int x = 0; x < width; x += sampleStep) {
            if (mask.data( x, y )) {
                float u = mf.u( x, y );
                float v = mf.v( x, y );
                float xDest = x + u;
                float yDest = y + v;
                if (mask.inBounds( xDest, yDest )) {
                    xVect.append( xDest );
                    yVect.append( yDest );
                    uVect.append( -u );
                    vVect.append( -v );
                }
            }
        }
    }
    return computeDenseMotion( xVect, yVect, uVect, vVect, width, height );
}


} // end namespace sbl

