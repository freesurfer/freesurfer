// Licensed under MIT license; see license.txt.

#include <sbl/image/ImageSeqUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
namespace sbl {


/// create a sequence of images
void initImageSeq( ImageGrayUSeq &seq, int width, int height, int length, bool clear, int clearValue ) {
    assertAlways( seq.count() == 0 );
    for (int z = 0; z < length; z++) {
        ImageGrayU *img = new ImageGrayU( width, height );
        if (clear)
            img->clear( clearValue );
        seq.append( img );
    }
}


/// create a sequence of images
void initImageSeq( ImageGrayFSeq &seq, int width, int height, int length, bool clear, float clearValue ) {
    assertAlways( seq.count() == 0 );
    for (int z = 0; z < length; z++) {
        ImageGrayF *img = new ImageGrayF( width, height );
        if (clear)
            img->clear( clearValue );
        seq.append( img );
    }
}


/// create a sequence of images
void initImageSeq( ImageGrayISeq &seq, int width, int height, int length, bool clear, int clearValue ) {
    assertAlways( seq.count() == 0 );
    for (int z = 0; z < length; z++) {
        ImageGrayI *img = new ImageGrayI( width, height );
        if (clear)
            img->clear( clearValue );
        seq.append( img );
    }
}


/// resize (in x, y) a sequence of images
void resizeSeqXY( ImageGrayUSeq &seq, int newWidth, int newHeight ) {
    for (int z = 0; z < seq.count(); z++) 
        seq.set( z, resize( seq[ z ], newWidth, newHeight, true ).release() );
}


/// resize (in x, y) a sequence of images
void resizeSeqXY( ImageGrayFSeq &seq, int newWidth, int newHeight ) {
    for (int z = 0; z < seq.count(); z++) 
        seq.set( z, resize( seq[ z ], newWidth, newHeight, true ).release() );
}


/// resize (in z) a sequence of images
void resizeSeqZ( ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq, int newLength ) {
    int oldLength = inSeq.count();
    assertAlways( oldLength );
    int width = inSeq[ 0 ].width(), height = inSeq[ 0 ].height();
    assertAlways( outSeq.count() == 0 );
    float zScale = (float) (oldLength - 1) / (float) (newLength - 1);
    for (int z = 0; z < newLength; z++) {
        ImageGrayU *out = new ImageGrayU( width, height );
        float zOld = z * zScale;
        int zOldInt = (int) zOld;
        if (zOldInt >= oldLength - 1) {
            const ImageGrayU &in = inSeq[ oldLength - 1 ];
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    out->data( x, y ) = in.data( x, y );
                }
            }
        } else {
            const ImageGrayU &in1 = inSeq[ zOldInt ];
            const ImageGrayU &in2 = inSeq[ zOldInt + 1 ];
            float frac = zOld - zOldInt;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    out->data( x, y ) = round( frac * (float) in2.data( x, y ) + (1.0f - frac) * (float) in1.data( x, y ) );
                }
            }
        }
        outSeq.append( out );
    }
}


/// resize (in z) a sequence of images
void resizeSeqZ( ImageGrayFSeq &inSeq, ImageGrayFSeq &outSeq, int newLength ) {
    int oldLength = inSeq.count();
    assertAlways( oldLength );
    int width = inSeq[ 0 ].width(), height = inSeq[ 0 ].height();
    assertAlways( outSeq.count() == 0 );
    float zScale = (float) (oldLength - 1) / (float) (newLength - 1);
    for (int z = 0; z < newLength; z++) {
        ImageGrayF *out = new ImageGrayF( width, height );
        float zOld = z * zScale;
        int zOldInt = (int) zOld;
        if (zOldInt >= oldLength - 1) {
            const ImageGrayF &in = inSeq[ oldLength - 1 ];
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    out->data( x, y ) = in.data( x, y );
                }
            }
        } else {
            const ImageGrayF &in1 = inSeq[ zOldInt ];
            const ImageGrayF &in2 = inSeq[ zOldInt + 1 ];
            float frac = zOld - zOldInt;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    out->data( x, y ) = frac * in2.data( x, y ) + (1.0f - frac) * in1.data( x, y );
                }
            }
        }
        outSeq.append( out );
    }
    assertAlways( outSeq.count() == newLength );
}


/// blur (in x, y) a sequence of images
void blurGaussSeqXY( ImageGrayUSeq &seq, float sigma ) {
    for (int z = 0; z < seq.count(); z++) 
        seq.set( z, blurGauss( seq[ z ], sigma ).release() );
}


/// blur (in x, y) a sequence of images
void blurGaussSeqXY( ImageGrayFSeq &seq, float sigma ) {
    for (int z = 0; z < seq.count(); z++) 
        seq.set( z, blurGauss( seq[ z ], sigma ).release() );
}


/// blur (in z) a sequences of images
void blurGaussSeqZ( const ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq, float sigma ) {
    int length = inSeq.count();
    int width = inSeq[ 0 ].width(), height = inSeq[ 0 ].height();

    // generate Gaussian table
    int tableRadius = (int) sigma * 3;
    double gFactor = gaussFactor( sigma );
    int tableSize = tableRadius * 2 + 1;
    double *table = new double[ tableSize ];
    for (int i = 0; i < tableSize; i++) {
        double x = (double) (i - tableRadius);
        table[ i ] = gauss( x * x, gFactor );
    }

    // compute filtered version of each image
    for (int z = 0; z < length; z++) {
        ImageGrayU *output = new ImageGrayU( width, height );
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {

                // fix(faster): try moving this loop outside the x,y loops (store sum in image)
                double sum = 0, sumWt = 0;
                for (int j = -tableRadius; j <= tableRadius; j++) {
                    int srcInd = z + j;
                    if (srcInd >= 0 && srcInd < length) {
                        double wt = table[ j + tableRadius ];
                        sum += wt * inSeq[ srcInd ].data( x, y );
                        sumWt += wt;
                    }
                }
                output->data( x, y ) = bound( sbl::round( sum / sumWt ), 0, 255 );
            }
        }
        outSeq.append( output );
    }
    delete [] table;
}


/// blur (in z) a sequences of images
void blurGaussSeqZ( const ImageGrayFSeq &inSeq, ImageGrayFSeq &outSeq, float sigma ) {
    int length = inSeq.count();
    int width = inSeq[ 0 ].width(), height = inSeq[ 0 ].height();

    // generate Gaussian table
    int tableRadius = (int) sigma * 3;
    double gFactor = gaussFactor( sigma );
    int tableSize = tableRadius * 2 + 1;
    double *table = new double[ tableSize ];
    for (int i = 0; i < tableSize; i++) {
        double x = (double) (i - tableRadius);
        table[ i ] = gauss( x * x, gFactor );
    }

    // compute filtered version of each image
    for (int z = 0; z < length; z++) {
        ImageGrayF *output = new ImageGrayF( width, height );
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {

                // fix(faster): try moving this loop outside the x,y loops (store sum in image)
                double sum = 0, sumWt = 0;
                for (int j = -tableRadius; j <= tableRadius; j++) {
                    int srcInd = z + j;
                    if (srcInd >= 0 && srcInd < length) {
                        double wt = table[ j + tableRadius ];
                        sum += wt * inSeq[ srcInd ].data( x, y );
                        sumWt += wt;
                    }
                }
                output->data( x, y ) = (float) (sum / sumWt);
            }
        }
        outSeq.append( output );
    }
    delete [] table;
}


/// perform bilinear interpolation to find value in image sequence
float interp( const Array<ImageGrayF> &seq, float x, float y, float z ) {
    int zInt = (int) z;
    float frac = z - (float) zInt;
    float v = 0;
    if (zInt == seq.count() - 1) {
        v = seq[ zInt ].interp( x, y );
    } else {
        assertAlways( zInt >= 0 && zInt < seq.count() - 1 );
        assertAlways( frac >= -0.00001 && frac <= 1.00001 );
        float v1 = seq[ zInt ].interp( x, y );
        float v2 = seq[ zInt + 1 ].interp( x, y );
        v = frac * v2 + (1.0f - frac) * v1; 
    } 
    return v;
}


/// perform bilinear interpolation to find value in image sequence
float interp( const Array<ImageGrayU> &seq, float x, float y, float z ) {
    int zInt = (int) z;
    float frac = z - (float) zInt;
    float v = 0;
    if (zInt == seq.count() - 1) {
        v = seq[ zInt ].interp( x, y );
    } else {
        assertAlways( zInt >= 0 && zInt < seq.count() - 1 );
        //assertAlways( frac >= -0.00001 && frac <= 1.00001 );
        if (frac < -0.00001 || frac > 1.00001)
           fatalError( "z: %.9f, zInt: %d, frac: %.9f", z, zInt, frac );
        float v1 = seq[ zInt ].interp( x, y );
        float v2 = seq[ zInt + 1 ].interp( x, y );
        v = frac * v2 + (1.0f - frac) * v1; 
    } 
    return v;
}


/// transpose the Y and Z axes in the image sequence
void transposeYZ( const ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq ) {
    int length = inSeq.count();
    assertAlways( length );
    int width = inSeq[ 0 ].width(), height = inSeq[ 0 ].height();
    assertAlways( outSeq.count() == 0 );
    for (int y = 0; y < height; y++) {
        ImageGrayU *out = new ImageGrayU( width, length );
        for (int z = 0; z < length; z++) {
            const ImageGrayU &in = inSeq[ z ];
            for (int x = 0; x < width; x++) {
                out->data( x, z ) = in.data( x, y );
            }
        }
        outSeq.append( out );
    }
}


// location to be processed for 3D flood fill
struct FillNode3D {
    int x;
    int y;
    int z;
};


/// fills starting at the given point, all values within the given range;
/// the fill value must *not* be within the fill range
int floodFillXYZ( ImageGrayUSeq &seq, int minRegionValue, int maxRegionValue, int fillValue, int x, int y, int z, 
                  float *xCentRet, float *yCentRet, float *zCentRet, 
                  int *xMinRet, int *xMaxRet, int *yMinRet, int *yMaxRet, int *zMinRet, int *zMaxRet ) {

    // alloc stack
    int width = seq[ 0 ].width(), height = seq[ 0 ].height(), length = seq.count();
    int stackSize = width * height * length * 3;
    FillNode3D *fillStack = new FillNode3D[ stackSize ];

    // region stats
    int xMin = width, xMax = 0;
    int yMin = height, yMax = 0;
    int zMin = length, zMax = 0;
    double xSum = 0, ySum = 0, zSum = 0;
    int count = 0;
    bool getStats = xCentRet || yCentRet || zCentRet || xMinRet || xMaxRet || yMinRet || yMaxRet || zMinRet || zMaxRet;

    // add first node to stack
    fillStack[ 0 ].x = x;
    fillStack[ 0 ].y = y;
    fillStack[ 0 ].z = z;
    int stackTop = 1;
    while (stackTop) {

        // pop node
        stackTop--;
        int cx = fillStack[ stackTop ].x;
        int cy = fillStack[ stackTop ].y;
        int cz = fillStack[ stackTop ].z;

        // check coords
        if (cx >= 0 && cx < width && cy >= 0 && cy < height && cz >= 0 && cz < length) {
            int current = seq[ cz ].data( cx, cy );
            if (current >= minRegionValue && current <= maxRegionValue) {
                seq[ cz ].data( cx, cy ) = fillValue;

                // accumulate stats
                count++;
                if (getStats) {
                    xSum += cx;
                    ySum += cy;
                    zSum += cz;
                    if (cx < xMin) xMin = cx;
                    if (cx > xMax) xMax = cx;
                    if (cy < yMin) yMin = cy;
                    if (cy > yMax) yMax = cy;
                    if (cz < zMin) zMin = cz;
                    if (cz > zMax) zMax = cz;
                }

                // recurse: push neighbors
                fillStack[ stackTop ].x = cx + 1;
                fillStack[ stackTop ].y = cy;
                fillStack[ stackTop ].z = cz;
                stackTop++;
                fillStack[ stackTop ].x = cx - 1;
                fillStack[ stackTop ].y = cy;
                fillStack[ stackTop ].z = cz;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy + 1;
                fillStack[ stackTop ].z = cz;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy - 1;
                fillStack[ stackTop ].z = cz;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy;
                fillStack[ stackTop ].z = cz + 1;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy;
                fillStack[ stackTop ].z = cz - 1;
                stackTop++;
                if (stackTop > stackSize) {
                    fatalError( "fillRegion stack overflow" );
                }
            }
        }
    }

    // clean up
    delete [] fillStack;

    // return results
    if (xCentRet) *xCentRet = (float) xSum / (float) count;
    if (yCentRet) *yCentRet = (float) ySum / (float) count;
    if (zCentRet) *zCentRet = (float) zSum / (float) count;
    if (xMinRet) *xMinRet = xMin;
    if (xMaxRet) *xMaxRet = xMax;
    if (yMinRet) *yMinRet = yMin;
    if (yMaxRet) *yMaxRet = yMax;
    if (zMinRet) *zMinRet = zMin;
    if (zMaxRet) *zMaxRet = zMax;
    return count;
}


} // end namespace sbl

