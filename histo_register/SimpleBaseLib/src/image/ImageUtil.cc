// Licensed under MIT license; see license.txt.

#include <cstddef> 
#include <sbl/image/ImageUtil.h>
#include <sbl/core/Command.h> // for filter registry
#include <sbl/math/MathUtil.h>
#include <sbl/math/MatrixUtil.h> // for mutual info
#include <sbl/image/Filter.h> // for filter registry
#ifdef USE_OPENCV
    #include <opencv/cv.h>
    #include <opencv/highgui.h>
#endif
#ifdef USE_GUI
    #include <sbl/gui/ImageSeqViewer.h>
#endif
namespace sbl {


//-------------------------------------------
// IMAGE CONVERSION
//-------------------------------------------


/// convert color image to gray image
aptr<ImageGrayU> toGray( const ImageColorU &input ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
#ifdef USE_OPENCV
    cvCvtColor( input.iplImage(), output->iplImage(), CV_BGR2GRAY );
#endif
    return output;
}


/// convert gray image to color image
aptr<ImageColorU> toColor( const ImageGrayU &input ) {
    int width = input.width(), height = input.height();
    aptr<ImageColorU> output( new ImageColorU( width, height ) );
#ifdef USE_OPENCV
    cvCvtColor( input.iplImage(), output->iplImage(), CV_GRAY2BGR );
#endif
    return output;
}


/// convert 8-bit image to float image
// fix(faster): use opencv?
aptr<ImageGrayF> toFloat( const ImageGrayU &input, float scaleFactor ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayF> output( new ImageGrayF( width, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            output->data( x, y ) = (float) input.data( x, y ) * scaleFactor;
    return output;
}


/// convert 8-bit image to float image
// fix(faster): use opencv?
aptr<ImageColorF> toFloat( const ImageColorU &input, float scaleFactor ) {
    int width = input.width(), height = input.height();
    aptr<ImageColorF> output( new ImageColorF( width, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            for (int c = 0; c < 3; c++) 
                output->data( x, y, c ) = (float) input.data( x, y, c ) * scaleFactor;
    return output;
}


/// convert float image to 8-bit image, automatically scaling values 
aptr<ImageGrayU> toUChar( const ImageGrayF &input ) {
    int width = input.width(), height = input.height();
    assertAlways( width && height );
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
    float min = input.data( 0, 0 );
    float max = input.data( 0, 0 );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float v = input.data( x, y );
            if (v < min)
                min = v;
            if (v > max)
                max = v;
        }
    }
    float factor = 255.0f / (max - min);
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            output->data( x, y ) = bound( round( (input.data( x, y ) - min) * factor ), 0, 255 );
    return output;
}


/// convert float image to 8-bit image, automatically scaling values 
aptr<ImageGrayU> toUChar( const ImageGrayF &input, float scaleFactor ) {
    int width = input.width(), height = input.height();
    assertAlways( width && height );
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            output->data( x, y ) = bound( round( input.data( x, y ) * scaleFactor ), 0, 255 );
    return output;
}


/// convert float image to 8-bit image, automatically scaling values 
aptr<ImageColorU> toUChar( const ImageColorF &input, float scaleFactor ) {
    int width = input.width(), height = input.height();
    assertAlways( width && height );
    aptr<ImageColorU> output( new ImageColorU( width, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
            for (int c = 0; c < 3; c++) 
                output->data( x, y, c ) = bound( round( input.data( x, y, c ) * scaleFactor ), 0, 255 );
    return output;
}


//-------------------------------------------
// IMAGE FILTERS
//-------------------------------------------


/// blur using box filter
template <typename ImageType> aptr<ImageType> blurBox( const ImageType &input, int boxSize ) {
    int width = input.width(), height = input.height();
    aptr<ImageType> output( new ImageType( width, height ) );
#ifdef USE_OPENCV
    cvSmooth( input.iplImage(), output->iplImage(), CV_BLUR, boxSize, boxSize );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU> blurBox( const ImageGrayU &input, int boxSize );
template aptr<ImageGrayF> blurBox( const ImageGrayF &input, int boxSize );
template aptr<ImageColorU> blurBox( const ImageColorU &input, int boxSize );
template aptr<ImageColorF> blurBox( const ImageColorF &input, int boxSize );


/// blur using Gaussian filter
template <typename ImageType> aptr<ImageType> blurGauss( const ImageType &input, float sigma ) {
    assertAlways( sigma > 0 );
    int width = input.width(), height = input.height();
    aptr<ImageType> output( new ImageType( width, height ) );
#ifdef USE_OPENCV
    cvSmooth( input.iplImage(), output->iplImage(), CV_GAUSSIAN, 0, 0, sigma, sigma );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU> blurGauss( const ImageGrayU &input, float sigma );
template aptr<ImageGrayF> blurGauss( const ImageGrayF &input, float sigma );
template aptr<ImageColorU> blurGauss( const ImageColorU &input, float sigma );
template aptr<ImageColorF> blurGauss( const ImageColorF &input, float sigma );


/// apply median filter (set each pixel to median of neighbors)
template <typename ImageType> aptr<ImageType> median( const ImageType &input, int apertureSize ) {
    int width = input.width(), height = input.height();
    aptr<ImageType> output( new ImageType( width, height ) );
#ifdef USE_OPENCV
    cvSmooth( input.iplImage(), output->iplImage(), CV_MEDIAN, apertureSize, apertureSize );
#endif
    return output;
}
template aptr<ImageGrayU> median( const ImageGrayU &input, int apertureSize );
template aptr<ImageGrayF> median( const ImageGrayF &input, int apertureSize );
template aptr<ImageColorU> median( const ImageColorU &input, int apertureSize );
template aptr<ImageColorF> median( const ImageColorF &input, int apertureSize );


/// image x gradient
template <typename ImageType> aptr<ImageGrayF> xGrad( const ImageType &input, int apertureSize ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayF> output( new ImageGrayF( width, height ) );
#ifdef USE_OPENCV
    cvSobel( input.iplImage(), output->iplImage(), 1, 0, apertureSize );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayF> xGrad( const ImageGrayU &input, int apertureSize );
template aptr<ImageGrayF> xGrad( const ImageGrayF &input, int apertureSize );


/// image y gradient
template <typename ImageType> aptr<ImageGrayF> yGrad( const ImageType &input, int apertureSize ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayF> output( new ImageGrayF( width, height ) );
#ifdef USE_OPENCV
    cvSobel( input.iplImage(), output->iplImage(), 0, 1, apertureSize );
#endif
    return output;
}
template aptr<ImageGrayF> yGrad( const ImageGrayU &input, int apertureSize );
template aptr<ImageGrayF> yGrad( const ImageGrayF &input, int apertureSize );


/// compute magnitude of each gradient vector
aptr<ImageGrayF> gradientMagnitude( const ImageGrayU &input, int apertureSize ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayF> gx = xGrad( input, apertureSize );
    aptr<ImageGrayF> gy = yGrad( input, apertureSize );
    aptr<ImageGrayF> mag( new ImageGrayF( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float gxi = gx->data( x, y );
            float gyi = gy->data( x, y );
            mag->data( x, y ) = sqrtf( gxi * gxi + gyi * gyi ); 
        }
    }
    return mag;
}


/// apply threshold to image values
template <typename ImageType> aptr<ImageType> threshold( const ImageType &input, float thresh, bool invert ) {
    int width = input.width(), height = input.height();
    aptr<ImageType> output( new ImageType( width, height ) );
#ifdef USE_OPENCV
    cvThreshold( input.iplImage(), output->iplImage(), thresh, 255, invert ? CV_THRESH_BINARY_INV : CV_THRESH_BINARY );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU> threshold( const ImageGrayU &input, float thresh, bool invert );
template aptr<ImageGrayF> threshold( const ImageGrayF &input, float thresh, bool invert );


/// invert a floating point image, assuming values in [0, 1]
// fix(faster): use openCV
aptr<ImageGrayF> invert( const ImageGrayF &input ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayF> output( new ImageGrayF( width, height ) );
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) 
            output->data( x, y ) = 1.0f - input.data( x, y );
    return output;
}


/// apply a box filter and then threshold
aptr<ImageGrayU> blurBoxAndThreshold( const ImageGrayU &input, int boxSize, int thresh ) {
    int width = input.width(), height = input.height();
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
#ifdef USE_OPENCV
    cvSmooth( input.iplImage(), output->iplImage(), CV_BLUR, boxSize, boxSize );
    cvThreshold( output->iplImage(), output->iplImage(), thresh, 255, CV_THRESH_BINARY );
#else
    fatalError( "not implemented" );
#endif
    return output;
}


/// apply a box filter and then threshold
aptr<ImageGrayU> blurBoxAndThreshold( const ImageGrayF &input, int boxSize, float thresh ) {
    int width = input.width(), height = input.height();
    ImageGrayF smooth( width, height );
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
#ifdef USE_OPENCV
    cvSmooth( input.iplImage(), smooth.iplImage(), CV_BLUR, boxSize, boxSize );
    cvThreshold( smooth.iplImage(), output->iplImage(), thresh, 255, CV_THRESH_BINARY );
#else
    fatalError( "not implemented" );
#endif
    return output;
}



/// multiply each pixel value by the specified factor
void multiply( const ImageGrayF &input, float factor, ImageGrayF &output ) {
    int width = input.width(), height = input.height();
    assertAlways( output.width() == width && output.height() == height );
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++) 
            output.data( x, y ) = input.data( x, y ) * factor;
}


/// perform local brightness normalization
aptr<ImageColorU> normalizeLocalBrightness( const ImageColorU &input, int windowSize ) {
    aptr<ImageGrayU> gray = toGray( input );
    aptr<ImageGrayU> blur = blurBox( *gray, windowSize );
    int width = input.width(), height = input.height();
    aptr<ImageColorU> output( new ImageColorU( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int offset = 128 - blur->data( x, y );
            for (int c = 0; c < 3; c++) {
                int val = input.data( x, y, c ) + offset;
                if (val < 0)
                    val = 0;
                if (val > 255)
                    val = 255;
                output->data( x, y, c ) = val;
            }
        }
    }
    return output;
}


/// perform local brightness normalization
aptr<ImageGrayU> normalizeLocalBrightness( const ImageGrayU &input, int windowSize ) {
    aptr<ImageGrayU> blur = blurBox( input, windowSize );
    int width = input.width(), height = input.height();
    aptr<ImageGrayU> output( new ImageGrayU( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int offset = 128 - blur->data( x, y );
            int val = input.data( x, y ) + offset;
            if (val < 0)
                val = 0;
            if (val > 255)
                val = 255;
            output->data( x, y ) = val;
        }
    }
    return output;
}


//-------------------------------------------
// IMAGE STATISTICS
//-------------------------------------------


/// compute mean pixel value
float mean( const ImageGrayU &img ) {
    int width = img.width(), height = img.height();
    double sum = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            sum += (double) img.data( x, y );
    return (float) (sum / (double) (width * height));
}


/// compute mean pixel value
float mean( const ImageGrayF &img ) {
    int width = img.width(), height = img.height();
    double sum = 0;
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            sum += img.data( x, y );
    return (float) (sum / (double) (width * height));
}


/// compute mean pixel value
float mean( const ImageColorU &img ) {
    int width = img.width(), height = img.height();
    int sum = 0;
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++)
            for (int c = 0; c < 3; c++) 
                sum += img.data( x, y, c );
    return (float) ((double) sum / (double) (width * height * 3));
}


/// compute mean pixel value of each color channel
void channelMean( const ImageColorU &img, float &rMean, float &gMean, float &bMean ) {
    int width = img.width(), height = img.height();
    int rSum = 0, gSum = 0, bSum = 0; // fix(later): could overflow
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            rSum += img.r( x, y );
            gSum += img.g( x, y );
            bSum += img.b( x, y );
        }
    }
    int size = width * height;
    if (size) {
        rMean = (float) (rSum / size);
        gMean = (float) (gSum / size);
        bMean = (float) (bSum / size);
    } else {
        rMean = 0;
        gMean = 0;
        bMean = 0;
    }
}


/// compute mean absolute difference between pixel values
// fix(faster): add OpenCV version
float meanAbsDiff( const ImageColorU &img1, const ImageColorU &img2 ) {
    int width = img1.width(), height = img1.height();
    assertAlways( img2.width() == width && img2.height() == height );
    int sum = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                int diff = img1.data( x, y, c ) - img2.data( x, y, c );
                if (diff < 0)
                    diff = -diff;
                sum += diff;
            }
        }
    }
    return (float) ((double) sum / (double) width * height * 3);
}


/// compute mean absolute difference between pixel values
// fix(faster): add OpenCV version
float meanAbsDiff( const ImageGrayU &img1, const ImageGrayU &img2, int xBorder, int yBorder ) {
    int width = img1.width(), height = img1.height();
    assertAlways( img2.width() == width && img2.height() == height );
    double sum = 0;
    for (int y = yBorder; y < height - yBorder; y++) {
        for (int x = xBorder; x < width - xBorder; x++) {
            int diff = img1.data( x, y ) - img2.data( x, y );
            if (diff < 0)
                diff = -diff;
            sum += diff;
        }
    }
    int size = (width - 2 * xBorder) * (height - 2 * yBorder);
    return (float) ((double) sum / (double) size);
}


/// compute mutual info between a pair of images
float mutualInfo( const ImageGrayU &img1, const ImageGrayU &img2, int xBorder, int yBorder, int bucketCount ) {
    int width = img1.width(), height = img1.height();
    assertAlways( img2.width() == width && img2.height() == height );

    // build histogram
    MatrixD histogram( bucketCount, bucketCount );
    histogram.clear( 0 );
    for (int y = yBorder; y < height - yBorder; y++) {
        for (int x = xBorder; x < width - xBorder; x++) {
            int i = bound( img1.data( x, y ) * (bucketCount - 1) / 255, 0, bucketCount - 1 );
            int j = bound( img2.data( x, y ) * (bucketCount - 1) / 255, 0, bucketCount - 1 );
            histogram.data( i, j )++;
        }
    }

    // normalize histogram
    double factor = 1.0 / histogram.sum();
    multiply( histogram, factor, histogram );

    // compute joint entropy
    double hxy = 0;
    for (int j = 0; j < bucketCount; j++) {
        for (int i = 0; i < bucketCount; i++) {
            double v = histogram.data( i, j );
            if (v > 1e-10)
                hxy -= v * log( v );
        }
    }

    // get marginal entropies
    double hx = entropy( rowSum( histogram ) );
    double hy = entropy( colSum( histogram ) );

    // return mutual information
    return (float) (hx + hy - hxy);
}


/// compute bounds of the non-zero mask region
void maskBounds( const ImageGrayU &mask, int &xMin, int &xMax, int &yMin, int &yMax ) {
    int width = mask.width(), height = mask.height();
    xMin = width;
    xMax = -1;
    yMin = height;
    yMax = -1;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y )) {
                if (x < xMin)
                    xMin = x;
                if (x > xMax)
                    xMax = x;
                if (y < yMin)
                    yMin = y;
                if (y > yMax)
                    yMax = y;
            }
        }
    }
}


/// count number of non-zero entries in mask
int maskCount( const ImageGrayU &mask ) {
    int width = mask.width(), height = mask.height();
    int count = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y ))
                count++;
        }
    }
    return count;
}



/// compute min/mean/max pixel values
void imageStats( const ImageGrayU &img, int &minRet, float &meanRet, int &maxRet ) {
    double sum = 0;
    int min = img.data( 0, 0 );
    int max = min;
    int width = img.width(), height = img.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int v = img.data( x, y );
            sum += (double) v;
            if (v < min)
                min = v;
            if (v > max)
                max = v;
        }
    }
    meanRet = (float) (sum / (double) (width * height));
    minRet = min;
    maxRet = max;
}


/// compute min/mean/max pixel values
void imageStats( const ImageGrayF &img, float &minRet, float &meanRet, float &maxRet ) {
    double sum = 0;
    float min = img.data( 0, 0 );
    float max = min;
    int width = img.width(), height = img.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float v = img.data( x, y );
            sum += v;
            if (v < min)
                min = v;
            if (v > max)
                max = v;
        }
    }
    meanRet = (float) (sum / (double) (width * height));
    minRet = min;
    maxRet = max;
}


/// compute histogram of image pixel values
VectorI imageHistogram( const ImageGrayU &image, int xMin, int xMax, int yMin, int yMax ) {
    VectorI hist( 256 );
    hist.clear( 0 );
    for (int y = yMin; y <= yMax; y++)
        for (int x = xMin; x <= xMax; x++)
            hist[ image.data( x, y ) ]++;
    return hist;
}


//-------------------------------------------
// IMAGE FILE I/O
//-------------------------------------------


/// save image to file (format determined by extension)
template <typename ImageType> void saveImage( const ImageType &img, const String &fileName ) {
#ifdef USE_OPENCV
    cvSaveImage( fileName.c_str(), img.iplImage() );
#else
    fatalError( "not implemented" );
#endif
}
template void saveImage( const ImageGrayU &img, const String &fileName );
template void saveImage( const ImageGrayF &img, const String &fileName );
template void saveImage( const ImageColorU &img, const String &fileName );
template void saveImage( const ImageColorF &img, const String &fileName );


/// load image from file (format determined by extension)
template <typename ImageType> aptr<ImageType> load( const String &fileName ) {
    aptr<ImageType> img;
#ifdef USE_OPENCV
    IplImage *iplImg = cvLoadImage( fileName.c_str(), CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR );
    if (iplImg) {
        img.reset( new ImageType( iplImg->width, iplImg->height ) );
        cvConvertImage( iplImg, img->iplImage(), CV_CVTIMG_FLIP ); // image is loaded using TL origin, but we assume BL origin
        cvReleaseImage( &iplImg );
    } else {
        warning( "unable to load image: %s", fileName.c_str() );
    }
#else
    fatalError( "not implemented" );
#endif
    return img;
}
template aptr<ImageGrayU> load<ImageGrayU>( const String &fileName );
template aptr<ImageGrayF> load<ImageGrayF>( const String &fileName );
template aptr<ImageColorU> load<ImageColorU>( const String &fileName );
template aptr<ImageColorF> load<ImageColorF>( const String &fileName );


//-------------------------------------------
// FLOOD FILL
//-------------------------------------------


// struct used to keep track of pending points for flood fill
struct FillNode {
    int x;
    int y;
};


/// fills starting at the given point, all values within the given range;
/// the fill value must *not* be within the fill range
int floodFill( ImageGrayU &img, int minRegionValue, int maxRegionValue, int fillValue, int x, int y, 
               float *xCentRet, float *yCentRet, int *xMinRet, int *xMaxRet, int *yMinRet, int *yMaxRet ) {

    // alloc stack
    int width = img.width(), height = img.height();
    int stackSize = width * height * 3;
    FillNode *fillStack = new FillNode[ stackSize ];

    // region stats
    int xMin = width, xMax = 0;
    int yMin = height, yMax = 0;
    double xSum = 0, ySum = 0;
    int count = 0;
    bool getStats = xCentRet || yCentRet || xMinRet || xMaxRet || yMinRet || yMaxRet;

    // add first node to stack
    fillStack[ 0 ].x = x;
    fillStack[ 0 ].y = y;
    int stackTop = 1;
    while (stackTop) {

        // pop node
        stackTop--;
        int cx = fillStack[ stackTop ].x, cy = fillStack[ stackTop ].y;

        // check coords
        if (cx >= 0 && cx < width && cy >= 0 && cy < height) {
            int current = img.data( cx, cy );
            if (current >= minRegionValue && current <= maxRegionValue) {
                img.data( cx, cy ) = fillValue;

                // accumulate stats
                count++;
                if (getStats) {
                    xSum += cx;
                    ySum += cy;
                    if (cx < xMin) xMin = cx;
                    if (cx > xMax) xMax = cx;
                    if (cy < yMin) yMin = cy;
                    if (cy > yMax) yMax = cy;
                }

                // recurse: push neighbors
                fillStack[ stackTop ].x = cx + 1;
                fillStack[ stackTop ].y = cy;
                stackTop++;
                fillStack[ stackTop ].x = cx - 1;
                fillStack[ stackTop ].y = cy;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy + 1;
                stackTop++;
                fillStack[ stackTop ].x = cx;
                fillStack[ stackTop ].y = cy - 1;
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
    if (xMinRet) *xMinRet = xMin;
    if (xMaxRet) *xMaxRet = xMax;
    if (yMinRet) *yMinRet = yMin;
    if (yMaxRet) *yMaxRet = yMax;
    return count;
}


/// keep only mask components that have a pixel count within the given range (assumes mask has only 0 and 255 values)
void filterMaskComponents( ImageGrayU &mask, int minSize, int maxSize ) {
    int width = mask.width(), height = mask.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y ) == 255) {
                int count = floodFill( mask, 255, 255, 200, x, y );
                if (count < minSize || count > maxSize) {
                    floodFill( mask, 200, 200, 0, x, y );
                }
            }
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y ) == 200)
                mask.data( x, y ) = 255;
        }
    }
}


/// fill holes in mask smaller than given size (assumes mask has only 0 and 255 values)
void fillMaskHoles( ImageGrayU &mask, int maxSize ) {
    int width = mask.width(), height = mask.height();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y ) == 0) {
                int count = floodFill( mask, 0, 0, 100, x, y );
                if (count < maxSize) {
                    floodFill( mask, 100, 100, 255, x, y );
                }
            }
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (mask.data( x, y ) == 100)
                mask.data( x, y ) = 0;
        }
    }
}


//-------------------------------------------
// OTHER IMAGE UTILS
//-------------------------------------------


/// returns true if specified line intersects non-zero mask values;
/// (assumes points are inside image bounds)
bool lineIntersects( const ImageGrayU &mask, float x1, float y1, float x2, float y2 ) {
    bool found = false;
    float xDiff = x2 - x1;
    float yDiff = y2 - y1;
    float len = sqrtf( xDiff * xDiff + yDiff * yDiff );
    float xStep = xDiff / len;
    float yStep = yDiff / len;
    for (float step = 0; step < len; step++) {
        int x = round( x1 + step * xStep );
        int y = round( y1 + step * yStep );
        if (mask.data( x, y )) {
            found = true;
            break;
        }
    }
    return found;
}


/// returns true if rectange intersects non-zero mask values;
/// (assumes points are inside image bounds)
bool rectIntersects( const ImageGrayU &mask, int x1, int y1, int x2, int y2 ) {
    assertAlways( y2 >= y1 );
    assertAlways( x2 >= x1 );
    for (int y = y1; y <= y2; y++) {
        for (int x = x1; x <= x2; x++) {
            if (mask.data( x, y ))
                return true;
        }
    }
    return false;
}


/// draw the an iso-contour of the mask on the given output image
void drawMaskBoundary( ImageColorU &output, const ImageGrayU &mask, int thresh, int r, int g, int b, int xBorder, int yBorder ) {
    int width = output.width(), height = output.height();
    assertAlways( mask.width() == width && mask.height() == height );
    for (int y = yBorder; y < height - 1 - yBorder; y++) {
        for (int x = xBorder; x < width - 1 - xBorder; x++) {
            int v = mask.data( x, y ) >= thresh;
            int vx = mask.data( x + 1, y ) >= thresh;
            int vy = mask.data( x, y + 1 ) >= thresh;
             if (v != vx) {
                output.setRGB( x, y, r, g, b );
                output.setRGB( x + 1, y, r, g, b );
             }
             if (v != vy) {
                output.setRGB( x, y, r, g, b );
                output.setRGB( x, y + 1, r, g, b );
             }
        }
    }
}


/// join images horizontally
aptr<ImageColorU> joinHoriz( const ImageColorU &img1, const ImageColorU &img2 ) {
    int height = img1.height();
    assertAlways( height == img2.height() );
    int width1 = img1.width();
    int width2 = img2.width();
    int newWidth = width1 + width2;
    aptr<ImageColorU> result( new ImageColorU( newWidth, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width1; x++)
            for (int c = 0; c < 3; c++) 
                result->data( x, y, c ) = img1.data( x, y, c );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width2; x++)
            for (int c = 0; c < 3; c++) 
                result->data( x + width1, y, c ) = img2.data( x, y, c );
    return result;
}


/// join images horizontally
aptr<ImageGrayU> joinHoriz( const ImageGrayU &img1, const ImageGrayU &img2 ) {
    int height = img1.height();
    assertAlways( height == img2.height() );
    int width1 = img1.width();
    int width2 = img2.width();
    int newWidth = width1 + width2;
    aptr<ImageGrayU> result( new ImageGrayU( newWidth, height ) );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width1; x++)
            result->data( x, y ) = img1.data( x, y );
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width2; x++)
            result->data( x + width1, y ) = img2.data( x, y );
    return result;
}


/// blend two images, apply alpha weight to each image (alpha values typically in [0, 1])
// fix(faster): add opencv version
aptr<ImageColorU> blend( const ImageColorU &image1, const ImageColorU &image2, float alpha1, float alpha2 ) {
    int width = image1.width(), height = image1.height();
    assertAlways( image2.width() == width && image2.height() == height );
    aptr<ImageColorU> output( new ImageColorU( width, height ) );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                int v = round( image1.data( x, y, c ) * alpha1 + image2.data( x, y, c ) * alpha2 );
                if (v < 0) v = 0;
                if (v > 255) v = 255;
                output->data( x, y, c ) = v;
            }
        }
    }
    return output;
}


/// show left half of image 1 with right half of image 2
aptr<ImageColorU> splitScreen( const ImageColorU &image1, const ImageColorU &image2 ) {
    int width = image1.width(), height = image1.height();
    assertAlways( image2.width() == width && image2.height() == height );
    aptr<ImageColorU> output( new ImageColorU( width, height ) );
    int xMid = width / 2;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < xMid; x++) {
            for (int c = 0; c < 3; c++) {
                output->data( x, y, c ) = image1.data( x, y, c );
            }
        }
        for (int x = xMid; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                output->data( x, y, c ) = image2.data( x, y, c );
            }
        }
    }
    return output;
}


/// display an image (if using GUI)
void dispImage( const ImageColorU &img, const String &label ) {
#ifdef USE_GUI
    if (ImageSeqViewer::instance()) {
        ImageSeqViewer::instance()->dispImage( img, label );
    }
#endif
}


/// display an image (if using GUI)
void dispImage( const ImageGrayU &img, const String &label ) {
    aptr<ImageColorU> color = toColor( img );
    dispImage( *color, label );
}


//-------------------------------------------
// FILTER REGISTRY
//-------------------------------------------


// generic filter version of blurBox function
aptr<ImageColorU> blurBox( const ImageColorU &input, Config &conf ) {
    int boxSize = conf.readInt( "boxSize" );
    return blurBox( input, boxSize );
}


// register commands, etc. defined in this module
void initImageUtil() {
    registerFilter( blurBox );
}


} // end namespace sbl

