#include "prep/HistoPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/other/Plot.h>
using namespace sbl;
namespace hb {


/// normalize intensity and mark background
void prepareHisto( ImageGrayU &image, int botThresh, int topThresh, int peakStart ) {
	int width = image.width(), height = image.height();
	int border = 100;

	// get smoothed image histogram
	if (botThresh == -1 || topThresh == -1) {
		VectorF histogram = toFloat( imageHistogram( image, border, width - border, border, height - border ) );
		histogram = gaussFilter( histogram, 3.0f );
		int peak = nearestMaxIndex( histogram, peakStart );
		topThresh = nearestMinIndex( histogram, peak - 3 );
		float sum = histogram.sum();
		float runSum = 0;
		botThresh = 0;
		for (int i = 0; i < 255; i++) {
			runSum += histogram[ i ];
			if (runSum > sum * 0.001) {
				botThresh = i;
				break;
			}
		}
		simplePlot( toDouble( histogram ) )->save( dataPath() + "histogram.svg" );
		disp( 1, "botThresh: %d (%f), topThresh: %d (%f), peak: %d (%f)", 
				botThresh, histogram[ botThresh ], topThresh, histogram[ topThresh ], peak, histogram[ peak ] );
	}

	// create mask
	int maskScale = 4;
	int maskWidth = width / maskScale;
	int maskHeight = height / maskScale;
	aptr<ImageGrayU> mask( new ImageGrayU( maskWidth, maskHeight ) );
	for (int y = 0; y < maskHeight; y++) {
		for (int x = 0; x < maskWidth; x++) {
			if (image.data( x * maskScale, y * maskScale ) <= topThresh) {
				mask->data( x, y ) = 255;
			} else {
				mask->data( x, y ) = 0;
			}
		}
	}

	// clean mask
	filterMaskComponents( *mask, 20 * 20, width * height / maskScale / maskScale );
	fillMaskHoles( *mask, 20 * 20 );
	aptr<ImageGrayU> maskBig = resize( *mask, width, height, true );
	maskBig = blurBoxAndThreshold( *maskBig, 5, 127 );
	saveImage( *mask, dataPath() + "mask.png" );

	// apply mapping
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (maskBig->data( x, y )) {
				int v = image.data( x, y );
				v = (v - botThresh) * 240 / (topThresh - botThresh);
				image.data( x, y ) = bound( v, 0, 240 );
			} else {
				image.data( x, y ) = 255;
			}
		}
	}
}


/// generate a test image for testing the fiber orientation estimation algorithm
void generateFiberTestImage( Config &conf ) {

	// get command parameters
	String outputFileName = addDataPath( conf.readString( "outputFileName" ) );
	int width = conf.readInt( "width", 6000 );
	int height = conf.readInt( "height", 6000 );
	int radius = conf.readInt( "radius", 2500 );

	// we will generate image at double resolution so we have less aliasing
	width *= 2;
	height *= 2;
	radius *= 2;

	// create white image with gray circle
	aptr<ImageGrayU> image( new ImageGrayU( width, height ) );
	image->clear( 255 );
	int xCent = width / 2;
	int yCent = height / 2;
	drawCircleFilled( *image, xCent, yCent, radius, 200 );

	// draw black lines at various angles
	for (float angle = 0; angle <= 2.0f * 3.13f; angle += 0.02f) {
		int xEnd = xCent + round( (float) radius * cosf( angle ) );
		int yEnd = yCent + round( (float) radius * sinf( angle ) );
		drawLine( *image, xCent, yCent, xEnd, yEnd, 0, true );
	}

	// blur and shrink the image
	image = blurGauss( *image, 2.0f );
	image = resize( *image, width / 2, height / 2, true );

	// save result
	saveImage( *image, outputFileName );
}


/// estimate local fiber orientation across the image
void computeFiberOrientation( Config &conf ) {

	// get command parameters
	String inputFileName = addDataPath( conf.readString( "inputFileName" ) );
	String outputPrefix = addDataPath( conf.readString( "outputPrefix" ) );
	float blurSigma = conf.readFloat( "blurSigma", 10.0f );

	// load input image
	aptr<ImageGrayU> inputImage = load<ImageGrayU>( inputFileName );
	int width = inputImage->width(), height = inputImage->height();

	// normalize intensity and mark background
	prepareHisto( *inputImage, -1, -1, 255 );

	// compute gradients
	aptr<ImageGrayF> xGradient = xGrad( *inputImage, 3 );
	aptr<ImageGrayF> yGradient = yGrad( *inputImage, 3 );

	// compute eroded mask
	aptr<ImageGrayU> mask = threshold( *inputImage, 250, true );
	mask = blurBoxAndThreshold( *mask, 7, 250 );

	// remove exterior boundary gradients
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask->data( x, y ) == 0) {
				xGradient->data( x, y ) = 0;
				yGradient->data( x, y ) = 0;
			}
		}
	}
	mask.release();

	// shrink everything
	width = width / 2;
	height = height / 2;
	xGradient = resize( *xGradient, width, height, true );
	yGradient = resize( *yGradient, width, height, true );
	inputImage = resize( *inputImage, width, height, true );

	// compute transform into 180-degree normalized psuedo-gradients
	aptr<ImageGrayF> px( new ImageGrayF( width, height ) );
	aptr<ImageGrayF> py( new ImageGrayF( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			float gxi = xGradient->data( x, y );
			float gyi = yGradient->data( x, y );
			float lenSqd = gxi * gxi + gyi * gyi;
			if (lenSqd > 0.01f * 0.01f) {
				float a = atan2( gyi, gxi );
				if (a < 0)
					a += 3.14159265f;
				a *= 2.0f;
				float len = sqrtf( lenSqd );
				px->data( x, y ) = cos( a ) * len;
				py->data( x, y ) = sin( a ) * len;
			} else {
				px->data( x, y ) = 0;
				py->data( x, y ) = 0;
			}
		}
	}

	// shrink again
	width = width / 2;
	height = height / 2;
	px = resize( *px, width, height, true );
	py = resize( *py, width, height, true );
	inputImage = resize( *inputImage, width, height, true );

	// blur the psuedo-gradients
	px = blurGauss( *px, blurSigma );
	py = blurGauss( *py, blurSigma );
	float pxMin = 0, pxMean = 0, pxMax = 0;
	imageStats( *px, pxMin, pxMean, pxMax );
	float pyMin = 0, pyMean = 0, pyMax = 0;
	imageStats( *py, pyMin, pyMean, pyMax );
	float valBound = max( max( -pxMin, pxMax ), max( -pyMin, pyMax ) );
	float valScale = 128.0f / (valBound - 5);

	// visualize
	aptr<ImageColorU> vis( new ImageColorU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (inputImage->data( x, y ) > 250) {
				vis->setRGB( x, y, 255, 255, 255 );
			} else {
				int r = 128;
				int g = 128 + round( px->data( x, y ) * valScale );
				int b = 128 + round( py->data( x, y ) * valScale );
				g = bound( g, 0, 255 );
				b = bound( b, 0, 255 );
				vis->setRGB( x, y, r, g, b );
			}
		}
	}

	// save visualization images
	aptr<ImageGrayU> xGradientVis = toUChar( *xGradient ); 
	aptr<ImageGrayU> yGradientVis = toUChar( *yGradient ); 
	saveImage( *xGradient, outputPrefix + ".xGrad.png" );
	saveImage( *yGradient, outputPrefix + ".yGrad.png" );
	saveImage( *vis, outputPrefix + ".vis.png" );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initHistoStats() {
	registerCommand( "hfiber", computeFiberOrientation );
	registerCommand( "hfibergen", generateFiberTestImage );
}


} // end namespace hb
