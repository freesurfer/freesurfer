#include "prep/BlockPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/TensorUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/image/Video.h>
#include <sbl/other/Plot.h>
using namespace sbl;
namespace hb {


/// takes a subset of the blockface images (to handle the case that multiple images were captured for each slice)
void selectBlockFaceImageSubset( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/extra" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/raw" ) );
	int imagesPerSlice = conf.readInt( "imagesPerSlice", 6 );
	int offset = conf.readInt( "offset", 5 );

	// create output dir
	createDir( outputPath );

	// loop over images
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	for (int inputIndex = offset; inputIndex < fileList.count(); inputIndex += imagesPerSlice) {
		disp( 1, "index: %d, file: %s", inputIndex, fileList[ inputIndex ].c_str() );

		// move file to output path
		moveFile( inputPath + "/" + fileList[ inputIndex ], outputPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


/// compute a rough mask of a block-face image
aptr<ImageGrayU> roughMask( const ImageColorU &image ) {

	// blur to remove outliers
	aptr<ImageColorU> blurImage = blurBox( image, 3 );
	int width = image.width(), height = image.height();

	// compute mask of candidate pixels
	int candCount = 0;
	aptr<ImageGrayU> mask( new ImageGrayU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int r = blurImage->r( x, y );
			int b = blurImage->b( x, y );
			if (r - b > 30) {
				mask->data( x, y ) = 255;
				candCount++;
			} else {
				mask->data( x, y ) = 0;
			}
		}
	}

	// clean up the mask
	mask = blurBoxAndThreshold( *mask, 5, 128 );

	// mark any component near the center
	int xMin = width / 2 - width / 8;
	int xMax = width / 2 + width / 8;
	int yMin = height / 2 - height / 8;
	int yMax = height / 2 + height / 8;
	for (int y = yMin; y <= yMax; y++) {
		for (int x = xMin; x <= xMax; x++) {
			if (mask->data( x, y ) == 255) {
				floodFill( *mask, 255, 255, 100, x, y );
			}
		}
	}

	// only keep marked components
	int finalCount = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask->data( x, y ) == 100) {
				mask->data( x, y ) = 255;
				finalCount++;
			} else if (mask->data( x, y )) {
				mask->data( x, y ) = 0;
			}
		}
	}
	disp( 1, "candCount: %d, finalCount: %d", candCount, finalCount );
	return mask;
}


/// estimate a transformation that registers the images before/after the microtome re-adjustment
aptr<ImageTransform> computeSplitTransform( const String &fileName1, const String &fileName2, const String &visPath ) {
	bool useInternal = true;

	// load boundary images
	aptr<ImageColorU> image1 = load<ImageColorU>( fileName1 );
	aptr<ImageColorU> image2 = load<ImageColorU>( fileName2 );

	// compute masks
	aptr<ImageGrayU> mask1 = roughMask( *image1 );
	aptr<ImageGrayU> mask2 = roughMask( *image2 );
	mask1 = blurBox( *mask1, 15 );
	mask2 = blurBox( *mask2, 15 );

	// perform initial registration using masks
	int paramCount = 6;
	int step = 2;
	int xBorder = 10;
	int yBorder = 10;
	float offsetBound = 200;
	aptr<ImageTransform> transform = registerUsingImageTransform( *mask1, *mask2, paramCount, step, xBorder, yBorder, offsetBound );
	disp( 1, "init transform:" );
	transform->display( 2 );

	// perform final registration using images inside masks
	if (useInternal) {
		aptr<ImageGrayU> regMask1 = threshold( *mask1, 250, false );
		aptr<ImageGrayU> regMask2 = threshold( *mask2, 250, false );
		aptr<ImageGrayU> regImage1 = blurBox( *toGray( *image1 ), 3 );
		aptr<ImageGrayU> regImage2 = blurBox( *toGray( *image2 ), 3 );
		transform = registerUsingImageTransform( *regImage1, *regImage2, paramCount, step, xBorder, yBorder, offsetBound, transform.get(), regMask1.get(), regMask2.get() );
		disp( 1, "final transform:" );
		transform->display( 2 );
	}

	// save diagnostic images
	int width = image1->width(), height = image1->height();
	aptr<ImageColorU> mapped = transform->mapForward( *image1, width, height, 200 );
	saveImage( *mapped, visPath + "/mapped.png" );
	saveImage( *image2, visPath + "/image.png" );
	saveImage( *mask1, visPath + "/mask1.png" );
	saveImage( *mask2, visPath + "/mask2.png" );
	return transform;
}


/// compute the approximate x-y-plane diameter of the MRI volume in pixels (assumes images are already segmented)
int computeMriSize( const String &mriPath ) {
	Array<String> fileList = dirFileList( mriPath, "", ".png" );
	int xMin = 100000, xMax = 0, yMin = 100000, yMax = 0;
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
		aptr<ImageGrayU> image = load<ImageGrayU>( mriPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (image->data( x, y )) {
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
	int xSize = xMax - xMin + 1;
	int ySize = yMax - yMin + 1;
	return (xSize + ySize) / 2;
}


/// crops, aligns, and resizes a set of block-face imges
void cropBlockFaceImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/raw" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/crop" ) );
	String mriPath = addDataPath( conf.readString( "mriPath", "mri/seg" ) );
	int visScale = 8;
	int pad = 100; // note: the pad will end up being larger than this because of the split transform

	// create output dir
	createDir( outputPath );

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// load first image to get dimensions
	aptr<ImageColorU> firstImage = load<ImageColorU>( inputPath + "/" + fileList[ 0 ] );
	int width = firstImage->width(), height = firstImage->height();
	firstImage.release();

	// open visualization video
	int visWidth = ((width / visScale) / 8) * 8; // make sure video size divisible by 8
	int visHeight = ((height / visScale) / 8) * 8;
	OutputVideo outputVideo( visPath + "/vis.avi", visWidth, visHeight );

	// compute MR bounds
	int mriSize = computeMriSize( mriPath );

	// first pass: compute bounds
	int xMin = 100000, xMax = 0, yMin = 100000, yMax = 0;
	VectorD countVect;
	int splitIndex = -1, maxDiff = 0;
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
	
		// load the input image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		assertAlways( image->width() == width && image->height() == height );

		// compute rough mask of tissue area
		aptr<ImageGrayU> mask = roughMask( *image );

		// find bounds of mask area
		int insideCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask->data( x, y ) == 255) {
					insideCount++;
					if (x < xMin) 
						xMin = x;
					if (x > xMax)
						xMax = x;
					if (y < yMin)
						yMin = y;
					if (y > yMax)
						yMax = y;
					image->setRGB( x, y, 255, 0, 0 );
				}
			}
		}
		
		// check for split
		int diff = 0;
		if (countVect.length())
			diff = iAbs( insideCount - (int) countVect.endValue() );
		if (diff > maxDiff) {
			maxDiff = diff;
			splitIndex = inputIndex - 1;
		}

		// store for counts
		disp( 1, "inputIndex: %d, file: %s, insideCount: %d, diff: %d", inputIndex, fileList[ inputIndex ].c_str(), insideCount, diff );
		countVect.append( (double) insideCount );

		// check for bad slice
		if (insideCount < 100000 || diff > 1000000) {
			warning( "bad slice: %s; stopping command", fileList[ inputIndex ].c_str() );
			return;
		}

		// create visualization image
		aptr<ImageColorU> visImage = resize( *image, visWidth, visHeight, true );
		outputVideo.append( *visImage );
//		saveImage( *visImage, visPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	disp( 1, "orig xMin: %d, xMax: %d, yMin: %d, yMax: %d", xMin, xMax, yMin, yMax );
	xMin = bound( xMin - pad, 0, width - 1 );
	xMax = bound( xMax + pad, 0, width - 1 );
	yMin = bound( yMin - pad, 0, width - 1 );
	yMax = bound( yMax + pad, 0, width - 1 );
	disp( 1, "crop xMin: %d, xMax: %d, yMin: %d, yMax: %d", xMin, xMax, yMin, yMax );

	// save plot of counts
	Plot plot;
	plot.add( countVect );
	plot.save( visPath + "/counts.svg" );

	// compute transformation
	String fileName1 = fileList[ splitIndex ];
	String fileName2 = fileList[ splitIndex + 1 ];
	disp( 1, "split index: %d, max diff: %d, file 1: %s, file 2: %s", splitIndex, maxDiff, fileName1.c_str(), fileName2.c_str() );
	aptr<ImageTransform> transform = computeSplitTransform( inputPath + "/" + fileName1, inputPath + "/" + fileName2, visPath );

	// compute final scale factor
	int xSize = xMax - xMin + 1;
	int ySize = yMax - yMin + 1;
	int blockSize = (xSize + ySize) / 2;
	int outputWidth = xSize * mriSize / blockSize;
	int outputHeight = ySize * mriSize / blockSize;
	outputWidth -= outputWidth & 7;
	outputHeight -= outputHeight & 7;
	disp( 1, "block size: %d, mri size: %d, output width: %d, output height: %d", blockSize, mriSize, outputWidth, outputHeight );

	// second pass: crop images
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
	
		// load the image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		disp( 1, "file: %s", fileList[ inputIndex ].c_str() );

		// if before split, apply transform
		if (inputIndex <= splitIndex) 
			image = transform->mapForward( *image, width, height, 255 );

		// crop and shrink the image
		image = crop( *image, xMin, xMax, yMin, yMax );
		image = resize( *image, outputWidth, outputHeight, true );

		// save to output path
		saveImage( *image, outputPath + "/" + fileList[ inputIndex ].leftOfLast( '.' ) + ".png" );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


// data type used for block-face histogram
typedef Tensor3F Histogram;


/// initialize a histogram with the specified number of bins in each dimension
void initHistogram( Histogram &histogram, int histogramSize ) {
	for (int i = 0; i < histogram.dimCount(); i++)
		histogram.setSize( i, histogramSize );
	histogram = 0;
}


/// compute a feature vector and transform to histogram indices
VectorI featureToIndices( int r, int g, int b, int histogramSize ) {
	VectorI featInd( 3 );
	featInd[ 0 ] = r * (histogramSize - 1) / 255;
	featInd[ 1 ] = g * (histogramSize - 1) / 255;
	featInd[ 2 ] = b * (histogramSize - 1) / 255;
	return featInd;
}


/// segment a set of block-face images
void segmentBlockFaceImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/crop" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/seg" ) );
	int histogramSize = conf.readInt( "histogramSize", 30 );
	int histogramBlurSize = conf.readInt( "histogramBlurSize", 5 );

	// create output dir
	createDir( outputPath );

	// if already run, stop here
//	if (dirFileList( outputPath, "", ".png" ).count()) {
//		disp( 1, "output path already has files; stopping" );
//		return;
//	}

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// initialize histograms for foreground and background distributions
	Histogram fgHistogramRaw, bgHistogramRaw;
	initHistogram( fgHistogramRaw, histogramSize );
	initHistogram( bgHistogramRaw, histogramSize );
	int fgCount = 0, bgCount = 0;

	// we'll skip over some pixels and some images for the sake of efficiency
	int pixelStep = 3;
	int imageStep = 3;

	// first pass: build models
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex += imageStep) {

		// load image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();

		// compute rough mask of tissue area
		aptr<ImageGrayU> mask = roughMask( *image );

		// blur the mask to obtain inside, outside, and uncertain areas
		mask = blurBox( *mask, 35 );

		// update histograms
		for (int y = 0; y < height; y += pixelStep) {
			for (int x = 0; x < width; x += pixelStep) {
				int r = image->r( x, y );
				int g = image->g( x, y );
				int b = image->b( x, y );
				VectorI featInd = featureToIndices( r, g, b, histogramSize );
				int m = mask->data( x, y );
				if (m == 0) {
					bgHistogramRaw.elem( featInd.dataPtr() )++;
					bgCount++;
				} else if (m == 255) {
					fgHistogramRaw.elem( featInd.dataPtr() )++;
					fgCount++;
				}
			}
		}
		status( "inputIndex: %d, fgCount: %d, bgCount: %d", inputIndex, fgCount, bgCount );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );
	disp( 1, "fgCount: %d, bgCount: %d", fgCount, bgCount );

	// create blurred histograms
	Histogram fgHistogram, bgHistogram;
	initHistogram( fgHistogram, histogramSize );
	initHistogram( bgHistogram, histogramSize );
	blurBox( fgHistogramRaw, fgHistogram, histogramBlurSize );
	blurBox( bgHistogramRaw, bgHistogram, histogramBlurSize );

	// normalize the histograms
	float fgNorm = 1.0f / fgCount;
	float bgNorm = 1.0f / bgCount;
	fgHistogram *= fgNorm;
	bgHistogram *= bgNorm;

	// second pass: perform (soft) segmentation
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {

		// load image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );
		int width = image->width(), height = image->height();

		// init mask
		aptr<ImageGrayU> mask( new ImageGrayU( width, height ) );

		// compute soft mask value using histograms
		int fgOutCount = 0, bgOutCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int r = image->r( x, y );
				int g = image->g( x, y );
				int b = image->b( x, y );
				VectorI featInd = featureToIndices( r, g, b, histogramSize );
				float fgVal = fgHistogram.elem( featInd.dataPtr() );
				float bgVal = bgHistogram.elem( featInd.dataPtr() );
				if (fgVal > 10 * bgVal) {
					mask->data( x, y ) = 255;
					fgOutCount++;
				} else {
					mask->data( x, y ) = 0;
					bgOutCount++;
				}
			}
		}
		status( "inputIndex: %d, fgCount: %d, bgCount: %d", inputIndex, fgOutCount, bgOutCount );

		// save the mask
		saveImage( *mask, outputPath + "/" + fileList[ inputIndex ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );
}


/// perform all steps needed to prepare a set of block-face images for registration
void prepareBlockFaceImages( Config &conf ) {
	cropBlockFaceImages( conf );
	execCommand( "vcross blockface/crop blockface/crop/vis/cross", false );
	segmentBlockFaceImages( conf );
	execCommand( "vcross blockface/seg blockface/seg/vis/cross", false );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initBlockPrep() {
	registerCommand( "bsubset", selectBlockFaceImageSubset );
	registerCommand( "bcrop", cropBlockFaceImages ); 
	registerCommand( "bseg", segmentBlockFaceImages ); 
	registerCommand( "bprep", prepareBlockFaceImages ); 
}


} // end namespace hb
