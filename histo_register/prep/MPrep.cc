#include "prep/MPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/TensorUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/other/Plot.h>
#include "prep/VolumeFile.h"
using namespace sbl;
namespace hb {


/// convert MRI volume files into sequences of images
void convertMData( Config &conf ) {
	conf.writeString( "fileName", "mri/flash20.mgz" );
	conf.writeString( "outputPath", "mri/rawFlash20" );
	convertMghFileToImages( conf );
	if (fileExists( dataPath() + "mri/T1.mgz" )) {
		conf.writeString( "fileName", "mri/T1.mgz" );
		conf.writeString( "outputPath", "mri/rawT1" );
		convertMghFileToImages( conf );
	}
}


/// threshold MRI data to obtain a foreground/background segmentation to be used for registration
void thresholdMData( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "mri/rawFlash20" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "mri/seg" ) );
	int thresh = conf.readInt( "thresh", 15 );

	// create output path
	createDir( outputPath );

	// loop over input images
	Array<String> fileList = dirFileList( inputPath, "", "png" );
	for (int i = 0; i < fileList.count(); i++) {

		// load input image
		String inputFileName = inputPath + "/" + fileList[ i ];
		disp( 1, "input: %s", inputFileName.c_str() );
		aptr<ImageGrayU> input = load<ImageGrayU>( inputFileName );
		int width = input->width(), height = input->height();

		// compute and process mask
		aptr<ImageGrayU> mask = threshold( *input, (float) thresh, false );
		filterMaskComponents( *mask, 20 * 20, width * height );
//		fillMaskHoles( *mask, 20 * 20 );

		// update input image based on mask
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask->data( x, y ) == 0)
					input->data( x, y ) = 0;
				else if (input->data( x, y ) < 5)
					input->data( x, y ) = 5;
			}
		}

		// save results
		saveImage( *input, outputPath + "/" + fileList[ i ] );

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
}


struct LabelPoint {
	LabelPoint( int xNew, int yNew, int zNew, int vNew ) : x( xNew ), y( yNew ), z( zNew ), v( vNew ) {}
	int x;
	int y;
	int z;
	int v;
};


int interpolate( int x, int y, int z, const Array<LabelPoint> &labels, float sigma ) {
	float factor = gaussFactor( sigma );
	float maxDist = sigma * 5.0f;
	float maxDistSqd = maxDist * maxDist;
	float vWtSum = 0, wtSum = 0;
	for (int i = 0; i < labels.count(); i++) {
		int xDiff = x - labels[ i ].x;
		int yDiff = y - labels[ i ].y;
		int zDiff = z - labels[ i ].z;
		float dSqd = (float) (xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
		if (dSqd < maxDistSqd) {
			float wt = gauss( dSqd, factor );
//			float wt = 1.0f - sqrtf( dSqd ) / maxDist; // compute weight from distance (linear)
			vWtSum += wt * labels[ i ].v;
			wtSum += wt;
		}
	}
	float vInterp = 128;
	if (wtSum)
		vInterp = vWtSum / wtSum;
	return round( vInterp );
}


/// segment MR data using PD/T1/flash values
void segmentMData( Config &conf ) {

	// get command parameters
	String inputPath1 = addDataPath( conf.readString( "inputPath1", "mri/rawFlash20" ) );
	String inputPath2 = addDataPath( conf.readString( "inputPath2", "mri/rawT1" ) );
	String labelPath = addDataPath( conf.readString( "labelPath", "mri/segLabel" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "mri/seg" ) );
	int blockSize = conf.readInt( "blockSize", 20 );
	float sigma = conf.readFloat( "sigma", 30 );

	// create output path
	createDir( outputPath );

	// prepare path for visualization output images
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// labelled points
	Array<LabelPoint> fgLabels;
	Array<LabelPoint> bgLabels;

	// load labelled images
	Array<String> labelFileList = dirFileList( labelPath, "", "png" );
	for (int i = 0; i < labelFileList.count(); i++) {

		// load input images
		String labelFileName = labelPath + "/" + labelFileList[ i ];
		String inputFileName = inputPath2 + "/" + labelFileList[ i ];
		aptr<ImageColorU> label = load<ImageColorU>( labelFileName );
		aptr<ImageGrayU> input2 = load<ImageGrayU>( inputFileName );
		if (label.get() == NULL) {
			warning( "unable to load: %s", labelFileName.c_str() );
			return;
		}
		if (input2.get() == NULL) {
			warning( "unable to load: %s", inputFileName.c_str() );
			return;
		}
		int width = label->width(), height = label->height();
		assertAlways( input2->width() == width && input2->height() == height );

		// blur the T1 image slightly
		input2 = blurGauss( *input2, 1.5 );

		// get z index from name
		int z = labelFileList[ i ].leftOfLast( '.' ).toInt();

		// find label points
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int r = label->r( x, y );
				int g = label->g( x, y );
				int b = label->b( x, y );
				if (r > 200 && g < 50 && b < 50) {
					int v = input2->data( x, y );
					fgLabels.appendCopy( LabelPoint( x, y, z, v ) );	
				} else if (b > 200 && g < 50 && r < 50) {
					int v = input2->data( x, y );
					bgLabels.appendCopy( LabelPoint( x, y, z, v ) );
				}
			}
		}
	}
	disp( 1, "fg label points: %d, bg label points: %d", fgLabels.count(), bgLabels.count() );
	if (fgLabels.count() == 0 || bgLabels.count() == 0) {
		warning( "requires foreground and background label points" );
		return;
	}

	// get list of input files
	Array<String> fileList = dirFileList( inputPath1, "", "png" );

	// store info for each slice
	VectorD insideCountVect, fgCountVect;

	// segment each image
	for (int z = 0; z < fileList.count(); z++) {

		// load input images
		aptr<ImageGrayU> input1 = load<ImageGrayU>( inputPath1 + "/" + fileList[ z ] );
		aptr<ImageGrayU> input2 = load<ImageGrayU>( inputPath2 + "/" + fileList[ z ] );
		assertAlways( input1.get() && input2.get() );
		int width = input1->width(), height = input1->height();
		assertAlways( input2->width() == width && input2->height() == height );

		// blur the T1 image slightly
		input2 = blurGauss( *input2, 1.5 );

		// compute a mask of the sample regions using the flash image
		aptr<ImageGrayU> segArea = threshold( *input1, 20, false );

		// erode the mask a bit
		segArea = blurBoxAndThreshold( *segArea, 7, 250 );

		// this will hold the computed mask 
		ImageGrayU mask( width, height );
		mask.clear( 0 );
		int insideCount = 0;
		VectorF thresholdVect;

		ImageGrayU fgValue( width, height );
		ImageGrayU bgValue( width, height );

		// loop over blocks, compute the mask for each block
		for (int yBlock = 0; yBlock < height; yBlock += blockSize) {
			for (int xBlock = 0; xBlock < width; xBlock += blockSize) {

				// get threshold from nearby label points
				int xCent = xBlock + blockSize / 2;
				int yCent = yBlock + blockSize / 2;
				int fgMean = interpolate( xCent, yCent, z, fgLabels, sigma );
				int bgMean = interpolate( xCent, yCent, z, bgLabels, sigma );
				int thresh = (fgMean + bgMean) / 2;
				thresholdVect.append( (float) thresh );

				// compute bounds of block
				int xMin = xBlock;
				int xMax = xBlock + blockSize;
				if (xMax > width - 1)
					xMax = width - 1;
				int yMin = yBlock;
				int yMax = yBlock + blockSize;
				if (yMax > height - 1)
					yMax = height - 1;

				// loop over the block computing mask for each pixel
				for (int y = yMin; y <= yMax; y++) {
					for (int x = xMin; x <= xMax; x++) {
						if (segArea->data( x, y )) {
							if (input2->data( x, y ) < thresh) {
								mask.data( x, y ) = 255;
							}
						}
						bgValue.data( x, y ) = bgMean;
						fgValue.data( x, y ) = fgMean;
					}
				}
			}
		}

		// clean mask
		filterMaskComponents( mask, 5 * 5, width * height );

		// apply mask: set bg to zero and make sure fg >= 5
		int fgCount = 0;
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				if (mask.data( x, y ) == 0) {
					input1->data( x, y ) = 0;
				} else {
					fgCount++;
					if (input1->data( x, y ) < 5)
						input1->data( x, y ) = 5;
				}
			}
		}
		disp( 1, "input: %s, threshold: %f/%f/%f, fg: %d", fileList[ z ].c_str(), thresholdVect.min(), thresholdVect.mean(), thresholdVect.max(), fgCount );

		// save results
		saveImage( *input1, outputPath + "/" + fileList[ z ] );

		// save diagnostic images
		saveImage( *input2, visPath + "/" + fileList[ z ] );
		saveImage( bgValue, visPath + "/" + fileList[ z ].leftOfLast( '.' ) + "_bg.png" );
		saveImage( fgValue, visPath + "/" + fileList[ z ].leftOfLast( '.' ) + "_fg.png" );

		// store counts for this slice
		insideCountVect.append( (double) insideCount );
		fgCountVect.append( (double) fgCount );

		// check for user cancel
		if (checkCommandEvents())
			return;
	}

	// save plots
	Plot areaPlot;
	areaPlot.setColor( 0, 0, 0 );
	areaPlot.add( insideCountVect );
	areaPlot.setColor( 0, 200, 0 );
	areaPlot.add( fgCountVect );
	areaPlot.save( visPath + "/area.svg" );
}


/// perform all steps needed to prepare MRI data for registration
void prepareMData( Config &conf ) {
	String labelPath = addDataPath( conf.readString( "labelPath", "mri/segLabel" ) );
	Array<String> labelFileList = dirFileList( labelPath, "", "png" );
	if (labelFileList.count()) {
		disp( 1, "found label files; fg/bg segmenting MR data" );
		segmentMData( conf );
	} else {
		disp( 1, "did not find label files; thresholding MR data" );
		thresholdMData( conf );
	}
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initMPrep() {
	registerCommand( "mconv", convertMData );
	registerCommand( "mthresh", thresholdMData );
	registerCommand( "mseg", segmentMData );
	registerCommand( "mprep", prepareMData );
}


} // end namespace hb
