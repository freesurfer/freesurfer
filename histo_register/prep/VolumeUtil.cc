#include "prep/VolumeFile.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
using namespace sbl;
namespace hb {


/// generate x and y cross sections from a volume represented as a sequence of images
void generateVolumeCrossSections( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputPrefix = addDataPath( conf.readString( "outputPrefix" ) );
	int zFactor = conf.readInt( "zFactor", 1 );

	// get input file list
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	if (fileList.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}

	// load first image to get dimensions
	aptr<ImageColorU> firstImage = load<ImageColorU>( inputPath + "/" + fileList[ 0 ] );
	int width = firstImage->width(), height = firstImage->height();
	firstImage.release();

	// allocate output images
	int depth = fileList.count();
	ImageColorU xCrossImage( depth * zFactor, height );
	ImageColorU yCrossImage( width, depth * zFactor );
	disp( 1, "width: %d, height: %d, depth: %d", width, height, depth );

	// compute output file names
	String xFileName = outputPrefix + ".x.png";
	String yFileName = outputPrefix + ".y.png";
	disp( 1, "xFileName: %s", xFileName.c_str() );
	disp( 1, "yFileName: %s", yFileName.c_str() );

	// loop over input images, filling output images
	for (int inputIndex = 0; inputIndex < fileList.count(); inputIndex++) {
		status( "%d", inputIndex );

		// load the input image
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ inputIndex ] );

		// add to slices
		int xCent = width / 2, yCent = height / 2;
		for (int x = 0; x < width; x++)
			for (int j = 0; j < zFactor; j++)
				for (int c = 0; c < 3; c++)
					yCrossImage.data( x, inputIndex * zFactor + j, c ) = image->data( x, yCent, c );
		for (int y = 0; y < height; y++)
			for (int j = 0; j < zFactor; j++)
				for (int c = 0; c < 3; c++)
					xCrossImage.data( inputIndex * zFactor + j, y, c ) = image->data( xCent, y, c );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );

	// save cross-section images
	saveImage( xCrossImage, xFileName );
	saveImage( yCrossImage, yFileName );
}


/// rotate a set of images 
void rotateVolume( Config &conf ) {

	// get command parameters
	// note: ok if inputPath is same as outputPath
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputPath = addDataPath( conf.readString( "outputPath" ) );
	int rotateAngleDegrees = conf.readInt( "rotateAngleDegrees" );
	int fillValue = conf.readInt( "fillValue", -1 );
	assertAlways( rotateAngleDegrees >= -180 && rotateAngleDegrees <= 180 );

	// make sure output path exists
	createDir( outputPath );

	// loop over input files
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
  if (fileList.count() == 0)
	  fileList = dirFileList( inputPath, "", ".PNG" );
  if (fileList.count() == 0)
	  fileList = dirFileList( inputPath, "", ".jpg" );
  if (fileList.count() == 0)
	  fileList = dirFileList( inputPath, "", ".JPG" );
    
	for (int i = 0; i < fileList.count(); i++) {

		// load input file
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ i ] );

		// rotate as requested
		if (rotateAngleDegrees == 90) {
			image = rotate90( *image );
		} else if (rotateAngleDegrees == -90) {
			image = rotate270( *image );
		} else if (rotateAngleDegrees == 180 || rotateAngleDegrees == -180) {
			image = rotate180( *image );
		} else {

			// if unspecified, get fill value from image corner
			if (fillValue == -1) 
				fillValue = image->data( 0, 0, 0 );

			// rotate by specified angle
			image = rotate( *image, (float) rotateAngleDegrees, fillValue );
		}

		// save output file
		saveImage( *image, outputPath + "/" + fileList[ i ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}

/// flop a set of images (horizontal flip)
void flopVolume( Config &conf ) {

	// get command parameters
	// note: ok if inputPath is same as outputPath
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputPath = addDataPath( conf.readString( "outputPath" ) );

	// make sure output path exists
	createDir( outputPath );

	// loop over input files
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	for (int i = 0; i < fileList.count(); i++) {

		// load input file
		aptr<ImageColorU> image = load<ImageColorU>( inputPath + "/" + fileList[ i ] );

	  // rotate by specified angle
		image = flipHoriz( *image );
		
		// save output file
		saveImage( *image, outputPath + "/" + fileList[ i ] );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}

//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initVolumeUtil() {
	registerCommand( "vcross", generateVolumeCrossSections );
	registerCommand( "vrotate", rotateVolume );
	registerCommand( "vflop", flopVolume );
}


} // end namespace hb
