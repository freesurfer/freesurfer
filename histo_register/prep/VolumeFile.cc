#include "prep/VolumeFile.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/ImageSeqUtil.h>
using namespace sbl;
namespace hb {


/// read a 32-bit integer from the file, swapping endianness
int readIntSwap( File &file ) {
	char data[ 4 ];
	file.readBlock( data, 4 );
	char swapData[] = { data[ 3 ], data[ 2 ], data[ 1 ], data[ 0 ] };
	int *swapDataInt = (int *) swapData;
	return *swapDataInt;
}


/// read a 32-bit float from the file, swapping endianness
float readFloatSwap( File &file ) {
	char data[ 4 ];
	file.readBlock( data, 4 );
	char swapData[] = { data[ 3 ], data[ 2 ], data[ 1 ], data[ 0 ] };
	float *swapDataFloat = (float *) swapData;
	return *swapDataFloat;
}


/// read a 16-bit int from the file, swapping endianness
short readShortSwap( File &file ) {
	char data[ 2 ];
	file.readBlock( data, 2 );
	char swapData[] = { data[ 1 ], data[ 0 ] };
	short *swapDataShort = (short *) swapData;
	return *swapDataShort;
}


/// write a 32-bit integer to the file, swapping endianness
void writeIntSwap( File &file, int val ) {
	char data[ 4 ];
	int *dataInt = (int *) data;
	*dataInt = val;
	char swapData[] = { data[ 3 ], data[ 2 ], data[ 1 ], data[ 0 ] };
	file.writeBlock( swapData, 4 );
}


/// write a 32-bit float to the file, swapping endianness
void writeFloatSwap( File &file, float val ) {
	char data[ 4 ];
	float *dataFloat = (float *) data;
	*dataFloat = val;
	char swapData[] = { data[ 3 ], data[ 2 ], data[ 1 ], data[ 0 ] };
	file.writeBlock( swapData, 4 );
}


/// write a 16-bit int to the file, swapping endianness
void writeShortSwap( File &file, short val ) {
	char data[ 2 ];
	short *dataShort = (short *) data;
	*dataShort = val;
	char swapData[] = { data[ 1 ], data[ 0 ] };
	file.writeBlock( swapData, 2 );
}


/// compute and save an axis flip/swap transformation matrix
void saveAxisTransform( const String &outputPath, bool xySwap, bool xzSwap, bool yzSwap, bool xFlip, bool yFlip, bool zFlip ) {
	int x[] = {1, 0, 0};
	int y[] = {0, 1, 0};
	int z[] = {0, 0, 1};
	if (xySwap) { 
		for (int i = 0; i < 3; i++) 
			swap( x[ i ], y[ i ] ); 
	}
	if (xzSwap) { 
		for (int i = 0; i < 3; i++) 
			swap( x[ i ], z[ i ] ); 
	}
	if (yzSwap) { 
		for (int i = 0; i < 3; i++) 
			swap( y[ i ], z[ i ] ); 
	}
	if (xFlip) {
		for (int i = 0; i < 3; i++) 
			x[ i ] = -x[ i ];
	}
	if (yFlip) {
		for (int i = 0; i < 3; i++) 
			y[ i ] = -y[ i ];
	}
	if (zFlip) {
		for (int i = 0; i < 3; i++) 
			z[ i ] = -z[ i ];
	}
	File axisFile( outputPath + "/volumeTransform.txt", FILE_WRITE, FILE_TEXT );
	if (axisFile.openSuccess()) {
		axisFile.writeF( "0, 0, 0, " ); // write zero translation offset so that this is a valid AffineTransform3 file
		axisFile.writeF( "%d, %d, %d, ", x[ 0 ], x[ 1 ], x[ 2 ] );
		axisFile.writeF( "%d, %d, %d, ", y[ 0 ], y[ 1 ], y[ 2 ] );
		axisFile.writeF( "%d, %d, %d\n", z[ 0 ], z[ 1 ], z[ 2 ] );
	} else {
		warning( "unable to write to output path: %s", outputPath.c_str() );
		return;
	}
}


/// convert an mgh/mgz file to a set of images
void convertMghFileToImages( const String &fileName, const String &outputPath, 
							 bool xySwap, bool xzSwap, bool yzSwap, 
							 bool xFlip, bool yFlip, bool zFlip, bool autoScaleValues ) {

	// open the file
	File file( fileName, FILE_READ, fileName.endsWith( "gz" ) ? FILE_GZIP_BINARY : FILE_BINARY );
	if (file.openSuccess() == false) {
		warning( "unable to open input file: %s", fileName.c_str() );
		return;
	}

	// make sure output path exists
	createDir( outputPath );

	// read header
	int version = readIntSwap( file );
	int width = readIntSwap( file );
	int height = readIntSwap( file );
	int depth = readIntSwap( file );
	int nFrames = readIntSwap( file );
	int type = readIntSwap( file );
	int dof = readIntSwap( file );
	int goodRASFlag = readShortSwap( file );
	float xSpacing = readFloatSwap( file );
	float ySpacing = readFloatSwap( file );
	float zSpacing = readFloatSwap( file );
	float xr = readFloatSwap( file );
	float xa = readFloatSwap( file );
	float xs = readFloatSwap( file );
	float yr = readFloatSwap( file );
	float ya = readFloatSwap( file );
	float ys = readFloatSwap( file );
	float zr = readFloatSwap( file );
	float za = readFloatSwap( file );
	float zs = readFloatSwap( file );
	float cr = readFloatSwap( file );
	float ca = readFloatSwap( file );
	float cs = readFloatSwap( file );

	// display header
	disp( 1, "filename: %s", fileName.c_str() );
	disp( 1, "version: %d", version );
	disp( 1, "width: %d", width );
	disp( 1, "height: %d", height );
	disp( 1, "depth: %d", depth );
	disp( 1, "nFrames: %d", nFrames );
	disp( 1, "type: %d", type );
	disp( 1, "dof: %d", dof );
	disp( 1, "goodRASFlag: %d", goodRASFlag );
	disp( 1, "spacing: %9.8f / %9.8f / %9.8f", xSpacing, ySpacing, zSpacing );
	disp( 1, "x: %9.8f / %9.8f / %9.8f", xr, xa, xs );
	disp( 1, "y: %9.8f / %9.8f / %9.8f", yr, ya, ys );
	disp( 1, "z: %9.8f / %9.8f / %9.8f", zr, za, zs );
	disp( 1, "c: %9.8f / %9.8f / %9.8f", cr, ca, cs );

	// save spacing
	File spacingFile( outputPath + "/spacing.txt", FILE_WRITE, FILE_TEXT );
	if (spacingFile.openSuccess()) {
		spacingFile.writeF( "%9.8f, %9.8f, %9.8f\n", xSpacing, ySpacing, zSpacing );
	} else {
		warning( "unable to write to output path: %s", outputPath.c_str() );
		return;
	}

	// save flip/swap params
	saveAxisTransform( outputPath, xySwap, xzSwap, yzSwap, xFlip, yFlip, zFlip );

	// sanity checks
	assertAlways( type == 1 || type == 3 || type == 4 );
	assertAlways( nFrames == 1 );

	// seek to data
	file.seek( 284, false );

	// compute output dimensions
	int xOutSize = width;
	int yOutSize = height;
	int zOutSize = depth;
	if (xySwap) swap( xOutSize, yOutSize );
	if (xzSwap) swap( xOutSize, zOutSize );
	if (yzSwap) swap( yOutSize, zOutSize );
	disp( 1, "output size: %d, %d, %d", xOutSize, yOutSize, zOutSize );

	// first pass: get value bounds
	float vMin = 0, vMax = 0;
	if (autoScaleValues) {
		disp( 1, "first pass" );
		VectorF vSample;
		int index = 0;
		for (int z = 0; z < depth; z++) {
			status( "        %d of %d", z + 1, depth );
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {

					// read a value
					float v = 0;
					if (type == 1) {
						v = (float) readIntSwap( file );
					} else if (type == 4) {
						v = (float) readShortSwap( file );
					} else {
						v = readFloatSwap( file );
					}

					// store a subset of the vluaes
					index++;
					if ((index % 37) == 0)
						vSample.append( v );
				}
			}

			// check for user cancel
			if (checkCommandEvents())
				break;
		}
		status( "\n" );

		// check for user cancel
		if (checkCommandEvents())
			return;

		// get bounds using percentiles
		int len = vSample.length();
		VectorI sortInd = sortIndex( vSample );
		vMin = vSample[ sortInd[ round( len * 0.001f ) ] ];
		vMax = vSample[ sortInd[ round( len * 0.999f ) ] ];
	//	float vMin = vSample.min();
	//	float vMax = vSample.max();
		disp( 1, "sample count: %d, vMin: %f, vMax: %f", len, vMin, vMax );
	} else {
		disp( 1, "auto scale values disabled; skipping first pass" );
	}

	// rewind to data
	file.seek( 284, false );

	// second pass: read the data into memory
	disp( 1, "second pass" );
	Array<ImageGrayU> imageSeq;
	initImageSeq( imageSeq, xOutSize, yOutSize, zOutSize, false, 0 );
	for (int z = 0; z < depth; z++) {
		status( "        %d of %d", z + 1, depth );
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int vInt = 0;

				// read and apply scaling determined in first pass
				if (autoScaleValues) {
					float v = 0;
					if (type == 1) {
						v = (float) readIntSwap( file );
					} else if (type == 4) {
						v = (float) readShortSwap( file );
					} else {
						v = readFloatSwap( file );
					}
					v = (v - vMin) / (vMax - vMin) * 255.0f;
					vInt = bound( round( v ), 0, 255 );
				} else {
					if (type == 1) {
						vInt = readIntSwap( file );
					} else if (type == 4) {
						vInt = readShortSwap( file );
					} else {
						vInt = sbl::round( readFloatSwap( file ) );
					}
					vInt = bound( vInt, 0, 255 );
				}

				// compute output coordinates
				int xOut = x;
				int yOut = y;
				int zOut = z;
				if (xySwap) swap( xOut, yOut );
				if (xzSwap) swap( xOut, zOut );
				if (yzSwap) swap( yOut, zOut );
				if (xFlip) xOut = xOutSize - xOut - 1;
				if (yFlip) yOut = yOutSize - yOut - 1;
				if (zFlip) zOut = zOutSize - zOut - 1;

				// store the value
				imageSeq[ zOut ].data( xOut, yOut ) = vInt;
			}
		}

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );

	// check for user cancel
	if (checkCommandEvents())
		return;

	// make sure output path exists
	createDir( outputPath );

	// third pass: save the data
	disp( 1, "saving" );
	for (int i = 0; i < imageSeq.count(); i++) {
		status( "        %d of %d", i + 1, imageSeq.count() );
		saveImage( imageSeq[ i ], outputPath + sprintF( "/%04d.png", i ) );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	status( "\n" );
}


/// convert a set of images to an mgz file
void convertImagesToMghFile( const String &inputPath, const String &outputFileName, float xSpacing, float ySpacing, float zSpacing, int width, int height, int valueScaleFactor ) {

	// open the output MGZ file
	File file( outputFileName, FILE_WRITE, FILE_GZIP_BINARY );
	if (file.openSuccess() == false) {
		warning( "unable to open input file: %s", outputFileName.c_str() );
		return;
	}

	// get input file list
	Array<String> inputFileNames = dirFileList( inputPath, "", ".png" );
	if (inputFileNames.count() == 0) { 
		warning( "no input files at %s", inputPath.c_str() );
		return;
	}
	int depth = inputFileNames.count();

	// if not specified, load first image to get dimensions
	if (width == 0 || height == 0) {
		aptr<ImageGrayU> firstImage = load<ImageGrayU>( inputPath + "/" + inputFileNames[ 0 ] );
		width = firstImage->width();
		height = firstImage->height();
		firstImage.release();
	}

	// write header
	writeIntSwap( file, 1 ); // version
	writeIntSwap( file, width ); // width
	writeIntSwap( file, height ); // height
	writeIntSwap( file, depth ); // depth
	writeIntSwap( file, 1 ); // nFrames
	writeIntSwap( file, 4 ); // type
	writeIntSwap( file, 0 ); // dof
	writeShortSwap( file, 1 ); // goodRASFlag
	writeFloatSwap( file, xSpacing ); // xSpacing
	writeFloatSwap( file, ySpacing ); // ySpacing
	writeFloatSwap( file, zSpacing ); // zSpacing
	writeFloatSwap( file, -1 ); // xr
	writeFloatSwap( file, 0 ); // xa
	writeFloatSwap( file, 0 ); // xs
	writeFloatSwap( file, 0 ); // yr
	writeFloatSwap( file, 0 ); // ya
	writeFloatSwap( file, -1 ); // ys
	writeFloatSwap( file, 0 ); // zr
	writeFloatSwap( file, -1 ); // za
	writeFloatSwap( file, 0 ); // zs
	writeFloatSwap( file, 0.0f ); // cr
	writeFloatSwap( file, 0.0f ); // ca
	writeFloatSwap( file, 0.0f ); // cs

	// seek to data
	file.seek( 284, false );

	// loop over files
	for (int i = 0; i < inputFileNames.count(); i++) {

		// load an image
		String name = inputFileNames[ i ];
		aptr<ImageGrayU> input = load<ImageGrayU>( inputPath + "/" + name );
		input = resize( *input, width, height, true );
		assertAlways( input->width() == width );
		assertAlways( input->height() == height );

		// write the data
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				short val = input->data( x, y ) * valueScaleFactor;
				writeShortSwap( file, val );
			}
		}
	}
}


//-------------------------------------------
// COMMANDS
//-------------------------------------------


/// convert an mgh/mgz file to a set of images
void convertMghFileToImages( Config &conf ) {
	String fileName = addDataPath( conf.readString( "fileName" ) );
	String outputPath = addDataPath( conf.readString( "outputPath" ) );
	bool xySwap = conf.readBool( "xySwap", true );
	bool xzSwap = conf.readBool( "xzSwap", false );
	bool yzSwap = conf.readBool( "yzSwap", true );
	bool xFlip = conf.readBool( "xFlip", false );
	bool yFlip = conf.readBool( "yFlip", true );
	bool zFlip = conf.readBool( "zFlip", false );
	bool autoScaleValues = conf.readBool( "autoScaleValues", true );
	convertMghFileToImages( fileName, outputPath, xySwap, xzSwap, yzSwap, xFlip, yFlip, zFlip, autoScaleValues );
}


/// convert a set of images to an mgz file
void convertImagesToMghFile( Config &conf ) {
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputFileName = addDataPath( conf.readString( "outputFileName" ) );
	float xSpacing = conf.readFloat( "xSpacing", 0.1f );
	float ySpacing = conf.readFloat( "ySpacing", 0.1f );
	float zSpacing = conf.readFloat( "zSpacing", 0.1f );
	int valueScaleFactor = conf.readInt( "valueScaleFactor", 16 );
	convertImagesToMghFile( inputPath, outputFileName, xSpacing, ySpacing, zSpacing, 0, 0, valueScaleFactor );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initVolumeFile() {
	registerCommand( "vread", convertMghFileToImages );
	registerCommand( "vwrite", convertImagesToMghFile );
}


} // end namespace hb
