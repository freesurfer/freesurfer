#include "prep/HistoPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/other/Plot.h>
#include "prep/VolumeFile.h"
using namespace sbl;
namespace hb {


/// extract mask of histology foreground regions
aptr<ImageGrayU> createHistologyMask( const ImageGrayU &image, const String &visPath, const String &fileName ) {
	int width = image.width(), height = image.height();
	int histogramBorder = 200;
	int initThresh = 200;
	int minSize = 200;
	bool savePlots = true;

	// compute a histogram
	VectorF hist = toFloat( imageHistogram( image, histogramBorder, width - histogramBorder, histogramBorder, height - histogramBorder ) );
	hist = gaussFilter( hist, 3 );
	int thresh = nearestMinIndex( hist, initThresh );
	disp( 1, "input: %s, thresh: %d, hist at thresh: %f", fileName.c_str(), thresh, hist[ thresh ] );

	// extract mask
	aptr<ImageGrayU> blur = blurGauss( image, 4 );
	aptr<ImageGrayU> mask = threshold( *blur, (float) thresh, true );
	filterMaskComponents( *mask, minSize * minSize, width * height );
	fillMaskHoles( *mask, 30 * 30 );

	// save histogram plot
	if (savePlots) {
		Plot plot;
		plot.add( toDouble( hist ) );
		plot.addVertLine( (double) thresh );
		String plotFileName = visPath + "/histogram_" + fileName.leftOfLast( '.' ) + ".svg";
		plot.save( plotFileName );
	}
	return mask;
}


/// find lower and upper thresholds in a histogram corresponding to upper and lower percentiles
void findHistogramBounds( const VectorF &histogram, float lowerFrac, float upperFrac, int &lowerThresh, int &upperThresh ) {
	VectorF histogramBlur = gaussFilter( histogram, 3 );
	float fullSum = histogramBlur.sum();
	float sum = 0;
	for (int i = 0; i < 256; i++) {
		sum += histogramBlur[ i ];
		if (sum > fullSum * lowerFrac && lowerThresh == 0)
			lowerThresh = i;
		if (sum > fullSum * upperFrac && upperThresh == 0)
			upperThresh = i;
	}
}


/// normalize the brightness of each channel separately
void normalizeChannels( ImageColorU &image ) {
	int width = image.width(), height = image.height();

	// loop over channels
	for (int c = 0; c < 3; c++) {

		// loop over image, creating histogram of non-background pixels
		VectorF histogram( 256 );
		histogram.clear( 0 );
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int v = image.data( x, y, c );
				if (v > 0 && v < 254) {
					histogram[ v - 1 ] += 0.5f;
					histogram[ v ] += 1.0f;
					histogram[ v + 1 ] += 0.5f;
				}
			}
		}

		// compute thresholds from histogram
		int lowerThresh = 0, upperThresh = 0;
		findHistogramBounds( histogram, 0.1f, 0.9f, lowerThresh, upperThresh );

		// apply thresholds
		if (lowerThresh != upperThresh ) {
			for (int y = 0; y < height; y++) {
				for (int x = 0; x < width; x++) {
					int v = image.data( x, y, c );
					if (v < 255) {
						v = 64 + (v - lowerThresh) * 128 / (upperThresh - lowerThresh);
						v = bound( v, 5, 250 );
						image.data( x, y, c ) = v;
					}
				}
			}
		}
	}
}


/// combine multiple histology wavelengths into a single image
aptr<ImageColorU> combineHistoWavelengths( const ImageGrayU &image700, const ImageGrayU &image800, const ImageGrayU &mask ) {
	int width = image700.width(), height = image700.height();
	assertAlways( image800.width() == width && image800.height() == height );

	// create combined image using each wavelength as a separate color channel
	aptr<ImageColorU> colorHisto( new ImageColorU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask.data( x, y )) {
				int r = image800.data( x, y );
				int g = image700.data( x, y );
				colorHisto->setRGB( x, y, r, g, 0 );
			} else {
				colorHisto->setRGB( x, y, 255, 255, 255 );
			}
		}
	}
	return colorHisto;
}


/// extract a foreground region (one tissue slice, hopefully) from the slide image
void extractRegion( const ImageColorU &input, ImageGrayU &mask, int xCheck, int yCheck, bool firstPass, int &outputWidth, int &outputHeight, const String &outputPath, int rotateAngle, int id ) {
	int border = 100;

	// find region nearest to check point
	int width = input.width(), height = input.height();
	int distSqdBest = width * width + height * height;
	int xBest = -1, yBest = -1;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			if (mask.data( x, y ) == 255) {
				int xDiff = x - xCheck;
				int yDiff = y - yCheck;
				int distSqd = xDiff * xDiff + yDiff * yDiff;
				if (distSqd < distSqdBest) {
					xBest = x;
					yBest = y;
					distSqdBest = distSqd;
				}
			}
		}
	}
	if (xBest == -1) {
		warning( "no points found" );
		return;
	} 

	// fill found region
	float xCent = 0, yCent = 0;
	int xMin = 0, xMax = 0, yMin = 0, yMax = 0;
	floodFill( mask, 255, 255, 100, xBest, yBest, &xCent, &yCent, &xMin, &xMax, &yMin, &yMax );

	// if first pass, compute desired output size and return
	if (firstPass) {

		// compute output width/height
		int curOutputWidth = xMax - xMin + 2 * border;
		int curOutputHeight = yMax - yMin + 2 * border;
		curOutputWidth -= (curOutputWidth & 7); // make divisible by 8
		curOutputHeight -= (curOutputHeight & 7); // make divisible by 8

		// update global output size
		if (curOutputWidth > outputWidth)
			outputWidth = curOutputWidth;
		if (curOutputHeight > outputHeight)
			outputHeight = curOutputHeight;
		return;
	}

	// compute border needed to center the image
	int xBorder = (outputWidth - (xMax - xMin)) / 2;
	int yBorder = (outputHeight - (yMax - yMin)) / 2;

	// create output image
	ImageColorU output( outputWidth, outputHeight );
	output.clear( 255, 255, 255 );
	for (int y = yMin; y <= yMax; y++) {
		for (int x = xMin; x <= xMax; x++) {
			if (mask.data( x, y ) == 100) {
				int xOut = x - xMin + xBorder;
				int yOut = y - yMin + yBorder;
				if (output.inBounds( xOut, yOut )) {
					for (int c = 0; c < 3; c++)
						output.data( xOut, yOut, c ) = input.data( x, y, c );
				}
				mask.data( x, y ) = 50;
			}
		}
	}

	// normalize image channels
	normalizeChannels( output );

	// save result
	String outputFileName = outputPath + "/" + sprintF( "%03d", id ) + ".png";
	saveImage( output, outputFileName );

	// save individual channels as grayscale images
	ImageGrayU output700( outputWidth, outputHeight );
	ImageGrayU output800( outputWidth, outputHeight );
	for (int y = 0; y < outputHeight; y++) {
		for (int x = 0; x < outputWidth; x++) {
			output700.data( x, y ) = output.g( x, y );
			output800.data( x, y ) = output.r( x, y );
		}
	}
	createDir( outputPath + "/gray700" );
	outputFileName = outputPath + "/gray700/" + sprintF( "%03d", id ) + "_700.png";
	saveImage( output700, outputFileName );
	createDir( outputPath + "/gray800" );
	outputFileName = outputPath + "/gray800/" + sprintF( "%03d", id ) + "_800.png";
	saveImage( output800, outputFileName );
}


/// prepare histology slide images: segment fg/bg, split multiple samples, normalize each sample
/// (note: we assume that slide photo images are a subset of the licor images)
void prepareHistologyImages( Config &conf ) {

	// get command parameters
	int rotateAngle = conf.readInt( "rotateAngle", 0 );
	String inputPath = addDataPath( conf.readString( "inputPath", "histo/raw" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "histo/split" ) );

	// create output directory
	createDir( outputPath );

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// these will hold the output size
	int outputWidth = 0, outputHeight = 0;

	// get list of input images
	Array<String> fileList = dirFileList( inputPath, "", "_700.png" );

	// do two passes:
	// 1. determine size bounds
	// 2. save images
	for (int pass = 0; pass < 2; pass++) {

		// loop over input images
		for (int i = 0; i < fileList.count(); i++) {

			// load input image
			String inputFileName = inputPath + "/" + fileList[ i ];
			aptr<ImageGrayU> image700 = load<ImageGrayU>( inputFileName );
			aptr<ImageGrayU> image800 = load<ImageGrayU>( inputFileName.leftOfLast( '_' ) + "_800.png" );
			int width = image700->width(), height = image700->height();

			// get slice numbers
			// we assume that the file name is of the form:
			//     x_sliceId1_sliceId2_y.png 
			// or:
			//     x_sliceId2_y.png 
			Array<String> split = fileList[ i ].split( "_" );
			if (split.count() < 3) {
				warning( "unexpected file name format (expected at least 2 underscores)" );
				return;
			}
			int sliceId1 = split[ split.count() - 3 ].toInt(); // will be zero if not numeric
			int sliceId2 = split[ split.count() - 2 ].toInt();
			disp( 2, "id1: %d, id2: %d", sliceId1, sliceId2 );

			// compute mask
			aptr<ImageGrayU> mask = createHistologyMask( *image700, visPath, fileList[ i ] );

			// combine histo wavelengths into single color image
			aptr<ImageColorU> image = combineHistoWavelengths( *image700, *image800, *mask );

			// if two IDs, extract two regions
			if (sliceId1) {
				extractRegion( *image, *mask, width / 4, height / 4, pass == 0, outputWidth, outputHeight, outputPath, rotateAngle, sliceId1 );
				extractRegion( *image, *mask, width * 3 / 4, height * 3 / 4, pass == 0, outputWidth, outputHeight, outputPath, rotateAngle, sliceId2 );

			// if one ID, extract one region
			} else {
				extractRegion( *image, *mask, width / 2, height / 2, pass == 0, outputWidth, outputHeight, outputPath, rotateAngle, sliceId2 );
			}

			// check for user cancel
			if (checkCommandEvents())
				return;
		}

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
}


// create an unregistered histology volume
void createHistologyVolume( Config &conf ) {
	String inputPath = addDataPath( conf.readString( "inputPath", "histo/split/gray800" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "histo/" ) );
	float xSpacing = conf.readFloat( "xSpacing", 1 );
	float ySpacing = conf.readFloat( "ySpacing", 1 );
	float zSpacing = conf.readFloat( "zSpacing", 2 );
    
	// create mapping file
	Array<String> fileList = dirFileList( inputPath, "", ".png" );
	String outputMapFileName = outputPath + "/splitVolMap.txt";
	File file( outputMapFileName, FILE_WRITE, FILE_TEXT );
	if (file.openSuccess()) {
		file.writeF( "mgz_index, histo_slice_index\n" );
		for (int i = 0; i < fileList.count(); i++) {
			int hSliceIndex = fileList[ i ].leftOfFirst( '.' ).toInt();
			file.writeF( "%d, %d\n", i, hSliceIndex );
		}
	}

	// create mgz file
	String outputVolumeFileName = outputPath + "/splitVol.mgz";
	convertImagesToMghFile( inputPath, outputVolumeFileName, xSpacing, ySpacing, zSpacing );
}


/// convert a volume generated by hvol back into images
/// (this command is the inverse of hvol)
void convertHistologyVolumeIntoImages( Config &conf ) {

	// get command parameters
	String volumeFileName = addDataPath( conf.readString( "volumeFileName", "splitVol.mgz" ) );
	String mapFileName = addDataPath( conf.readString( "mapFileName", "splitVolMap.txt" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "." ) );
	bool autoScaleValues = conf.readBool( "autoScaleValues", false );

	// split the volume into individual images
	convertMghFileToImages( volumeFileName, outputPath, false, false, false, false, false, false, autoScaleValues );

	// open file containing histo index to mgz index mapping
	File mapFile( mapFileName, FILE_READ, FILE_TEXT );
	if (mapFile.openSuccess() == false) {
		warning( "unable to open map file: %s", mapFileName.c_str() );
		return;
	}

	// rename files according to map
	mapFile.readLine(); // skip header
	while (mapFile.endOfFile() == false) {
		String line = mapFile.readLine();
		Array<String> parts = line.split( "," );
		if (parts.count() == 2) {
			int mgzIndex = parts[ 0 ].strip().toInt();
			int histoIndex = parts[ 1 ].strip().toInt();
			String oldFileName = outputPath + sprintF( "/%04d.png", mgzIndex );
			String newFileName = outputPath + sprintF( "/_%04d.png", histoIndex );
			disp( 1, "rename: %d -> %d", mgzIndex, histoIndex /*oldFileName.c_str(), newFileName.c_str()*/ );
			moveFile( oldFileName, newFileName );
		}
	}

	// check whether any files were not renamed
	Array<String> fileList = dirFileList( outputPath, "", ".png" );
	for (int i = 0; i < fileList.count(); i++) {
		String fileName = fileList[ i ];
		if (fileName.startsWith( "_" ) == false) {
			warning( "deleting unexpected file: %s", fileName.c_str() );
			deleteFile( outputPath + "/" + fileName );
		}
	}

	// rename files again to remove underscores
	fileList = dirFileList( outputPath, "_", ".png" );
	for (int i = 0; i < fileList.count(); i++) {
		String fileName = fileList[ i ];
		moveFile( outputPath + "/" + fileName, outputPath + "/" + fileName.rightOf( 0 ) );
	}
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initHistoPrep() {
	registerCommand( "hprep", prepareHistologyImages );
	registerCommand( "hvol", createHistologyVolume );
	registerCommand( "hslice", convertHistologyVolumeIntoImages );
}


} // end namespace hb
