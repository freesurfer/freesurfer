#include "prep/Polarization.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/other/Plot.h>
using namespace sbl;
namespace hb {


/// read a set of polarization test images and extract summary information
void readPolarizationImages( Config &conf ) {
	
	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputFileName = addDataPath( conf.readString( "outputFileName", "polarization.csv" ) );
	int angleCount = conf.readInt( "angleCount", 18 );
	int angleStep = conf.readInt( "angleStep", 10 );
	int polarizationStateCount = conf.readInt( "polarizationStateCount", 16 );

	// get list of images
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	if (angleCount * polarizationStateCount != fileList.count()) {
		warning( "found images: %d, expected images: %d", fileList.count(), angleCount * polarizationStateCount );
		return;
	}

	// open output file
	File outputFile( outputFileName, FILE_WRITE, FILE_TEXT );
	if (outputFile.openSuccess() == false) {
		warning( "unable to open output file: %s", outputFileName.c_str() );
		return;
	}

	// loop over images
	int index = 0;
	for (int angleIndex = 0; angleIndex < angleCount; angleIndex++) {
		int angle = angleIndex * angleStep;
		outputFile.writeF( "%d", angle );

		// loop over images for this angle
		for (int i = 0; i < polarizationStateCount; i++) {

			// load the image
			assertAlways( index < fileList.count() );
			aptr<ImageGrayU> input = load<ImageGrayU>( inputPath + "/" + fileList[ index ] );

			// compute bounds of area of interest
			int width = input->width(), height = input->height();
			int xMin = width / 2 - 100;
			int xMax = width / 2 + 100;
			int yMin = height / 2 - 100;
			int yMax = height / 2 + 100;

			// get brightness
			double sum = 0; 
			for (int y = yMin; y <= yMax; y++) {
				for (int x = xMin; x <= xMax; x++) {
					sum += input->data( x, y );
				}
			}
			double value = sum / (double) ((xMax - xMin + 1) * (yMax - yMin + 1));
			disp( 1, "index: %d, angle: %d, file: %s, value: %f", index, angle, fileList[ index ].c_str(), value );

			// add to output file
			outputFile.writeF( ", %f", value );

			// move to next file
			index++;

			// check for user cancel
			if (checkCommandEvents())
				break;
		}

		// start new line in output file
		outputFile.writeF( "\n" );
		outputFile.flush();

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


/// use information from polarization test images to predict the specimen polarization 
void predictFromPolarizationData( Config &conf ) {

	// get command parameters
	String inputFileName = addDataPath( conf.readString( "inputFileName", "polarization.csv" ) );

	// open summary data file
	File inputFile( inputFileName, FILE_READ, FILE_TEXT );
	if (inputFile.openSuccess() == false) {
		warning( "unable to open input file: %s", inputFileName.c_str() );
		return;
	}

	// will contain summary data
	VectorI angles;
	Array<VectorF> values;

	// read summary data
	while (inputFile.endOfFile() == false) {
		String line = inputFile.readLine();
		Array<String> split = line.split( "," );
		if (split.count() > 1) {
			angles.append( split[ 0 ].toInt() );
			VectorF value;
			for (int i = 1; i < split.count(); i++) {
				float v = split[ i ].toFloat();
				value.append( v );
			}
			values.appendCopy( value );
		}
	}

	// error for each prediction
	VectorF error;

	// predict angle from values
	for (int i = 1; i < angles.length() - 1; i++) {
		const VectorF &v = values[ i ];
		const VectorF &v1 = values[ i - 1 ];
		const VectorF &v2 = values[ i + 1 ];
		float aBest = -1, distSqdBest = 0;
		for (float a = 0; a <= 1; a += 0.0001f) {

			// compute interpolated vector
			VectorF interp = v1;
			for (int j = 0; j < v1.length(); j++) 
				interp[ j ] += a * (v2[ j ] - v1[ j ]);	

			// compute distance to observation
			float dSqd = distSqd( interp, v );

			// keep best
			if (dSqd < distSqdBest || aBest == -1) {
				distSqdBest = dSqd;
				aBest = a;
			}
		}
		float a1 = (float) angles[ i - 1 ];
		float a2 = (float) angles[ i + 1 ];
		float estAngle = a1 + aBest * (a2 - a1);
		error.append( fAbs( (float) angles[ i ] - estAngle ) );
		disp( 1, "angle: %d, est angle: %f", angles[ i ], estAngle );
	}
	disp( 1, "error: %f / %f / %f", error.min(), error.mean(), error.max() );
}


/// shrink polarization calibration images so that they can be loaded more quickly
void preparePolarizationImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath", "blockface/polarization/calib" ) );
	int calibScaleFactor = conf.readInt( "calibScaleFactor", 16 );

	// loop over input files, processing each one
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	for (int i = 0; i < fileList.count(); i++) {
		disp( 1, "index: %d, file: %s", i, fileList[ i ].c_str() );
		aptr<ImageGrayU> input = load<ImageGrayU>( inputPath + "/" + fileList[ i ] );

		// blur, shrink, and blur again
		input = blurGauss( *input, 10 );
		int newWidth = input->width() / calibScaleFactor;
		int newHeight = input->height() / calibScaleFactor;
		input = resize( *input, newWidth, newHeight, true );
		input = blurGauss( *input, 5 );

		// save for future use
		String outputFileName = inputPath + "/" + fileList[ i ].leftOfLast( '.' ) + ".png";
		saveImage( *input, outputFileName );

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
}


/// compute polarization orientation colorization
/// (assumes angle in [0,180] and magnitude in [0,1])
void orientationColor( float angle, float magnitude, int &r, int &g, int &b ) {
	float theta = 2 * angle * 3.14159f / 180.0f;
	r = 128 + round( 128.0f * cosf( theta ) * magnitude * magnitude );
	g = 128 + round( 128.0f * sinf( theta ) * magnitude * magnitude );
	b = 128;
	r = bound( r, 0, 255 );
	g = bound( g, 0, 255 );
}


/// estimate orientation for a single tissue slice
void processPolarizationSlice( const Array< Array<ImageGrayU> > &calibImages, 
							   const String &samplePath, const Array<String> &sampleFileList, const String &outputPath, 
							   int sliceIndex, int polarizationStateCount, int angleCount, int angleStep, 
							   int sampleScaleFactor, float angleInterpSigma ) {

	// enable this for more diagnostic plots
	bool savePixelPlots = true;

	// normalization mode
	bool normBounds = false;
	bool normCirc = false;

	// enable angle interpolation
	bool interpolate = false;

	// load sample images
	Array<ImageGrayU> sampleImages;
	for (int i = sliceIndex * polarizationStateCount; i < (sliceIndex + 1) * polarizationStateCount; i++) {
		assertAlways( i < sampleFileList.count() );
		aptr<ImageGrayU> input = load<ImageGrayU>( samplePath + "/" + sampleFileList[ i ] );

		// shrink
		int newWidth = input->width() / sampleScaleFactor;
		int newHeight = input->height() / sampleScaleFactor;
		input = resize( *input, newWidth, newHeight, true );

		// store in collection of sample images
		sampleImages.append( input.release() );

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
	disp( 1, "done loading sample images, slice: %d", sliceIndex );

	// this will be the output
	int width = sampleImages[ 0 ].width();
	int height = sampleImages[ 0 ].height();
	ImageColorU outputAngle( width, height );
	ImageColorU outputAngleStrength( width, height );
	ImageGrayU outputBright( width, height );
	ImageGrayF outputStDev( width, height );
	outputAngle.clear( 255, 255, 255 );
	outputAngleStrength.clear( 255, 255, 255 );
	outputStDev.clear( 0 );
	ImageColorU outputMarker( width, height );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int v = sampleImages[ 0 ].data( x, y );
			outputMarker.setRGB( x, y, v, v, v );
		}
	}

	// gather stats
	VectorF stDevVect, brightVect, angleVect;

	// compute factor used for interpolating angles
	float factor = gaussFactor( angleInterpSigma );

	// compute output image
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			// compute feature vector
			VectorF sampleFeat( polarizationStateCount );
			for (int i = 0; i < polarizationStateCount; i++) 
				sampleFeat[ i ] = sampleImages[ i ].data( x, y );
			float brightness = sampleFeat.mean();
			VectorF sampleFeatPreNorm = sampleFeat;
			if (normBounds)
				normalizeBounds( sampleFeat ); 
			if (normCirc)
				normalize( sampleFeat );
			add( sampleFeat, -sampleFeat.mean(), sampleFeat );

			// compute position in calibration images
			// note: this may be off a little bit, but shouldn't matter
			int calibWidth = calibImages[ 0 ][ 0 ].width();
			int calibHeight = calibImages[ 0 ][ 0 ].height();
			int xCalib = x * calibWidth / sampleImages[ 0 ].width();
			int yCalib = y * calibHeight / sampleImages[ 0 ].height();
			if (xCalib > calibWidth - 1)
				xCalib = calibWidth - 1;
			if (yCalib > calibHeight - 1)
				yCalib = calibHeight - 1;

			// for diagnostics
			Array<VectorF> calibFeats;

			// compute feature vector distance to each calibration angle
			VectorF distSqd( angleCount );
			float distSqdBest = 0;
			int angleIndexBest = -1;
			for (int angleIndex = 0; angleIndex < angleCount; angleIndex++) {

				// compute feature vector from calibration images 
				VectorF calibFeat( polarizationStateCount );
				for (int i = 0; i < polarizationStateCount; i++)
					calibFeat[ i ] = calibImages[ angleIndex ][ i ].data( xCalib, yCalib );
				if (normBounds)
					normalizeBounds( calibFeat );
				if (normCirc)
					normalize( calibFeat );
				add( calibFeat, -calibFeat.mean(), calibFeat );
				calibFeats.appendCopy( calibFeat );

				// compute distance 
				float dSqd = sbl::distSqd( sampleFeat, calibFeat );
				distSqd[ angleIndex ] = dSqd;
				if (dSqd < distSqdBest || angleIndexBest == -1) {
					distSqdBest = dSqd;
					angleIndexBest = angleIndex;
				}
			}
			float angle = angleIndexBest * (float) angleStep;

			// take weighted average (by distance) of angles near the best
			if (interpolate) {
				int minAngleIndex = angleIndexBest - 5;
				int maxAngleIndex = angleIndexBest + 5;
				double wtSum = 0, wtAngleSum = 0;
				for (int angleIndex = minAngleIndex; angleIndex <= maxAngleIndex; angleIndex++) {

					// wrap to [0, angleCount - 1]
					int angleIndexWrapped = angleIndex;
					if (angleIndexWrapped < 0)
						angleIndexWrapped += angleCount;
					if (angleIndexWrapped >= angleCount)
						angleIndexWrapped -= angleCount;
					
					// update average using weight computed from distance
					double wt = gauss( distSqd[ angleIndexWrapped ], factor );
					double a = angleIndex * (double) angleStep;
					wtSum += wt;
					wtAngleSum += wt * a;
				}
				if (wtSum > 0.0001)
					angle = (float) (wtAngleSum / wtSum);
			}

			// compute magnitude from stdev of distances
			float stDev = sbl::stDev( distSqd, distSqd.mean() );
			float magnitude = stDev;
			magnitude = bound( magnitude, 0.0f, 1.0f );

			// store stats
			stDevVect.append( stDev );
			brightVect.append( brightness );
			angleVect.append( angle );

			// store colorized angle
			if (brightness > 10) {
				int r = 0, g = 0, b = 0;
				orientationColor( angle, 1.0f, r, g, b );
				outputAngle.setRGB( x, y, r, g, b );
				orientationColor( angle, magnitude, r, g, b );
				outputAngleStrength.setRGB( x, y, r, g, b );
				outputStDev.data( x, y ) = stDev * 150;
			}

			// compute brightness image
			int v = round( brightness ) * 2;
			v = bound( v, 0, 255 );
			outputBright.data( x, y ) = v;

			// save diagnostic plots for select pixels
			if (savePixelPlots) {
				if (((y == height / 2 && (x % 120) == 0) || (x == width / 2 && (y % 120) == 0))
					&& x > 50 && x < width - 50 && y > 50 && y < height - 50) {
					Plot featPlot;
					featPlot.setAxisLabels( "polarization state", "brightness" );
					VectorF calibFeatMean( calibFeats[ 0 ].length() );
					calibFeatMean.clear( 0 );
					featPlot.setColor( 230, 230, 230 );
					for (int i = 0; i < calibFeats.count(); i++) {
						featPlot.add( toDouble( calibFeats[ i ] ) );
						for (int j = 0; j < calibFeatMean.length(); j++)
							calibFeatMean[ j ] += (float) calibFeats[ i ][ j ] / (float) calibFeats.count();
					}
					featPlot.setColor( 200, 200, 0 );
					featPlot.add( toDouble( calibFeatMean ) );
					featPlot.setColor( 200, 0, 0 );
					featPlot.add( toDouble( sampleFeat ) );
					featPlot.save( outputPath + sprintF( "/feat_%04d_%d_%d_feat.svg", sliceIndex, x, y ) );
					aptr<Plot> distPlot = simplePlot( toDouble( distSqd ) );
					featPlot.setAxisLabels( "polarization state", "distance" );
					distPlot->setColor( 255, 0, 0 );
					distPlot->addVertLine( (angle / (float) angleStep) + 1.0f );
					distPlot->save( outputPath + sprintF( "/dist_%04d_%d_%d.svg", sliceIndex, x, y ) );
					disp( 2, "pre-norm: %f / %f / %f", sampleFeatPreNorm.min(), sampleFeatPreNorm.mean(), sampleFeatPreNorm.max() );
					drawCross( outputMarker, x, y, 10, 255, 0, 0, false );
					String label1 = sprintF( "%d %d", x, y );
					drawText( outputMarker, label1, x + 2, y, 255, 0, 0 );
					String label2 = sprintF( "%3.1f", angle );
					drawText( outputMarker, label2, x + 2, y - 14, 255, 255, 0 );
				}
			}
		}

		// check for user cancel
		if (checkCommandEvents())
			return;
	}

	// display statistics
	disp( 2, "st. dev: %f / %f / %f", stDevVect.min(), stDevVect.mean(), stDevVect.max() );
	disp( 2, "bright: %f / %f / %f", brightVect.min(), brightVect.mean(), brightVect.max() );
	disp( 2, "angle: %f / %f / %f", angleVect.min(), angleVect.mean(), angleVect.max() );

	// save results
	saveImage( outputAngle, outputPath + sprintF( "/angle_%04d.jpg", sliceIndex ) );
	saveImage( outputAngleStrength, outputPath + sprintF( "/angleStrength_%04d.jpg", sliceIndex ) );
	saveImage( outputBright, outputPath + sprintF( "/brightness_%04d.jpg", sliceIndex ) );
	saveImage( outputStDev, outputPath + sprintF( "/strength_%04d.jpg", sliceIndex ) );
	aptr<Plot> plot = histogramPlot( toDouble( angleVect ), 100 );
	plot->save( outputPath + sprintF( "/angleHistogram_%04d.svg", sliceIndex ) );
	saveImage( outputMarker, outputPath + sprintF( "/marker_%04d.png", sliceIndex ) );
}


/// compute polarization angles on a sample given a set of calibration images
void computeFullPolarizationImage( Config &conf ) {

	// get command parameters
	String calibPath = addDataPath( conf.readString( "calibPath", "blockface/polarization/calib" ) );
	String samplePath = addDataPath( conf.readString( "samplePath", "blockface/polarization/input" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "blockface/polarization/output" ) );
	int angleCount = conf.readInt( "angleCount", 18 );
	int angleStep = conf.readInt( "angleStep", 10 );
	int polarizationStateCount = conf.readInt( "polarizationStateCount", 16 );
	int sampleScaleFactor = conf.readInt( "sampleScaleFactor", 8 );
	float angleInterpSigma = conf.readFloat( "angleInterpSigma", 1.0f );

	// make sure output path exists
	createDir( outputPath );

	// get list of calibration images (small versions)
	Array<String> calibFileList = dirFileList( calibPath, "", ".png" );
	if (angleCount * polarizationStateCount != calibFileList.count()) {
		warning( "found calib images: %d, expected images: %d", calibFileList.count(), angleCount * polarizationStateCount );
		return;
	}

	// get list of sample images (full size versions)
	Array<String> sampleFileList = dirFileList( samplePath, "", ".JPG" );
	int sliceCount = sampleFileList.count() / polarizationStateCount;
	if (sliceCount == 0 || sliceCount * polarizationStateCount != sampleFileList.count()) {
		warning( "found sample images: %d, expected a (positive) multiple of %d images", sampleFileList.count(), polarizationStateCount );
		return;
	}

	// this will the calibration images (small versions)
	Array< Array<ImageGrayU> > calibImages;

	// loop over calib images
	int index = 0;
	for (int angleIndex = 0; angleIndex < angleCount; angleIndex++) {
		calibImages.append( new Array<ImageGrayU>() );

		// loop over images for this angle
		for (int i = 0; i < polarizationStateCount; i++) {

			// load and store the image
			assertAlways( index < calibFileList.count() );
			aptr<ImageGrayU> input = load<ImageGrayU>( calibPath + "/" + calibFileList[ index ] );
			calibImages[ angleIndex ].append( input.release() );
			index++;

			// check for user cancel
			if (checkCommandEvents())
				return;
		}

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
	disp( 1, "done loading calib images" );

	// loop over slices, estimating orientations for each one
	for (int sliceIndex = 0; sliceIndex < sliceCount; sliceIndex++) {
		processPolarizationSlice( calibImages, samplePath, sampleFileList, outputPath, sliceIndex, polarizationStateCount, angleCount, angleStep, sampleScaleFactor, angleInterpSigma );

		// check for user cancel
		if (checkCommandEvents())
			return;
	}
}


/// create key for reading polarization orientation colorization
void createPolarizationKey( Config &conf ) {
	int size = 500;
	int radius = 200;
	ImageColorU output( size, size );
	output.clear( 255, 255, 255 );
	for (float angle = 0; angle < 180; angle += 0.1f) {
		for (float dist = 0; dist < radius; dist += 0.5f) {
			float magnitude = (float) dist / (float) radius;
			int r = 0, g = 0, b = 0;
			orientationColor( angle, magnitude, r, g, b );
			float theta = angle * 3.14159f / 180.0f;
			int c = round( cosf( theta ) * (float) dist );
			int s = round( sinf( theta ) * (float) dist );
			int x = size / 2 + c;
			int y = size / 2 + s;
			output.setRGB( x, y, r, g, b );
			x = size / 2 - c;
			y = size / 2 - s;
			output.setRGB( x, y, r, g, b );
		}
	}
	saveImage( output, dataPath() + "angleKey.jpg" );
}


void checkLinear( Config &conf ) {
	
	// get command parameters
	String inputFileName1 = addDataPath( conf.readString( "inputFileName1" ) );
	String inputFileName2 = addDataPath( conf.readString( "inputFileName2" ) );

	// load input images
//	aptr<ImageColorU> image1 = load<ImageColorU>( inputFileName1 );
//	aptr<ImageColorS> image2 = load<ImageColorS>( inputFileName2 );
//	disp( 1, "image 1: %d x %d", image1->width(), image1->height() );
//	disp( 1, "image 2: %d x %d", image2->width(), image2->height() );

	// generate plot
	//aptr<Plot> plot( new Plot( "linear vs. jpg" ) );
	

	// save plot
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initPolarization() {
	registerCommand( "pread", readPolarizationImages );
	registerCommand( "ppred", predictFromPolarizationData );
	registerCommand( "pprep", preparePolarizationImages );
	registerCommand( "pfull", computeFullPolarizationImage );
	registerCommand( "pkey", createPolarizationKey );
	registerCommand( "plin", checkLinear );
}


} // end namespace hb
