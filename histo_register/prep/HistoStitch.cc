#include "prep/HistoStitch.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/Optimizer.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/system/FileSystem.h>
#include <sbl/other/Plot.h>
#include "prep/StitchNode.h"
using namespace sbl;
namespace hb {


//-------------------------------------------
// SHARED OPTIMIZATION UTILITIES
//-------------------------------------------


/// compute the pixel value error between a pair of registered nodes
void computeError( const StitchNode &node1, const StitchNode &node2, double &sumError, int &countError, bool mapBrightness, bool verbose ) {
	if (verbose) {
		disp( 3, "node1: %d, %d, node2: %d, %d", node1.xImageIndex(), node1.yImageIndex(), node2.xImageIndex(), node2.yImageIndex() );
		disp( 3, "bright1: %f, %f, bright2: %f, %f", node1.minBrightness(), node1.maxBrightness(), node2.minBrightness(), node2.maxBrightness() );
	}
//	disp( 3, "[%f,%f]-[%f,%f]", node1.transform().xOffset(), node1.transform().yOffset(), node2.transform().xOffset(), node2.transform().yOffset() );
	const ImageGrayU &image1 = node1.image();
	const ImageGrayU &image2 = node2.image();
	int width = image1.width(), height = image1.height();
	double vSum1 = 0, vSum2 = 0, mSum1 = 0, mSum2 = 0;
	int count = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			// compute coordinate in output image (global coord)
			Point2 outPoint = node1.transform().mapForward( Point2( x, y ) );

			// compute coordinate in second node (local coord)
			Point2 inPoint = node2.transform().mapBackward( outPoint );

			// if in bounds, accumulate error
			int xIn = sbl::round( inPoint.x );
			int yIn = sbl::round( inPoint.y );
			if (image2.inBounds( xIn, yIn )) {
				double err = 0;
				if (verbose) {
					vSum1 += image1.data( x, y );
					vSum2 += image2.data( xIn, yIn );
					mSum1 += node1.mapBrightness( image1.data( x, y ) );
					mSum2 += node2.mapBrightness( image2.data( xIn, yIn ) );
					count++;
				}
				if (mapBrightness)
					err = node1.mapBrightness( image1.data( x, y ) ) - node2.mapBrightness( image2.data( xIn, yIn ) );
				else
					err = (double) (image1.data( x, y ) - image2.data( xIn, yIn ));
				if (err < 0)
					err = -err;
				sumError += err;
				countError++;
			}
		}
	}
	if (verbose) {
		if (count) {
			vSum1 /= (double) count;
			vSum2 /= (double) count;
			mSum1 /= (double) count;
			mSum2 /= (double) count;
			disp( 3, "v1: %f, v2: %f, m1: %f, m2: %f", vSum1, vSum2, mSum1, mSum2 );
		}
	}
}


/// update the position of each node using the specified global transformation
void updateTransformsGlobal( StitchNodeSet &nodeSet, float xStep, float yStep, float xSkew, float ySkew, double scaleFactor ) {
	for (int xImageIndex = 0; xImageIndex < nodeSet.xCount(); xImageIndex++) {
		for (int yImageIndex = 0; yImageIndex < nodeSet.yCount(); yImageIndex++) {
			StitchNode &node = nodeSet.node( xImageIndex, yImageIndex );
			float xOffset = (float) xImageIndex * xStep + (float) yImageIndex * xSkew; 
			float yOffset = (float) yImageIndex * yStep + (float) xImageIndex * ySkew;
			node.transform().setOffset( xOffset * (float) scaleFactor, yOffset * (float) scaleFactor );
		}
	}
}


/// count number of pixels inside mask on a given horizontal line
int horizCount( const ImageGrayU &mask, int y ) {
	int width = mask.width();
	int count = 0;
	for (int x = 0; x < width; x++)
		if (mask.data( x, y ))
			count++;
	return count;
}


/// count number of pixels inside mask on a given vertical line
int vertCount( const ImageGrayU &mask, int x ) {
	int height = mask.height();
	int count = 0;
	for (int y = 0; y < height; y++)
		if (mask.data( x, y ))
			count++;
	return count;
}


/// remove black borders from around a slide
aptr<ImageGrayU> removeBorders( const ImageGrayU &input ) {

	// we'll work at a smaller scale just for efficiency
	int scaleFactor = 8;
	int width = input.width() / scaleFactor;
	int height = input.height() / scaleFactor;
	aptr<ImageGrayU> small = resize( input, width, height, false );

	// compute mask of foreground areas
	aptr<ImageGrayU> mask = blurBoxAndThreshold( *small, 15, 5 );

	// find top bound (first line above center with too few fg pixels)
	int yMax = height - 1;
	for (int y = height / 2; y < height; y++) {
		if (horizCount( *mask, y ) < width / 4) {
			yMax = y;
			break;
		}
	}

	// find bottom bound (first line below center with too few fg pixels)
	int yMin = 0;
	for (int y = height / 2; y >= 0; y--) {
		if (horizCount( *mask, y ) < width / 4) {
			yMin = y;
			break;
		}
	}

	// find right bound (first line right of center with too few fg pixels)
	int xMax = width - 1;
	for (int x = width / 2; x < width; x++) {
		if (vertCount( *mask, x ) < height / 4) {
			xMax = x;
			break;
		}
	}

	// find left bound (first line left of center with too few fg pixels)
	int xMin = 0;
	for (int x = width / 2; x >= 0; x--) {
		if (vertCount( *mask, x ) < height / 4) {
			xMin = x;
			break;
		}
	}

	// move boundary inward a bit
	xMin += 10;
	xMax -= 10;
	yMin += 10;
	yMax -= 10;

	// crop the image
	xMin *= scaleFactor;
	xMax *= scaleFactor;
	yMin *= scaleFactor;
	yMax *= scaleFactor;
	disp( 2, "crop xMin: %d, xMax: %d, width: %d", xMin, xMax, input.width() );
	disp( 2, "crop yMin: %d, yMax: %d, height: %d", yMin, yMax, input.height() );
	return crop( input, xMin, xMax, yMin, yMax );
}


//-------------------------------------------
// GLOBAL OPTIMIZATION
//-------------------------------------------


/// The GlobalStitchObjective class is used to evaluate global stitching parameters (offset and skew between images).
class GlobalStitchObjective : public Objective {
public:

	// basic constructor
	GlobalStitchObjective( StitchNodeSet &nodeSet ) : m_nodeSet( nodeSet ) {}

	/// evaluate objective function at given point
	double eval( const VectorD &point ) {
		float xStep = (float) point[ 0 ];
		float yStep = (float) point[ 1 ];
		float xSkew = 0;
		float ySkew = 0;
		if (point.length() >= 4) {
			xSkew = (float) point[ 2 ];
			ySkew = (float) point[ 3 ];
		}

		// update transforms
		updateTransformsGlobal( m_nodeSet, xStep, yStep, xSkew, ySkew, 1.0 );

		// perform evaluation
		double sumError = 0;
		int countError = 0;
		for (int xImageIndex = 0; xImageIndex < m_nodeSet.xCount(); xImageIndex++) {
			for (int yImageIndex = 0; yImageIndex < m_nodeSet.yCount(); yImageIndex++) {
				StitchNode &node = m_nodeSet.node( xImageIndex, yImageIndex );

				// if node to right, compute error with it
				if (xImageIndex + 1 < m_nodeSet.xCount()) {
					StitchNode &otherNode = m_nodeSet.node( xImageIndex + 1, yImageIndex );
					computeError( node, otherNode, sumError, countError, false, false );
				}

				// if node above, compute error with it
				if (yImageIndex + 1 < m_nodeSet.yCount()) {
					StitchNode &otherNode = m_nodeSet.node( xImageIndex, yImageIndex + 1 );
					computeError( node, otherNode, sumError, countError, false, false );
				}
			}
		}
		double meanError = 0;
		if (countError)
			meanError = sumError / (double) countError;
		disp( 3, "xStep: %4.2f, yStep: %4.2f, xSkew: %6.4f, ySkew: %6.4f, countError: %d meanError: %f", xStep, yStep, xSkew, ySkew, countError, meanError );
		return meanError;
	}

private:

	// set of nodes being optimized
	StitchNodeSet &m_nodeSet;
};


/// optimize the global stitching transformation (offset and skew) parameteres
void optimizeGlobalTransform( StitchNodeSet &nodeSet, int &xStep, int &yStep, float &xSkew, float &ySkew, int scaleFactor ) {
	bool optSkew = true;

	// prepare data for optimizer
	GlobalStitchObjective objective( nodeSet );
	int paramCount = optSkew ? 4 : 2;
	VectorD start( paramCount ), lBound( paramCount ), uBound( paramCount );
	start[ 0 ] = (double) xStep / (double) scaleFactor;
	start[ 1 ] = (double) yStep / (double) scaleFactor;
	if (paramCount >= 4) {
		start[ 2 ] = (double) xSkew / (double) scaleFactor;
		start[ 3 ] = (double) ySkew / (double) scaleFactor;
	}
	lBound[ 0 ] = start[ 0 ] - 50;
	lBound[ 1 ] = start[ 1 ] - 50;
	uBound[ 0 ] = start[ 0 ] + 50;
	uBound[ 1 ] = start[ 1 ] + 50;
	if (paramCount >= 4) {
		lBound[ 2 ] = start[ 2 ] - 50;
		lBound[ 3 ] = start[ 3 ] - 50;
		uBound[ 2 ] = start[ 2 ] + 50;
		uBound[ 3 ] = start[ 3 ] + 50;
	}

	// initialize the optimizer
	SimplexOptimizer optimizer( objective );
	optimizer.setStart( start );
	optimizer.setBounds( lBound, uBound );

	// run optimizer
	VectorD result = optimizer.run();

	// return results
	xStep = sbl::round( result[ 0 ] * scaleFactor );
	yStep = sbl::round( result[ 1 ] * scaleFactor );
	if (result.length() >= 4) {
		xSkew = (float) result[ 2 ] * scaleFactor;
		ySkew = (float) result[ 3 ] * scaleFactor;
	}
}


//-------------------------------------------
// PER-NODE OPTIMIZATION
//-------------------------------------------


/// The NodeObjective class is used to optimize the parameters of a single node within a set of nodes.
class NodeObjective : public Objective {
public:

	// basic constructor
	NodeObjective( StitchNodeSet &nodeSet, StitchNode &node ) : m_nodeSet( nodeSet ), m_node( node ) {
		m_verbose = false;
	}

	/// evaluate objective function at given point
	double eval( const VectorD &point ) {
		int xIndex = m_node.xImageIndex();
		int yIndex = m_node.yImageIndex();
		double sumError = 0;
		int countError = 0;

		m_node.setBrightnessBounds( (float) point[ 0 ], (float) point[ 1 ] );

		// compute error with neighboring nodes
		if (xIndex - 1 >= 0) {
			StitchNode &otherNode = m_nodeSet.node( xIndex - 1, yIndex );
			computeError( m_node, otherNode, sumError, countError, true, m_verbose );
		}
		if (yIndex - 1 >= 0) {
			StitchNode &otherNode = m_nodeSet.node( xIndex, yIndex - 1 );
			computeError( m_node, otherNode, sumError, countError, true, m_verbose );
		}
		double meanError = 0;
		if (countError)
			meanError = sumError / (double) countError;
		if (m_verbose)
			disp( 3, "min: %6.4f, max: %6.4f, countError: %d meanError: %f", point[ 0 ], point[ 1 ], countError, meanError );
		return meanError;
	}

	/// enable diagnostic output
	inline void setVerbose() { m_verbose = true; }

private:

	// internal data
	StitchNodeSet &m_nodeSet;
	StitchNode &m_node;
	bool m_verbose;
};


/// optimize the parameters for a single node, holding fixed the parameters of the other nodes in the set
void optimizeNode( StitchNodeSet &nodeSet, StitchNode &node ) {

	// prepare data for optimizer
	NodeObjective objective( nodeSet, node );
	VectorD start( 2 ), lBound( 2 ), uBound( 2 );
	start[ 0 ] = node.minBrightness();
	start[ 1 ] = node.maxBrightness();
	lBound[ 0 ] = start[ 0 ] - 50;
	lBound[ 1 ] = start[ 1 ] - 50;
	uBound[ 0 ] = start[ 0 ] + 50;
	uBound[ 1 ] = start[ 1 ] + 50;

	// initialize the optimizer
	SimplexOptimizer optimizer( objective );
	optimizer.setStart( start );
	optimizer.setBounds( lBound, uBound );

	// run optimizer
	VectorD result = optimizer.run();

	// store/display results
	objective.setVerbose();
	objective.eval( result );
	node.setBrightnessBounds( (float) result[ 0 ], (float) result[ 1 ] );
	disp( 1, "node: %d, %d, bright: %f, %f", node.xImageIndex(), node.yImageIndex(), result[ 0 ], result[ 1 ] );
}


/// stitch together the images for a single slide
void stitchSingleSlide( const String &inputPath, const Array<String> &fileList, int &imageIndex, 
					    const String &bgFileName, const String &outputPath, const String &outputFileName ) {

	// optimization parameters
	int inputScaleFactor = 1;
	int optScaleFactor = 4;
	int brightnessBase = 220;
	bool runGlobalOpt = true;
	bool runPerNodeOpt = true;

	// load acquisition parameters
	String fileName = inputPath + "/scan.conf";
	Config scanConf;
	scanConf.load( fileName );
	int xSizeMil = scanConf.readInt( "xSizeMil" );
	int ySizeMil = scanConf.readInt( "ySizeMil" );
	int xOverlapMil = scanConf.readInt( "xOverlapMil" );
	int yOverlapMil = scanConf.readInt( "yOverlapMil" );
	int xSubCount = scanConf.readInt( "xSubCount" );
	int ySubCount = scanConf.readInt( "ySubCount" );

	// load background image
	aptr<ImageGrayU> bgImage = load<ImageGrayU>( bgFileName );
	int origWidth = bgImage->width();
	int origHeight = bgImage->height();

	// compute output size
	int xStepMil = xSizeMil - xOverlapMil;
	int yStepMil = ySizeMil - yOverlapMil;
	int fullWidth = origWidth / inputScaleFactor;
	int fullHeight = origHeight / inputScaleFactor;
	int outputWidthMil = xStepMil * (xSubCount - 1 ) + xSizeMil;
	int outputHeightMil = yStepMil * (ySubCount - 1 ) + ySizeMil;
	int outputWidth = outputWidthMil * fullWidth / xSizeMil;
	int outputHeight = outputHeightMil * fullHeight / ySizeMil;
	disp( 1, "output size (mils): %d, %d, output size (pixels): %d, %d", 
		outputWidthMil, outputHeightMil, outputWidth, outputHeight );

	// resize bg image if needed
	if (inputScaleFactor > 1)
		bgImage = resize( *bgImage, fullWidth, fullHeight, true );
	aptr<ImageGrayU> bgImageSmall = resize( *bgImage, bgImage->width() / optScaleFactor, bgImage->height() / optScaleFactor, true );

	// prepare optimization parameters
	int xStep = xStepMil * fullWidth / xSizeMil;
	int yStep = yStepMil * fullHeight / ySizeMil;
	float xSkew = 0;
	float ySkew = 0;

	// collection of images for this slide
	StitchNodeSet nodeSet( xSubCount, ySubCount );

	// first pass: load stitch images
	for (int xImageIndex = 0; xImageIndex < xSubCount; xImageIndex++) {
		for (int yImageIndex = 0; yImageIndex < ySubCount; yImageIndex++) {

			// load image file (rotate 180 on load since camera is rotated)
			String imageName = fileList[ imageIndex ];
			String fileName = inputPath + "/" + imageName;
			aptr<ImageGrayU> input = load<ImageGrayU>( fileName );
			input = rotate180( *input );
			if (inputScaleFactor > 1) 
				input = resize( *input, fullWidth, fullHeight, true );
			aptr<ImageGrayU> inputSmall = resize( *input, input->width() / optScaleFactor, input->height() / optScaleFactor, true );
			disp( 1, "loaded: %s", fileName.c_str() );
			assertAlways( input->width() == fullWidth );
			assertAlways( input->height() == fullHeight );

			// compute normalized image
			int smallWidth = inputSmall->width();
			int smallHeight = inputSmall->height();
			for (int y = 0; y < smallHeight; y++) {
				for (int x = 0; x < smallWidth; x++) {
					int v = inputSmall->data( x, y ) + brightnessBase - bgImageSmall->data( x, y );
					inputSmall->data( x, y ) = bound( v, 0, 255 );
				}
			}

			// blur the small image a bit
			inputSmall = blurGauss( *inputSmall, 6 );

			// add to collection of nodes
			nodeSet.append( new StitchNode( inputSmall, xImageIndex, yImageIndex, fileName ) );

			// done with this image
			imageIndex++;

			// check for user cancel
			if (checkCommandEvents())
				break;
		}

		// check for user cancel
		if (checkCommandEvents())
			break;
	}

	// perform global optimization
	if (runGlobalOpt) {
		optimizeGlobalTransform( nodeSet, xStep, yStep, xSkew, ySkew, optScaleFactor );
		double xStepMilNew = (double) xStep * (double) xSizeMil / (double) fullWidth;
		double yStepMilNew = (double) yStep * (double) ySizeMil / (double) fullHeight;
		disp( 2, "new step (full pix): %d, %d, new step (mils): %f, %f", 
			xStep, yStep, xStepMilNew, yStepMilNew );
		disp( 2, "new size (mils): %d, %d", sbl::round( xStepMilNew + xOverlapMil ), sbl::round( yStepMilNew + yOverlapMil ) );
		disp( 2, "new skew: %f, %f", xSkew, ySkew );
	}
	updateTransformsGlobal( nodeSet, (float) xStep, (float) yStep, xSkew, ySkew, 1.0 / optScaleFactor );

	// optimize each node (except first)
	if (runPerNodeOpt) {
		for (int i = 1; i < nodeSet.count(); i++) {
			optimizeNode( nodeSet, nodeSet.node( i ) );

			// check for user cancel
			if (checkCommandEvents())
				break;
		}
	}

	// fix(later): replace with scale all if do local affine opt
	updateTransformsGlobal( nodeSet, (float) xStep, (float) yStep, xSkew, ySkew, 1.0 );

	// second pass: generate output image
	aptr<ImageGrayU> output( new ImageGrayU( outputWidth, outputHeight ) );
	output->clear( 255 );
	int nodeIndex = 0;
	for (int xImageIndex = 0; xImageIndex < xSubCount; xImageIndex++) {
		for (int yImageIndex = 0; yImageIndex < ySubCount; yImageIndex++) {
			StitchNode &node = nodeSet.node( xImageIndex, yImageIndex );
			nodeIndex++;

			// load image file
			String fileName = node.fileName();
			aptr<ImageGrayU> input = load<ImageGrayU>( fileName );
			input = rotate180( *input );
			if (inputScaleFactor > 1) 
				input = resize( *input, fullWidth, fullHeight, true );
			disp( 1, "loaded: %s", fileName.c_str() );
			assertAlways( input->width() == fullWidth );
			assertAlways( input->height() == fullHeight );

			// add to output
			for (int y = 0; y < fullHeight; y++) {
				for (int x = 0; x < fullWidth; x++) {
					Point2 outPoint = node.transform().mapForward( Point2( x, y ) );
					int xOut = sbl::round( outPoint.x );
					int yOut = sbl::round( outPoint.y );
					if (output->inBounds( xOut, yOut )) {
						int v = input->data( x, y ) + brightnessBase - bgImage->data( x, y );
						output->data( xOut, yOut ) = bound( sbl::round( node.mapBrightness( v ) ), 0, 255 );
					}
				}
			}

			// check for user cancel
			if (checkCommandEvents())
				return;
		}
	}

	// save the result for this slide
	output = removeBorders( *output );
	saveImage( *output, outputPath + "/" + outputFileName );

	// save small version
	output = resize( *output, output->width() / 16, output->height() / 16, true );
	saveImage( *output, outputPath + "/small/" + outputFileName );
}


//-------------------------------------------
// COMMANDS
//-------------------------------------------


/// for each histology slice, stitches together a set of sub-images into a single large image
void stitchHistoImages( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String bgFileName = addDataPath( conf.readString( "bgFileName" ) );
	String outputPath = addDataPath( conf.readString( "outputPath" ) );
	String outputPrefix = conf.readString( "outputPrefix", "" );

	// load acquisition parameters
	String fileName = inputPath + "/scan.conf";
	Config scanConf;
	scanConf.load( fileName );
	if (scanConf.entryCount() == 0) {
		warning( "unable to load scanning configuration: %s", fileName.c_str() );
		return;
	}
	int xSubCount = scanConf.readInt( "xSubCount" );
	int ySubCount = scanConf.readInt( "ySubCount" );

	// open slide ID file
	fileName = inputPath + "/slides.txt";
	File file( fileName, FILE_READ, FILE_TEXT );
	if (file.openSuccess() == false) {
		warning( "unable to open slide ID file: %s", fileName.c_str() );
		return;
	}

	// read slide ID file
	Array< Array<String> > slideIds;
	int xSlideCount = 0;
	while (file.endOfFile() == false) {
		Array<String> lineSplit = file.readLine().split( "," );
		if (lineSplit.count()) {
			Array<String> ids;
			for (int i = 0; i < lineSplit.count(); i++) {
				String id = lineSplit[ i ].strip();
				if (id.length())
					ids.appendCopy( id );
			}
			slideIds.appendCopy( ids );
			if (xSlideCount == 0)
				xSlideCount = ids.count();
			if (xSlideCount != ids.count()) {
				warning( "must have same number of slides on each line" );
				return;
			}
		}
	}
	int ySlideCount = slideIds.count();
	disp( 1, "xSlideCount: %d, ySlideCount: %d", xSlideCount, ySlideCount );

	// get list of input files
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );
	int imageCount = xSlideCount * ySlideCount * xSubCount * ySubCount;
	if (fileList.count() != imageCount) {
		warning( "found %d images; expected %d images", fileList.count(), imageCount );
		return;
	}

	// make sure output path exists
	createDir( outputPath );
	createDir( outputPath + "/small" );

	// loop over slides
	int imageIndex = 0;
	for (int x = 0; x < xSlideCount; x++) {
		for (int y = 0; y < ySlideCount; y++) {
			String outputFileName;
			if (outputPrefix.length())
				outputFileName += outputPrefix + "_";
			outputFileName += slideIds[ ySlideCount - y - 1 ][ x ] + ".png";
			stitchSingleSlide( inputPath, fileList, imageIndex, bgFileName, outputPath, outputFileName );

			// check for user cancel
			if (checkCommandEvents())
				return;
		}
	}
}


/// compute a background image from set of images of empty glass
void computeHistoStitchBackground( Config &conf ) {

	// get command parameters
	String inputPath = addDataPath( conf.readString( "inputPath" ) );
	String outputFileName = addDataPath( conf.readString( "outputFileName" ) );
	int scaleFactor = 8;

	// get list of input files
	Array<String> fileList = dirFileList( inputPath, "", ".JPG" );

	// sum of inputs
	aptr<ImageGrayI> inputSum;
	int inputSumCount = 0;
	int fullWidth = 0, fullHeight = 0;

	// loop over images
	for (int imageIndex = 0; imageIndex < fileList.count(); imageIndex++) {

		// load input image
		String fileName = inputPath + "/" + fileList[ imageIndex ];
		disp( 1, "%s", fileName.c_str() );
		aptr<ImageGrayU> input = load<ImageGrayU>( fileName );

		// shrink and blur
		// note: we use box blur because gauss blur makes image brighter
		fullWidth = input->width();
		fullHeight = input->height();
		int width = input->width() / scaleFactor;
		int height = input->height() / scaleFactor;
		input = blurBox( *input, 11 );
		input = resize( *input, width, height, true );
		input = blurBox( *input, 7 );

		// allocate sum image if needed
		if (inputSum.get() == NULL) {
			inputSum.reset( new ImageGrayI( width, height ) );
			inputSum->clear( 0 );
		}

		// update sum
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				inputSum->data( x, y ) += input->data( x, y );
			}
		}
		inputSumCount++;

		// check for user cancel
		if (checkCommandEvents())
			return;
	}

	// generate result
	if (inputSum.get()) {

		// compute mean
		int width = inputSum->width(), height = inputSum->height();
		aptr<ImageGrayU> meanImage( new ImageGrayU( width, height ) );
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
				meanImage->data( x, y ) = inputSum->data( x, y ) / inputSumCount;

		// expand and blur image
		// note: we use box blur because gauss blur makes image brighter
		meanImage = blurBox( *meanImage, 5 );
		meanImage = resize( *meanImage, fullWidth, fullHeight, true );
		meanImage = blurBox( *meanImage, 11 );

		// save image
		saveImage( *meanImage, outputFileName );
	}
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initHistoStitch() {
	registerCommand( "hstitch", stitchHistoImages );
	registerCommand( "hstitchbg", computeHistoStitchBackground );
}


} // end namespace hb
