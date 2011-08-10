#include "prep/HistoPrep.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#include <sbl/image/ImageSeqUtil.h>
#include "registration/VarCorres3D.h"
#include "registration/VarCorres3DUtil.h"
using namespace sbl;
namespace hb {


// transpose Y axis of correspondence field with Z axis of corresponce field
void transposeYZ( const Array<CorresField3D> &cfSeqInput, Array<CorresField3D> &cfSeqOutput ) {
	int length = cfSeqInput.count();
	assertAlways( length );
	int width = cfSeqInput[ 0 ].width(), height = cfSeqInput[ 0 ].height();
	assertAlways( cfSeqOutput.count() == 0 );
	for (int y = 0; y < height; y++) {
		CorresField3D *cfOut = new CorresField3D( width, length );
		for (int z = 0; z < length; z++) {
			const CorresField3D &cfIn = cfSeqInput[ z ];
			for (int x = 0; x < width; x++) {
				cfOut->u( x, z ) = cfIn.u( x, y );
				cfOut->v( x, z ) = cfIn.w( x, y ); // note: transposing w and v
				cfOut->w( x, z ) = cfIn.v( x, y );
			}
		}
		cfSeqOutput.append( cfOut );
	}
}


/// create a synthetic image volume for testing the 3D correspondence estimation algorithm
void initTestSequence( Array<ImageGrayU> &seq, int width, int height, int length ) {
	bool flat = (width <= 2 || height <= 2 || length <= 2);
	float smallBlurSigma = 1;
	float largeBlurSigma = 2;
	int splatCount = flat ? 10 : 30;
	int border = 2;
	int radius = 2;
	int smallFactor = 2;
	int smallWidth = width > 2 ? width / smallFactor : width;
	int smallHeight = height > 2 ? height / smallFactor : height;
	int smallLength = length > 2 ? length / smallFactor : length;
	Array<ImageGrayU> seqSmall;
	initImageSeq( seqSmall, smallWidth, smallHeight, smallLength, true, 0 );
	for (int i = 0; i < splatCount; i++) {
		int xCent = width > 2 ? randomInt( radius + border, smallWidth - radius - border - 1 ) : 0;
		int yCent = height > 2 ? randomInt( radius + border, smallHeight - radius - border - 1 ) : 0;
		int zCent = length > 2 ? randomInt( radius + border, smallLength - radius - border - 1 ) : 0;
		int xMin = bound( xCent - radius, 0, width - 1 );
		int xMax = bound( xCent + radius, 0, width - 1 );
		int yMin = bound( yCent - radius, 0, height - 1 );
		int yMax = bound( yCent + radius, 0, height - 1 );
		int zMin = bound( zCent - radius, 0, length - 1 );
		int zMax = bound( zCent + radius, 0, length - 1 );
		for (int z = zMin; z <= zMax; z++) {
			for (int y = yMin; y <= yMax; y++) {
				for (int x = xMin; x <= xMax; x++) {
					seqSmall[ z ].data( x, y ) = 255;
				}
			}
		}
	}
	Array<ImageGrayU> blurSeqSmall, seqLarge;
	blurGaussSeqXY( seqSmall, smallBlurSigma );
	blurGaussSeqZ( seqSmall, blurSeqSmall, smallBlurSigma );
	resizeSeqXY( blurSeqSmall, width, height );
	resizeSeqZ( blurSeqSmall, seqLarge, length );
	blurGaussSeqXY( seqLarge, largeBlurSigma );
	blurGaussSeqZ( seqLarge, seq, largeBlurSigma );
}


/// create component (u, v, or w) of a synthetic 3D corresponence volume
void initCorresComponent( Array<ImageGrayF> &seq, int width, int height, int length ) {
	bool flat = (width <= 2 || height <= 2 || length <= 2);
	float smallBlurSigma = 3;
	float largeBlurSigma = 4;
	float maxMag = 6.0f;
	int splatCount = flat ? 10 : 30;
	int border = 2;
	int radius = 3;
	int smallFactor = 2;
	int smallWidth = width > 2 ? width / smallFactor : width;
	int smallHeight = height > 2 ? height / smallFactor : height;
	int smallLength = length > 2 ? length / smallFactor : length;
	Array<ImageGrayF> seqSmall;
	initImageSeq( seqSmall, smallWidth, smallHeight, smallLength, true, 0 );
	for (int i = 0; i < splatCount; i++) {
		int xCent = width > 2 ? randomInt( radius + border, smallWidth - radius - border - 1 ) : 0;
		int yCent = height > 2 ? randomInt( radius + border, smallHeight - radius - border - 1 ) : 0;
		int zCent = length > 2 ? randomInt( radius + border, smallLength - radius - border - 1 ) : 0;
		float v = randomFloat( -maxMag, maxMag );
		int xMin = bound( xCent - radius, 0, width - 1 );
		int xMax = bound( xCent + radius, 0, width - 1 );
		int yMin = bound( yCent - radius, 0, height - 1 );
		int yMax = bound( yCent + radius, 0, height - 1 );
		int zMin = bound( zCent - radius, 0, length - 1 );
		int zMax = bound( zCent + radius, 0, length - 1 );
		for (int z = zMin; z <= zMax; z++) {
			for (int y = yMin; y <= yMax; y++) {
				for (int x = xMin; x <= xMax; x++) {
					seqSmall[ z ].data( x, y ) = v;
				}
			}
		}
	}
	Array<ImageGrayF> blurSeqSmall, seqLarge;
	blurGaussSeqXY( seqSmall, smallBlurSigma );
	blurGaussSeqZ( seqSmall, blurSeqSmall, smallBlurSigma );
	resizeSeqXY( blurSeqSmall, width, height );
	resizeSeqZ( blurSeqSmall, seqLarge, length );
	blurGaussSeqXY( seqLarge, largeBlurSigma );
	blurGaussSeqZ( seqLarge, seq, largeBlurSigma );
}


/// create a synthetic 3D corresponence volume for testing the 3D correspondence estimation algorithm
void initCorres( Array<CorresField3D> &cfSeq, int width, int height, int length ) {
	Array<ImageGrayF> uSeq, vSeq, wSeq;
	initCorresComponent( uSeq, width, height, length );
	initCorresComponent( vSeq, width, height, length );
	initCorresComponent( wSeq, width, height, length );
	for (int z = 0; z < length; z++) {
		CorresField3D *cf = new CorresField3D( width, height );
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				cf->u( x, y ) = width > 2 ? uSeq[ z ].data( x, y ) : 0;
				cf->v( x, y ) = height > 2 ? vSeq[ z ].data( x, y ) : 0;
				cf->w( x, y ) = length > 2 ? wSeq[ z ].data( x, y ) : 0;
			}
		}
		cfSeq.append( cf );
	}
}


/// test 3D correspondence estimation algorithm using synthetic data 
void testCorres3D( Config &conf ) {
	int width = 60;
	int height = 60;
	int length = 60;
	bool transpose = true;

	// create test image sequence
	Array<ImageGrayU> destSeq;
	initTestSequence( destSeq, width, height, length );
	assertAlways( destSeq.count() == length && destSeq[ 0 ].width() == width && destSeq[ 0 ].height() == height );

	// create test correspondence sequence
	Array<CorresField3D> cfSeq;
	initCorres( cfSeq, width, height, length );

	// create src sequence by projecting dest seq backward according to correspondences
	Array<ImageGrayU> srcSeq;
	initImageSeq( srcSeq, width, height, length, false, 0 );
	mapBack( cfSeq, destSeq, srcSeq );

	// estimate correspondences
	Array<CorresField3D> cfSeqEst;
	if (transpose) {
		Array<ImageGrayU> srcSeqTransposed, destSeqTransposed;
		transposeYZ( srcSeq, srcSeqTransposed );
		transposeYZ( destSeq, destSeqTransposed );
		Array<CorresField3D> cfSeqEstTransposed;
		varCorres3D( srcSeqTransposed, destSeqTransposed, cfSeqEstTransposed, 1 );
		transposeYZ( cfSeqEstTransposed, cfSeqEst );
	} else {
		varCorres3D( srcSeq, destSeq, cfSeqEst, 1 );
	}

	// map back according to estiamted correspondences
	Array<ImageGrayU> destSeqMapped;
	initImageSeq( destSeqMapped, width, height, length, false, 0 );
	mapBack( cfSeqEst, destSeq, destSeqMapped );

	// evaluate/visualize results
	String outputPath = dataPath() + "testCorres3D";
	double corresDiffSum = 0, imageDiffSum = 0;
	for (int z = 0; z < length; z++) {

		// plot estimated and actual correspondences
		aptr<ImageColorU> cfVis = cfSeq[ z ].colorize();
		aptr<ImageColorU> cfEstVis = cfSeqEst[ z ].colorize();
		cfVis = resize( *cfVis, width * 4, height * 4, true ); // make bigger just so easier to see
		cfEstVis = resize( *cfEstVis, width * 4, height * 4, true );
		saveImage( *cfVis, outputPath + sprintF( "/cf-%04d.png", z ) );
		saveImage( *cfEstVis, outputPath + sprintF( "/cf-%04d-e.png", z ) );

		// plot source vs dest projected via estimated correspondences
		ImageColorU vis( width, height );
		vis.clear( 255, 255, 255 );
		drawMaskBoundary( vis, srcSeq[ z ], 128, 0, 0, 255 );
		drawMaskBoundary( vis, destSeqMapped[ z ], 128, 255, 0, 0 );
		aptr<ImageColorU> visLarge = resize( vis, width * 4, height * 4, true ); // make bigger just so easier to see
		saveImage( *visLarge, outputPath + sprintF( "/proj-%04d.png", z ) );
		aptr<ImageGrayU> destLarge = resize( destSeq[ z ], width * 4, height * 4, true ); // make bigger just so easier to see
		saveImage( *destLarge, outputPath + sprintF( "/dest-%04d.png", z ) );

		// compute and display stats
		double corresDiff = meanAbsDiff( cfSeq[ z ], cfSeqEst[ z ] );
		double imageDiff = meanAbsDiff( srcSeq[ z ], destSeqMapped[ z ], 0, 0 );
//		disp( 1, "corres diff: %f, image diff: %f", corresDiff, imageDiff );
		corresDiffSum += corresDiff;
		imageDiffSum += imageDiff;
	}
	disp( 1, "mean corres diff: %f, mean image diff: %f",
		corresDiffSum / (double) length, imageDiffSum / (double) length );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initTestCorres3D() {
	registerCommand( "cortest", testCorres3D );
}


} // end namespacs hb
