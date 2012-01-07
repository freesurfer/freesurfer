#include "registration/BlockRegister.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageDraw.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/image/ImageSeqUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/Optimizer.h>
#include <sbl/math/Geometry.h>
#include "prep/VolumeFile.h"
#include "registration/VarCorres3D.h"
#include "registration/Normalization.h"
using namespace sbl;
namespace hb {


//-------------------------------------------
// LINEAR REGISTRATION
//-------------------------------------------


// a class used to evaluate a block-face registration
class BlockRegisterData : public Objective {
public:

	// basic constructor
	BlockRegisterData( const String &outputPath, float blockSlicePerMrSlice, int xCrop, int yCrop ) {
		m_saveResults = false;
		m_outputPath = outputPath;
		m_xCrop = xCrop;
		m_yCrop = yCrop;
		m_blockSlicePerMrSlice = blockSlicePerMrSlice;
	}

	/// an objective function for evaluating a linear block-face-to-MRI registration
	double eval( const VectorD &params );

	// access registration data
	Array<ImageGrayU> &blockNormImages() { return m_blockNormImages; }
	Array<ImageGrayU> &blockOrigImages() { return m_blockOrigImages; }
	Array<ImageGrayU> &mrNormImages() { return m_mrNormImages; }
	Array<ImageGrayU> &mrOrigImages() { return m_mrOrigImages; }
	const Array<ImageGrayU> &blockNormImages() const { return m_blockNormImages; }
	const Array<ImageGrayU> &blockOrigImages() const { return m_blockOrigImages; }
	const Array<ImageGrayU> &mrNormImages() const { return m_mrNormImages; }
	const Array<ImageGrayU> &mrOrigImages() const { return m_mrOrigImages; }
	inline const String &outputPath() const { return m_outputPath; }
	inline void setSaveResults( bool saveResults ) { m_saveResults = saveResults; }

private:

	// registration data (public for simplicity)
	Array<ImageGrayU> m_blockNormImages;
	Array<ImageGrayU> m_blockOrigImages;
	Array<ImageGrayU> m_mrNormImages;
	Array<ImageGrayU> m_mrOrigImages;
	bool m_saveResults;
	String m_outputPath;
	float m_blockSlicePerMrSlice;
	int m_xCrop;
	int m_yCrop;
};


/// an objective function for evaluating a linear block-face-to-MRI registration
double BlockRegisterData::eval( const VectorD &params ) {
	int step = 2;
	int xCrop = m_xCrop;
	int yCrop = m_yCrop;
	double zBoundPenalty = 100.0;
	AffineTransform3 transform( params );
	if (params.length() == 3)
		transform.setDiag( 1.0, 1.0, 1.0 / m_blockSlicePerMrSlice );
	double sumDiff = 0;
	int count = 0;
	if (m_saveResults) {
		step = 1;
		disp( 1, "saving results..." );
	}
	
	// loop over block face images
	for (int bz = 0; bz < m_blockNormImages.count(); bz++) {

		// get block-face slice for quick reference
		ImageGrayU &bNorm = m_blockNormImages[ bz ];
		ImageGrayU &bOrig = m_blockOrigImages[ bz ];
		int bWidth = bNorm.width(), bHeight = bNorm.height();
	
		// prepare visualization images if needed
		aptr<ImageGrayU> mNormProj, mOrigProj;
		if (m_saveResults) {
			mNormProj.reset( new ImageGrayU( bWidth, bHeight ) );
			mNormProj->clear( 0 );
			mOrigProj.reset( new ImageGrayU( bWidth, bHeight ) );
			mOrigProj->clear( 0 );
		}

		// loop over block-face slice
		for (int by = yCrop; by < bHeight - yCrop; by += step) {
			for (int bx = xCrop; bx < bWidth - xCrop; bx += step) {
				Point3 bPt( bx, by, bz );
				Point3 mPt = transform.transform( bPt );

				// check whether z coordinate falls within MR volume
				int mz = sbl::round( mPt.z );
				if (mz < 0 || mz >= m_mrNormImages.count()) {
					sumDiff += zBoundPenalty;
					count++;
				} else {
					ImageGrayU &mNorm = m_mrNormImages[ mz ];
					int mx = sbl::round( mPt.x );
					int my = sbl::round( mPt.y );
					if (mNorm.inBounds( mx, my )) {
						int diff = bNorm.data( bx, by ) - mNorm.data( mx, my );
						if (diff < 0)
							diff = -diff;
						sumDiff += diff;
						count++;

						// if requested, create a visualization of the registration
						if (m_saveResults) {
							if (mPt.z > 0) {
								int vNorm = round( interp( m_mrNormImages, (float) mPt.x, (float) mPt.y, (float) mPt.z ) );
								int vOrig = round( interp( m_mrOrigImages, (float) mPt.x, (float) mPt.y, (float) mPt.z ) );
								if (bx + step <= bWidth && by + step <= bHeight) {
									for (int y = by; y < by + step; y++) {
										for (int x = bx; x < bx + step; x++) {
											mNormProj->data( x, y ) = vNorm;
											mOrigProj->data( x, y ) = vOrig;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		// if requested, save a visualization of the registration
		if (m_saveResults) {
			ImageColorU vis( bWidth, bHeight );
			vis.clear( 255, 255, 255 );
			drawMaskBoundary( vis, bNorm, 128, 255, 0, 0 );
			drawMaskBoundary( vis, *mNormProj, 128, 0, 0, 255 );
			saveImage( vis, m_outputPath + "/vis/" + sprintF( "%04d.png", bz ) );
			saveImage( *mOrigProj, m_outputPath + "/" + sprintF( "%04d.png", bz ) );
			aptr<ImageGrayU> both = joinHoriz( bOrig, *mOrigProj );
			saveImage( *both, m_outputPath + "/vis/" + sprintF( "both_%04d.png", bz ) );
		}
	}
	double score = sumDiff / (double) count;
	disp( 1, "t: %3.1f, %3.1f, %3.1f, s: %4.2f, %4.2f, %4.2f, score: %f",
		transform.b.x, transform.b.y, transform.b.z, 
		transform.a.data[ 0 ][ 0 ], transform.a.data[ 1 ][ 1 ], transform.a.data[ 2 ][ 2 ], 
		score );
	if (m_saveResults) {
		disp( 1, "a[ 0 ] = { %4.2f, %4.2f, %4.2f }", 
			transform.a.data[ 0 ][ 0 ], transform.a.data[ 0 ][ 1 ], transform.a.data[ 0 ][ 2 ] );
		disp( 1, "a[ 1 ] = { %4.2f, %4.2f, %4.2f }", 
			transform.a.data[ 1 ][ 0 ], transform.a.data[ 1 ][ 1 ], transform.a.data[ 1 ][ 2 ] );
		disp( 1, "a[ 2 ] = { %4.2f, %4.2f, %4.2f }", 
			transform.a.data[ 2 ][ 0 ], transform.a.data[ 2 ][ 1 ], transform.a.data[ 2 ][ 2 ] );
		File logFile( m_outputPath + "/vis/log.txt", FILE_APPEND, FILE_TEXT );
		if (logFile.openSuccess()) {
			logFile.writeF( "t: %3.1f, %3.1f, %3.1f, a: [%4.2f, %4.2f, %4.2f], [%4.2f, %4.2f, %4.2f], [%4.2f, %4.2f, %4.2f], score: %f\n",
				transform.b.x, transform.b.y, transform.b.z, 
				transform.a.data[ 0 ][ 0 ], transform.a.data[ 0 ][ 1 ], transform.a.data[ 0 ][ 2 ],
				transform.a.data[ 1 ][ 0 ], transform.a.data[ 1 ][ 1 ], transform.a.data[ 1 ][ 2 ],
				transform.a.data[ 2 ][ 0 ], transform.a.data[ 2 ][ 1 ], transform.a.data[ 2 ][ 2 ],
				score );
		}
		File paramFile( m_outputPath + "/linearTransform.txt", FILE_WRITE, FILE_TEXT );
		if (paramFile.openSuccess()) 
			saveTransform( paramFile, transform );
	}
	return score;
}


/// compute a linear registration of a block-face volume with an MR volume
void registerBlockToMrLinear( BlockRegisterData &brd ) {
	
	// parameters
	double translateBound = 400.0f;
	double affineBound = 0.3;
	bool runOpt = true;
	int iterCount = 4;

	// create optimization bounds and start
	AffineTransform3 start;
	start.setTranslation( 0, 0, 0 );
	start.setDiag( 1.0, 1.0, 0.5 ); // assuming 2 block slices for each one MR slice
	VectorD params = start.toVector();

	// perform additional affine iterations
	for (int iter = 0; iter < iterCount; iter++) {
		VectorD lBound = params;
		VectorD uBound = params;
		for (int i = 0; i < 3; i++) {
			lBound[ i ] -= translateBound;
			uBound[ i ] += translateBound;
		}
		for (int i = 3; i < 12; i++) {
			lBound[ i ] -= affineBound;
			uBound[ i ] += affineBound;
		}

		// prepare optimizer
		SimplexOptimizer opt( brd );
		opt.setPenaltyFactor( 1000.0 );
		opt.setStart( params );
		opt.setBounds( lBound, uBound );
		opt.setFinalTolerance( 0.0001 );

		// run optimizer
		if (runOpt) 
			params = opt.run();

		// save/visualize registration parameters and quality
		brd.setSaveResults( true );
		brd.eval( params );
		brd.setSaveResults( false );
	}
}


//-------------------------------------------
// NON-LINEAR REGISTRATION
//-------------------------------------------


/// estimate a non-linear registration of a block-face volume with an MR volume, assuming a linear registration has already been completed
void registerBlockToMrNonLinear( const BlockRegisterData &brd ) {

	// estimate correspondences
	const Array<ImageGrayU> &srcNormSeq = brd.blockNormImages();
	const Array<ImageGrayU> &destNormSeq = brd.mrNormImages();
	const Array<ImageGrayU> &srcOrigSeq = brd.blockOrigImages();
	const Array<ImageGrayU> &destOrigSeq = brd.mrOrigImages();
	Array<CorresField3D> cfSeq;
	disp( 1, "srcNormSeq: %d, destNormSeq: %d, srcOrigSeq: %d, destOrigSeq: %d",
		srcNormSeq.count(), destNormSeq.count(), srcOrigSeq.count(), destOrigSeq.count() );
	varCorres3D( srcNormSeq, destNormSeq, cfSeq, 5 );

	// use correspondences to map destination sequance (MR) back to src sequence (blockface)
	int depth = cfSeq.count();
	int width = cfSeq[ 0 ].width(), height = cfSeq[ 0 ].height();
	Array<ImageGrayU> destNormSeqMapped;
	initImageSeq( destNormSeqMapped, width, height, depth, false, 0 );
	mapBack( cfSeq, destNormSeq, destNormSeqMapped );
	Array<ImageGrayU> destOrigSeqMapped;
	initImageSeq( destOrigSeqMapped, width, height, depth, false, 0 );
	mapBack( cfSeq, destOrigSeq, destOrigSeqMapped );

	// save correspondence sequence
	saveCorresSeq( cfSeq, brd.outputPath() + "/cfSeq.cfs" );

	// evaluate and visualize results
	double sumDiff = 0;
	for (int z = 0; z < depth; z++) {

		// plot source vs dest projected via estimated correspondences
		ImageColorU vis( width, height );
		vis.clear( 255, 255, 255 );
		drawMaskBoundary( vis, srcNormSeq[ z ], 128, 0, 0, 255 );
		drawMaskBoundary( vis, destNormSeqMapped[ z ], 128, 255, 0, 0 );
		saveImage( vis, brd.outputPath() + sprintF( "/vis/%04d.png", z ) );
		saveImage( destOrigSeqMapped[ z ], brd.outputPath() + sprintF( "/%04d.png", z ) );

		// save some sample images
		if (z == depth / 2 ) {
			saveImage( srcNormSeq[ z ], brd.outputPath() + "/vis/srcNormSeq.png" );
			saveImage( destNormSeq[ z ], brd.outputPath() + "/vis/destNormSeq.png" );
			saveImage( destNormSeqMapped[ z ], brd.outputPath() + "/vis/destNormSeqMapped.png" );
		}

		// compute matching quality
		double diff = meanAbsDiff( srcNormSeq[ z ], destNormSeqMapped[ z ], 0, 0 );
		sumDiff += diff;
	}
	double meanDiff = sumDiff / (double) cfSeq.count();
	disp( 1, "mean diff: %f", meanDiff );
}


/// blur and normalize an image volume
void normalizeVolume( const Array<ImageGrayU> &origImages, const Array<ImageGrayU> &masks, Array<ImageGrayU> &normImages, const String &visPrefix ) {

	// check inputs
	assertAlways( origImages.count() );
	assertAlways( origImages.count() == masks.count() );

	// blur the masks
	Array<ImageGrayU> blurredMasks;
	blurGaussSeqZ( masks, blurredMasks, 5.0f );
	blurGaussSeqXY( blurredMasks, 5.0f );

	// blur the images
	Array<ImageGrayU> blurredOrigImages;
	blurGaussSeqZ( origImages, blurredOrigImages, 1.0f );
	blurGaussSeqXY( blurredOrigImages, 1.0f );

	// normalize each slice
	for (int i = 0; i < origImages.count(); i++) {

		// compute normalized image
		aptr<ImageGrayU> normImage = normalize( blurredOrigImages[ i ], masks[ i ], blurredMasks[ i ] );
		normImages.append( normImage.release() );
		
		// save diagnostics
		if (i == origImages.count() / 2) {
			saveImage( origImages[ i ], visPrefix + ".orig.png" );
			saveImage( blurredOrigImages[ i ], visPrefix + ".origBlur.png" );
			saveImage( masks[ i ], visPrefix + ".mask.png" );
			saveImage( blurredMasks[ i ], visPrefix + ".maskBlur.png" );
			saveImage( normImages[ i ], visPrefix + ".norm.png" );
		}

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


//-------------------------------------------
// COMMANDS
//-------------------------------------------


/// create mgz files for the MR and blockface volumes in order to validate the registration
void createRegistrationVolumes( Config &conf ) {

	// get command parameters
	String blockCropPath = addDataPath( conf.readString( "blockCropPath", "blockface/crop" ) );
	String blockSegPath = addDataPath( conf.readString( "blockSegPath", "blockface/seg" ) );
	String mrRawPath = addDataPath( conf.readString( "mrRawPath", "mri/rawFlash20" ) );
	String mrRegLinPath = addDataPath( conf.readString( "mrRegLinPath", "mri/regLin" ) );
	String mrRegPath = addDataPath( conf.readString( "mrRegPath", "mri/reg" ) );
	String histoRegPath = addDataPath( conf.readString( "histoRegPath", "histo/reg" ) );

	// load spacing from MR file
	String spacingFileName = mrRawPath + "/spacing.txt";
	File spacingFile( spacingFileName, FILE_READ, FILE_TEXT );
	if (spacingFile.openSuccess() == false) {
		warning( "unable to open spacing file: %s", spacingFileName.c_str() );
		return;
	}
	String line = spacingFile.readLine();
	Array<String> lineSplit = line.split( "," );
	if (lineSplit.count() != 3) {
		warning( "invalid spacing file: %s", spacingFileName.c_str() );
		return;
	}
	float xSpacing = lineSplit[ 0 ].strip().toFloat();
	float ySpacing = lineSplit[ 1 ].strip().toFloat();
	float zSpacing = lineSplit[ 2 ].strip().toFloat();
	if (xSpacing <= 0 || ySpacing <= 0 || zSpacing <= 0) {
		warning( "invalid spacing file: %s", spacingFileName.c_str() );
		return;
	}

	// load linear MR/block transformation
	String transformFileName = mrRegLinPath + "/linearTransform.txt";
	bool useTransform = fileExists( transformFileName );
	if (useTransform) {
		AffineTransform3 transform;
		File transformFile( transformFileName, FILE_READ, FILE_TEXT );
		if (transformFile.openSuccess()) {
			transform = loadTransform( transformFile );
		} else {
			warning( "invalid transform file: %s", transformFileName.c_str() );
			return;
		}

		// compute new spacing
		if (dAbs( transform.a.data[ 0 ][ 0 ] ) < 1e-6 || dAbs( transform.a.data[ 1 ][ 1 ] ) < 1e-6 || dAbs( transform.a.data[ 2 ][ 2 ] ) < 1e-6) {
			warning( "invalid transform: %f, %f, %f", transform.a.data[ 0 ][ 0 ], transform.a.data[ 1 ][ 1 ], transform.a.data[ 2 ][ 2 ] );
			return;
		}
		xSpacing /= (float) transform.a.data[ 0 ][ 0 ];
		ySpacing /= (float) transform.a.data[ 1 ][ 1 ];
		zSpacing /= (float) transform.a.data[ 2 ][ 2 ];
	}

	// convert files
	if (useTransform) {
		convertImagesToMghFile( blockCropPath, dataPath() + "blockface/blockCrop.mgz", xSpacing, ySpacing, zSpacing );
		convertImagesToMghFile( blockSegPath, dataPath() + "blockface/blockSeg.mgz", xSpacing, ySpacing, zSpacing );
		convertImagesToMghFile( mrRegLinPath, dataPath() + "blockface/mrRegLin.mgz", xSpacing, ySpacing, zSpacing );
		convertImagesToMghFile( mrRegPath, dataPath() + "blockface/mrReg.mgz", xSpacing, ySpacing, zSpacing );
	} else {
		convertImagesToMghFile( blockCropPath, dataPath() + "blockface/blockCropNoTransform.mgz", xSpacing, ySpacing, zSpacing );
		convertImagesToMghFile( blockSegPath, dataPath() + "blockface/blockSegNoTransform.mgz", xSpacing, ySpacing, zSpacing );
	}

	// we don't handle this case currently; need to add code to handle missing histo slices
	//convertImagesToMghFile( histoRegPath, dataPath() + "blockface/histoReg.mgz", xSpacing, ySpacing, zSpacing, width, height );
}


/// perform a linear or non-linear registration of a block-face volume with an MR volume
void registerBlockToMrOnePass( Config &conf ) {

	// get command parameters
	float blockSlicePerMrSlice = conf.readFloat( "blockSlicePerMrSlice", 0.5f );
	String blockCropPath = addDataPath( conf.readString( "blockCropPath", "blockface/crop" ) );
	String blockSegPath = addDataPath( conf.readString( "blockSegPath", "blockface/seg" ) );
	String mrPath = addDataPath( conf.readString( "mrPath" ) );
	String outputPath = addDataPath( conf.readString( "outputPath" ) );
	bool linear = conf.readBool( "linear" );

	// make sure output path exists
	createDir( outputPath );
	createDir( outputPath + "/vis" );

	// make sure output directories are empty
	disp( 1, "removing existing images from output path: %s", outputPath.c_str() );
	String command = sprintF( "rm -f %s/*.png", outputPath.c_str() );
	system( command.c_str() );
	command = sprintF( "rm -f %s/vis/*.png", outputPath.c_str() );
	system( command.c_str() );

	// data used for registration
	BlockRegisterData brd( outputPath, blockSlicePerMrSlice, 0, 0 );

	// load block face images
	// we assume that bprep has made the blockface images roughly the same size as the MR images
	Array<String> blockFileList = dirFileList( blockCropPath, "", ".png" );
	Array<ImageGrayU> blockMaskImages;
	disp( 1, "loading %d blockface images", blockFileList.count() );
	for (int i = 0; i < blockFileList.count(); i++) {
		aptr<ImageGrayU> origInput = load<ImageGrayU>( blockCropPath + "/" + blockFileList[ i ] );
		aptr<ImageGrayU> segInput = load<ImageGrayU>( blockSegPath + "/" + blockFileList[ i ] );
		aptr<ImageGrayU> segMask = threshold( *segInput, 128, false );
		blockMaskImages.append( segMask.release() );
		brd.blockOrigImages().append( origInput.release() );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	if (checkCommandEvents())
		return;

	// load segmented MR images
	Array<String> mrFileList = dirFileList( mrPath, "", ".png" );
	Array<ImageGrayU> mrMaskImages;
	disp( 1, "loading %d MR images", mrFileList.count() );
	for (int i = 0; i < mrFileList.count(); i++) {
		aptr<ImageGrayU> origInput = load<ImageGrayU>( mrPath + "/" + mrFileList[ i ] );
		aptr<ImageGrayU> segMask = threshold( *origInput, 1, false );
		mrMaskImages.append( segMask.release() );
		brd.mrOrigImages().append( origInput.release() );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
	if (checkCommandEvents())
		return;

	// display diagnostics
	disp( 1, "block images: %d, MR images: %d", brd.blockOrigImages().count(), brd.mrOrigImages().count() ); 
	if (brd.blockOrigImages().count() == 0 || brd.mrOrigImages().count() == 0)
		return;

	// prepare the volumes
	normalizeVolume( brd.blockOrigImages(), blockMaskImages, brd.blockNormImages(), outputPath + "/vis/block" );
	if (checkCommandEvents())
		return;
	normalizeVolume( brd.mrOrigImages(), mrMaskImages, brd.mrNormImages(), outputPath + "/vis/mr" );
	if (checkCommandEvents())
		return;

	// perform registration
	if (linear) {
		registerBlockToMrLinear( brd );
	} else {
		registerBlockToMrNonLinear( brd );
	}
}


/// perform a linear then non-linear registration of a block-face volume with an MR volume
void registerBlockToMr( Config &conf ) {

	// first pass: linear registration
	conf.writeString( "mrPath", "mri/seg" );
	conf.writeString( "outputPath", "mri/regLin" );
	conf.writeBool( "linear", true );
	registerBlockToMrOnePass( conf );

	// second pass: non-linear registration
	conf.writeString( "mrPath", "mri/regLin" );
	conf.writeString( "outputPath", "mri/reg" );
	conf.writeBool( "linear", false );
	registerBlockToMrOnePass( conf );
}


/// create a table showing the mapping between histology images and blockface images
void createHistologyBlockIndexTable( Config &conf ) {

	// get command parameters
	// blockOffset is index of blockface image corresponding to slice #1 (not necessarily the first histology image)
	int blockOffset = conf.readInt( "blockOffset" );
	String blockPath = addDataPath( conf.readString( "blockPath", "blockface/raw" ) );
	String outputFileName = addDataPath( conf.readString( "outputFileName", "blockface/mapping.txt" ) );

	// get input file lists
	Array<String> bFileList = dirFileList( blockPath, "", ".JPG" );
	disp( 1, "blockface file count: %d", bFileList.count() );

	// open output file
	File outFile( outputFileName, FILE_WRITE, FILE_TEXT );
	if (outFile.openSuccess() == false) {
		warning( "unable to open output file: %s", outputFileName.c_str() );
		return;
	}
	outFile.writeF( "histo_slice_number, blockface_image\n" );

	// output histology slice index for each block face file
	for (int bIndex = 0; bIndex < bFileList.count(); bIndex++) {
		int hIndex = bIndex - blockOffset;
		outFile.writeF( "%d, %s\n", hIndex, bFileList[ bIndex ].c_str() );
	}
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initBlockRegister() {
	registerCommand( "rvol", createRegistrationVolumes );
	registerCommand( "bregpass", registerBlockToMrOnePass );
	registerCommand( "breg", registerBlockToMr );
	registerCommand( "btable", createHistologyBlockIndexTable );
}


} // end namespace hb
