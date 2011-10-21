#include "registration/HistoRegister.h"
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/Geometry.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageRegister.h>
#include <sbl/image/MotionFieldUtil.h>
#include <sbl/image/Video.h> // for eval visualization
#include <sbl/system/FileSystem.h>
#include <sbl/other/Plot.h>
#include <pvl/VarMotion.h>
#include <pvl/VarMotionUtil.h>
#include "registration/Normalization.h"
#include "registration/CorresField3D.h"
using namespace sbl;
using namespace pvl;
namespace hb {


/// create a color visualization of the difference between a pair of grayscale images
aptr<ImageColorU> colorizeDifference( const ImageGrayU &input1, const ImageGrayU &input2 ) {
	int width = input1.width(), height = input1.height();
	assertAlways( input2.width() == width && input2.height() == height );
	aptr<ImageColorU> output( new ImageColorU( width, height ) );
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int r = 255, g = 255, b = 255;
			int diff = input1.data( x, y ) - input2.data( x, y );
			if (diff < 0) {
				g = 255 + diff;
				b = 255 + diff;
			} else {
				r = 255 - diff;
				g = 255 - diff;
			}
			output->setRGB( x, y, r, g, b );
		}
	}
	return output;
}


/// register a single histology image with an MR slice
void registerHistologySlice( int hIndex, int blockOffset,
							 const String &hPath, const String &bPath, const String &mPath, String &outputPath,
							 const Array<String> &hFileList, const Array<String> &bFileList, const Array<String> &mFileList, 
							 OutputVideo &outputVideo, 
							 VectorD &mutInfo, VectorD &normDiff, VectorD &flowGrad ) {

	// parameters
	int linearParamCount = 6;
	float maskBlurSigma = 5.0f;
	int motionScaleFactor = 1;
	bool computeFlow = true;
	bool drawLinearBoundary = false;
	bool saveFlow = true;
	String varMotionConfigFileName = "varMotion.conf"; // note: using current directory

	// load histo input file
	String hFileName = hPath + "/" + hFileList[ hIndex ];
  // here load histo as gray image (probably averaging channels)
	//aptr<ImageGrayU> hImage = load<ImageGrayU>( hFileName );
	//int hWidth = hImage->width(), hHeight = hImage->height();
  // here only use one channel (green):
  aptr<ImageColorU> hImageColor = load<ImageColorU>( hFileName );
  int hWidth = hImageColor->width(), hHeight = hImageColor->height();
  aptr<ImageGrayU> hImage( new ImageGrayU( hWidth, hHeight ) );
  for (int y = 0; y < hHeight; y++) {
    for (int x = 0; x < hWidth; x++) {
      hImage->data( x, y ) = hImageColor->g( x, y );
    }
  }
  
	aptr<ImageGrayU> hMask = threshold( *hImage, 254, true );
	int hSliceIndex = hFileList[ hIndex ].leftOfLast( '.' ).toInt();

	// compute block/MR slice index (they are in correspondence at this point)
	int mIndex = hSliceIndex + blockOffset;
	if (mIndex < 0 || mIndex >= mFileList.count()) 
		return;

	// get corresponding MR slice
	String mFileName = mPath + "/" + mFileList[ mIndex ];
	aptr<ImageGrayU> mImage = load<ImageGrayU>( mFileName );
	int mWidth = mImage->width(), mHeight = mImage->height();
	aptr<ImageGrayU> mMask = threshold( *mImage, 1, false );

	// get block file name
	String bFileName = bFileList[ mIndex ];

	// display info
	disp( 1, "histo: %s, %d, %d", hFileName.c_str(), hWidth, hHeight );
	disp( 1, "block: %s", bFileName.c_str() );
	disp( 1, "mr: %s, %d, %d", mFileName.c_str(), mWidth, mHeight );

	// resize the histo image
	// fix(later): preserve aspect ratio: resize then pad/crop
	// note: if change this, also need to change HistoTransform
	hWidth = mWidth; // hWidth / scaleFactor;
	hHeight = mHeight; // hHeight / scaleFactor;
	hImage = resize( *hImage, hWidth, hHeight, true );
	hMask = resize( *hMask, hWidth, hHeight, true );
	aptr<ImageGrayU> hMaskBlurred = blurGauss( *hMask, maskBlurSigma );
	aptr<ImageGrayU> mMaskBlurred = blurGauss( *mMask, maskBlurSigma ); // note: we could use a volume blur for the MR data, but want it amount of blur to match histology (where we can't do a volume blur prior to alignment)

	// compute registration images
	aptr<ImageGrayU> hNorm = normalize( *hImage, *hMask, *hMaskBlurred );
	aptr<ImageGrayU> mNorm = normalize( *mImage, *mMask, *mMaskBlurred );

	// do linear alignment
	aptr<ImageTransform> transform = registerUsingImageTransform( *mNorm, *hNorm, linearParamCount, 1, 10, 10, 100 );
	aptr<ImageGrayU> hMaskLinear = transform->mapBackward( *hMask, mWidth, mHeight, 0 );
	aptr<ImageGrayU> hNormLinear = transform->mapBackward( *hNorm, mWidth, mHeight, 0 );
	aptr<ImageGrayU> hImageLinear = transform->mapBackward( *hImage, mWidth, mHeight, 255 );
	String transformFileName = outputPath + sprintF( "/transform_%d.dat", hSliceIndex );
	File transformFile( transformFileName, FILE_WRITE, FILE_BINARY );
	if (transformFile.openSuccess()) {
		transform->save( transformFile );
	}

	// compute optical flow
	aptr<MotionField> mf;
	aptr<ImageGrayU> hMaskFlow, hNormFlow, hImageFlow;
	if (computeFlow) {
		mf = varMotion( *mNorm, *hNormLinear, varMotionConfigFileName, NULL, NULL, motionScaleFactor );
		hMaskFlow = mf->mapBackward( *hMaskLinear, 255, 0.5f );
		hNormFlow = mf->mapBackward( *hNormLinear, 255, 0.5f );
		hImageFlow = mf->mapBackward( *hImageLinear, 255, 0.5f );
		if (saveFlow) {
			String outputFileName = outputPath + sprintF( "/flow_%d.mf", hSliceIndex );
			saveMotionField( *mf, outputFileName );
		}
	}

	// save diagnostic images
	ImageColorU vis( mWidth, mHeight );
	vis.clear( 255, 255, 255 );
	drawMaskBoundary( vis, *mMask, 128, 0, 0, 255, 10, 10 );
	if (drawLinearBoundary)
		drawMaskBoundary( vis, *hMaskLinear, 128, 0, 255, 0 );
	if (computeFlow)
		drawMaskBoundary( vis, *hMaskFlow, 128, 255, 0, 0, 10, 10 );
	String name = mFileList[ mIndex ].leftOfLast( '.' );
	saveImage( vis, outputPath + "/vis/reg_" + name + "_boundaries.png" );
	if (computeFlow) {
		saveImage( *hImageFlow, outputPath + "/vis/reg_" + name + "_hOrigFlow.png" );
		saveImage( *hNormFlow, outputPath + "/vis/reg_" + name + "_hNormFlow.png" );
	}
	saveImage( *mNorm, outputPath + "/vis/reg_" + name + "_mNorm.png" );
	saveImage( *hNorm, outputPath + "/vis/reg_" + name + "_hNorm.png" );
	saveImage( *mImage, outputPath + "/vis/reg_" + name + "_mOrig.png" );
	saveImage( *hImage, outputPath + "/vis/reg_" + name + "_hOrig.png" );

	// only do eval if have flow
	if (computeFlow) {

		// update eval video
		aptr<ImageColorU> mfVis = colorizeMotion( *mf );
		aptr<ImageColorU> diff = colorizeDifference( *hNormFlow, *mNorm );
		aptr<ImageColorU> vis = joinHoriz( *mfVis, *diff );
		int visWidth = (vis->width() / 8) * 8;
		int visHeight = vis->height();
		vis = resize( *vis, visWidth, visHeight, false );
		if (outputVideo.openSuccess() == false) 
			outputVideo.open( outputPath + "/vis/eval.avi", vis->width(), vis->height(), 6 );
		outputVideo.append( *vis );

		// store eval stats
		int xBorder = 20, yBorder = 20;
		int bucketCount = 50;
		mutInfo.append( mutualInfo( *hImageFlow, *mImage, xBorder, yBorder, bucketCount ) );
		normDiff.append( meanAbsDiff( *hNormFlow, *mNorm, xBorder, yBorder ) );
		flowGrad.append( meanGradMag( *mf ) );
	}
}


/// register histology with MRI, assuming MRI already registered with block-face images 
/// (so histology registration is handle on a slice-by-slice basis)
void registerHistology( Config &conf ) {

	// get command parameters
	// blockOffset is index of blockface image corresponding to slice #1 (not necessarily the first histology image)
	int blockOffset = conf.readInt( "blockOffset" );
	int histoStep = conf.readInt( "histoStep", 1 );
	String histoPath = addDataPath( conf.readString( "histoPath", "histo/split" ) );
	String blockPath = addDataPath( conf.readString( "blockPath", "blockface/crop" ) );
	String mrPath = addDataPath( conf.readString( "mrPath", "mri/reg" ) );
	String outputPath = addDataPath( conf.readString( "outputPath", "histo/reg" ) );

	// read file names from input paths
	Array<String> hFileList = dirFileList( histoPath, "", ".png" );
	Array<String> bFileList = dirFileList( blockPath, "", ".png" );
	Array<String> mFileList = dirFileList( mrPath, "", ".png" );

	// check input lists
	disp( 1, "histo image count: %d, block image count: %d, MR image count: %d", hFileList.count(), bFileList.count(), mFileList.count() );
	if (hFileList.count() == 0) {
		warning( "no histo images found" );
		return;
	}
	if (bFileList.count() != mFileList.count()) {
		warning( "expected same number of block images and MR images" );
		return;
	}

	// create output directory
	createDir( outputPath );

	// prepare visualization path
	String visPath = outputPath + "/vis";
	createDir( visPath );

	// store block offset
	File blockOffsetFile( outputPath + "/blockOffset.txt", FILE_WRITE, FILE_TEXT );
	if (blockOffsetFile.openSuccess()) {
		blockOffsetFile.writeF( "%d\n", blockOffset );
	} else {
		warning( "unable to write to output path: %s", outputPath.c_str() );
		return;
	}

	// diagnostic data
	OutputVideo outputVideo;
	VectorD mutInfo, normDiff, flowGrad;

	// loop over histo input files
	for (int hIndex = 0; hIndex < hFileList.count(); hIndex += histoStep) {

		// register this slice
		registerHistologySlice( hIndex, blockOffset, histoPath, blockPath, mrPath, outputPath, 
								hFileList, bFileList, mFileList, outputVideo, mutInfo, normDiff, flowGrad );

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initHistoRegister() {
	registerCommand( "hreg", registerHistology );
}


} // end namespace hb
