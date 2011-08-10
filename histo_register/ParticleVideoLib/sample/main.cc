#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <sbl/image/Video.h>
#include <sbl/image/MotionFieldUtil.h>
#include <sbl/system/FileSystem.h>
#include <pvl/VarMotion.h>
#include <pvl/SimpleParticleBuild.h>
using namespace sbl;
using namespace pvl;


// estimate optical flow between frames of an input video
void estimateOpticalFlow( int argc, char *argv[] ) {

	// get command parameters
	String inputVideoFileName = addDataPath( argv[ 2 ] );
	int frameCount = String( argv[ 3 ] ).toInt();
	bool saveText = true;
	String varMotionConfFileName = dataPath() + "varMotion.conf";

	// note: processing frameCount frames will produce frameCount - 1 optical flow fields

	// open input video
	InputVideo inputVideo( inputVideoFileName );
	if (inputVideo.openSuccess() == false) {
		warning( "unable to open input video: %s", inputVideoFileName.c_str() );
		return;
	}
	if (frameCount > inputVideo.length())
		frameCount = inputVideo.length();

	// open output video
	String outputVideoFileName = inputVideoFileName.leftOfLast( '.' ) + "_mfVis.avi";
	OutputVideo outputVideo( outputVideoFileName, inputVideo.width(), inputVideo.height(), 6 );
	if (outputVideo.openSuccess() == false) {
		warning( "unable to open output video: %s", outputVideoFileName.c_str() );
		return;
	}

	// create output path
	String outputPath = inputVideoFileName.leftOfLast( '.' ) + "_motion";
	createDir( outputPath );

	// keep track of the previous video frame
	aptr<ImageColorU> lastFrame;

	// loop over input video
	for (int frameIndex = 0; frameIndex < frameCount; frameIndex++) {

		// load video frame
		aptr<ImageColorU> frame = inputVideo.frame( frameIndex );

		// estimate optical flow
		if (frame.get() && lastFrame.get()) {
			aptr<MotionField> mf = varMotion( *lastFrame, *frame, varMotionConfFileName, NULL, NULL );

			// save estimated motion field
			String outputFileName = outputPath + sprintF( "/motion.%05d.%05d.mf", frameIndex - 1, frameIndex );
			disp( 1, "output: %s", outputFileName.c_str() );
			saveMotionField( *mf, outputFileName );
			if (saveText) {
				String textOutputFileName = outputPath + sprintF( "/motion.%05d.%05d.txt", frameIndex - 1, frameIndex );
				saveMotionFieldText( *mf, textOutputFileName );
			}

			// visualize motion
			aptr<ImageColorU> vis = colorizeMotion( *mf );
			outputVideo.append( *vis );

			// display diagnostic info
			dispMotionStats( 1, *mf );
		}

		// keep frame for next iter
		lastFrame = frame;

		// check for user cancel
		if (checkCommandEvents())
			break;
	}
}


// construct particle video using frames of an input video
void buildParticleVideo( int argc, char *argv[] ) {

	// get command parameters
	String inputVideoFileName = addDataPath( argv[ 2 ] );
	int frameCount = String( argv[ 3 ] ).toInt();
	bool doClustering = false;

	// load particle config
	String spConfFileName = dataPath() + "simpleParticle.conf";
	Config spConf;
	spConf.load( spConfFileName );
	if (spConf.entryCount() == 0) {
		warning( "unable to load simple particle config: %s", spConfFileName.c_str() );
		return;
	}

	// open input video
	InputVideo inputVideo( inputVideoFileName );
	if (inputVideo.openSuccess() == false) {
		warning( "unable to open input video: %s", inputVideoFileName.c_str() );
		return;
	}
	if (frameCount > inputVideo.length())
		frameCount = inputVideo.length();

	// open output video
	String outputVideoFileName = inputVideoFileName.leftOfLast( '.' ) + "_spVis.avi";
	OutputVideo outputVideo( outputVideoFileName, inputVideo.width(), inputVideo.height(), 6 );
	if (outputVideo.openSuccess() == false) {
		warning( "unable to open output video: %s", outputVideoFileName.c_str() );
		return;
	}

	// create particle video object
	SimpleParticleSet particleSet;

	// path containing optical flow data
	String motionPath = inputVideoFileName.leftOfLast( '.' ) + "_motion";

	// loop over input frames
	for (int frameIndex = 0; frameIndex < frameCount; frameIndex++) {

		// get video frame
		aptr<ImageColorU> frame = inputVideo.frame( frameIndex );

		// load motion fields
		String motionFileName = motionPath + sprintF( "/motion.%05d.%05d.mf", frameIndex, frameIndex + 1 );
		aptr<MotionField> mf = loadMotionField( motionFileName );
		if (mf.get() == NULL) {
			warning( "unable to load motion field: %s", motionFileName.c_str() );
			return;
		}
		motionFileName = motionPath + sprintF( "/motion.%05d.%05d.mf", frameIndex - 1, frameIndex );
		aptr<MotionField> mfPrev = loadMotionField( motionFileName );

		// update particles
		if (frame.get()) {
			simpleParticleFrame( particleSet, frameIndex, *frame, mfPrev.get(), *mf, spConf );

			// draw particles on current frame
			aptr<ImageColorU> vis = particleSet.draw( frameIndex, *frame );
			outputVideo.append( *vis );
		}

		// check for user cancel
		if (checkCommandEvents())
			break;
	}

	// if requested, perform clustering
	if (doClustering) {
		buildParticleClusters( particleSet, inputVideo.width(), inputVideo.height(), spConf, true );
	}

	// save particle set
	String particleSetFileName = inputVideoFileName.leftOfLast( '.' ) + ".sps";
	particleSet.saveSimple( particleSetFileName );

	// display diagnostics
	particleSet.dispStats( 1 );
}


// display program argument usage
void dispUsage( const char *programName ) {
	disp( 0, "usage:" );
	disp( 1, "%s flow [videoFileName] [frameCount]", programName );
	disp( 1, "%s pv [videoFileName] [frameCount]", programName );
	disp( 0, "notes:" );
	disp( 1, "- program path should include path.conf that specifies dataPath" );
	disp( 1, "- input files should be in dataPath" );
	disp( 1, "- simpleParticle.conf and varMotion.conf should be in dataPath" );
	disp( 1, "- pv command expects flow to be estimated first" );
}


// program entry point
int main( int argc, char *argv[] ) {

	// check args
	if (argc != 4) {
		dispUsage( argv[ 0 ] );
		return 1;
	}

	// start timing
	Timer timer;
	timer.start();

	// execute desired command
	String commandName = argv[ 1 ];
	if (commandName == "flow")
		estimateOpticalFlow( argc, argv );
	else if (commandName == "pv") 
		buildParticleVideo( argc, argv );
	else
		dispUsage( argv[ 0 ] );

	// display total run time
	timer.stop();
	disp( 1, "total time: %f seconds", timer.timeSum() );

	// clean up
	runCleanUp();
	return 0;
}
