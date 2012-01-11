#include <sbl/core/Command.h>
#include <sbl/core/Init.h>
#include <sbl/core/PathConfig.h>
#include "prep/BlockPrep.h"
#include "prep/HistoStats.h"
#include "prep/VolumeFile.h"
#include "prep/VolumeUtil.h"
#include "prep/MPrep.h"
#include "prep/HistoPrep.h"
#include "prep/HistoStitch.h"
#include "prep/Polarization.h"
#include "registration/BlockRegister.h"
#include "registration/HistoRegister.h"
#include "registration/HistoTransform.h"
#include "registration/TestCorres3D.h"
using namespace sbl;
using namespace hb;


// initialize modules specific to this program
void initHistoRegisterModules() {

	// init preparation modules
	initBlockPrep();
	initHistoPrep();
	initHistoStats();
	initVolumeFile();
	initVolumeUtil();
	initMPrep();
	initHistoStitch();
	initPolarization();

	// init registration modules
	initBlockRegister();
	initHistoRegister();
	initHistoTransform();
	initTestCorres3D();
}


// program entry point
int main( int argc, char *argv[] ) {

	// register commands, etc.
	sbl::initModules();
	initHistoRegisterModules();

	// check for path as first arg
	int firstArgPos = 1;
	if (argc > firstArgPos) {
		String firstArg = argv[ firstArgPos ];
		if (firstArg.contains( '/' ) || firstArg.contains( '\\' ) || firstArg.contains( '.' )) {
			setDataPath( firstArg );
			firstArgPos++;
		}
	}

	// display current data path (from command line or path.conf)
	disp( 0, "dataPath: %s", dataPath().c_str() );

	// if command-line command, run it
	if (argc > firstArgPos) {
		String command = argv[ firstArgPos ];
		for (int i = firstArgPos + 1; i < argc; i++) {
			command += " ";
			command += argv[ i ];
		}
		execCommand( command, false ); 

	// otherwise run default start-up script
	} else {
    execCommand( "cmdlist", false );
		//runScript( dataPath() + "histo_register_script.txt" );
	}

	// clean up
	runCleanUp();
	return 0;
}
