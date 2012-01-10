#ifndef _VOLUME_FILE_H_
#define _VOLUME_FILE_H_
#include <sbl/core/Config.h>
using namespace sbl;
namespace hb {


/*! \file VolumeFile.cc
	\brief The VolumeFile module provides functions for reading and writing .mgz files.
*/


// register commands, etc. defined in this module
void initVolumeFile();


/// convert an mgh/mgz file to a set of images
void convertMghFileToImages( const String &fileName, const String &outputPath, 
							 bool xySwap, bool xzSwap, bool yzSwap, 
							 bool xFlip, bool yFlip, bool zFlip, bool autoScaleValues );


/// convert a set of images to an mgz file
void convertImagesToMghFile( const String &inputPath, const String &outputFileName, 
							 float xSpacing, float ySpacing, float zSpacing, 
							 int width = 0, int height = 0, int valueScaleFactor = 0 );


/// convert an mgh/mgz file to a set of images
// fix(clean): remove this version from header file
void convertMghFileToImages( Config &conf );


} // end namespace hb
#endif // _VOLUME_FILE_H_
