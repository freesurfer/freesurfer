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
// fix(clean): create non-command version and expose it here
void convertMghFileToImages( Config &conf );


/// convert a set of images to an mgz file
void convertImagesToMghFile( const String &inputPath, const String &outputFileName, float xSpacing, float ySpacing, float zSpacing, int width = 0, int height = 0 );


} // end namespace hb
#endif // _VOLUME_FILE_H_
