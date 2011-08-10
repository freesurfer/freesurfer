#ifndef _SBL_FILTER_H_
#define _SBL_FILTER_H_
#include <sbl/core/Config.h>
#include <sbl/core/Pointer.h>
#include <sbl/image/Image.h>
namespace sbl {


/*! \file Filter.h
    \brief The Filter module represents generalized image filters (with a uniform function
    specification that allows run-time selection of filters).
*/


// register commands, etc. defined in this module
void initFilter();


// type definitions for generalized image filters
typedef aptr<ImageColorU> (*ColorImageFilterCallback)( const ImageColorU &input, Config &conf );
typedef aptr<ImageGrayU> (*GrayImageFilterCallback)( const ImageGrayU &input, Config &conf );


/// register an image processing filter 
#define registerFilter( filter ) registerFilterInternal( #filter, filter ) 
void registerFilterInternal( const String &name, ColorImageFilterCallback callback );
void registerFilterInternal( const String &name, GrayImageFilterCallback callback );


} // end namespace sbl
#endif // _SBL_FILTER_H_

