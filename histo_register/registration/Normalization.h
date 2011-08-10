#ifndef _NORMALIZATION_H_
#define _NORMALIZATION_H_
#include <sbl/core/Pointer.h>
#include <sbl/image/Image.h>
using namespace sbl;
namespace hb {


/*! \file Normalization.cc
    \brief The Normalization module is used to generate transform different imaging 
	modalities into a common form.
*/


/// compute a entropy image in which each pixel is set according to the amount of local entropy;
/// see: C. Wachinger and N. Navab, "Structural image representation for image registration," CVPRW 2010
aptr<ImageGrayU> entropyImage( const ImageGrayU &image );


/// compute a normalized image by combining a blurred mask with an entropy image
aptr<ImageGrayU> normalize( const ImageGrayU &image, const ImageGrayU &mask, const ImageGrayU &maskBlurred );


} // end namespace hb
#endif // _NORMALIZATION_H_
