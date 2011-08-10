#ifndef _IMAGE_SET_SEQ_H_
#define _IMAGE_SET_SEQ_H_
#include <sbl/core/Array.h>
#include <sbl/image/Image.h>
using namespace sbl;
namespace hb {


/// The ImageSetSeq class represents a asequence of image sets (e.g. a 3D volume of multi-channel images)
class ImageSetSeq {
public:

	/// create an empty sequence of image sets
	ImageSetSeq() {}

	/// the number of image sets in the sequence
	inline int count() const { return m_images.count(); }

	/// access an element (image set) in the sequence by index
	const Array<ImageGrayF> &operator[]( int index ) const { return m_images[ index ]; }

	/// append a set of images (e.g. a multi-channel image) to sequence
	void append( Array<ImageGrayF> *imgSet ) { m_images.append( imgSet ); }

	/// compute the value of channel k at the given coordinates using linear interpolation
	float interp( int k, float x, float y, float z ) const;

private:

	// the images in the sequence of sets
	Array<Array<ImageGrayF> > m_images;

	// disable copy constructor and assignment operator
	ImageSetSeq( const ImageSetSeq &x );
	ImageSetSeq &operator=( const ImageSetSeq &x );
};


} // end namespace hb
#endif // _IMAGE_SET_SEQ_H_
