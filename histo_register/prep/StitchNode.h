#ifndef _STITCH_NODE_H_
#define _STITCH_NODE_H_
#include <sbl/core/Pointer.h>
#include <sbl/core/Array.h>
#include <sbl/image/Image.h>
#include <sbl/image/ImageTransform.h>
using namespace sbl;
namespace hb {


//-------------------------------------------
// STITCH NODE CLASS
//-------------------------------------------


/// The StitchNode class represents a single input image to be stitched together into a larger image
class StitchNode {
public:

	// basic constructor; takes ownership of image
	StitchNode( aptr<ImageGrayU> image, int xImageIndex, int yImageIndex, const String &fileName );

	/// access the (small) image for this node
	inline const ImageGrayU &image() const { return *m_image; }

	/// access the current transform for this node (in small image coords)
	inline ImageTransform &transform() { return *m_transform; }
	inline const ImageTransform &transform() const { return *m_transform; }

	/// the x and y indices within this slide
	inline int xImageIndex() const { return m_xImageIndex; }
	inline int yImageIndex() const { return m_yImageIndex; }

	/// the file name for the large 
	inline const String &fileName() const { return m_fileName; }
	
	/// the brightness mapping currently specified for this image
	inline float minBrightness() const { return m_minBrightness; }
	inline float maxBrightness() const { return m_maxBrightness; }

	/// the brightness mapping for this image
	inline void setBrightnessBounds( float minBrightness, float maxBrightness ) { m_minBrightness = minBrightness; m_maxBrightness = maxBrightness; }

	/// map a pixel value according to the brightness bounds
	inline float mapBrightness( int v ) const { return ((float) v - m_minBrightness) * 255.0f / (m_maxBrightness - m_minBrightness); }

private:

	// internal data
	aptr<ImageGrayU> m_image;
	aptr<ImageTransform> m_transform;
	int m_xImageIndex;
	int m_yImageIndex;
	String m_fileName;
	float m_minBrightness;
	float m_maxBrightness;
};


//-------------------------------------------
// STITCH NODE SET CLASS
//-------------------------------------------


/// The StitchNodeSet class represents a set of images to be stitched together into a single larger image
class StitchNodeSet {
public:

	// basic constructor
	StitchNodeSet( int xCount, int yCount ) {
		m_xCount = xCount;
		m_yCount = yCount;
	}

	/// append node to the set
	inline void append( StitchNode *node ) { m_nodes.append( node ); }

	/// the number of noxes in the x direction
	inline int xCount() const { return m_xCount; }

	/// the number of nodes in the y direction
	inline int yCount() const { return m_yCount; }

	/// the total number of nodes
	inline int count() const { return m_nodes.count(); }

	/// access a node by index
	inline StitchNode &node( int index ) { return m_nodes[ index ]; }

	/// find a node given a node's x, y indices
	StitchNode &node( int xImageIndex, int yImageIndex );

private:

	// internal data
	Array<StitchNode> m_nodes;
	int m_xCount;
	int m_yCount;
};


} // end namespace hb
#endif // _STITCH_NODE_H_
