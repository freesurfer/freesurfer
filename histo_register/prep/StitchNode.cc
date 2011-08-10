#include "prep/StitchNode.h"
namespace hb {


//-------------------------------------------
// STITCH NODE CLASS
//-------------------------------------------


// basic constructor; takes ownership of image
StitchNode::StitchNode( aptr<ImageGrayU> image, int xImageIndex, int yImageIndex, const String &fileName ) {
	m_image = image;
	m_xImageIndex = xImageIndex;
	m_yImageIndex = yImageIndex;
	m_fileName = fileName;
	m_transform.reset( new ImageTransform( 0, 0 ) );
	m_minBrightness = 0;
	m_maxBrightness = 255;
}


//-------------------------------------------
// STITCH NODE SET CLASS
//-------------------------------------------


/// find a node given a node's x, y indices
StitchNode &StitchNodeSet::node( int xImageIndex, int yImageIndex ) {
	int nodeIndex = xImageIndex * m_yCount + yImageIndex;
	assertAlways( nodeIndex < m_nodes.count() );
	StitchNode &node = m_nodes[ nodeIndex ];
	assertAlways( node.xImageIndex() == xImageIndex && node.yImageIndex() == yImageIndex );
	return node;
}


} // end namespace hb
