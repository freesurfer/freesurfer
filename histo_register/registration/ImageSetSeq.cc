#include "registration/ImageSetSeq.h"
#include <sbl/math/MathUtil.h>
namespace hb {


/// compute the value of channel k at the given coordinates using linear interpolation
float ImageSetSeq::interp( int k, float x, float y, float z ) const {
	int zInt = (int) z;
	float frac = z - (float) zInt;
	assertAlways( zInt >= 0 && zInt < m_images.count() - 1 );
	assertAlways( frac >= -0.00001 && frac <= 1.00001 );
	float v1 = m_images[ zInt ][ k ].interp( x, y );
	float v2 = m_images[ zInt + 1 ][ k ].interp( x, y );
	return frac * v2 + (1.0f - frac) * v1; 
}


} // end namespace hb
