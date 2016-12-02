#ifndef kvlTetrahedronInteriorIterator_hxx
#define kvlTetrahedronInteriorIterator_hxx

#include "kvlTetrahedronInteriorIterator.h"

namespace kvl
{


//
//
//  
template< typename TPixel >
TetrahedronInteriorIterator< TPixel >
::TetrahedronInteriorIterator( ImageType *ptr,
                               const PointType& p0, 
                               const PointType& p1, 
                               const PointType& p2, 
                               const PointType& p3 )
: TetrahedronInteriorConstIterator< TPixel >( ptr, p0, p1, p2, p3 )
{
  
}


} // end namespace kvl

#endif
