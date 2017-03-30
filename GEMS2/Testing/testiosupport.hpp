#pragma once

#include <iostream>

#include "dimensioncuda.hpp"

namespace std {
   template<unsigned char nDims, typename IndexType>
   ostream& operator<<(ostream& os,
		       const kvl::cuda::Dimension<nDims,IndexType> d) {
     os << "( ";
     
     os << d[0];
     
     for( unsigned char i=1; i<nDims; i++ ) {
       os << ", " << d[i];
     }
     os << " )";
     return os;
   }
}
