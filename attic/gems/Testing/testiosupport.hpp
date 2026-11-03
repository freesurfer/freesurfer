#pragma once

#include <iostream>
#include <iomanip>

#ifdef CUDA_FOUND
#include "dimensioncuda.hpp"
#endif
#include "stopwatch.hpp"

namespace std {
#ifdef CUDA_FOUND
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
#endif

  static ostream& operator<<(ostream& os,
			     const kvl::Stopwatch s) {
    os << "Elapsed ";
    os << setprecision(2) << setw(8) << fixed << s.GetElapsedTime() << " ms";
    os << " Avg ";
    os << setprecision(2) << setw(8) << fixed << s.GetAverageElapsedTime() << " ms";

    return os;
  }
}
