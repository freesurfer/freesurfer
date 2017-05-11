#pragma once

#include <stdexcept>

#include "cudaimage.hpp"

namespace kvl {
  namespace cuda {
    template<typename T,typename Internal>
    void RunVisitCounterSimpleCUDA( CudaImage<int,3,unsigned short>& d_output,
				    const CudaImage<T,3,size_t>& d_tetrahedra ) {
      throw std::runtime_error("Must call RunVisitCounterSimpleCUDA with float or double");
    }

    template<>
    void RunVisitCounterSimpleCUDA<float,float>( CudaImage<int,3,unsigned short>& d_output,
						 const CudaImage<float,3,size_t>& d_tetrahedra );

    template<>
    void RunVisitCounterSimpleCUDA<double,double>( CudaImage<int,3,unsigned short>& d_output,
						   const CudaImage<double,3,size_t>& d_tetrahedra );

    template<>
    void RunVisitCounterSimpleCUDA<float,double>( CudaImage<int,3,unsigned short>& d_output,
						  const CudaImage<float,3,size_t>& d_tetrahedra );
  }
}
