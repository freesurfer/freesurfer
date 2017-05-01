#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"



namespace kvl {
  namespace cuda {
    template<>
    void RunVisitCounterSimpleCUDA( CudaImage<int,3,unsigned short>& d_output,
				    const CudaImage<float,3,size_t>& d_tetrahedra ) {
      throw std::runtime_error("RunVisitCounterSimpleCUDA not yet implemented for float");
    }
  }
}
