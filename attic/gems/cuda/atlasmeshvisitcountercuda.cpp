#include <stdexcept>

#include "atlasmeshvisitcountercuda.hpp"

namespace kvl {
  namespace cuda {
    void AtlasMeshVisitCounterCUDA::SetRegions( const kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType& region ) {
      throw std::runtime_error("Not implemented");
    }

    void AtlasMeshVisitCounterCUDA::VisitCount( const kvl::AtlasMesh* mesh ) {
      throw std::runtime_error("Not implemented");
    }

    const AtlasMeshVisitCounterCUDA::ImageType* AtlasMeshVisitCounterCUDA::GetImage() const {
      throw std::runtime_error("Not implemented");
    }
  }
}
