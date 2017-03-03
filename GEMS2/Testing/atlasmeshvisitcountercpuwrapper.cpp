#include "atlasmeshvisitcountercpuwrapper.hpp"

namespace kvl {
  void AtlasMeshVisitCounterCPUWrapper::SetRegions( const ImageType::RegionType&  region ) {
    this->impl->SetRegions(region);
  }

  void AtlasMeshVisitCounterCPUWrapper::VisitCount( const AtlasMesh* mesh ) {
    this->impl->Rasterize(mesh);
  }

  const AtlasMeshVisitCounterCPUWrapper::ImageType* AtlasMeshVisitCounterCPUWrapper::GetImage() const {
    return this->impl->GetImage();
  }
}
