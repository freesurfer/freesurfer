#include "atlasmeshvisitcountercpuwrapper.hpp"

namespace kvl {
  void AtlasMeshVisitCounterCPUWrapper::SetRegions( const ImageType::RegionType&  region ) {
    this->tSetRegions.Start();
    this->impl->SetRegions(region);
    this->tSetRegions.Stop();
  }

  void AtlasMeshVisitCounterCPUWrapper::VisitCount( const AtlasMesh* mesh ) {
    this->tVisitCount.Start();
    this->impl->Rasterize(mesh);
    this->tVisitCount.Stop();
  }

  const AtlasMeshVisitCounterCPUWrapper::ImageType* AtlasMeshVisitCounterCPUWrapper::GetImage() const {
    return this->impl->GetImage();
  }
}
