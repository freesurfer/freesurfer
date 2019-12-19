#include "atlasmeshalphadrawercpuwrapper.hpp"

namespace kvl {
  void AtlasMeshAlphaDrawerCPUWrapper::SetRegions( const ImageType::RegionType&  region ) {
    this->tSetRegions.Start();
    this->impl->SetRegions(region);
    this->tSetRegions.Stop();
  }

  void AtlasMeshAlphaDrawerCPUWrapper::Interpolate( const kvl::AtlasMesh* mesh ) {
    this->tInterpolate.Start();
    this->impl->Rasterize(mesh);
    this->tInterpolate.Stop();
  }

  const AtlasMeshAlphaDrawerCPUWrapper::ImageType* AtlasMeshAlphaDrawerCPUWrapper::GetImage() const {
    return this->impl->GetImage();
  }
  
  void AtlasMeshAlphaDrawerCPUWrapper::SetClassNumber( const int classNumber ) {
    this->impl->SetClassNumber( classNumber );
  }
}
